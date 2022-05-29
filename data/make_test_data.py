#!/usr/bin/env python3
import argparse
import itertools
import multiprocessing
import os
import pathlib
import re
import json
import subprocess
import tempfile
import numpy as np
from dataclasses import dataclass, field
from random import choices

import cyvcf2
import tskit


def _get_metadata(tablerow):
    md_string = tablerow.metadata.decode()
    if md_string:
        md = json.loads(tablerow.metadata.decode())
    else:
        md = {}
    return md



def partition_individuals(individuals: list):
    """Partition individuals into reference individual and samples"""

    def _filter_reference(ind):
        return ind.metadata["is_reference"]

    return list(filter(_filter_reference, individuals)), list(
        itertools.filterfalse(_filter_reference, individuals)
    )


def write_variants(ref_ind, individuals, vcf_tmp_fn, outdir, prefix):
    """Write variants to output vcf.

    If reference vcf is provided flip alleles where necessary
    """
    vcf_ref = None
    if len(ref_ind) > 0:
        vcf_ref = cyvcf2.VCF(vcf_tmp_fn, samples=ref_ind[0].metadata["tskit_id"])

    vcf_out = outdir / f"{prefix}.vcf.gz"
    print(f"Writing output vcf {vcf_out}")
    samples = [ind.metadata["tskit_id"] for ind in individuals]
    vcf_tmp = cyvcf2.VCF(vcf_tmp_fn, samples=samples)
    vcfwriter = cyvcf2.cyvcf2.Writer(str(vcf_out), vcf_tmp, "wz")
    for ind in individuals:
        vcfwriter.add_to_header(ind.metadata["vcfheader"])
    if vcf_ref is None:
        vcfwriter.set_samples(vcf_tmp.samples)
        for variant in vcf_tmp:
            vcfwriter.write_record(variant)
    else:
        # in cyvcf2: Pick one individual *haplotype* as reference: this
        # individual should have only 0's, so all calls at a site with a
        # derived allele should be flipped for all individuals.
        for gt, variant in zip(vcf_ref, vcf_tmp):
            ref = gt.genotypes[0][0]
            if ref == 1:
                # Flip states for all other genotypes
                for i in range(len(variant.genotypes)):
                    alleles = variant.genotypes[i]
                    variant.genotypes[i] = [1 - alleles[0], 1 - alleles[1], alleles[2]]
                variant.genotypes = variant.genotypes
                variant.REF = gt.ALT[0]
                variant.ALT = [gt.REF]
            vcfwriter.write_record(variant)
    vcfwriter.close()
    subprocess.run(["tabix", vcf_out], check=True)




@dataclass
class SeqRecord:
    haplotype: int = -1
    node: int = -1
    population: str = None

    def __init__(self, *, seq, name, id):
        self.seq = seq
        self.name = name
        self.id = id

    @property
    def seqid(self):
        return f"{self.population}-{self.id}-{self.haplotype}"

    @property
    def description(self):
        return (
            f"id:{self.id}, node:{self.node}, name:{self.name}, "
            f"haplotype:{self.haplotype}, "
            f"population:{self.population}"
        )

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        s = f">{self.id} {self.description}\n"
        n = 72
        chunks = [self.seq[i : i + n] for i in range(0, len(self), n)]
        s += "\n".join(chunks)
        s += "\n"
        return s

    def __repr__(self):
        return str(self)




# The node is what makes this unique
def make_sequence_from_tree_sequence(ts, individual, node, haplotype, reference):
    """Create DNA sequence with variants at sites"""
    seqid = f"tsk_{individual.metadata['population'].metadata['name']}-{individual.id}-{haplotype}"
    for site, variant in zip(ts.sites(), list(ts.haplotypes())[node]):
        i = int(site.position) - 1
        reference[i] = variant
    record = SeqRecord(
        seq="".join(reference),
        name=individual.metadata["name"],
        id=seqid,
    )
    record.population = individual.metadata['population'].metadata['name']
    record.haplotype = haplotype

    return record


def make_reads(outdir, **kw):
    kwargs = []
    for k, v in kw.items():
        kwargs.append(f"--{k}")
        if isinstance(v, list):
            for item in v:
                kwargs.append(str(item))
        else:
            kwargs.append(str(v))
    args = [
        "iss",
        "generate",
        "-z",
        "--model",
        "HiSeq",
        "--output",
        str(outdir / "reads"),
    ] + kwargs
    print(args)
    subprocess.run(args, check=True)


def make_ooa_old(args):
    TARGET = "ooa/ooa.reference.fasta"
    if os.path.exists(TARGET):
        if not args.force:
            print(f"Target {TARGET} exists; skipping. Force remake with --force flag")
            return
    outdir = pathlib.Path("ooa")
    prefix = "ooa"
    SEQLENGTH = 1e5
    RECOMBINATION_RATE = 1e-8
    MUTATION_RATE = 1e-9
    RANDOM_SEED = 10
    REFERENCE_CHROMOSOME = "CHB:0"

    ts = tskit.load(outdir / "ooa.mut.ts")
    individuals = update_individual_metadata(ts)
    individuals = set_reference_individual(REFERENCE_CHROMOSOME, individuals)
    ref_ind, individuals = partition_individuals(individuals)

    vcf_tmp_fn = tempfile.mkstemp(suffix=".vcf")[1]
    print(f"Writing temporary ts vcf {vcf_tmp_fn}")
    with open(vcf_tmp_fn, "w") as fh:
        ts.write_vcf(fh)

    write_variants(ref_ind, individuals, vcf_tmp_fn, outdir, "ooa_1_PASS")

    print("Writing fasta sequences")
    # Make fasta sequences
    reference = None
    dna = ["A", "C", "G", "T"]
    reference = choices(dna, k=int(ts.sequence_length))
    # Reference individual; only print first node
    if len(ref_ind) > 0:
        ind = ref_ind[0]
        node = ind.nodes[0]
        haplotype = node % len(ind.nodes)
        rec = make_sequence(ts, ind, node, haplotype, reference)
        ref_out = outdir / f"{prefix}.reference.fasta"
        with open(ref_out, "w") as fh:
            fh.write(repr(rec))
    for ind in individuals:
        indseq = []
        print(ind.metadata)
        seq_out = outdir / f"{ind.id}.fasta"
        nodes = ind.nodes
        for node in nodes:
            haplotype = node % len(nodes)
            rec = make_sequence(ts, ind, node, haplotype, reference)
            seq_out = outdir / f"{rec.id}.fasta"
            indseq.append(seq_out)
            with open(seq_out, "w") as fh:
                fh.write(repr(rec))
        make_reads(outdir, genomes=indseq, cpus=args.cpus, n_reads=100)
        for seq in indseq:
            os.unlink(seq)

@dataclass
class Individual:
    uniqueid: int
    id: int
    population: int
    population_name: str = None
    metadata: dict = None
    paternal_chromosome: SeqRecord = None
    maternal_chromosome: SeqRecord = None
    is_reference: bool = False

    def __str__(self):
        return f"{self.population_name}-{self.id}"

    def __repr__(self):
        return (f"Individual(uniqueid={self.uniqueid}, id={self.id}, "
                f"population={self.population}, population_name={self.population_name}, "
                f"metadata={self.metadata})"
                )

    def haplotype(self, which=0):
        if which == 0:
            return self.maternal_chromosome
        elif which == 1:
            return self.paternal_chromosome
        else:
            print("no such haplotype {which}")

    @property
    def tskit_id(self):
        return f"tsk_{self.uniqueid}"

    def simulate_reads(self, coverage):
        pass


@dataclass
class Population:
    id: int
    name: str
    metadata: dict = field(default_factory = dict)
    individuals: list[Individual] = field(default_factory = list)

    @property
    def size(self):
        return len(self.individuals)

    def get_individual(self, index: int):
        return self.individuals[index]

    def set_reference(self, index: int):
        self.individuals[index].is_reference = True


@dataclass
class DemesModel:
    name: str
    path: pathlib.Path
    demesfile: str
    ts: tskit.trees.TreeSequence = None
    populations: list[Population] = field(default_factory=list)
    individuals: list[Individual] = field(default_factory=list)

    def __post_init__(self):
        self._census = 0

    @property
    def census(self):
        return self._census

    @census.setter
    def census(self, census):
        self._census = census

    def add_individuals(self, **kw):
        for pop, count in kw.items():
            if pop in self.popdict.keys():
                population = self.get_population(pop)
            else:
                population = Population(id=len(self.populations), name=pop)
                self.populations.append(population)
            popindex = population.id
            popsize = len(self.populations[popindex].individuals)
            for i in np.arange(popsize, popsize + count):
                ind = Individual(uniqueid=self.census, id=i, population=popindex, population_name=pop)
                self.populations[popindex].individuals.append(ind)
                self.individuals.append(ind)
                self.census += 1

    @property
    def reference(self):
        # loop samples
        for ind in self.individuals:
            if ind.is_reference:
                return ind

    def set_reference(self, individual: int, population: str):
        for pop in self.populations:
            if pop.name == population:
                pop.set_reference(individual)

    def get_population(self, name):
        for pop in self.populations:
            if pop.name == name:
                return pop

    @property
    def popdict(self):
        return dict((p.name, p.size) for p in self.populations)

    def make_tree_sequence(self):
        pass

    def update_ts_metadata(self):
        pass

    # FIXME: should be agnostic to msprime / SLiM
    def simulate(self, *, seqlength, recombination_rate,
                 mutation_rate, anc_seed=42, mut_seed=10):
        tsfile = self._ts_ancestry(seqlength, recombination_rate, anc_seed)
        tsfile = self._ts_mutate(tsfile, mutation_rate, mut_seed)
        self.ts = tskit.load(tsfile)

    def _ts_ancestry(self, seqlength, recombination_rate, seed):
        populations = [f"{k}:{v}" for k, v in self.popdict.items()]
        tsfile = self.path / f"{self.name}.ts"
        args = ["msp", "ancestry",
                "--demography", self.demesfile,
                "--length", str(seqlength),
                "--recombination-rate", str(recombination_rate),
                "--random-seed", str(seed),
                "-o", tsfile,
            ] + populations

        # Run msp ancestry
        subprocess.run(args, check=True)
        return tsfile

    def _ts_mutate(self, tsfile, mutation_rate, seed):
        tsmutfile = self.path / f"{self.name}.mut.ts"
        args = ["msp", "mutations",
                "-o", tsmutfile,
                "--random-seed", str(seed),
                str(mutation_rate), tsfile]
        subprocess.run(args, check=True)
        os.unlink(tsfile)
        return tsmutfile


    # FIXME: should sync back to Individual?
    def sync_metadata(self):
        """Sync tree sequence and DemesModel metadata"""
        tc = self.ts.dump_tables()
        indrows = list(getattr(tc, "individuals"))
        getattr(tc, "individuals").reset()
        oldname = None
        for i, row in enumerate(indrows):
            ind = self.ts.individuals()[i]
            popindex = self.ts.nodes()[ind.nodes[0]].population
            pop = self.ts.population(popindex)
            popname = pop.metadata["name"]
            if oldname is not None and oldname != popname:
                i = 0
            oldname = popname
            name = f"tsk_{ind.id}_{popname}_{i}"
            md = {
                "id": ind.id,
                "tskit_id": f"tsk_{ind.id}",
                "name": name,
                "popname": popname,
                "popindex": popindex,
                "ind_index": i,
                "is_reference": self.individuals[i].is_reference,
            }
            md["vcfheader"] = (
                f"##<SAMPLE=<ID={md['tskit_id']},Name={md['name']},"
                f"Index={md['id']},Population={md['popname']},"
                f"Description=\"{pop.metadata['description']}\">"
            )
            row = row.replace(metadata=json.dumps(md).encode())
            getattr(tc, "individuals").append(row)
        self.ts = tc.tree_sequence()


    def write_variants(self):
        """Write variants to output vcf.

        If reference vcf is provided flip alleles where necessary
        """
        vcf_tmp_fn = tempfile.mkstemp(suffix=".vcf")[1]
        print(f"Writing ts to vcf {vcf_tmp_fn}")
        with open(vcf_tmp_fn, "w") as fh:
            self.ts.write_vcf(fh)

        vcf_ref = None
        if self.reference is not None:
            vcf_ref = cyvcf2.VCF(vcf_tmp_fn, samples=self.reference.tskit_id)

        vcf_out = self.path / f"{self.name}.vcf.gz"
        print(f"Writing output vcf {vcf_out}")
        samples = [ind.tskit_id for ind in self.individuals]
        vcf_tmp = cyvcf2.VCF(vcf_tmp_fn, samples=samples)
        vcfwriter = cyvcf2.cyvcf2.Writer(str(vcf_out), vcf_tmp, "wz")
        for ind in self.individuals:
            tsk_ind = self.ts.individuals()[ind.uniqueid]
            vcfwriter.add_to_header(_get_metadata(tsk_ind)["vcfheader"])
        if vcf_ref is None:
            vcfwriter.set_samples(vcf_tmp.samples)
            for variant in vcf_tmp:
                vcfwriter.write_record(variant)
        else:
            # in cyvcf2: Pick one individual *haplotype* as reference: this
            # individual should have only 0's, so all calls at a site with a
            # derived allele should be flipped for all individuals.
            for gt, variant in zip(vcf_ref, vcf_tmp):
                ref = gt.genotypes[0][0]
                if ref == 1:
                    # Flip states for all other genotypes
                    for i in range(len(variant.genotypes)):
                        alleles = variant.genotypes[i]
                        variant.genotypes[i] = [1 - alleles[0], 1 - alleles[1], alleles[2]]
                    variant.genotypes = variant.genotypes
                    variant.REF = gt.ALT[0]
                    variant.ALT = [gt.REF]
                vcfwriter.write_record(variant)
        vcfwriter.close()
        subprocess.run(["tabix", vcf_out], check=True)


def make_ooa(args):
    model = DemesModel(name="ooa", path=pathlib.Path("ooa"), demesfile="ooa/ooa_with_outgroup.demes.yaml")
    popdict = {'CHB': 6, 'CEU': 7, 'YRI': 6, 'gorilla': 1,
               'chimpanzee': 1, 'orangutan': 1}
    model.add_individuals(**popdict)
    model.set_reference(population="CEU", individual=6)
    model.simulate(seqlength=1e5, recombination_rate=1e-8, mutation_rate=1e-9)
    model.sync_metadata()
    model.write_variants()


def main():
    parser = argparse.ArgumentParser(description="Make test data")
    parser.add_argument(
        "--force",
        "-F",
        help=("force dataset creation"),
        action="store_true",
        default=False,
    )

    subparsers = parser.add_subparsers(dest="subcommand")
    subparsers.required = True

    ooaparser = subparsers.add_parser("ooa", help="Make ooa dataset (neutral)")
    ooaparser.set_defaults(runner=make_ooa)
    ooaparser.add_argument(
        "--cpus", "-p", help="Set number of cpus", default=1, dest="cpus"
    )

    args = parser.parse_args()
    args.runner(args)


if __name__ == "__main__":
    main()
