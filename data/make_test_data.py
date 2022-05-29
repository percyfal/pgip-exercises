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
    haplotype: list[SeqRecord] = field(default_factory=list)
    is_reference: bool = False

    def __str__(self):
        return f"{self.population_name}-{self.id}"

    def __repr__(self):
        return (f"Individual(uniqueid={self.uniqueid}, id={self.id}, "
                f"population={self.population}, population_name={self.population_name}, "
                f"metadata={self.metadata})")

    @property
    def tskit_id(self):
        return f"tsk_{self.uniqueid}"

    def simulate_reads(self, coverage):
        pass

    @property
    def reference(self):
        try:
            return self.haplotype[0]
        except Excetption as e:
            print(e)
            pass

    def make_sequence(self, ts, reference=None):
        """Create DNA sequence with variants at sites"""
        if reference is None:
            return
        for haplotype, node in enumerate(ts.individual(self.uniqueid).nodes):
            if self.is_reference and haplotype == 0:
                next
            seqid = f"{self.tskit_id}-{self.population_name}-{self.id}-{haplotype}"
            for site, variant in zip(ts.sites(), list(ts.haplotypes())[node]):
                i = int(site.position) - 1
                reference[i] = variant
            record = SeqRecord(
                seq="".join(reference),
                name=seqid,
                id=seqid,
            )
            record.population = self.population_name
            record.haplotype = haplotype
            self.haplotype.append(record)

    def write_haplotypes(self, path):
        with open(path, "w") as fh:
            for hap in self.haplotype:
                fh.write(str(hap))

    def simulate_reads(self, path, *, coverage=10, readlength=125, **kw):
        nreads = int((len(self.haplotype[0]) * coverage) / (readlength * 2))
        tmp = tempfile.mkstemp(suffix=".fasta")[1]
        self.write_haplotypes(tmp)
        args = ["iss", "generate", "--output", path / f"{str(self)}",
                "--genomes", tmp, "-z", "--model", "HiSeq", "-n", str(nreads)]
        for k, v in kw.items():
            args += [f"--{k}"]
            args += [v]
        subprocess.run(args, check=True)
        os.unlink(tmp)



@dataclass
class Population:
    id: int
    name: str
    metadata: dict = field(default_factory = dict)
    individuals: list[Individual] = field(default_factory = list)

    @property
    def size(self):
        return len(self.individuals)



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
                pop.individuals[individual].is_reference = True

    def make_reference_sequence(self):
        if self.ts is None:
            return
        dna = ["A", "C", "G", "T"]
        if len(self.reference.haplotype) == 0:
            self.reference.haplotype.append(choices(dna, k=int(self.ts.sequence_length)))
        else:
            self.reference.haplotype[0] = choices(dna, k=int(self.ts.sequence_length))

    def get_population(self, name):
        for pop in self.populations:
            if pop.name == name:
                return pop

    @property
    def popdict(self):
        return dict((p.name, p.size) for p in self.populations)

    # FIXME: should be agnostic to msprime / SLiM
    def simulate(self, *, seqlength, recombination_rate,
                 mutation_rate, anc_seed=42, mut_seed=10):
        tsfile = self._ts_ancestry(seqlength, recombination_rate, anc_seed)
        tsfile = self._ts_mutate(tsfile, mutation_rate, mut_seed)
        self.ts = tskit.load(tsfile)
        os.unlink(tsfile)

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

    # FIXME: dump more than ts?
    def dump(self, **kw):
        path = self.path / f"{self.name}.ts"
        print("Dumping tree sequences to", str(path))
        self.ts.dump(path, **kw)

    def write_haplotypes(self):
        for ind in self.individuals:
            path = self.path / f"{str(ind)}"
            ind.write_haplotypes(path)


def make_ooa(args):
    model = DemesModel(name="ooa", path=pathlib.Path("ooa"), demesfile="ooa/ooa_with_outgroup.demes.yaml")
    popdict = {'CHB': 6, 'CEU': 7, 'YRI': 6, 'gorilla': 1,
               'chimpanzee': 1, 'orangutan': 1}
    model.add_individuals(**popdict)
    model.set_reference(population="CEU", individual=6)
    model.simulate(seqlength=1e5, recombination_rate=1e-8, mutation_rate=1e-9)
    model.sync_metadata()
    model.write_variants()
    model.make_reference_sequence()

    for ind in model.individuals:
        ind.make_sequence(model.ts, model.reference.haplotype[0])
        if ind.is_reference:
            next
        ind.simulate_reads(pathlib.Path("foo"), cpus=args.cpus)

    model.dump()


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
