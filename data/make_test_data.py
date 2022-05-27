#!/usr/bin/env python3
import argparse
import itertools
import multiprocessing
import os
import pathlib
import re
import subprocess
import tempfile
from dataclasses import dataclass
from random import choices

import cyvcf2
import tskit


def update_individual_metadata(ts):
    """Update individual metadata in ts object.

    Returns new list of individuals with updated metadata
    """
    individuals = []
    oldname = None
    i = 0
    for ind in ts.individuals():
        popindex = ts.nodes()[ind.nodes[0]].population
        popname = ts.population(popindex).metadata["name"]
        if oldname is not None and oldname != popname:
            i = 0
        oldname = popname
        name = f"tsk_{ind.id}_{popname}_{i}"
        md = {
            "id": ind.id,
            "tskit_id": f"tsk_{ind.id}",
            "name": name,
            "popname": popname,
            "ind_index": i,
            "population": ts.population(popindex),
            "description": (
                f"tskit_individual:tsk_{ind.id}, name:{name}, "
                f"population:{ts.population(popindex)}"
            ),
            "is_reference": False,
        }
        md["vcfheader"] = (
            f"##<SAMPLE=<ID={md['tskit_id']},Name={md['name']},"
            f"Index={md['id']},Population={md['popname']},"
            f"Description=\"{md['population'].metadata['description']}\">"
        )
        newind = ind.replace(metadata=md)
        individuals.append(newind)
        i = i + 1
    return individuals


def set_reference_individual(ref, individuals):
    """Add metadata tag to indicate reference individual. States will
    be flipped wrt to the reference indivual's states."""
    if ref is None:
        return individuals
    values = ref.split(":")
    refpop = int(values[0]) if re.match(r"\d+", values[0]) else values[0]
    refind = int(values[1])
    for i, ind in enumerate(individuals):
        if ind.metadata["ind_index"] == refind:
            if (
                ind.metadata["population"].id == refpop
                or ind.metadata["popname"] == refpop
            ):
                ind.metadata["is_reference"] = True
                individuals[i] = ind
                break
    return individuals


def partition_individuals(individuals):
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


class SeqRecord:
    def __init__(self, *, seq, name, id, description):
        self.seq = seq
        self.name = name
        self.id = id
        self.desc = description

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        s = f">{self.id} {self.desc}\n"
        n = 72
        chunks = [self.seq[i : i + n] for i in range(0, len(self), n)]
        s += "\n".join(chunks)
        s += "\n"
        return s

    def __repr__(self):
        return str(self)


@dataclass
class Population:
    id: str
    metadata: dict


@dataclass
class Sample:
    id: str
    population: Population = None
    metadata: dict = None
    paternal_chromosome: str = None
    maternal_chromosome: str = None

    @property
    def sampleid(self):
        return f"{self.population}-{self.id}"

    def haplotype(self, which=0):
        desc = (
            f"id:{individual.id}, node:{node}, name:{individual.metadata['name']}, "
            f"haplotype:{haplotype}, "
            f"population:{individual.metadata['population']}"
        )
        if which == 0:
            return self.maternal_chromosome
        return self._haplotype[which]

    @hapotype.setter
    def haplotype(self, sec):
        self._haplotype

    def simulate_reads(self, coverage):
        pass


def make_sequence(ts, individual, node, haplotype, reference):
    """Create DNA sequence with variants at sites"""
    seqid = f"tsk_{individual.metadata['population'].metadata['name']}-{individual.id}-{haplotype}"
    for site, variant in zip(ts.sites(), list(ts.haplotypes())[node]):
        i = int(site.position) - 1
        reference[i] = variant
    desc = (
        f"id:{individual.id}, node:{node}, name:{individual.metadata['name']}, "
        f"haplotype:{haplotype}, "
        f"population:{individual.metadata['population']}"
    )
    record = SeqRecord(
        seq="".join(reference),
        name=individual.metadata["name"],
        id=seqid,
        description=desc,
    )
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


def make_ooa(args):
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

    # Run msp ancestry
    subprocess.run(
        [
            "msp",
            "ancestry",
            "--demography",
            outdir / "ooa_with_outgroup.demes.yaml",
            "CHB:7",
            "CEU:6",
            "YRI:6",
            "chimpanzee:1",
            "gorilla:1",
            "orangutan:1",
            "--length",
            str(SEQLENGTH),
            "--recombination-rate",
            str(RECOMBINATION_RATE),
            "--random-seed",
            str(RANDOM_SEED),
            "-o",
            outdir / "ooa.ts",
        ],
        check=True,
    )

    RANDOM_SEED = 42
    # Run msp mutations
    subprocess.run(
        [
            "msp",
            "mutations",
            "-o",
            outdir / "ooa.mut.ts",
            "--random-seed",
            str(RANDOM_SEED),
            str(MUTATION_RATE),
            outdir / "ooa.ts",
        ],
        check=True,
    )
    os.unlink(outdir / "ooa.ts")

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
        # s = Sample(ind.metadata["name"])
        seq_out = outdir / f"{ind.id}.fasta"
        nodes = ind.nodes
        for node in nodes:
            haplotype = node % len(nodes)
            rec = make_sequence(ts, ind, node, haplotype, reference)
            seq_out = outdir / f"{rec.id}.fasta"
            indseq.append(seq_out)
            with open(seq_out, "w") as fh:
                fh.write(repr(rec))
        make_reads(outdir, genomes=indseq, cpus=args.cpus)
        for seq in indseq:
            os.unlink(seq)


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
