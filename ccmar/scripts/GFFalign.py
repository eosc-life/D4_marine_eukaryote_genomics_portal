#!/usr/bin/python

import argparse
import os

import gffutils
from Bio import SeqIO


class GeneComp:
    """ This class aims to discover the gene position in the alignment and extract information on his 
        characteristics """
    def __init__(self, qg_gene, qg_aln_start, tg_gene, tg_aln_start, extract, tol, gattr, output=None):
        # gene are the gene from maf-tab file
        # start, end are the coordinates
        self.tg_db_result_start = int(tg_gene.start) - tg_aln_start
        self.tg_db_result_end = int(tg_gene.end) - tg_aln_start
        self.qg_db_result_start = int(qg_gene.start) - qg_aln_start
        self.qg_db_result_end = int(qg_gene.end) - qg_aln_start
        self.qg_gene_length = self.qg_db_result_end - self.qg_db_result_start  # query gene length
        self.rg_gene_length = self.tg_db_result_end - self.tg_db_result_start  # target gene length
        self.percentage_diff_rg_tg_gene = (self.qg_gene_length / self.rg_gene_length) * 100  # how much longer or 
        # shorter the gene are, in percentage 
        
        self.tol = tol  # how much tolerance in nucleotides
        self.qg_region_name = qg_gene.chrom
        self.rg_region_name = tg_gene.chrom
        self.q_outstart = qg_aln_start + self.tg_db_result_start
        self.q_outend = qg_aln_start + self.tg_db_result_end
        self.genelist = []
        self.qg_aln_start = qg_aln_start
        if extract:
            self.fastafile = extract[0]
            self.fastaout = extract[1]
        self.output = output
        self.out = ""
        self.gattr = f"New_annotation='{gattr}'"

    def is_equal(self):
        if (self.tg_db_result_start - self.tol) <= self.qg_db_result_start <= (self.tg_db_result_start + self.tol) and (
                self.tg_db_result_end - self.tol) <= self.qg_db_result_end <= (self.tg_db_result_end + self.tol):
            self.out = f"{self.qg_region_name}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t\
            {self.qg_aln_start}\t.\tID={self.qg_region_name};{self.gattr};Note=confirmed"
            return self.out

    def is_shorter(self):
        under_length = (self.rg_gene_length - self.qg_gene_length) / self.rg_gene_length * 100

        condition1 = (self.tg_db_result_start + self.tol) < self.qg_db_result_start and \
                     (self.tg_db_result_end - self.tol) > self.qg_db_result_end
        condition2 = (self.tg_db_result_start - self.tol) <= self.qg_db_result_start <= \
                     (self.tg_db_result_start + self.tol) and (
                             self.tg_db_result_end - self.tol) > self.qg_db_result_end
        condition3 = (self.tg_db_result_start + self.tol) < self.qg_db_result_start and \
                     (self.tg_db_result_end - self.tol) <= self.qg_db_result_end <= (
                             self.tg_db_result_end + self.tol)

        if condition1 or condition2 or condition3:
            self.out = f"{self.qg_region_name}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qg_aln_start}\t.\tID={self.qg_region_name};{self.gattr};\
                    Note=shorter;Size={under_length:,.2f}% shorter"
            return self.out

    def is_longer(self):
        over_length = (self.qg_gene_length - self.rg_gene_length) / self.qg_gene_length * 100

        condition1 = (self.tg_db_result_start - self.tol) > self.qg_db_result_start and (self.tg_db_result_end + self.tol) < self.qg_db_result_end
        condition2 = (self.tg_db_result_start - self.tol) <= self.qg_db_result_start <= (self.tg_db_result_start + self.tol) and (
                self.tg_db_result_end + self.tol) < self.qg_db_result_end
        condition3 = (self.tg_db_result_start - self.tol) > self.qg_db_result_start and (self.tg_db_result_end - self.tol) <= self.qg_db_result_end <= (
                self.tg_db_result_end + self.tol)

        if condition1 or condition2 or condition3:
            self.out = f"{self.qg_region_name}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qg_aln_start}\t.\tID={self.qg_region_name};{self.gattr};\
                    Note=longer;Size={over_length:,.2f}% longer"
            return self.out

    def is_offset(self):
        note = None

        if (self.tg_db_result_start + self.tol) < self.qg_db_result_start < (self.tg_db_result_end - self.tol) and (self.tg_db_result_end + self.tol) < self.qg_db_result_end:
            note = ';Note=offset_right'

        if (self.tg_db_result_start - self.tol) > self.qg_db_result_start and (self.tg_db_result_start - self.tol) < self.qg_db_result_end < self.tg_db_result_end + self.tol:
            self.out = f"{self.qg_region_name}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qg_aln_start}\t.\tID={self.qg_region_name};{self.gattr};Note=offset_left"
            note = ';Note=offset_left'

        self.out = f"{self.qg_region_name}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qg_aln_start}\t.\tID={self.qg_region_name};{self.gattr} {note}"

        return note

    def is_different(self):
        if self.tg_db_result_start - self.tol > self.qg_db_result_end or self.tg_db_result_end + self.tol < self.tg_db_result_start:
            self.out = f"{self.qg_region_name}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qg_aln_start}\t.\tID={self.qg_region_name};{self.gattr};Note=new"
            return self.out

    def extract_fasta(self):
        with open(self.fastaout, "a") as fileout:
            with open(self.fastafile) as filefasta:
                for record in SeqIO.parse(filefasta, "fasta"):
                    if record.id == self.qg_region_name:
                        fileout.write(
                            f">{self.qg_region_name}_{self.q_outstart}:{self.q_outend}\n{record.seq[self.q_outstart:self.q_outend]}\n")

    def fout(self):
        try:
            if self.__class__.fout.called:
                with open(self.output, "a") as fileout:
                    fileout.write(self.out + "\n")
        except AttributeError:
            try:
                if os.path.isfile(self.output):
                    os.remove(self.output)
                with open(self.output, "w") as fileout:
                    # fileout.write("##gff-version 3\n")
                    fileout.write(f"##gff-version 3\n")
                self.__class__.fout.called = True
                self.__class__.fout(self)
            except TypeError:
                print("Problems reading the file...")

def diff_gene(qg_genes, tg_genes, tg_aln_start, qg_aln_start, query_db):
    # are the  two genes the same?
    extract = args.extract or None
    tol = args.tolerance or 30

    for target_gene in tg_genes:
        # gattr is used to suggest a functional annotation based on target gene
        gattr = str(target_gene).split("\t")[-1]
        for query_gene in qg_genes:
            algene = GeneComp(qg_gene=query_gene, qg_aln_start=qg_aln_start, tg_gene=target_gene,
                              tg_aln_start=tg_aln_start,
                              extract=extract, tol=tol, gattr=gattr)
            if "new" in args.verbosity or "all" in args.verbosity and algene.is_different():
                algene_out = algene.out
                algene_name = algene.qg_region_name
                algene_start = int(algene_out.split("\t")[3])
                algene_end = int(algene_out.split("\t")[4])

                if not list(
                        query_db.region(region=(algene_name, algene_start, algene_end), completely_within=False)):
                    if args.output:
                        algene.output = args.output
                        algene.fout()
                    else:
                        # this print allow to show on the console with the user does not want to create a file
                        print(algene.is_different())
                        if args.extract:
                            algene.extract_fasta()

            elif "shorter" in args.verbosity or "all" in args.verbosity and algene.is_shorter():
                if args.output:
                    algene.output = args.output
                    algene.fout()
                else:
                    print(algene.is_shorter())
                    if args.extract:
                        algene.extract_fasta()

            elif "longer" in args.verbosity or "all" in args.verbosity and algene.is_longer():
                if args.output:
                    algene.output = args.output
                    algene.fout()
                else:
                    print(algene.is_longer())
                    if args.extract:
                        algene.extract_fasta()

            elif "offset" in args.verbosity or "all" in args.verbosity and algene.is_offset():
                if args.output:
                    algene.output = args.output
                    algene.fout()
                else:
                    print(algene.is_offset())
                    if args.extract:
                        algene.extract_fasta()

            elif "confirmed" in args.verbosity and algene.is_equal():
                if args.output:
                    algene.output = args.output
                    algene.fout()
                else:
                    print(algene.is_equal())
                    if args.extract:
                        algene.extract_fasta()


def gff_gene_check(GffName, db_name, memory=0):
    """ The IDs on the GFFs must be unique. Only the gene information are needed, all the other information must be removed
    from the GFFs. """
    tempgff = ""
    for line in open(GffName):
        if not line.startswith("#") and line.split()[2] == 'gene':
            tempgff += line

    if memory:
        # Write the db in memory and return it as variable, so it can be used as subclass of _DBCreator
        dbout = gffutils.create_db(tempgff, ":memory:", from_string=True)
        return dbout
    else:
        gffutils.create_db(tempgff, db_name, from_string=True)


def get_coordinates_by_strand_type(aln_start, aln_length, strand):
    if strand == "+":
        end = aln_start + aln_length
    else:
        end = aln_start - aln_length
        aln_start, end = end, aln_start
    return aln_start, end


def check_position(line, query_db, target_db):
    # check if there is a gene in the aligned area lets' start with the coordinates characteristics
    line_content_list = line.split('\t')

    # tg - target genome
    # qg - query genome
    # aln - alignment

    rg_aln_region_name = line_content_list[1]
    rg_aln_start = int(line_content_list[2])
    rg_aln_length = int(line_content_list[3])
    rg_aln_strand = line_content_list[4]
    qg_aln_region_name = line_content_list[6]
    qg_aln_start = int(line_content_list[7])
    qg_aln_length = int(line_content_list[8])
    qg_aln_strand = line_content_list[9]

    # check the strand, if - reverse the start and the end
    qg_aln_start, qg_aln_end = get_coordinates_by_strand_type(qg_aln_start, qg_aln_length, qg_aln_strand)
    rg_aln_start, rg_aln_end = get_coordinates_by_strand_type(rg_aln_start, rg_aln_length, rg_aln_strand)

    # lists of gene within the coordinates
    rg_genes = [gene for gene in list(target_db.region(region=(rg_aln_region_name, rg_aln_start, rg_aln_end),
                                                       completely_within=True))]
    qg_genes = [gene for gene in list(query_db.region(region=(qg_aln_region_name, qg_aln_start, qg_aln_end),
                                                      completely_within=True))]
    if len(rg_genes):
        diff_gene(qg_genes=qg_genes, tg_genes=rg_genes, tg_aln_start=rg_aln_start, qg_aln_start=qg_aln_start,
                  query_db=query_db)

    return 1


def main(args):
    # import the GFF library and create (has to be done only once) and the GFF databases

    target_db = query_db = None
    db_query_name = args.queryGff + "_db"
    db_target_name = args.targetGff + "_db"

    if args.force_database:
        os.remove(db_query_name)
        os.remove(db_target_name)

    if args.use_query_database:
        # use the db passed by the user
        query_db = gffutils.FeatureDB(args.use_query_database)
    else:
        # create the db
        query_genome_DB = gff_gene_check(args.queryGff, db_query_name, args.memory)

        if args.memory:
            query_db = gffutils.FeatureDB(query_genome_DB.dbfn)
        else:
            query_db = gffutils.FeatureDB(db_query_name)

    if args.use_target_database:
        target_db = gffutils.FeatureDB(args.use_target_database)
    else:
        target_genome_DB = gff_gene_check(args.targetGff, db_target_name, args.memory)

        if args.memory:
            target_db = gffutils.FeatureDB(target_genome_DB.dbfn)
        else:
            target_db = gffutils.FeatureDB(db_target_name)

    fcounter = 0

    # put the content of the DB in objects
    with open(args.aln) as maftab:
        for line in maftab:
            if not line[0].startswith('#'):
                fcounter += check_position(line, query_db, target_db)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Tool to extract genes coordinates from a whole genome alignent.
                        This script needs an alignment in TAB format and two gff files''')
    parser.add_argument('aln', help='''Alignment file in TAB format. The suggested way to obtain it is to run Last and 
                        than convert the file from MAF to TAB with maf-convert''')
    parser.add_argument('queryGff',
                        help='''Gff file of the query organism. The gene IDs in the GFF must be unique.''')
    parser.add_argument('targetGff',
                        help='''Gff file of the "target" organism. The gene IDs in the GFF must be unique.''')
    parser.add_argument("-uq", "--use-query-database",
                        help='''Use this parament if you already have a query gffutils formatted database or
                                if it\'s not the first time you run this script''', type=str)
    parser.add_argument("-ut", "--use-target-database",
                        help='''Use this parament if you already have a target gffutils formatted database or
                                if it\'s not the first time you run this script''', type=str)
    parser.add_argument("-fd", "--force-database", help='''delete old gffutils databases and create new ones''',
                        action='store_true')
    parser.add_argument("-m", "--memory", help='''create an in-memory database. This option can't be used with the 
                        other DB options. Probably usefully in Galaxy integration''', action='store_true')
    parser.add_argument("-e", "--extract",
                        help='''Extract the fasta sequence of the new suggested gene. It takes two argument: the fasta 
                        file of the genome and the name of the output file. This will slow down the process A LOT.''',
                        nargs=2, type=str)
    parser.add_argument("-o", "--output", help='''Name of the output file. Default output is stout''', type=str)
    parser.add_argument("-t", "--tolerance",
                        help='''Interval, in nucleotide, within a gene is considered in the same position. 
                        Default 30''', default=30, type=int)
    parser.add_argument("-v", "--verbosity",
                        help='''Output options. If not specify the software shows only the genes that are in the exact 
                        position of the genes in the target. It\'s possible to show annotated genes that are in aligned 
                        regions but that have different lengths or in slightly different positions. It's possible to 
                        select multiple, space separated, values.''',
                        choices=["all", "shorter", "longer", "offset", "new", "confirmed"], nargs='*', default='new',
                        type=str)
    parser.add_argument('--version', action='version', version='0.1.0')
    args = parser.parse_args()
    main(args)
