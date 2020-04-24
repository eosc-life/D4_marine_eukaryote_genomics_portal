#!/usr/bin/python3 

import gffutils
from Bio import SeqIO
from collections import defaultdict
import argparse
import os
#import sys, getopt

### TODO: add strand control


class GeneComp:
    # This cass aims to discover the position of the genes in the alignment
    # and to extract informations on the genes characteristics

    def __init__(self, ccgene, qstart, qend, otgene, dstart, dend, extract, tol, output=None):
        # d* / ot* are the db result, q*  / cc the query
        # gene are the gene from maf-tab file 
        # start, end are the coordinates  '''
        self.otbeg = int(otgene.start) - dstart
        self.otend = int(otgene.end) - dstart
        self.ccbeg = int(ccgene.start) - qstart
        self.ccend = int(ccgene.end) - qstart
        self.tol = tol # how much tolerance in nucleotides
        self.ccname = ccgene.chrom
        self.otname = otgene.chrom
        self.q_outstart = qstart+self.otbeg
        self.q_outend = qstart+self.otend
        self.genelist = []
        self.qstart = qstart
        if extract:
            self.fastafile = extract[0]
            self.fastaout = extract[1]
        self.output = output
        self.out = ""

    def is_equal(self):
        if (self.otbeg - self.tol) <= self.ccbeg <= (self.otbeg + self.tol) and (self.otend - self.tol) <= self.ccend <= (self.otend + self.tol):
            self.out = f"{self.ccname}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qstart}\t.\tID={self.ccname};Note=confirmed"
            return self.out

    def is_shorter(self):
        if (self.otbeg + self.tol) < self.ccbeg and (self.otend - self.tol) > self.ccend:
            self.out =  f"{self.ccname}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qstart}\t.\tID={self.ccname};Note=shorter"
            return self.out
        elif (self.otbeg - self.tol) <= self.ccbeg <= (self.otbeg + self.tol) and (self.otend - self.tol) > self.ccend:
            self.out =  f"{self.ccname}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qstart}\t.\tID={self.ccname};Note=shorter_right"
            return self.out
        elif (self.otbeg + self.tol) < self.ccbeg and (self.otend - self.tol) <= self.ccend <= (self.otend + self.tol):
            self.out =  f"{self.ccname}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qstart}\t.\tID={self.ccname};Note=shorter_left"
            return self.out

    def is_longer(self):
        if (self.otbeg - self.tol) > self.ccbeg and (self.otend + self.tol) < self.ccend:
            self.out =  f"{self.ccname}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qstart}\t.\tID={self.ccname};Note=longer"
            return self.out
        elif (self.otbeg - self.tol) <= self.ccbeg <= (self.otbeg + self.tol) and (self.otend + self.tol) < self.ccend:
            self.out =  f"{self.ccname}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qstart}\t.\tID={self.ccname};Note=longer_right"
            return self.out
        elif (self.otbeg - self.tol) > self.ccbeg and (self.otend - self.tol) <= self.ccend <= (self.otend + self.tol):
            self.out =  f"{self.ccname}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qstart}\t.\tID={self.ccname};Note=longer_left"
            return self.out

    def is_offset(self):
        if (self.otbeg + self.tol) < self.ccbeg < (self.otend - self.tol) and (self.otend + self.tol) < self.ccend:
            self.out =  f"{self.ccname}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qstart}\t.\tID={self.ccname};Note=offset_right"
            return self.out
        if (self.otbeg - self.tol) > self.ccbeg and (self.otbeg - self.tol) < self.ccend < self.otend + self.tol:
            self.out =  f"{self.ccname}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qstart}\t.\tID={self.ccname};Note=offset_left"
            return self.out

    def is_different(self):
        if self.otbeg - self.tol > self.ccend or self.otend + self.tol < self.otbeg:
            self.out =  f"{self.ccname}\tprediction\tgene\t{self.q_outstart}\t{self.q_outend}\t.\t{self.qstart}\t.\tID={self.ccname};Note=different"
            return self.out

    def extract_fasta(self):
        with open(self.fastaout,"a") as fileout:
            with open(self.fastafile) as filefasta:
                for record in SeqIO.parse(filefasta,"fasta"):
                    if record.id == self.ccname:
                        fileout.write(f">{self.ccname}_{self.q_outstart}:{self.q_outend}\n{record.seq[self.q_outstart:self.q_outend]}\n")

    # def extract_stout(self):
    #     with open(self.fastafile) as filefasta:
    #             for record in SeqIO.parse(filefasta,"fasta"):
    #                 if record.id == self.ccname:
    #                     print(f">{self.ccname}_{self.q_outstart}:{self.q_outend}\n{record.seq[self.q_outstart:self.q_outend]}")


    #def __str__(self):
    #    return f"{self.ccname}\t{self.q_outstart}\t{self.q_outend}"


    def fout(self):
        try:
            if self.__class__.fout.called:
                with open(self.output,"a") as fileout:
                    fileout.write(self.out+"\n")
        except AttributeError:
            try:
                if os.path.isfile(self.output):
                    os.remove(self.output)
                with open(self.output,"a") as fileout:
                        #fileout.write("##gff-version 3\n")
                    fileout.write(f"##gff-version 3\n{self.out}\n")
                self.__class__.fout.called = True
                self.__class__.fout(self)
            except TypeError:
                pass


def diff_gene(query_genes, target_genes, dstart, dend, qstart, qend, query_db):
    # are the  two genes the same?
    if args.extract:
        extract = args.extract
    else:
        extract = None

    if args.tollerance:
        tol = args.tollerance
    else:
        tol = 30

    for otgene in target_genes:
    #print(args.use_query_database)
       # otbeg = int(otgene.start) - dstart
       # otend = int(otgene.end) - dstart
        for ccgene in query_genes:
            algene = GeneComp(ccgene, qstart, qend, otgene, dstart, dend, extract, tol)
            if "new" in args.verbosity or "all" in args.verbosity:
                if algene.is_different():
                    algene_out = algene.is_different()
                    algene_name = algene_out.split("\t")[0]
                    algene_start = int(algene_out.split("\t")[3])
                    algene_end = int(algene_out.split("\t")[4])
                    #print(al_name, al_start, al_end)
                    #print(list(target_db.region(region=(al_name, al_start, al_end), completely_within=True)))
                    if not list(query_db.region(region=(algene_name, algene_start, algene_end), completely_within=False)):
                        if args.output:
                            algene.output=args.output
                            algene.fout()
                        else:
                            print(algene.is_different())
                            if args.extract:
                                algene.extract_fasta()

            if "shorter" in args.verbosity or "all" in args.verbosity:
                if algene.is_shorter():
                    if args.output:
                        algene.output=args.output
                        algene.fout()
                    else:
                        print(algene.is_shorter())
                        if args.extract:
                            algene.extract_fasta()
            if "longer" in args.verbosity or "all" in args.verbosity:
                if algene.is_longer():
                    if args.output:
                        algene.output=args.output
                        algene.fout()
                    else:
                        print(algene.is_longer())
                        if args.extract:
                                algene.extract_fasta()
            if "offset" in args.verbosity or "all" in args.verbosity:
                if algene.is_offset():
                    if args.output:
                        algene.output=args.output
                        algene.fout()
                    else:
                        print(algene.is_offset())
                        if args.extract:
                                algene.extract_fasta()
            if "confirmed" in args.verbosity or "all" in args.verbosity:
                if algene.is_equal():
                    algene.output=args.output
                    algene.fout()
                else:
                    print(algene.is_equal())
                    if args.extract:
                        algene.extract_fasta()




def check_strand(start,leng,strand):
    if strand == "+":
        end = start+leng
    else:
        end = start-leng
        start,end = end,start
    return(start,end)

def check_position(line, query_db,target_db):
    # check if there is a gene in the aligned area
    # lets' start with the coordinatescharacteristics
    elsp = line.split('\t')
    dname = elsp[1]
    dstart = int(elsp[2])
    dlen = int(elsp[3])
    dstrand = elsp[4]
    qname = elsp[6]
    qstart = int(elsp[7])
    qlen = int(elsp[8])
    qstrand = elsp[9]

    # check the strand, if - reverse the start and the end
    qstart,qend = check_strand(qstart,qlen,qstrand)
    dstart,dend = check_strand(dstart,dlen,dstrand)

    #counter
    d_counter = 0

    # lists of gene within the coordinates
    target_genes = [ gene for gene in list(target_db.region(region=(dname, dstart, dend), completely_within=True))]
    query_genes = [ gene for gene in list(query_db.region(region=(qname, qstart, qend), completely_within=True))]



    if len(target_genes):
    # and len(query_genes):
        if len(target_genes) >= len(query_genes): # the number of genes in the aligned area must be bigger in the target (?)
            diff_gene(query_genes, target_genes, dstart, dend, qstart, qend, query_db)


    if d_counter>=1:
        return(d_counter)
    else:
        return(0)


def main(args):
    fcounter = 0
    # import the GFF library and create (has to be done only once)
    # the GFF databases
    target_db = query_db = None
    db_query_name=args.queryGff + "_db"
    db_target_name=args.targetGff + "_db"

    if args.force_database:
        os.remove(db_query_name)
        os.remove(db_target_name)

    if args.use_query_database:
        query_db = gffutils.FeatureDB(args.use_query_database)
    else:
        gffutils.create_db(args.queryGff, db_query_name)
        query_db = gffutils.FeatureDB(db_query_name)
        # except ValueError:
        #     print("Please check your GFF file")
        # except:
        #     raise Exception(f"It seems you already have a db called {db_query_name}. Use the -qu if you want to use it or delete the")


    if args.use_target_database:
        target_db = gffutils.FeatureDB(args.use_target_database)
    else:
        gffutils.create_db(args.targetGff, db_target_name)
        target_db = gffutils.FeatureDB(db_target_name)
        # except ValueError:
        #     print("Please check your GFF file")if "all" in args.verbosity:
        # except:
        #     raise Exception(f"It seems you already have a db called {db_target_name}. Use the -du if you want to use it or delete the")

    # put the content of the DB in objects
    with open(args.aln) as maftab:
        for line in maftab:
            fcounter += check_position(line, query_db, target_db)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = '''Tool to extract genes coordinates from a whole genome alignent.
                        This script needs an alignement in TAB format and two gff files''')
    parser.add_argument('aln', help='alignment file in TAB format. The suggested way to obtain it is to run Last and\
                                     than convert the file from MAF to TAB with maf-convert')
    parser.add_argument('queryGff',
                        help='''Gff file of the query organism. The gene IDs in the GFF must be unique. 
                        To solve the problem please extract only the "gene" lines. Try to format
                        the file with AWK: awk '{if ($3==\"gene\") print $0}' GFFfile''')
    parser.add_argument('targetGff',
                        help='''Gff file of the "target" organism. The gene IDs in the GFF must be unique.
                        To solve the problem please extract only the "gene" lines as explained in queryGff''')
    parser.add_argument("-uq","--use-query-database",
                        help='''Use this parament if you already have a query gffutils formatted database or
                                if it\'s not the first time you run this script''', type=str)
    parser.add_argument("-ut","--use-target-database",
                        help='''Use this parament if you already have a target gffutils formatted database or
                                if it\'s not the first time you run this script''', type=str)
    parser.add_argument("-fd", "--force-database", help="delete old gffutils databases and create new ones", action='store_true')
    parser.add_argument("-e","--extract",
                        help='''Extract the fasta sequence of the new suggested gene. It takes two argument: the fasta file
                        of the genome and the name of the output file. This will slow down the process A LOT.''', nargs=2, type=str)
    #parser.add_argument("-es", "--extract_stout", help='Like -e but it will print the result in the standard output. FASTER than -e. It need the fasta file', type=str)
    parser.add_argument("-o", "--output", help='Name of the output file. Default output is stout', type=str)
    parser.add_argument("-t", "--tollerance", help='Interval, in nucleotide, within a gene is considered in the same position. Default 30', default=30, type=int)
    parser.add_argument("-v", "--verbosity",
                        help='''Output options. If not specify the software shows only the genes that are in the  exact position of the
                        genes in the target. It\'s possible to show annotated genes that are in aligned regions but that have different lengths or in slightly 
                        different positions. It's possible to select multiple, space separated, values.''',
                        choices=["all","shorter", "longer", "offset", "new", "confirmed"], nargs='*', default='new', type=str)
                        #metavar=('all','shorter','longer','shifted','new'), 
                        #action='append', 
                        #choices=["all","shorter", "longer", "shifted", "new"], default="new")
    parser.add_argument('--version', action='version', version='0.1.0')
    args = parser.parse_args()
    main(args)
