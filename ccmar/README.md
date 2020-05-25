# Demonstrator Workflow

This is an example of the steps that can be part of the Galaxy Workflow.

To run this example you need two genomes in fasta format and their corresponding GFF files.
Note:
  - The query genome is the genome you want to "annotate"
  - The target genome is the genome you want to use to annotate your genome
  - This is just a small and fast example. Last has A LOT of different that can be setted to improve the results
  - The computational time is just to have an idea on how long it should take. It can vary a lot.


## Workflow Steps
 
### 1. Genome Masking

The target genome must be masked. If the target genome is not pre-masked it's possible to use different software. 
I choose WindowMasker because is fast and gave a decent results. Repeatmasker it's another possible choice.
Run WindowMasker on the Genome you want to user as target (2 steps ~15 minutes (5+10)):
 ```bash
    windowmasker -mk_counts -in target_genome.fasta > target_genome.wm
    windowmasker -ustat target_genome.wm -outfmt fasta -in target_genome.fasta > target_genome-wm.fasta
```
### 2. Target genome indexing

Target genome indexing using lastdb, using the masked genome (~ minutes 3 with 56 threads):
 ```bash
	lastdb -P {threads} -uNEAR -R11 -c target_db target_genome-wm.fasta
```

### 3. Last-train

Finds the rates of insertion, deletion, and substitutions between the genomes (~1  minute with 56 threads):
 ```bash
	last-train -P {threads} target_db query_genome.fasta > train.out
```
### 4. Alignment (~12 minutes with 56 threads):

The genomes alignment. In Galaxy, probably, this command has to be splitted.
 ```bash
	lastal -P {threads} -p train.out target_cgenome query_genome.fasta | last-split -m1 > genome_aligned.maf
```
- to increase the sensitivity (but the computational time will be about 1 hour) use: `-m100 -l10`
				( -l is the minimal initial match, reduce the sensitivity and increase the speed)

### 5. From MAF to TXT (~1 minute):

Covert the maf output in txt.
 ```bash
	maf-convert -j1e5 tab genome_aligned.maf > genome_aligned.tab 
```
Again, 1e5 it's only an example

### 6. GFF cleaning

The IDs on the GFFs must be unique. Moreover, because only the gene information are needed, all the other information must
be removed from the GFFs.
 ```bash
	awk -i inplace'{if ($3=="gene") print $0}' *.GFF
```
In some cases you want to run awk without the __-i inplace__ parameter and redirect the output in a new file.

### 7. Writing the GFF

Run the script called __Gff_sel_out.py__ that is in the scripts folder:

```bash
    Gff_sel_out.py -o suggested_genes.gff genome_aligned.tab query_genome.gff target_genome.gff 
```

Also in this case there are different parameters that can be selected.
