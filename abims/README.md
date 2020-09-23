# Workflow

In accordance with the workflow set up by CCMAR, ABiMS has put online the various tools necessary for its execution on the [galaxy.sb-roscoff.fr](https://galaxy.sb-roscoff.fr/) website.

## Workflow steps

1. Repeatmasker

[Repeatmasker](http://www.repeatmasker.org/) is accessible to mask repetitive regions of a genome.


2. Last

5 commands of [Last](http://last.cbrc.jp/) package have been wrapped as Galaxy tools :

- [lastdb](http://last.cbrc.jp/doc/lastdb.html) : for genome indexing.

- [last-train](http://last.cbrc.jp/doc/last-train.html) : to find the rates of insertion, deletion, and substitutions between the genomes.

- [lastal](http://last.cbrc.jp/doc/lastal.html) : for genome alignments.

- [last-split](http://last.cbrc.jp/doc/last-split.html) : to find split or spliced alignments.

- [maf-convert](http://last.cbrc.jp/doc/maf-convert.html) : to convert lastal maf output file in a tabular file.


3. Synteny GFF align

To complete the workflow, we also wrapped the CCMAR tool that transfer annotations between related species using genome synteny relationship.


## Docker 

A docker Galaxy stable containing the workflow tools is also available. You can find it on https://github.com/abims-sbr/D4_marine_eukaryote_genomics_portal_docker.


# Expectations : GGA in the cloud

We plan to provide the annotation pipeline and the GGA environmentavailable through the European Open Science Cloud. We are currently developing a script to automatically deploy the GGA containers and load data into them.

Our approach : https://docs.google.com/presentation/d/1xY6DrAX-3eAjke7emOS8BGr_HQF9CsoBoSa5_O2d3Mc/edit?usp=sharing