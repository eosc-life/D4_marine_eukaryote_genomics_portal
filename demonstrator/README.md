# Demonstrator 4 - Marine Eukaryote Genomics Portal

The demonstrator D6 is addressing the transfer of annotations among marine biology species using [ORCAE](../orcae/) technology. They will develop a workflow for these annotation transfer and the annotations themselves will be accessible though a portal.

## Questions for demonstrator D4
### What is the scientific problem you want to solve?
Community-driven genome annotation portals (e.g. Orcae) do not currently have a means of transfering annotation between closely related species. As genomes are annotated differently and at different rates, it would be convenient to have a means of transfering annotations between closely related species. We intend to provide a tool that calculates synteny blocks between genomes and suggests reciprocal annotations which would require the intervention of the researcher to then update a genome. We intend to provide this tool as a Galaxy workflow (maybe also documented in CWL) and standalone tool (containerized) that is able to inject annotations into a Orcae annotation platform.

### What is the workflow (in terms of steps/tasks to be executed)?
The project aims to build a computational workflow to transfer genome annotations between closely related species – as a test case, pelagic fishes - using genome synteny relationships
The workflow will consist of 3 parts:
Alignment of the genomes and extraction of synteny relationships
Visualization of the synteny blocks and selection of possible missing annotation
Injection of the new potential annotation on the ORCAE portal


### How do you currently run your workflow? e.g.
It’s not currently a documented workflow, but we run the pipeline on a computational cluster

### What are the known bottlenecks/limitations of your workflows? Any information that you can provide based on past experience
Synteny computations can be quite demanding if run fully, but we are only looking for the blocks of highest similarity so it is hoped we can reduce the computational burden

### How do you document your workflows?
It is not currently documented. We intend to provide a Galaxy workflow, and/or CWL workflow.

### How do you develop and manage your workflows?
Currently, we don’t.

### How do you disseminate your workflows?
We don’t have any workflows to disseminate.

### Please specify if there is already a working prototype of the workflow

No there is not.

**If yes at the previous question: Please specify if the entire workflow or some tools have been used previously in a cloud computing environment**
N/A

### How does the workflow access data?

Through Galaxy.

### Does the workflow require human interaction?

Yes.

### Does the workflow require use of public/private/access to controlled access datasets?

No.

