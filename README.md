# NextFlow Data Integration Pipelines
Nextflow is a reactive workflow framework and a programming DSL that eases the writing of data-intensive computational pipelines. 

Here we provide a NextFlow 'automation glue script' which allow users running MOFA and WGCNA in tandem. The glue script first collect user data entry and then triggers the execution of the WCGNA (https://github.com/cran/WGCNA; https://www.ncbi.nlm.nih.gov/pubmed/19114008) and MOFA (https://github.com/bioFAM/MOFA) pipelines.

To execute the NextFlow glue script: 
1) first, make sure you have NextFlow installed (https://www.nextflow.io/docs/latest/getstarted.html).

2) download the ```WGCNA_MOFA_NextFlowPipeline.nf``` script and the 'GO' folder to your computer. The template GO folder provided here contains a file mouse gene IDs which are mapped to GO terms (to be used in the functional annotation steps in WGCNA and MOFA). Then navigate to the the location on your computer which contains the nf script and GO folder. 

3) Create a folder named 'data' which is where you will then place your input (omics) data sets. 

4) Type the command in the Terminal: 
```
nextflow run WGCNA_MOFA_NextFlowPipeline.nf 
```
A typical 2 data sets run will take approximitly 10 minutes to complete. 
