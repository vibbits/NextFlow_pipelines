# NextFlow Data Integration Pipelines
Nextflow is a reactive workflow framework and a programming DSL that eases the writing of data-intensive computational pipelines. 

In order to extends this approach, here we provide automation script to allow users running MOFA and WGCNA in tandem, using a glue script to collect user data entry and trigger the execution of the WCGNA and MOFA scripts.

So, to execute your pipeline something like this should work:
First, download the script and the GO folder to your computer.
Second, make sure you create a folder named "data" which is where you will then place your input (omics) data sets
Next, make sure you have NextFlow installed (https://www.nextflow.io/docs/latest/getstarted.html)
Then navigate to the the location on your computer whoich contains the nf script and GO folder and type the command in the Terminal: 
```
nextflow run WGCNA_MOFA_NextFlowPipeline.nf 
```
A typical 2Xdata sets will run for approximitly 10 minutes. 
