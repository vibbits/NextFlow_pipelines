# NextFlow Data Integration Pipelines
Nextflow is a reactive workflow framework and a programming DSL that eases the writing of data-intensive computational pipelines. NextFlow provides many simple but powerful command-line and scripting tools that, when chained together, facilitate complex data manipulations.

In order to extends this approach, adding the ability to define complex program interactions and a high-level parallel computational environment based on the dataflow programming model.

Nextflow is designed for automation of MOFA and WGCNA pipelines, using a glue script to collect user data entry and trigger the execution of the WCGNA and MOFA scripts.

So, to execute your pipeline something like this should work:
```
nextflow run WGCNA_MOFA_NextFlowPipeline.nf 
```

