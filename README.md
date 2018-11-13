# NextFlow_pipelines
Tricks and tips for designing a NextFlow pipeline

So I was wondering if we can allow stops in a NextFlow pipeline at some key steps where the user need to examine the output before continue the pipeline. 

the use case where you ask the user from within the pipeline for input is useful in case in the second part of the pipeline, you still want to be able to access intermediate results from the last part.
That wouldn't be possible with separate sub-pipeline as they don't share a cache.

So we have two questions:
- ask for user input during pipeline
- is cache invalidated when changing the config

OPTION 1

Seems essential that there's a single point in the pipeline where you want to make a break.
So, I think we can solve it, by collecting all output items at some point. That way the `leftOrRight` choice process only triggers when all input values are processed.

(I guess `collect()` could also be replaced by other operators that trigger only when all items are processed. Like `last()`
but then the best I can come up with now, is give the choice to continue (maybe choice of different paths) or stop (and then manually resume). Which... boils down to splitting up your pipeline in multiple pipelines.

OPTION 2

Group your parameters in the config file per process (see ```cachedParams.nf``` and ```nextflow.config```) and call them like in the test above, when changing one of the values, only the scripts involved will be re-run on a `-resume`

OPTION 3
you could also add a command line parameter to specifiy whether to execute the second part of the pipeline or not. You can then use either an if statement or the when directive to specify which processes belong to the second part.
So, something like this should work:
```
params.runPart2 = false

process p1 {
    ...
}

process p2 {
    ...
}

process p3 {
    ...
    when:
    runPart2 == true
    ...
}
```
And then execute your pipeline like this: 
```
nextflow run ...
# Check the intermediate output of process p1 and p2
nextflow run ... --runPart2 -resume
```
OPTION 4 (suggested by nf developer (https://github.com/nextflow-io/nextflow/issues/923#issuecomment-437887855)
Nextflow is designed for automation therefore interactive data entry is not a supported scenario. I would try to split the pipeline in two, using a glue Bash script to collect user entry and trigger the execution of the second script. Cache is unaffected by the actual name of the script executed. As long as you specify the `-resume` option it should work as expected."
