# NextFlow_pipelines
Tricks and tips for designing a NextFlow pipeline


So I was wondering if we can allow stops in a NextFlow pipeline at some key steps where the user need to examine the output before continue the pipeline. 

So we have two questions:
- ask for user input during pipeline
- is cache invalidated when changing the config

OPTION 1

1) Seems essential that there's a single point in the pipeline where you want to make a break.
So, I think we can solve it, by collecting all output items at some point. That way the `leftOrRight` choice process only triggers when all input values are processed.

(I guess `collect()` could also be replaced by other operators that trigger only when all items are processed. Like `last()`
but then the best I can come up with now, is give the choice to continue (maybe choice of different paths) or stop (and then manually resume). Which... boils down to splitting up your pipeline in multiple pipelines.

OPTION 2

2) Group your parameters in the config file per process (see ```cachedParams.nf``` and ```nextflow.config```) and call them like in the test above, when changing one of the values, only the scripts involved will be re-run on a `-resume`
