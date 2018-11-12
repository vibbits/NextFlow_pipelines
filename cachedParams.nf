//Create a simple channel
inputChannel = Channel.from([1,2,3,4,5])

process init{
    tag "$x"

    input:
	val x from inputChannel

    output:
	val y into nextChannel

    exec:
	y = 2*x
}


//Do something for every value
//Process parameters are set in nextflow.config
//Changing these will re-run this process.
process doSomething{
    tag "$x"
    
    input:
	val x from nextChannel

    output:
	set val(x), file(outputFile) into outputChannel
    
    script:
	"""
        echo "$params.tweak.a $params.tweak.b $params.tweak.c   : $x"> outputFile
        """
}

process nextStep{
    tag "$x"
    
    input:
	set val(x), file(myFile) from outputChannel

    script:
	"""
        """
}
