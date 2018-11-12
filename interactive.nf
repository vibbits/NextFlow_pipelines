
//Create a simple channel
inputChannel = Channel.from([1,2,3,4,5])

//Do something for every value
process doSomething{
    tag "$x"
    
    input:
	val x from inputChannel

    output:
	val x into passChannel
    
    script:
	'''
        ##Do nothing special
        '''
}

//Collect all output values
passChannel
    .collect()
    .set{ allDone }


//When all values are collected: ask for user input
process leftOrRight{
    input:
	val resultList from allDone

    output:
	val answer into goOn

    exec:
	print resultList
        println ""
        print "What now ([L]eft/[R]ight/[Q]uit)?   "
        answer = new InputStreamReader(System.in).readLine()
}

goOn.into{ goOnLeft ; goOnRight ; stopPipeline }

//Take action based on user input
process left{
    input:
	val answer from goOnLeft

    when:
	answer=='L'

    exec:
	print "Run process LEFT"
}

process right{
    input:
	val answer from goOnRight

    when:
	answer=='R'

    exec:
	print "Run process RIGHT"
}

process stop{
    input:
	val answer from stopPipeline

    when:
	answer=='Q'

    exec:
	print "Bye!"
}
