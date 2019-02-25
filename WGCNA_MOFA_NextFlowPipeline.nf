#!/usr/bin/env nextflow
myDir = file("./")
myData = file("./data")


process MOFA {

	storeDir myData

	output:
	stdout mofa

	"""
    #!/usr/bin/env Rscript

	library(Hmisc)
	library(WGCNA)
	library(pheatmap)
	enableWGCNAThreads()

	setwd(paste0("$myDir","/data"))
	getwd()

	files <- list.files(pattern = "*.csv")
	print(files)

    i <- 1
    setnames <- c()
    superset <- list()


	for (file in files) {

    	setnames[i] <- gsub(".csv", "", file)
    	data <- read.csv(file,row.names=1, comment.char="", colClasses=c("character",rep("numeric",9)), strip.white=FALSE)
    	dataExpr = as.data.frame(t(data))

    	gsg = goodSamplesGenes(dataExpr, verbose = 3);
    	if (!gsg\$allOK) {
		  # Optionally, print the gene and sample names that were removed:
		  if (sum(!gsg\$goodGenes) > 0)
		    printFlush(paste("Removing genes:", paste(names(dataExpr)[!gsg\$goodGenes], collapse = ", ")));
		  if (sum(!gsg\$goodSamples) > 0)
		    printFlush(paste("Removing samples:", paste(rownames(dataExpr)[!gsg\$goodSamples], collapse = ", ")));
		  # Remove the offending genes and samples from the data:
		  dataExpr = dataExpr[gsg\$goodSamples, gsg\$goodGenes]
		  data = as.data.frame(t(dataExpr))
		}

		superset[[i]] <- data

    	i <- i + 1
    }



    names(superset) = setnames


	library(MOFAtools)
	library(magrittr)
	library(Rmisc)


	setwd("$myDir")	#storeDir or publishDir won't work
	dir.create("MOFA")
	setwd("MOFA")
	MOFAobject <- createMOFAobject(superset)	#the name "superset" is used when loading the data
	MOFAobject

	pdf("TilesData.pdf")
	plotTilesData(MOFAobject)
	dev.off()

	DataOptions <- getDefaultDataOptions()
	
	ModelOptions <- getDefaultModelOptions(MOFAobject)
	ModelOptions\$numFactors <- 25


	TrainOptions <- getDefaultTrainOptions()
	#TrainOptions\$maxiter <- 5		# optional parameter
	TrainOptions\$DropFactorThreshold <- 0.02
	TrainOptions\$seed <- 2017

	MOFAobject <- prepareMOFA(
  		MOFAobject, 
  		DataOptions = DataOptions,
  		ModelOptions = ModelOptions,
  		TrainOptions = TrainOptions
  		)

  	MOFAobject

    MOFAobject <- runMOFA(MOFAobject, outfile = "./MOFAobject")

    MOFAobject
    r2 <- calculateVarianceExplained(MOFAobject)
	r2
	

	for (object in viewNames(MOFAobject)) {

		pdf(object)
		plotWeightsHeatmap(MOFAobject, object, factors="all", show_colnames=F)
		dev.off()
		
		plot2 <- plotWeights(MOFAobject, object, factor = 1)
		plot3 <- plotTopWeights(MOFAobject, object, 1)
		
		png(paste0(object,"2"))
		multiplot(plot2,plot3)
		dev.off()
	}

	if ("mRNA" %in% viewNames(MOFAobject)) {
		data("reactomeGS")


		pdf("enrichment")
		# perfrom enrichment analysis
		##################################### CHECKPOINT!!!
		fsea.out <- runEnrichmentAnalysis(MOFAobject, "mRNA", reactomeGS, alpha = 0.01)

		plotEnrichmentBars(fsea.out, alpha=0.01)
		interestingFactors <- 4:5
		for (factor in interestingFactors) {
		  lineplot <- plotEnrichment(MOFAobject, fsea.out, factor, alpha=0.01, max.pathways=10)
		  print(lineplot)
		}
		####### NEED TO CHANGE HERE FROM SPECIFIC NAMING TO GENERIC ONES:
		#Ordination of samples by factors to reveal clusters and graadients in the sample space
		plotFactorScatter(MOFAobject, factors = 1:2, color_by = "IGHV", shape_by = "trisomy12")
		plotFactorScatters(MOFAobject, factors = 1:4, color_by = "IGHV")
		#A single factor can be visualised using the ‘FactorBeeswarmPlot’ function
		plotFactorBeeswarm(MOFAobject, factors = 1, color_by = "IGHV")
		##Customized analysis
		## For customized exploration of weights and factors, you can directly fetch the variables from the model using ‘get’ functions: ‘getWeights’, ‘getFactors’ and ‘getTrainData’:
		MOFAweights <- getWeights(MOFAobject, views="all", factors="all", as.data.frame = T)
		head(MOFAweights)
		MOFAfactors <- getFactors(MOFAobject, factors=c(1,2), as.data.frame = F)
		head(MOFAfactors)
		dev.off()
		}

	"""
}


process setter {
	output:
	stdout files

	"""
	#!/usr/bin/env python

	import os
	import glob
	os.chdir("{}/data/".format("$myDir"))
	result = [i for i in glob.glob('*.csv')]
	maxit = len(result)
	mySet = []

	for x in range(maxit):
		for y in range(x + 1, maxit):
			print (result[x], result[y])

	"""
}


process wgcna {
	input:
	each x from files.splitText()

	output:
	stdout test

	"""
	#!/usr/bin/env Rscript
	dir.create(paste0("$myDir","/WGCNA"))
	setwd(paste0("$myDir","/WGCNA"))

	#####################################
	#get filenames
	#####################################
	library(tools)
	
	omicsset <- strsplit(trimws("$x"), " ")

	omics1 <- omicsset[[1]][1]
	omics2 <- omicsset[[1]][2]

	print(omics1)
	print(omics2)

	headDir = paste0(file_path_sans_ext(omics1, compression = FALSE),"_-_",file_path_sans_ext(omics2, compression = FALSE))


	############################################################################################################################################################################################
	# First part to read and make the data format in the correct format
	############################################################################################################################################################################################
	# WGCNA tutorial website https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/index.html
	
	#install.packages("Hmisc")
	#install.packages("BiocManager")
	#BiocManager::install("WGCNA")

	library(Hmisc)
	library(WGCNA)
	library(pheatmap)
	#enableWGCNAThreads()
	#allowWGCNAThreads() 
	options(stringsAsFactors = FALSE)
	dir.create(headDir)
	setwd(headDir)
	dir.create("Figures")
	dir.create("Info")
	dir.create("ModuleMembership_vs_GeneSignificance")
	dir.create("GeneInfo")
	data = read.csv((paste0("$myData","/",omics1)), header = T);  ## GE matrix ####
	dim(data);
	names(data);

	# The first row of your transposed matrix was rownames of the original matrix, you should first remove them from data,
	dataEhsan = data[,-1]
	# And then add them as rownames to the new data frame
	rownames(dataEhsan) = as.character(data[,1])
	# And finally transpose it
	datExpr0 = as.data.frame(t(dataEhsan))

	gsg = goodSamplesGenes(datExpr0, verbose = 3);
	gsg\$allOK # returns FALSE so we clean all bad rows

	if (!gsg\$allOK) {
	  # Optionally, print the gene and sample names that were removed:
	  if (sum(!gsg\$goodGenes) > 0)
	    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg\$goodGenes], collapse = ", ")));
	  if (sum(!gsg\$goodSamples) > 0)
	    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg\$goodSamples], collapse = ", ")));
	  # Remove the offending genes and samples from the data:
	  datExpr0 = datExpr0[gsg\$goodSamples, gsg\$goodGenes]
	}
	# Next we cluster the samples (in contrast to clustering genes that will come later) to see for obvious outliers.  

	sampleTree = hclust(dist(datExpr0), method = "average");

	sizeGrWindow(12,9)
	pdf(file = "Figures/sampleClustering.pdf", width = 12, height = 9);
	par(cex = 0.6);
	par(mar = c(0,4,2,0))
	plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
	     cex.axis = 1.5, cex.main = 2)
	# Plot a line to show the cut
	abline(h = 45000, col = "red");
	dev.off()

	###if removing an outlier sample, do it with a height cut-off
	#clust = cutreeStatic(sampleTree, cutHeight = 45000, minSize = 10)
	#table(clust)
	### clust 1 contains the samples we want to keep.
	#keepSamples = (clust==1)
	#datExpr = datExpr0[keepSamples, ]
	# IF NOT REMOVING ANY SAMPLE: 
	datExpr = datExpr0
	nGenes = ncol(datExpr)
	nSamples = nrow(datExpr)
	# datExpr; # The variable datExpr now contains the expression data ready for network analysis.

	# Save and return final data frame ready for analysis
	save(datExpr, file = "saveToFile.RData")

	######################################################
	# Use this to parse phenotype data that you want to correlate with the expression data
	# Can be metabolite data, individual traits, etc.
	traitData = read.csv((paste0("$myData","/",omics2)), header = TRUE);
	dim(traitData)
	names(traitData)
	rownames(traitData) = traitData[,1];
	traitData = traitData [,-1]
	allTraits = t(traitData);

	MLsamples = rownames(datExpr);

	traitRows = match(MLsamples, rownames(allTraits));

	datTraits = allTraits[traitRows,];
	rownames(datTraits) = rownames(allTraits[traitRows,]);
	collectGarbage();
	  
	sampleTree2 = hclust(dist(datExpr), method = "average")
	# Convert traits to a color representation: white means low, red means high, grey means missing entry
	traitColors = numbers2colors(datTraits, signed = FALSE);
	plotDendroAndColors(sampleTree2, traitColors,
	                    groupLabels = names(datTraits),
	                    main = "Sample dendrogram and trait heatmap")

	save(datExpr, datTraits, file = "01-dataInput.RData")




	############################################################################################################################################################################################
	# Automatic, one-step network construction and module detection
	############################################################################################################################################################################################

	#=====================================================================================
	#
	#  Code chunk 1
	#
	#=====================================================================================

	# Display the current working directory
	#getwd();
	# If necessary, change the path below to the directory where the data files are stored. 
	# "." means current directory. On Windows use a forward slash / instead of the usual .
	#workingDir = ".";
	#setwd(headDir); 
	# Load the WGCNA package
	library(WGCNA)
	# The following setting is important, do not omit.
	options(stringsAsFactors = FALSE);
	# Allow multi-threading within WGCNA. 
	# Caution: skip this line if you run RStudio or other third-party R environments.
	# See note above.
	enableWGCNAThreads()
	# Load the data saved in the first part
	lnames = load(file = "01-dataInput.RData");
	# The variable lnames contains the names of loaded variables.
	lnames
	# Get the number of sets in the multiExpr structure.
	# nSets = checkSets(multiExpr)\$nSets
	nSets = 1
	multiExpr = datExpr

	print("1 done")

	#=====================================================================================
	#
	#  Code chunk 2
	#
	#=====================================================================================

	# Choose a set of soft-thresholding powers
	powers = c(seq(4,10,by=1), seq(12,20, by=2));
	# Initialize a list to hold the results of scale-free analysis
	powerTables = vector(mode = "list", length = nSets);
	# Call the network topology analysis function for each set in turn
	for (set in 1:nSets)
	  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr, powerVector=powers,
	                                                     verbose = 2)[[2]]);
	collectGarbage();
	# Plot the results:
	colors = c("black", "red")
	# Will plot these columns of the returned scale free analysis tables
	plotCols = c(2,5,6,7)
	colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
	"Max connectivity");
	# Get the minima and maxima of the plotted points
	ylim = matrix(NA, nrow = 2, ncol = 4);
	for (set in 1:nSets){
	  for (col in 1:length(plotCols))
	  {
	    ylim[1, col] = min(ylim[1, col], powerTables[[set]]\$data[, plotCols[col]], na.rm = TRUE);
	    ylim[2, col] = max(ylim[2, col], powerTables[[set]]\$data[, plotCols[col]], na.rm = TRUE);
	  }
	}
	# Plot the quantities in the chosen columns vs. the soft thresholding power
	sizeGrWindow(8, 6)
	pdf(file = "Figures/scaleFreeAnalysis.pdf", wi = 8, he = 6)
	par(mfcol = c(2,2));
	par(mar = c(4.2, 4.2 , 2.2, 0.5))
	cex1 = 0.7;
	for (col in 1:length(plotCols)) for (set in 1:nSets)
	{
	  if (set==1)
	  {
	    plot(powerTables[[set]]\$data[,1], -sign(powerTables[[set]]\$data[,3])*powerTables[[set]]\$data[,2],
	         xlab = "Soft Threshold (power)", ylab = colNames[col],type = "n", ylim = ylim[, col],
	         main = colNames[col]);
	    addGrid();
	  }
	  if (col==1)
	  {
	    text(powerTables[[set]]\$data[,1], -sign(powerTables[[set]]\$data[,3])*powerTables[[set]]\$data[,2],
	         labels = powers, cex = cex1, col = colors[set]);
	  } else
	    text(powerTables[[set]]\$data[,1], powerTables[[set]]\$data[,plotCols[col]],
	         labels = powers, cex = cex1, col = colors[set]);
	  # if (col==1)
	  # {
	  #   legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
	  # } else
	  #   legend("topright", legend = setLabels, col = colors, pch = 20) ;
	}
	dev.off();

	print("2 done")

	#=====================================================================================
	#
	#  Code chunk 3
	#
	#=====================================================================================

	# Form multi-set expression data: columns starting from 9 contain actual expression data.
	multiExpr = vector(mode = "list", length = nSets)
	multiExpr[[1]] = list(data = datExpr);
	names(multiExpr[[1]]\$data) = colnames(datExpr);
	rownames(multiExpr[[1]]\$data) = rownames(datExpr);


	# Check that the data has the correct format for many functions operating on multiple sets:
	exprSize = checkSets(multiExpr)

	for (i in 1:ncol(multiExpr[[1]]\$data)) {
	  multiExpr[[1]]\$data[,i] <- as.numeric(as.character(multiExpr[[1]]\$data[,i]))
	  }

	softTreshold <- 8


	net <- blockwiseConsensusModules(
	        multiExpr, power = softTreshold , minModuleSize = 20, deepSplit = 2,
	        pamRespectsDendro = FALSE, 
	        mergeCutHeight = 0.25, numericLabels = TRUE,
	        minKMEtoStay = 0,
	        saveTOMs = TRUE, verbose = 5)

	print("3 done")

	#=====================================================================================
	#
	#  Code chunk 4
	#
	#=====================================================================================

	consMEs = net\$multiMEs;
	moduleLabels = net\$colors;
	# Convert the numeric labels to color labels
	moduleColors = labels2colors(moduleLabels)
	consTree = net\$dendrograms[[1]]; 
	# saving modules and genes in a table
	anno = read.csv(file = "../../GO/Mouse2GO.csv", sep = ",");
	anno = anno[,c(1:2)]
	# Save annotation and data expression files to .csv file
	write.csv(anno, file = "Info/anno.csv", row.names = FALSE)
	write.csv(data\$X, file = "Info/datExpr.csv", row.names = FALSE)

	t <- colnames(multiExpr[[1]]\$data)
	t <- as.matrix(t)
	head(t)
	colnames(t) <- "GENE"

	# Right merge of annotation and data datasets
	#genesAnno <- merge(x = anno, y = t, by= "GENE" ,all.y = TRUE). # THIS ONE TRIGGERS AN ERROR (OT)
	genesAnno <- merge(x = anno, y = t, by.x= 'gene_id', by.y = 'GENE'  ,all.y = TRUE)

	genesAnno <- genesAnno[,c(1:2)]
	colnames(genesAnno) <- c("GENES", "Annotation")

	library(dplyr)

	test <- genesAnno %>% group_by(GENES) %>% summarise(Annotations = paste(Annotation, collapse = ", "))
	colnames(test)

	# Write text file with genes and corresponding modules and annotations
	mod <- data.frame("Gene" = t, "Modules" = labels2colors(moduleLabels), "Annotation" = test\$Annotations)
	mod <- mod[order(mod\$Modules),]
	write.table(mod, file = "Info/ModuleMembers.txt", row.names = F, sep = "\t", quote = F)

	save(consMEs, moduleLabels, moduleColors, consTree, file = "02-networkConstruction-auto.RData")

	print("4 done")

	############################################################################################################################################################################################
	# 3 - Relating modules to external clinical traits and identifying important genes
	############################################################################################################################################################################################

	print("start part 3")

	#=====================================================================================
	#
	#  Code chunk 1
	#
	#=====================================================================================

	# Display the current working directory
	#getwd();
	# If necessary, change the path below to the directory where the data files are stored. 
	# "." means current directory. On Windows use a forward slash / instead of the usual .
	workingDir = ".";
	setwd(workingDir); 
	# Load the WGCNA package
	library(WGCNA)
	# The following setting is important, do not omit.
	options(stringsAsFactors = FALSE);
	# Load the expression and trait data saved in the first part
	lnames = load(file = "01-dataInput.RData");
	#The variable lnames contains the names of loaded variables.
	lnames
	# Load network data saved in the second part.
	lnames = load(file = "02-networkConstruction-auto.RData");
	lnames

	print("3.1 done")

	#=====================================================================================
	#
	#  Code chunk 2
	#
	#=====================================================================================

	# Define numbers of genes and samples
	nGenes = ncol(datExpr);
	nSamples = nrow(datExpr);
	# Recalculate MEs with color labels
	MEs0 = moduleEigengenes(datExpr, moduleColors)\$eigengenes
	MEs = orderMEs(MEs0)
	moduleTraitCor = cor(MEs, datTraits, use = "p");
	moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

	write.table(moduleTraitCor, file = "Info/metabolites_transcriptomics_corr_matrix.txt", row.names = T, sep = "\t", quote = F)

	print("3.2 done")

	#=====================================================================================
	#
	#  Code chunk 3
	#
	#=====================================================================================

	sizeGrWindow(10,6)
	# Will display correlations and their p-values
	textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
	                           signif(moduleTraitPvalue, 1), ")", sep = "");
	dim(textMatrix) = dim(moduleTraitCor)
	# Display the correlation values within a heatmap plot
	pdf(file = "Figures/Heatmap_metabolites_transcriptomics_corr_matrix.pdf", wi = 10, he = 6)
	par(mar = c(3, 6, 7, 3));
	labeledHeatmap(Matrix = moduleTraitCor,
	               xLabels = colnames(moduleTraitCor),
	               xLabelsPosition = "top",
	               xLabelsAngle = -90,
	               yLabels = names(MEs),
	               ySymbols = names(MEs),
	               colorLabels = FALSE,
	               colors = greenWhiteRed(50),
	               setStdMargins = FALSE,
	               cex.text = 0.5,
	               cex.lab.x = 0.45,
	               cex.lab.y = 0.3,
	               zlim = c(-1,1))
	title("Module-trait relationships", line = 5)
	dev.off()

	print("3.3 done")

	#=====================================================================================
	#
	#  Code chunk 4
	#
	#=====================================================================================

	# 
	metabo = list()
	for (i in colnames(moduleTraitPvalue)) {
	  # Define variable meta containing the metabolite columns of datTrait
	  meta <- data.frame(datTraits[,i]);
	  meta\$i <- i
	  metabo[[i]] <- meta
	}

	metabol = do.call(cbind, metabo)
	metabol <- metabol[c(T,F)]
	colnames(metabol) = gsub(".datTraits...i.", "", colnames(metabol))
	  
	# Names (colors) of the modules
	modNames = substring(names(MEs), 3)
	geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
	MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples), append = T);

	names(geneModuleMembership) = paste("MM", modNames, sep="");
	names(MMPvalue) = paste("p.MM", modNames, sep="");
	  
	geneTraitSignificance = as.data.frame(cor(datExpr, metabol, use = "p"));
	GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
	  
	names(geneTraitSignificance) = paste("GS.", names(metabol), sep="");
	names(GSPvalue) = paste("p.GS.", names(metabol), sep="");

	print("3.4 done")


	#=====================================================================================
	#
	#  Code chunk 5 
	#
	#=====================================================================================

	# Calculate the module membership values (aka. module eigengene based connectivity kME):
	datKME <- signedKME(datExpr, MEs)  # equals geneModuleMembership
	colnames(datKME) <- sub("kME", "MM.", colnames(datKME))

	colorOfColumn <- substring(names(datKME),4)

	metabolites <- colnames(moduleTraitPvalue)
	for (metabolite in metabolites) {
	  cat(paste("\nmetabolite:", metabolite, "\n"))
	  pdf(paste0("ModuleMembership_vs_GeneSignificance/", metabolite, "_ModuleMembershipVsgeneSig.pdf"))
	  par(mar=c(6, 8, 4, 4) + 0.1)
	  modNames = substring(names(MEs), 3)
	  for (module in modNames) {
	    column <- match(module,colorOfColumn)
	    moduleGenes <- moduleColors == module;
	    GS = paste("GS.", metabolite, sep = "")
	    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
	                       abs(geneTraitSignificance[moduleGenes, GS]),
	                       xlab = paste("Module membership in", module, "module"),
	                       ylab = paste("Gene significance for", metabolite),
	                       xlim = c(0,1),
	                       ylim = c(0,1),
	                       main = paste("Module membership vs. gene significance\n"),
	                       cex.main = 1.3,
	                       cex.lab = 1.1,
	                       cex.axis = 1,
	                       pch = 21, 
	                       col = "dark gray", 
	                       bg = module)
	  }
	  
	  dev.off()
	} # end of loop

	print("3.5 done")


	#=====================================================================================
	#
	#  Code chunk 6
	#
	#=====================================================================================


	# Define top genes for every metabolite in every module 
	topGene <- list()
	topGene_ <- list()
	for(metabolite in metabolites) {
	  for(module in modNames) {
	    GS = paste("GS.", metabolite, sep = "")
	    MM = paste("MM", module, sep = "")
	    topGenes = abs(geneModuleMembership[,MM]) > 0.8 & abs(geneTraitSignificance[,GS]) > 0.8
	    topGene[[module]] <- topGenes
	  }
	  topGene_[[metabolite]] <- topGene
	} # end of loop

	# Link with gene names
	topGene2name <- list()
	topGene2name_ <- list()
	for(metabolite in metabolites) {
	  for(module in modNames) {
	    topGenes2names <- dimnames(datExpr)[[2]][topGene_[[metabolite]][[module]]]
	    topGene2name[[module]] <- topGenes2names
	  }
	  topGene2name_[[metabolite]] <- topGene2name
	} # end of loop

	# Filter empty list 
	filterGene1_ <- list()
	for(metabolite in names(topGene2name_)) {
	  filter = topGene2name_[[metabolite]][lapply(topGene2name_[[metabolite]],length)>0]
	  filterGene1_[[metabolite]] <- filter
	} # end of loop


	# Filter empty list 
	filterGene_ <- list()
	for (items in filterGene1_){
	  filter2 = filterGene1_[lapply(filterGene1_,length)>0]
	  filterGene_ <- filter2
	}

	print("3.6 done")



	#=====================================================================================
	#
	#  Code chunk 7
	#
	#=====================================================================================

	# Write output to files
	metabolites2 <-names(filterGene_)
	for(metabolite in metabolites2) {
	  modules <-names(filterGene_[[metabolite]])
	  GTS <- paste("GS.", metabolite, sep = "")
	  
	  for(module in modules) {
	    cat(paste("\nCreating output for top genes in", module, "for", metabolite, "\n"))
	    genes <- filterGene_[[metabolite]][[module]]
	    GMM = paste("MM", module, sep = "")
	    # Create the starting data frame
	    infoGenes = list(genes = genes,
	                           moduleColor = moduleColors[topGene_[[metabolite]][[module]]],
	                           annotation = test\$Annotations[topGene_[[metabolite]][[module]]],
	                           geneTraitSignificance = geneTraitSignificance[filterGene_[[metabolite]][[module]], GTS],
	                           geneModuleMembership = geneModuleMembership[filterGene_ [[metabolite]][[module]], GMM])
	    

	    # Order the genes in the infoTopGenes variable first by module color, then by geneTraitSignificance
	    #geneOrder <- order(infoGenes\$moduleColor, -abs(infoGenes[,4]))
	    #infoGenes = infoGenes[geneOrder, ]
	    
	    write.csv(infoGenes, file = paste0("GeneInfo/", module, "_", metabolite, "_GeneInfo.csv"), row.names = FALSE)
	  }
	} # end of loop

	print("3.7 done")



	"""
}
