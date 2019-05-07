
cat("\n Started! \n")
#sessionInfo()


# Send R errors to stderr
options(show.error.messages = F, error = function(){cat(geterrmessage(), file = stderr()); q("no", 1, F)})

# Avoid crashing Galaxy with an UTF8 error on German LC settings
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")
suppressPackageStartupMessages({
library(HTqPCR)
library(base)
library(Biobase)
library(utils)
library(stats)
library(graphics)
library(grDevices)
library(RColorBrewer)
library(limma)
library(RankProd)
library(methods)
library(impute)
library(BBmisc)
library(affy)
library(psych)
#library(gmp)
library(zoo)
library(nondetects)
library(Hmisc)
#library(missForest)
#library(mice)
})

cat("\n R libraries...loaded!\n")

args = commandArgs(trailingOnly=TRUE)
dpfiles<-basename(args[1])
path000<-dirname(args[1])
format<-args[2]
nfeatures<-args[3]
rawout<-args[4]
path <- args[5]
dcCtmin<-args[6]
dcCtmax<-args[7]
dcflag<-args[8]
x<-args[9]
normalizationMethod<-args[10]
if (normalizationMethod=="deltaCt") {
    normalizers<-args[11]
    outputNorm<-args[12]
    outputECDF<-args[13]
    percentofnastoremove<-args[14]
    outputRemaining<-args[15]
    imputeMethod<-args[16]
    if (imputeMethod=="knn") {
      kappa<- args[17]
      macsp<-args[18]
      outputIMP<-args[19]

      DEAMethod<-args[20]
      if (DEAMethod=="ttest") {
          alternative<- args[21]
          paired<-args[22]
          replicates<- args[23]
          sort<-args[24]
          stringent<- args[25]
          padjust<-args[26]
          outputDEA<-args[27]
          filtnames<-args[28]
        } else {
          outputDEA<-args[21]
          filtnames<-args[22]
        }
    } else {
      #mean, median, nondetects, cubic 
        outputIMP<-args[17]
        DEAMethod<-args[18]
        if (DEAMethod=="ttest") {
            alternative<- args[19]
            paired<-args[20]
            replicates<- args[21]
            sort<-args[22]
            stringent<- args[23]
            padjust<-args[24]
            outputDEA<-args[25]
            filtnames<-args[26]
          } else {
            outputDEA<-args[19]
            filtnames<-args[20]
          }
    }
  }else {
    outputNorm<-args[11]
    outputECDF<-args[12]
    percentofnastoremove<-args[13]
    outputRemaining<-args[14]
    imputeMethod<-args[15]
   
    if (imputeMethod=="knn") {
      kappa<- args[16]
      macsp<-args[17]
      outputIMP<-args[18]

      DEAMethod<-args[19]
     
      if (DEAMethod=="ttest") {
          alternative<- args[20]
          paired<-args[21]
          replicates<- args[22]
          sort<-args[23]
          stringent<- args[24]
          padjust<-args[25]
          outputDEA<-args[26]
          filtnames<-args[27]
        } else {
          outputDEA<-args[20]
          filtnames<-args[21]
        }
     } else {
       #mean, median, nondetects, cubic 
        outputIMP<-args[16]
        DEAMethod<-args[17]
          if (DEAMethod=="ttest") {
            alternative<- args[18]
            paired<-args[19]
            replicates<- args[20]
            sort<-args[21]
            stringent<- args[22]
            padjust<-args[23]
            outputDEA<-args[24]
            filtnames<-args[25]
           
          } else {
            outputDEA<-args[18]
            filtnames<-args[19]
            
          }
      }
  }
cat("\n Initialization completed! \n")

.readCtEDS	<-
function(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...)
{
	# Scan through beginning of file, max 100 lines
	file.header <- readLines(con=readfile, n=100)
	n.header	<- grep("^Well", file.header)
	if (length(n.header)==0)
		n.header	<- 0
	# Read data, skip the required lines
	out	<- read.delim(file=readfile, header=TRUE, colClasses="character", nrows=nspots*n.data[i], skip=n.header-1, strip.white=TRUE, ...)
	out
} # .readCtEDS


.readCtPlain	<- 
function(readfile=readfile, header=header, n.features=n.features, n.data=n.data, i=i, ...)
{
	# A check for file dimensions. Read a single file.
	sample	<- read.delim(file=readfile, header=header, ...)
	nspots	<- nrow(sample)
	if (nspots != n.features*n.data[1])
		warning(paste(n.features, "gene names (rows) expected, got", nspots))
	# Read in the required file
	out	<- read.delim(file=readfile, header=header, colClasses="character", nrows=nspots*n.data[i], ...)
	# Return
	out
} # .readCtPlain

.readCtSDS	<- 
function(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...)
{
	# Scan through beginning of file, max 100 lines
	file.header <- readLines(con=readfile, n=100)
	n.header	<- grep("^#", file.header)
	if (length(n.header)==0) 
		n.header	<- 0
	# Read data, skip the required lines
	out	<- read.delim(file=readfile, header=FALSE, colClasses="character", nrows=nspots*n.data[i], skip=n.header, strip.white=TRUE, ...)
	# Return
	out
} # .readCtSDS

.readCtLightCycler	<- 
function(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...)
{
	# Read data, skip the required lines
	out	<- read.delim(file=readfile, header=TRUE, as.is=TRUE, nrows=nspots*n.data[i], skip=1, strip.white=TRUE, ...)
	# Return
	out
} # .readCtLightCycler

.readCtCFX	<- function(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...)
{
	# Read data, skip the required lines
	out	<- read.csv(file=readfile, header=TRUE, as.is=TRUE, nrows=nspots*n.data[i], strip.white=TRUE, ...)
	# Return
	out
} # .readCtCFX

.readCtOpenArray	<- 
function(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...)
{
	# Read data
	out	<- read.csv(file=readfile, header=TRUE, as.is=TRUE, nrows=nspots*n.data[i], strip.white=TRUE, ...)
	# Regard those marked as outliers as "Unreliable"
	out$ThroughHole.Outlier[out$ThroughHole.Outlier=="False"]	<- "OK"
	out$ThroughHole.Outlier[out$ThroughHole.Outlier=="True"]	<- "Unreliable"	
	# Return
	out
} # .readCtOpenArray

.readCtBioMark	<- 
function(readfile=readfile, n.data=n.data, i=i, nspots=nspots, ...)
{
	# Scan through beginning of file, max 100 lines
	file.header <- readLines(con=readfile, n=100)
	n.header	<- grep("^ID", file.header)-1
	if (length(n.header)==0) 
		n.header	<- 0
	# Read data, skip the required lines
	out	<- read.csv(file=readfile, header=TRUE, as.is=TRUE, nrows=nspots*n.data[i], skip=n.header, strip.white=TRUE, ...)
	# Convert the calls into flags
	out$Call[out$Call=="Pass"] <- "OK"
	out$Call[out$Call=="Fail"] <- "Undetermined"	
	# Return
	out
} # .readCtBioMark



readCtDataDav<-
function (files, path = NULL, n.features = 384, format = "plain",
    column.info, flag, feature, type, position, Ct, header = FALSE,
    SDS = FALSE, n.data = 1, samples, na.value = 40, sample.info,
    ...)
{
    if (missing(files))
        stop("No input files specified")
    if (length(files) != length(n.data) & length(n.data) != 1)
        stop("n.data must either be a single integer, or same length as number of files")
    if (length(n.data) == 1)
        n.data <- rep(n.data, length(files))
    nsamples <- sum(n.data)
    ncum <- cumsum(n.data)
    s.names <- NULL
    nspots <- n.features
    if (SDS) {
        warning("Please use format='SDS'. The SDS' parameter is retained for backward compatibility only.")
        format <- "SDS"
    }
    if (!missing(flag) | !missing(feature) | !missing(type) |
        !missing(position) | !missing(Ct)) {
        warning("Please use 'column.info' for providing a list of column numbers containing particular information. The use of 'flag', 'feature', 'type', 'position' and 'Ct' is deprecated and will be removed in future versions.")
    }
    if (missing(column.info)) {
        column.info <- switch(format, EDS = list(flag="EXPFAIL", feature="Target.Name",  position="Well.Position", Ct="CT"),
          plain = list(flag = 4, feature = 6, type = 7, position = 3, Ct = 8), 
            SDS = list(flag = 4,feature = 6, type = 7, position = 3, Ct = 8), 
            #SDS = list(flag = "Omit",feature = "Detector", type = "Task", position = "Pos", Ct = "Avg.Ct"), 
            LightCycler = list(feature = "Name",
            position = "Pos", Ct = "Cp"), CFX = list(feature = "Content",
            position = "Well", Ct = "Cq.Mean"), OpenArray = list(flag = "ThroughHole.Outlier",
            feature = "Assay.Assay.ID", type = "Assay.Assay.Type",
            position = "ThroughHole.Address", Ct = "ThroughHole.Ct"),
            BioMark = list(flag = "Call", feature = "Name.1",
                position = "ID", Ct = "Value"))
    }
    X <- matrix(0, nspots, nsamples)
    X.flags <- as.data.frame(X)
    X.cat <- data.frame(matrix("OK", ncol = nsamples, nrow = nspots),
        stringsAsFactors = FALSE)
    for (i in seq_along(files)) {
        if (i == 1) {
            cols <- 1:ncum[i]
        }
        else {
            cols <- (ncum[i - 1] + 1):ncum[i]
        }
        readfile <- ifelse(is.null(path), files[i], file.path(path,
            files[i]))
        sample <- switch(format, EDS =.readCtEDS(readfile = readfile,
        n.data = n.data, i = i, nspots = nspots, ...), plain = .readCtPlain(readfile = readfile,
            header = header, n.features = n.features, n.data = n.data,
            i = i, ...), SDS = .readCtSDS(readfile = readfile,
            n.data = n.data, i = i, nspots = nspots, ...), LightCycler = .readCtLightCycler(readfile = readfile,
            n.data = n.data, i = i, nspots = nspots, ...), CFX = .readCtCFX(readfile = readfile,
            n.data = n.data, i = i, nspots = nspots, ...), OpenArray = .readCtOpenArray(readfile = readfile,
            n.data = n.data, i = i, nspots = nspots, ...), BioMark = .readCtBioMark(readfile = readfile,
            n.data = n.data, i = i, nspots = nspots, ...))
        #if (format == "SDS") {
        #    if("Avg Ct" %in% colnames(n.data)){
        #        data <- matrix(sample[, column.info[["Avg.Ct"]]], ncol = n.data[i])
        #    } elseif {
        #      cat("\n Unsupported SDS format! ")
        #      }
        #}else{
      data <- matrix(sample[, column.info[["Ct"]]], ncol = n.data[i])
       # }
        undeter <- apply(data, 2, function(x) x %in% c("Undetermined",
            "No Ct"))
        X.cat[, cols][undeter] <- "Undetermined"
        nas <- c("Undetermined", "No Ct", "999", "N/A")
        if (is.null(na.value)) {
            data[data %in% nas | data == 0] <- NA
        }
        else {
            data[data %in% nas | is.na(data) | data == 0] <- na.value
        }
        X[, cols] <- apply(data, 2, function(x) as.numeric(as.character(x)))
        if ("flag" %in% names(column.info)) {
            flags <- matrix(sample[, column.info[["flag"]]],
                ncol = n.data[i])
            flags[flags == "-"] <- "Failed"
            flags[flags == "+"] <- "Passed"
            X.flags[, cols] <- flags
        }
        else {
            X.flags[, cols] <- "Passed"
        }
        if (format == "OpenArray") {
            s.names <- c(s.names, unique(sample$SampleInfo.SampleID))
        }
        else if (format %in% c("EDS","plain", "SDS")) {
            s.names <- c(s.names, unique(sample[, 2]))
        }
        else {
            s.names <- s.names
        }
        if (i == 1) {
            featPos <- paste("feature", 1:nspots, sep = "")
            if ("position" %in% names(column.info))
                featPos <- as.character(sample[1:nspots, column.info[["position"]]])
            featType <- factor(rep("Target", nspots))
            if ("type" %in% names(column.info))
                featType <- sample[1:nspots, column.info[["type"]]]
            featName <- paste("feature", 1:nspots, sep = "")
            if ("feature" %in% names(column.info))
                featName <- as.character(sample[1:nspots, column.info[["feature"]]])
            df <- data.frame(featureNames = featName, featureType = as.factor(featType),
                featurePos = featPos)
            metaData <- data.frame(labelDescription = c("Name of the qPCR feature (gene)",
                "Type pf feature", "Position on assay"))
            featData <- AnnotatedDataFrame(data = df, varMetadata = metaData)
        }
    }
    if (!missing(samples)) {
        if (length(samples) < nsamples) {
            warning("Not enough sample names provided; using Sample1, Sample2, ... instead\n")
            samples <- paste("Sample", 1:nsamples, sep = "")
        }
        else if (length(samples) == nsamples) {
            samples <- samples
        }
    }
    else if (missing(samples)) {
        if (length(files) == nsamples) {
            samples <- gsub("(.+)\\..+", "\\1", files)
        }
        else if (length(s.names) == nsamples) {
            samples <- s.names
        }
        else {
            samples <- paste("Sample", 1:nsamples, sep = "")
        }
    }
    samples <- make.unique(samples)
    if (any(is.na(X)))
        warning("One or more samples contain NAs. Consider replacing these with e.g. Ct=40 now.")
    if (missing(sample.info)) {
        pdata <- data.frame(sample = 1:length(samples), row.names = samples)
        sample.info <- new("AnnotatedDataFrame", data = pdata,
            varMetadata = data.frame(labelDescription = "Sample numbering",
                row.names = "Sample names"))
    }
    X.hist <- data.frame(history = capture.output(match.call(readCtData)),
        stringsAsFactors = FALSE)
    out <- new("qPCRset", exprs = X, phenoData = sample.info,
        featureData = featData, featureCategory = X.cat, flag = X.flags,
        CtHistory = X.hist)
    out
}


head(read.delim(file.path(path000, dpfiles), sep="\t"))
files <- read.delim(file.path(path000, dpfiles), sep="\t")
switch(format,
  "EDS"={
      columns<- list(flag="EXPFAIL", feature="Target.Name",  position="Well.Position", Ct="CT")
      metadata <- data.frame(labelDescription = c("sampleName", "Treatment"),  row.names = c("sampleName", "Treatment"))
      phenoData <- new("AnnotatedDataFrame", data = files, varMetadata = metadata)
      rownames(phenoData)=as.vector(files$sampleName)
      raw<- readCtDataDav(files = as.vector(files$sampleName), header=TRUE,  format="EDS", column.info=columns, path = path,sample.info=phenoData,n.features = as.numeric(nfeatures))
    },
    "plain"={
      metadata <- data.frame(labelDescription = c("sampleName", "Treatment"),  row.names = c("sampleName", "Treatment"))
      phenoData <- new("AnnotatedDataFrame", data = files, varMetadata = metadata)
      rownames(phenoData)=as.vector(files$sampleName)
      raw<- readCtDataDav(files = as.vector(files$sampleName), header=FALSE,  format="plain", path = path, sample.info=phenoData,n.features = as.numeric(nfeatures))
    },
    "SDS"={
      #columns<- list(feature=3, Ct=6, flag=11)
      columns <-list(flag = "Omit",feature = "Detector", type = "Task", position = "Wells", Ct = "Avg.Ct")
      metadata <- data.frame(labelDescription = c("sampleName", "Treatment"),  row.names = c("sampleName", "Treatment"))
      phenoData <- new("AnnotatedDataFrame", data = files, varMetadata = metadata)
      rownames(phenoData)=as.vector(files$sampleName)
      raw<- readCtDataDav(files = files$sampleName, format="SDS",column.info=columns, path = path, sample.info=phenoData, n.features=as.numeric(nfeatures))
    },
    "LightCycler"={
      metadata <- data.frame(labelDescription = c("sampleName", "Treatment"),  row.names = c("sampleName", "Treatment"))
      phenoData <- new("AnnotatedDataFrame", data = files, varMetadata = metadata)
      rownames(phenoData)=as.vector(files$sampleName)
      raw <- readCtDataDav(files = files$sampleName, path = path, format = "LightCycler", sample.info=phenoData,n.features = as.numeric(nfeatures))
    },
    "CFX"={
      metadata <- data.frame(labelDescription = c("sampleName", "Treatment"),  row.names = c("sampleName", "Treatment"))
      phenoData <- new("AnnotatedDataFrame", data = files, varMetadata = metadata)
      rownames(phenoData)=as.vector(files$sampleName)
      raw <- readCtDataDav(files = files$sampleName, path = path, format = "CFX", sample.info=phenoData,n.features = as.numeric(nfeatures))
    },
    "OpenArray"={
      metadata <- data.frame(labelDescription = c("sampleName", "Treatment"),  row.names = c("sampleName", "Treatment"))
      phenoData <- new("AnnotatedDataFrame", data = files, varMetadata = metadata)
      rownames(phenoData)=as.vector(files$sampleName)
      raw <- readCtDataDav(files = files$sampleName, path = path, format = "OpenArray", sample.info=phenoData,n.features = as.numeric(nfeatures))
    },
    "BioMark"={
      metadata <- data.frame(labelDescription = c("sampleName", "Treatment"),  row.names = c("sampleName", "Treatment"))
      phenoData <- new("AnnotatedDataFrame", data = files, varMetadata = metadata)
      rownames(phenoData)=as.vector(files$sampleName)
      raw <- readCtDataDav(files = files$sampleName, path = path, format = "BioMark", sample.info=phenoData, n.features = as.numeric(nfeatures))
    },
    stop("Enter something that switches me!")
)
cat("\n Files read correctly! ")

write.table(exprs(raw), file=rawout, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")
####################################################################################################################
#Set a new categories for the values meeting two criterions
switch(format,
    "EDS"={
      unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="Y", quantile=NULL)
    },
    "plain"={
      unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="Flagged", quantile=NULL)
    },
    "SDS"={
      unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="TRUE", quantile=NULL)
    },
    "LightCycler"={
      unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="Y", quantile=NULL)
    },
    "CFX"={
      unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="Y", quantile=NULL)
    },
    "OpenArray"={
      unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="TRUE ", quantile=NULL)
    },
    "BioMark"={
      unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="Fail", quantile=NULL)
    },
    stop("Enter something that switches me!")
)
#unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="Y", quantile=NULL)

####################################################################################################################
#Filter out the values of the new category
xFilter <- filterCategory(unreliable)

cat("\n Categorization completed! ")

png(x,    # create PNG for the heat map
  width = 10*300,        # 5 x 300 pixels
  height = 10*300,
  res = 300,            # 300 pixels per inch
  pointsize = 8)
  plotCtBoxes(xFilter, stratify=NULL, xlab = "Samples", ylab="Ct", names=as.character(seq(1, ncol(xFilter), 1)))       # smaller font size
dev.off()

#write.table(exprs(xFilter), file=x, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")

####################################################################################################################
#NORMALIZATION
#method version 3.5.1 + Global mean
normalizeCtDataDav <-
function(q,
	norm	= "deltaCt",
	deltaCt.genes	= NULL,
	scale.rank.samples,
	rank.type	= "pseudo.median",
	Ct.max	= dcCtmax,
	geo.mean.ref,
	verbose	= TRUE)
{
	# Extract the data
	data	<- exprs(q)
	data.norm	<- data
	# Get the normalisation method
	method	<- match.arg(norm, c("quantile", "scale.rankinvariant", "norm.rankinvariant", "deltaCt", "geometric.mean", "globalmean"))
	# Some general stuff that will be used by both rank.invariant methods
	if (method %in% c("scale.rankinvariant", "norm.rankinvariant")) {
		# Index to use for too high Ct values
		#Ct.index	<- data>Ct.max
		data.Ctmax	<- data
    Ct.index	<- is.na(data)
		#data.Ctmax[Ct.index]	<- NA
    data<-na.spline(data)
		# Define what to rank against
		if (rank.type=="pseudo.median") {
			ref.data	<- apply(data.Ctmax, 1, median, na.rm=TRUE)
		} else if (rank.type=="pseudo.mean") {
			ref.data	<- apply(data.Ctmax, 1, mean, na.rm=TRUE)
		}
		# Mark + replace NA values with something temporary
		na.index	<- is.na(ref.data)
		ref.data[na.index] <- 30
		# Run the rank.invariant function
		data.rankinvar	<- apply(data, 2, normalize.invariantset, ref=ref.data)
	}
	# The actual normalisation
	switch(method,
		quantile = {
			# Use an internal limma function
			data.norm	<- normalizeQuantiles(data)
		},
		scale.rankinvariant = {
			# Get all the rank invariant genes
			ri.genes	<- sapply(data.rankinvar, "[[", "i.set")
			# Remove those with too high Ct values
			ri.genes[Ct.index]	<- FALSE
			# Remove those that were all NA for potentially other reasons
			ri.genes[na.index,]	<- FALSE
			# Select those to use here
			ri.count	<- rowSums(ri.genes)
			if (missing(scale.rank.samples))
				scale.rank.samples	<- ncol(data)-1
			ri.index	<- ri.count >= scale.rank.samples
			if (sum(ri.index)==0)
				stop(paste("No rank invariant genes were found across", scale.rank.samples, "samples"))
			# Extract the corresponding Ct values; average
			ri.mean	<- colMeans(data[ri.index,,drop=FALSE])
			ri.scale	<- ri.mean/ri.mean[1]
			# Correct the data
			data.norm	<- t(t(data)*ri.scale)
			# Print info
			if (verbose) {
				cat(c("Scaling Ct values\n\tUsing rank invariant genes:", paste(featureNames(q)[ri.index], collapse=" "), "\n"))
				cat(c("\tScaling factors:", format(ri.scale, digits=3), "\n"))
			}
		},
		norm.rankinvariant = {
			# Print info
			if (verbose)
				cat("Normalizing Ct values\n\tUsing rank invariant genes:\n")
			# Correct the data based on the calculations above
			for (i in 1:ncol(data)) {
				# Check if there are any rank invariant genes
				ri.sub	<- data.rankinvar[[i]]
				ri.genes	<- ri.sub[["i.set"]]
				# Remove those that don't pass the Ct.max criteria
				ri.genes[Ct.index[,i]]	<- FALSE
				# Remove those that are NA for other reasons
				ri.genes[na.index]	<- FALSE
				if (sum(ri.genes)==0) {
	           		warning(paste("\tNo rank invariant genes were found for sample ", sampleNames(q)[i], "; sample not normalized\n", sep=""))
	           		next
	           	}
	           	# If verbose, print some info
	           	if (verbose)
	           		cat(paste("\t", sampleNames(q)[i], ": ", sum(ri.genes), " rank invariant genes\n", sep=""))
	           # The actual correction
	           data.norm[,i]	<- as.numeric(approx(ri.sub$n.curve$y, ri.sub$n.curve$x, xout=data[,i], rule=2)$y)
	       }
 		},
		deltaCt	= {
			# Which are the reference genes (endogenous controls)
			if (is.null(deltaCt.genes))
				deltaCt.genes	<- unique(featureNames(q)[featureType(q)=="Endogenous Control"])
			c.index	<- featureNames(q) %in% deltaCt.genes
			if (verbose) {
				cat(c("Calculating deltaCt values\n\tUsing control gene(s):", paste(deltaCt.genes, collapse=" "), "\n"))
			}
			# Run though all cards; perform internal normalisation
			for (c in 1:ncol(data)) {
				# Calculate the control genes
				refCt	<- mean(data[c.index,c], na.rm=TRUE)
				refsd	<- sd(data[c.index,c], na.rm=TRUE)
				# Difference for target and controls
				data.norm[,c] <- data[,c]-refCt
				# Print results
				if (verbose)
					cat(paste("\tCard ", c, ":\tMean=", format(refCt, dig=4), "\tStdev=", format(refsd, dig=3), "\n", sep=""))
			}
		},
		geometric.mean = {
			# For each column, calculate the geometric mean of Ct values<Ct.max
			#geo.mean	<- apply(data, 2, function(x) {
			#						xx <- log2(subset(x, x<Ct.max))
			#						2^mean(xx)})
      geo.mean	<- apply(data, 2, function(x) {
									xx <- subset(x, x<=Ct.max)
									geometric.mean(xx)})
			# Which sample to scale to
			#if (missing(geo.mean.ref))
			#	geo.mean.ref <- 1
			# Calculate the scaling factor
			#geo.scale	<- geo.mean/geo.mean[geo.mean.ref]
			# Adjust the data accordingly
			data.norm <- t(t(data) - geo.mean)
			if (verbose) {
				cat(c("Scaling Ct values\n\tUsing geometric mean within each sample\n"))
				#cat(c("\tScaling factors:", format(geo.scale, digits=3), "\n"))
			}
		} # switch
    ,globalmean = {
        glo <- apply(data, 2, function(x) {
            xx <- subset(x, x <= Ct.max)
            mean(xx)
        })
        data.norm <- t(t(data) - glo)
    }
	)
	# Replace with the normalised Ct exprs
	exprs(q)	<- data.norm
	# Add to the history of the object
	if (nrow(getCtHistory(q))==0)
		setCtHistory(q)	<- data.frame(history="Manually created qPCRset object.", stringsAsFactors=FALSE)
	setCtHistory(q)	<- rbind(getCtHistory(q), capture.output(match.call(normalizeCtData)))
	# Return the normalised object
	q
}
#library(NormqPCR)

#delete.na <- function(DF, n=0) {
 # DF[rowSums(is.na(DF)) <= n,]
#}

#user_number=5
#genorm <- selectHKs(t(delete.na(as.matrix(exprs(xGlico)),0)), method = "geNorm", Symbols = rownames(as.matrix(delete.na(exprs(xGlico),0))), minNrHK = as.numeric(user_number), log = TRUE)
#genorm
#normfinder <- selectHKs(as.matrix(t(delete.na(exprs(xGlico),0))), group= files$Treatment , method = "NormFinder", Symbols =rownames(as.matrix(delete.na(exprs(xGlico),0))), minNrHK = as.numeric(user_number), log = TRUE)
#normfinder
#intersection= intersect(normfinder$ranking, genorm$ranking[1:as.numeric(user_number)])

#cat("\n GeNorm and NormFinder transcripts selected as housekeeping for normalization! \n")
#intersection
#dnorm <- normalizeCtData(xGlico , norm="deltaCt", deltaCt.genes=as.vector(intersection)) 

switch(normalizationMethod,
    "deltaCt"={
      normalizedDataset <- normalizeCtDataDav(xFilter, norm="deltaCt",  deltaCt.genes =explode(normalizers, sep = ","))
    },
    "quantile"={
      normalizedDataset <- normalizeCtDataDav(xFilter, norm=normalizationMethod)
    },
    "scale.rankinvariant"={
      normalizedDataset <- normalizeCtDataDav(xFilter, norm=normalizationMethod)
    },
    "norm.rankinvariant"={
       normalizedDataset <- normalizeCtDataDav(xFilter, norm=normalizationMethod)
    },
    "geometric.mean"={
      normalizedDataset <- normalizeCtDataDav(xFilter, norm=normalizationMethod)
    },
    "globalmean"={
      normalizedDataset <- normalizeCtDataDav(xFilter, norm=normalizationMethod)
    },
    stop("Enter something that switches me!")
)

  #if (normalizationMethod=="deltaCt") {
#normalize CT data

#normalizedDataset <- normalizeCtDataDav(xFilter, norm="deltaCt",  deltaCt.genes =explode(normalizers, sep = ","))
#} else {
#normalizedDataset <- normalizeCtDataDav(xFilter, norm=normalizationMethod)

#}
cat("\n Data normalized correctly! \n")
write.table(exprs(normalizedDataset), file=outputNorm, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")


#normalizedDataset
####################################################################################################################
#Check noise reduction by empirical cumulative distribution

#X = rnorm(100) # X is a sample of 100 normally distributed random variables
# P = ecdf(X)    # P is a function giving the empirical CDF of X
#Y = rnorm(1000) # X is a sample of 100 normally distributed random variables
# PY = ecdf(Y)
#plot℗

#lines(PY)
png(outputECDF,    # create PNG for the heat map
  width = 10*300,        # 5 x 300 pixels
  height = 10*300,
  res = 300,            # 300 pixels per inch
  pointsize = 8)        # smaller font size
vec=c()
for (i in 1:nrow(exprs(xFilter))){
   CVX<-(sd(2^-exprs(xFilter)[i,], na.rm = TRUE)/mean(2^-exprs(xFilter)[i,], na.rm = TRUE))*100
vec=c(vec, c(CVX))
}
vec<-na.omit(vec)
P = ecdf(vec)
gm=c()
for (i in 1:nrow(exprs(normalizedDataset))){
   CVGM<-(sd(2^-exprs(normalizedDataset)[i,], na.rm = TRUE)/mean(2^-exprs(normalizedDataset)[i,], na.rm = TRUE))*100
gm=c(gm, c(CVGM))
}
gm<-na.omit(gm)

PY = ecdf(gm)

plot_colors <- c(rgb(r=0.0,g=0.0,b=0.9), "red", "forestgreen",rgb(r=0.0,g=0.0,b=0.0),rgb(r=0.5,g=0.0,b=0.3),rgb(r=0.0,g=0.4,b=0.4))

plot(P,col=plot_colors[1],xlim=c(0.0,600),  ylim=c(0.0,1),xaxp = c(0.0, 600, 6),yaxp = c(0.0, 1, 10), cex=1.3, lwd=5, main=NULL,xlab="CV(%)",ylab="Empirical Cumulative Distribution")
lines(PY, lwd=5, col=plot_colors[6],cex=1.3)
legend("bottomright", c("not normalized", "normalized"), cex=1.3, col=c(plot_colors[1],plot_colors[6]), lwd=c(5,5));
dev.off()

#Two-sample Kolmogorov-Smirnov
ks.test(vec,gm)


#write.table(exprs(qFiltNAs), file=outputIMP, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")
png(outputIMP,    # create PNG for the heat map
  width = 10*300,        # 5 x 300 pixels
  height = 10*300,
  res = 300,            # 300 pixels per inch
  pointsize = 8)
  plotCtBoxes(normalizedDataset, stratify=NULL, xlab = "Samples", ylab="DeltaCt", names=as.character(seq(1, ncol(normalizedDataset), 1)))       # smaller font size
dev.off()

################################################## Filtering based on number of NAs##################################################

#FILTERING on the basis of NAs
#qPCRset.R
setMethod("exprs", signature(object="qPCRset"), definition =
            function (object) {x <- assayDataElement(object, "exprs"); rownames(x) <- featureNames(object); colnames(x) <- sampleNames(phenoData(object));x}
)
ncatn <- as.integer(n.samples(normalizedDataset))*as.integer(percentofnastoremove)/100

qFiltNAs <- filterCtData(normalizedDataset, remove.category=c("Undetermined","Unreliable"), n.category=as.integer(ncatn),remove.name=explode(filtnames, sep = ","))

cat("\n Data filtered correctly! \n")

if (anyNA(exprs(qFiltNAs))){
switch(imputeMethod,
    "knn"={
     imp<-impute.knn(exprs(qFiltNAs) ,k = as.integer(kappa), maxp = as.integer(macsp), rng.seed=362436069)
      exprs(qFiltNAs)=imp$data
    },
    "mestdagh"={
     #Mesdagh
      #sostituisce a NA -1000
      for (i in 1:nrow(exprs(qFiltNAs))){
        for(j in 1:ncol(exprs(qFiltNAs))){
          if(is.na(exprs(qFiltNAs)[i,j])>0)exprs(qFiltNAs)[i,j]<- -1000
        }
      }
      temp<-exprs(qFiltNAs)
      for (i in 1:nrow(exprs(qFiltNAs))){
        for(j in 1:ncol(exprs(qFiltNAs))){
          if(exprs(qFiltNAs)[i,j]<(-100))exprs(qFiltNAs)[i,j]<- max(temp[i,])+1
        }
      }
    },
    "cubic"={
     exprs(qFiltNAs) <- na.spline(exprs(qFiltNAs))
    },
    "mean"={
      exprs(qFiltNAs)<-impute(exprs(qFiltNAs),mean)
    },
    "median"={
      exprs(qFiltNAs)<-impute(exprs(qFiltNAs),median)
    },
    "nondetects"={
      qFiltNAs <- qpcrImpute(qFiltNAs, outform=c("Single"),linkglm = c("logit"))
    },
    stop("Enter something that switches me!")
)

cat("\n Imputation completed! \n")
}else{
  cat("\n Nothing to impute! \n")
}

write.table(exprs(qFiltNAs), file=outputRemaining, quote=FALSE,  row.names=TRUE, col.names=TRUE, sep = "\t")

if (DEAMethod=="ttest") {
 #Differential expression analysis (paired t test+BH). Returns Fold change in linear scale.
 DEG<-ttestCtData(qFiltNAs, groups = files$Treatment, alternative = alternative, paired = ifelse(paired=="TRUE", TRUE, FALSE), replicates =ifelse(replicates=="TRUE", TRUE, FALSE), sort=ifelse(sort=="TRUE", TRUE, FALSE), stringent=ifelse(stringent=="TRUE", TRUE, FALSE), p.adjust=padjust)
 write.table(DEG, file=outputDEA, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")
} 
if (DEAMethod=="rp") {
 DEG<-RP(exprs(qFiltNAs), as.numeric(pData(qFiltNAs)$Treatment)-1, num.perm = 1000,logged = TRUE, gene.names = featureNames(qFiltNAs), huge=TRUE, plot = FALSE, rand = 123)
 write.table(DEG[1:5], file=outputDEA, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")
}
cat("\n Differential expression analysis completed correctly! \n")
cat("\n Workflow ended correctly! \n")
