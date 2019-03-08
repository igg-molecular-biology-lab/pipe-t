
cat("\n Started! ")
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
#library(gmp)
})

cat("\n R libraries...loaded!")

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
#cat("\n format ", format)
#cat("\n nfeatures ", nfeatures)
#cat("\n dpfiles ", dpfiles)
#cat("\n path000 ", path000)
#cat("\n rawout ", rawout)
#cat("\n path ", path)
#cat("\n dcCtmin ", dcCtmin)
#cat("\n dcCtmax ", dcCtmax)
#cat("\n dcflag ", dcflag)
#cat("\n x ", x)
#cat("\n normalizationMethod ", normalizationMethod)
if (normalizationMethod=="deltaCt") {
    normalizers<-args[11]
    outputNorm<-args[12]
    outputECDF<-args[13]
    percentofnastoremove<-args[14]
    outputRemaining<-args[15]
    imputeMethod<-args[16]
    #cat("\n normalizers ", normalizers)
    #cat("\n outputNorm ", outputNorm)
    #cat("\n outputECDF ", outputECDF)
    #cat("\n percentofnastoremove ", percentofnastoremove)
    #cat("\n outputRemaining ", outputRemaining)
    #cat("\n imputeMethod ", imputeMethod)
    if (imputeMethod=="knn") {
      kappa<- args[17]
      macsp<-args[18]
      outputIMP<-args[19]

      DEAMethod<-args[20]
      #cat("\n kappa ", kappa)
      #cat("\n macsp ", macsp)
      #cat("\n outputIMP ", outputIMP)
      #cat("\n DEAMethod ", DEAMethod)
      #cat("\n Dea ", DEAMethod)

        if (DEAMethod=="ttest") {
          alternative<- args[21]
          paired<-args[22]
          replicates<- args[23]
          sort<-args[24]
          stringent<- args[25]
          padjust<-args[26]
          outputDEA<-args[27]
          filtnames<-args[28]
          #cat("\n alternative ", alternative)
          #cat("\n paired ", paired)
          #cat("\n replicates ", replicates)
          #cat("\n sort ", sort)
          #cat("\n stringent ", stringent)
          #cat("\n padjust ", padjust)
          #cat("\n outputDEA ", outputDEA)
          #cat("\n filtnames ", filtnames)
        } else {
          outputDEA<-args[21]
          filtnames<-args[22]
          #cat("\n outputDEA ", outputDEA)
          #cat("\n filtnames ", filtnames)
        }
    } else {
      #globalmean, modified globalmean, etc
        outputIMP<-args[17]
        DEAMethod<-args[18]
        #cat("\n outputIMP ", outputIMP)
        #cat("\n DEAMethod ", DEAMethod)
          if (DEAMethod=="ttest") {
            alternative<- args[19]
            paired<-args[20]
            replicates<- args[21]
            sort<-args[22]
            stringent<- args[23]
            padjust<-args[24]
            outputDEA<-args[25]
            filtnames<-args[26]
            #cat("\n alternative ", alternative)
            #cat("\n paired ", paired)
            #cat("\n replicates ", replicates)
            #cat("\n sort ", sort)
            #cat("\n stringent ", stringent)
            #cat("\n padjust ", padjust)
            #cat("\n outputDEA ", outputDEA)
            #cat("\n filtnames ", filtnames)
          } else {
            outputDEA<-args[19]
            filtnames<-args[20]
            #cat("\n outputDEA ", outputDEA)
            #cat("\n filtnames ", filtnames)
          }
    }
  }else {
    outputNorm<-args[11]
    outputECDF<-args[12]
    percentofnastoremove<-args[13]
    outputRemaining<-args[14]
    imputeMethod<-args[15]
    #cat("\n outputNorm ", outputNorm)
    #cat("\n outputECDF ", outputECDF)
    #cat("\n percentofnastoremove ", percentofnastoremove)
    #cat("\n outputRemaining ", outputRemaining)
    #cat("\n imputeMethod ", imputeMethod)

    if (imputeMethod=="knn") {
      kappa<- args[16]
      macsp<-args[17]
      outputIMP<-args[18]

      DEAMethod<-args[19]
      #cat("\n kappa ", kappa)
      #cat("\n macsp ", macsp)
      #cat("\n outputIMP ", outputIMP)
      #cat("\n DEAMethod ", DEAMethod)
      if (DEAMethod=="ttest") {
          alternative<- args[20]
          paired<-args[21]
          replicates<- args[22]
          sort<-args[23]
          stringent<- args[24]
          padjust<-args[25]
          outputDEA<-args[26]
          filtnames<-args[27]
          #cat("\n alternative ", alternative)
          #cat("\n paired ", paired)
          #cat("\n replicates ", replicates)
          #cat("\n sort ", sort)
          #cat("\n stringent ", stringent)
          #cat("\n padjust ", padjust)
          #cat("\n outputDEA ", outputDEA)

        } else {
          outputDEA<-args[20]
          filtnames<-args[21]
          #cat("\n outputDEA ", outputDEA)
        }
     } else {
        outputIMP<-args[16]
        DEAMethod<-args[17]
        #cat("\n outputIMP ", outputIMP)
        #cat("\n DEAMethod ", DEAMethod)
          if (DEAMethod=="ttest") {
            alternative<- args[18]
            paired<-args[19]
            replicates<- args[20]
            sort<-args[21]
            stringent<- args[22]
            padjust<-args[23]
            outputDEA<-args[24]
            filtnames<-args[25]
            #cat("\n alternative ", alternative)
            #cat("\n paired ", paired)
            #cat("\n replicates ", replicates)
            #cat("\n sort ", sort)
            #cat("\n stringent ", stringent)
            #cat("\n padjust ", padjust)
            #cat("\n outputDEA ", outputDEA)
            #cat("\n filtnames ", filtnames)
          } else {
            outputDEA<-args[18]
            filtnames<-args[19]
            #cat("\n outputDEA ", outputDEA)
            #cat("\n filtnames ", filtnames)
          }
      }
  }
  cat("\n Initialization completed! ")
#cat("\n Norma ", as.character(normalizers))
#cat("\n base name ", dpfiles)
#cat("\n raw ", rawout)
#cat("\n paired ", paired)
#cat("\n rawout ", rawout)
#cat("\n min ", dcCtmin)
#cat("\n max ", dcCtmax)
#cat("\n flag ", dcflag)
#cat("\n x ", x)
#cat("\n dirname ", dirname(args[1]))
#<validator type="metadata" check="sampleName,Treatment" message="We accept only sampleName and Treatment column names at the moment">Wrong format</validator>
#  <validator type="regex" message="Should be of form name1,name2,name3 e.g U6 snRNA-001973 or U6 snRNA-001973,hsa-miR-328-000543">^[a-z][A-Z][0-9](,[a-z][A-Z][0-9])*$</validator>

head(read.delim(file.path(path000, dpfiles), sep="\t"))
files <- read.delim(file.path(path000, dpfiles), sep="\t")
switch(format,
    "plain"={
      columns<- list(flag="EXPFAIL", feature="Target.Name",  position="Well.Position", Ct="CT")
      metadata <- data.frame(labelDescription = c("sampleName", "Treatment"),  row.names = c("sampleName", "Treatment"))
      phenoData <- new("AnnotatedDataFrame", data = files, varMetadata = metadata)
      rownames(phenoData)=as.vector(files$sampleName)
      raw<- readCtData(files = as.vector(files$sampleName), header=TRUE,  format="plain", column.info=columns, path = path,sample.info=phenoData)
    },
    "SDS"={
      columns<- list(feature=3, Ct=6, flag=11)
      metadata <- data.frame(labelDescription = c("sampleName", "Treatment"),  row.names = c("sampleName", "Treatment"))
      phenoData <- new("AnnotatedDataFrame", data = files, varMetadata = metadata)
      rownames(phenoData)=as.vector(files$sampleName)
      raw<- readCtData(files = files$sampleName, format="SDS",column.info=columns, path = path, sample.info=phenoData, n.features=as.numeric(nfeatures))
    },
    "LightCycler"={
      raw <- readCtData(files = files$sampleName,, path = path, format = "LightCycler", n.features = as.numeric(nfeatures))
    },
    "CFX"={
      raw <- readCtData(files = files$sampleName,, path = path, format = "CFX", n.features = as.numeric(nfeatures))
    },
    "OpenArray"={
      raw <- readCtData(files = files$sampleName,, path = path, format = "OpenArray", n.features = as.numeric(nfeatures))
    },
    "BioMark"={
      raw <- readCtData(files = files$sampleName,, path = path, format = "BioMark", n.features = as.numeric(nfeatures))
    },
    stop("Enter something that switches me!")
)
cat("\n Files read correctly! ")

write.table(exprs(raw), file=rawout, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")
####################################################################################################################
#Set a new categories for the values meeting two criterions
switch(format,
    "plain"={
      unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="Y", quantile=NULL)
    },
    "SDS"={
      unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="FALSE", quantile=NULL)
    },
    "LightCycler"={
      unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="Y", quantile=NULL)
    },
    "CFX"={
      unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="Y", quantile=NULL)
    },
    "OpenArray"={
      unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="FALSE", quantile=NULL)
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
		Ct.index	<- data>Ct.max
		data.Ctmax	<- data
		data.Ctmax[Ct.index]	<- NA
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
			geo.mean	<- apply(data, 2, function(x) {
									xx <- log2(subset(x, x<Ct.max))
									2^mean(xx)})
			# Which sample to scale to
			if (missing(geo.mean.ref))
				geo.mean.ref <- 1
			# Calculate the scaling factor
			geo.scale	<- geo.mean/geo.mean[geo.mean.ref]
			# Adjust the data accordingly
			data.norm <- t(t(data) * geo.scale)
			if (verbose) {
				cat(c("Scaling Ct values\n\tUsing geometric mean within each sample\n"))
				cat(c("\tScaling factors:", format(geo.scale, digits=3), "\n"))
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

if (normalizationMethod=="deltaCt") {
#normalize CT data

normalizedDataset <- normalizeCtDataDav(xFilter, norm="deltaCt",  deltaCt.genes =explode(normalizers, sep = ","))
} else {
normalizedDataset <- normalizeCtDataDav(xFilter, norm=normalizationMethod)

}
cat("\n Data normalized correctly! ")
write.table(exprs(normalizedDataset), file=outputNorm, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")
 #cat("\n arg 7 ", args[7])
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
#cat("\n Significance of noise reduction ",  ks.test(vec,gm)$p.value)

################################################## 6 NAs##################################################

#FILTERING on the basis of NAs
#qPCRset.R
setMethod("exprs", signature(object="qPCRset"), definition =
            function (object) {x <- assayDataElement(object, "exprs"); rownames(x) <- featureNames(object); colnames(x) <- sampleNames(phenoData(object));x}
)
ncatn <- as.integer(n.samples(normalizedDataset))*as.integer(percentofnastoremove)/100

qFiltNAs <- filterCtData(normalizedDataset, remove.category=c("Undetermined","Unreliable"), n.category=as.integer(ncatn),remove.name=explode(filtnames, sep = ","))

cat("\n Data filtered correctly! ")

if (imputeMethod=="knn") {
 imp<-impute.knn(exprs(qFiltNAs) ,k = as.integer(kappa), maxp = as.integer(macsp), rng.seed=362436069)
 exprs(qFiltNAs)=imp$data
} else {
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

}
cat("\n Imputation completed! ")

#write.table(exprs(qFiltNAs), file=outputIMP, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")
png(outputIMP,    # create PNG for the heat map
  width = 10*300,        # 5 x 300 pixels
  height = 10*300,
  res = 300,            # 300 pixels per inch
  pointsize = 8)
  plotCtBoxes(qFiltNAs, stratify=NULL, xlab = "Samples", ylab="DeltaCt", names=as.character(seq(1, ncol(qFiltNAs), 1)))       # smaller font size
dev.off()

write.table(2^-exprs(qFiltNAs), file=outputRemaining, quote=FALSE,  row.names=TRUE, col.names=TRUE, sep = "\t")

if (DEAMethod=="ttest") {
 #Differential expression analysis (paired t test+BH). Returns Fold change in linear scale.
 DEG<-ttestCtData(qFiltNAs, groups = files$Treatment, alternative = alternative, paired = ifelse(paired=="TRUE", TRUE, FALSE), replicates =ifelse(replicates=="TRUE", TRUE, FALSE), sort=ifelse(sort=="TRUE", TRUE, FALSE), stringent=ifelse(stringent=="TRUE", TRUE, FALSE), p.adjust=padjust)
 write.table(DEG, file=outputDEA, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")
} else {
 DEG<-RP(exprs(qFiltNAs), as.numeric(pData(qFiltNAs)$Treatment)-1, num.perm = 1000,logged = TRUE, gene.names = featureNames(qFiltNAs), plot = FALSE, rand = 123)
 write.table(DEG[1:6], file=outputDEA, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")
}
cat("\n Differential expression analysis completed correctly! ")
cat("\n Workflow ended correctly! ")
