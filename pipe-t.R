#conda install -c r r=3.1.2

#cat("\n Entrato in IGG_RTqPCR_file_parser.R! ")

#dovremmo settore la WD di R nella job_working_directory unica creata per ogni sessione

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
#sessionInfo()
args = commandArgs(trailingOnly=TRUE)

dpfiles<-basename(args[1])
path000<-dirname(args[1])
rawout<-args[2]
path <- args[3]
dcCtmin<-args[4]
dcCtmax<-args[5]
dcflag<-args[6]
x<-args[7]
normalizationMethod<-args[8]
cat("\n dpfiles ", dpfiles)
cat("\n path000 ", path000)
cat("\n rawout ", rawout)
cat("\n path ", path)
cat("\n dcCtmin ", dcCtmin)
cat("\n dcCtmax ", dcCtmax)
cat("\n dcflag ", dcflag)
cat("\n x ", x)
cat("\n normalizationMethod ", normalizationMethod)
if (normalizationMethod=="globalmean") {
  outputNorm<-args[9]
  outputECDF<-args[10]
  percentofnastoremove<-args[11]
  outputRemaining<-args[12]
  imputeMethod<-args[13]
  cat("\n outputNorm ", outputNorm)
  cat("\n outputECDF ", outputECDF)
  cat("\n percentofnastoremove ", percentofnastoremove)
  cat("\n outputRemaining ", outputRemaining)
  cat("\n imputeMethod ", imputeMethod)

  if (imputeMethod=="knn") {
    kappa<- args[14]
    macsp<-args[15]
    outputIMP<-args[16]

    DEAMethod<-args[17]
    cat("\n kappa ", kappa)
    cat("\n macsp ", macsp)
    cat("\n outputIMP ", outputIMP)
    cat("\n DEAMethod ", DEAMethod)
    if (DEAMethod=="ttest") {
        alternative<- args[18]
        paired<-args[19]
        replicates<- args[20]
        sort<-args[21]
        stringent<- args[22]
        padjust<-args[23]
        outputDEA<-args[24]
        filtnames<-args[25]
        cat("\n alternative ", alternative)
        cat("\n paired ", paired)
        cat("\n replicates ", replicates)
        cat("\n sort ", sort)
        cat("\n stringent ", stringent)
        cat("\n padjust ", padjust)
        cat("\n outputDEA ", outputDEA)

      } else {
        outputDEA<-args[18]
        filtnames<-args[19]
        cat("\n outputDEA ", outputDEA)
      }
   } else {
      outputIMP<-args[14]
      DEAMethod<-args[15]
      cat("\n outputIMP ", outputIMP)
      cat("\n DEAMethod ", DEAMethod)
        if (DEAMethod=="ttest") {
          alternative<- args[16]
          paired<-args[17]
          replicates<- args[18]
          sort<-args[19]
          stringent<- args[20]
          padjust<-args[21]
          outputDEA<-args[22]
          filtnames<-args[23]
          cat("\n alternative ", alternative)
          cat("\n paired ", paired)
          cat("\n replicates ", replicates)
          cat("\n sort ", sort)
          cat("\n stringent ", stringent)
          cat("\n padjust ", padjust)
          cat("\n outputDEA ", outputDEA)
          cat("\n filtnames ", filtnames)
        } else {
          outputDEA<-args[16]
          filtnames<-args[17]
          cat("\n outputDEA ", outputDEA)
          cat("\n filtnames ", filtnames)
        }
    }
} else {
    normalizers<-args[9]
    outputNorm<-args[10]
    outputECDF<-args[11]
    percentofnastoremove<-args[12]
    outputRemaining<-args[13]
    imputeMethod<-args[14]
    cat("\n normalizers ", normalizers)
    cat("\n outputNorm ", outputNorm)
    cat("\n outputECDF ", outputECDF)
    cat("\n percentofnastoremove ", percentofnastoremove)
    cat("\n outputRemaining ", outputRemaining)
    cat("\n imputeMethod ", imputeMethod)
    if (imputeMethod=="knn") {
      kappa<- args[15]
      macsp<-args[16]
      outputIMP<-args[17]

      DEAMethod<-args[18]
      cat("\n kappa ", kappa)
      cat("\n macsp ", macsp)
      cat("\n outputIMP ", outputIMP)
      cat("\n DEAMethod ", DEAMethod)
      #cat("\n Dea ", DEAMethod)

        if (DEAMethod=="ttest") {
          alternative<- args[19]
          paired<-args[20]
          replicates<- args[21]
          sort<-args[22]
          stringent<- args[23]
          padjust<-args[24]
          outputDEA<-args[25]
          filtnames<-args[26]
          cat("\n alternative ", alternative)
          cat("\n paired ", paired)
          cat("\n replicates ", replicates)
          cat("\n sort ", sort)
          cat("\n stringent ", stringent)
          cat("\n padjust ", padjust)
          cat("\n outputDEA ", outputDEA)
          cat("\n filtnames ", filtnames)
        } else {
          outputDEA<-args[19]
          filtnames<-args[20]
          cat("\n outputDEA ", outputDEA)
          cat("\n filtnames ", filtnames)
        }
    } else {
        outputIMP<-args[15]
        DEAMethod<-args[16]
        cat("\n outputIMP ", outputIMP)
        cat("\n DEAMethod ", DEAMethod)
          if (DEAMethod=="ttest") {
            alternative<- args[17]
            paired<-args[18]
            replicates<- args[19]
            sort<-args[20]
            stringent<- args[21]
            padjust<-args[22]
            outputDEA<-args[23]
            filtnames<-args[24]
            cat("\n alternative ", alternative)
            cat("\n paired ", paired)
            cat("\n replicates ", replicates)
            cat("\n sort ", sort)
            cat("\n stringent ", stringent)
            cat("\n padjust ", padjust)
            cat("\n outputDEA ", outputDEA)
            cat("\n filtnames ", filtnames)
          } else {
            outputDEA<-args[17]
            filtnames<-args[18]
            cat("\n outputDEA ", outputDEA)
            cat("\n filtnames ", filtnames)
          }
    }
  }
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
#head(read.delim(file.path(path,"files_4d.txt"), sep="\t"))
#files <- read.delim(file.path(path, "files_4d.txt"), sep="\t")
columns<- list(flag="EXPFAIL", feature="Target.Name",  position="Well.Position", Ct="CT")
metadata <- data.frame(labelDescription = c("sampleName", "Treatment"),  row.names = c("sampleName", "Treatment"))
phenoData <- new("AnnotatedDataFrame", data = files, varMetadata = metadata)
rownames(phenoData)=as.vector(files$sampleName)
#cat("\n path ", path)
#cat("\n files$sampleName ", files$sampleName)
#cat("\n columns ", columns)
#cat("\n phenoData ", phenoData)
#write.table(files$sampleName, file=rawout, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")
raw<- readCtData(files = as.vector(files$sampleName), header=TRUE,  format="plain", column.info=columns, path = path,sample.info=phenoData)
write.table(exprs(raw), file=rawout, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")
####################################################################################################################
#Set a new categories for the values meeting two criterions
#raw
unreliable<-setCategory(raw, Ct.max=dcCtmax, Ct.min=dcCtmin,replicates=FALSE,  flag=dcflag, flag.out="Y", quantile=NULL)

####################################################################################################################
#Filter out the values of the new category

xFilter <- filterCategory(unreliable)
#xFilter
write.table(exprs(xFilter), file=x, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")

####################################################################################################################
#NORMALIZATION
#redefine the globalmean normalization function
normalizeCtDataDav<-function (q, norm = "deltaCt", deltaCt.genes = NULL, scale.rank.samples,
    rank.type = "pseudo.median", Ct.max = dcCtmax, geo.mean.ref, verbose = TRUE)
{
    data <- exprs(q)
    data.norm <- data
    method <- match.arg(norm, c("quantile", "scale.rankinvariant",
        "norm.rankinvariant", "deltaCt", "geometric.mean", "globalmean"))
    if (method %in% c("scale.rankinvariant", "norm.rankinvariant")) {
        Ct.index <- data > Ct.max
        data.Ctmax <- data
        data.Ctmax[Ct.index] <- NA
        if (rank.type == "pseudo.median") {
            ref.data <- apply(data.Ctmax, 1, median, na.rm = TRUE)
        }
        else if (rank.type == "pseudo.mean") {
            ref.data <- apply(data.Ctmax, 1, mean, na.rm = TRUE)
        }
        na.index <- is.na(ref.data)
        ref.data[na.index] <- 30
        data.rankinvar <- apply(data, 2, normalize.invariantset,
            ref = ref.data)
    }
    switch(method, quantile = {
        data.norm <- normalizeQuantiles(data)
    }, scale.rankinvariant = {
        ri.genes <- sapply(data.rankinvar, "[[", "i.set")
        ri.genes[Ct.index] <- FALSE
        ri.genes[na.index, ] <- FALSE
        ri.count <- rowSums(ri.genes)
        if (missing(scale.rank.samples)) scale.rank.samples <- ncol(data) -
            1
        ri.index <- ri.count >= scale.rank.samples
        if (sum(ri.index) == 0) stop(paste("No rank invariant genes were found across",
            scale.rank.samples, "samples"))
        ri.mean <- colMeans(data[ri.index, , drop = FALSE])
        ri.scale <- ri.mean/ri.mean[1]
        data.norm <- t(t(data) * ri.scale)
        if (verbose) {
            cat(c("Scaling Ct values\n\tUsing rank invariant genes:",
                paste(featureNames(q)[ri.index], collapse = " "),
                "\n"))
            cat(c("\tScaling factors:", format(ri.scale, digits = 3),
                "\n"))
        }
    }, norm.rankinvariant = {
        if (verbose) cat("Normalizing Ct values\n\tUsing rank invariant genes:\n")
        for (i in 1:ncol(data)) {
            ri.sub <- data.rankinvar[[i]]
            ri.genes <- ri.sub[["i.set"]]
            ri.genes[Ct.index[, i]] <- FALSE
            ri.genes[na.index] <- FALSE
            if (sum(ri.genes) == 0) {
                warning(paste("\tNo rank invariant genes were found for sample ",
                  sampleNames(q)[i], "; sample not normalized\n",
                  sep = ""))
                next
            }
            if (verbose) cat(paste("\t", sampleNames(q)[i], ": ",
                sum(ri.genes), " rank invariant genes\n", sep = ""))
            data.norm[, i] <- as.numeric(approx(ri.sub$n.curve$y,
                ri.sub$n.curve$x, xout = data[, i], rule = 2)$y)
        }
    }, deltaCt = {
        if (is.null(deltaCt.genes)) deltaCt.genes <- unique(featureNames(q)[featureType(q) ==
            "Endogenous Control"])
        c.index <- featureNames(q) %in% deltaCt.genes
        if (verbose) {
            cat(c("Calculating deltaCt values\n\tUsing control gene(s):",
                paste(deltaCt.genes, collapse = " "), "\n"))
        }
        for (c in 1:ncol(data)) {
            refCt <- mean(data[c.index, c], na.rm = TRUE)
            refsd <- sd(data[c.index, c], na.rm = TRUE)
            data.norm[, c] <- data[, c] - refCt
            if (verbose) cat(paste("\tCard ", c, ":\tMean=",
                format(refCt, dig = 4), "\tStdev=", format(refsd,
                  dig = 3), "\n", sep = ""))
        }
    }, geometric.mean = {
        geo.mean <- apply(data, 2, function(x) {
            xx <- log2(subset(x, x < Ct.max))
            2^mean(xx)
        })
        if (missing(geo.mean.ref)) geo.mean.ref <- 1
        geo.scale <- geo.mean/geo.mean[geo.mean.ref]
        data.norm <- t(t(data) * geo.scale)
        if (verbose) {
            cat(c("Scaling Ct values\n\tUsing geometric mean within each sample\n"))
            cat(c("\tScaling factors:", format(geo.scale, digits = 3),
                "\n"))
        }
    },globalmean = {
        glo <- apply(data, 2, function(x) {
            xx <- subset(x, x <= Ct.max)
            mean(xx)
        })
        data.norm <- t(t(data) - glo)

    })
    exprs(q) <- data.norm
    if (nrow(getCtHistory(q)) == 0)
        setCtHistory(q) <- data.frame(history = "Manually created qPCRset object.",
            stringsAsFactors = FALSE)
    setCtHistory(q) <- rbind(getCtHistory(q), capture.output(match.call(normalizeCtData)))
    q
}

if (normalizationMethod=="globalmean") {
#normalize CT data
normalizedDataset <- normalizeCtDataDav(xFilter, norm=normalizationMethod)

} else {
normalizedDataset <- normalizeCtDataDav(xFilter, norm="deltaCt",  deltaCt.genes =explode(normalizers, sep = ","))

}
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

ncatn <- as.integer(n.samples(normalizedDataset))*as.integer(percentofnastoremove)/100

qFiltNAs <- filterCtData(normalizedDataset, remove.category=c("Undetermined","Unreliable"), n.category=as.integer(ncatn),remove.name=explode(filtnames, sep = ","))

#sampleNames(qFiltNAs)=sampleNames(normalizedDataset)

#write.table(exprs(qFiltNAs), file=outputRemaining, quote=FALSE,  row.names=TRUE, col.names=TRUE, sep = "\t")
#cat("\n names ",  colnames(qFiltNAs))
#cat("\n names ",  colnames(normalizedDataset))
#cat("\n impute method ", imputeMethod)
#qFiltNAs

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
#pData(qFiltNAs)
#sampleNames(qFiltNAs)=sampleNames(normalizedDataset)
write.table(exprs(qFiltNAs), file=outputIMP, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")

write.table(2^-exprs(qFiltNAs), file=outputRemaining, quote=FALSE,  row.names=TRUE, col.names=TRUE, sep = "\t")
#cat("\n ncol", ncol(qFiltNAs))
#cat("\n gloval", ncol(normalizedDataset))
if (DEAMethod=="ttest") {
 #Differential expression analysis (paired t test+BH). Returns Fold change in linear scale.
 DEG<-ttestCtData(qFiltNAs, groups = pData(qFiltNAs)$Treatment, alternative = alternative, paired = ifelse(paired=="TRUE", TRUE, FALSE), replicates =replicates, sort=sort, stringent=stringent, p.adjust=padjust)

} else {
 #library(RankProd)

 #ATTENZIONE: QUESTO DIPENDE DA QUANTI CAMPIONI HO IN CIASCUN GRUPPO (1 è trattato e 0 è controllo )dobbiamo dirgli le classi controllo e trattato
 #cl<-c(rep(1,5),rep(0,10))

 DEG<-RP(exprs(qFiltNAs), as.numeric(pData(qFiltNAs)$Treatment)-1, num.perm = 1000,logged = TRUE, gene.names = featureNames(qFiltNAs), plot = FALSE, rand = 123)
}

write.table(DEG, file=outputDEA, quote=FALSE,  row.names=TRUE, col.names=TRUE,sep = "\t")
