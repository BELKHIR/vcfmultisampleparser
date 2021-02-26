# This version will  get  statistics over multi samples from a pre-processed vcf file
# pre-processing must be maid by parseMultisampleVCF.sh

# usage starting with a vcf file in 2 steps
# get pre-processed file from vcf file optionnaly gziped
# bash ./parseMultiSamplesVCF.sh myfile.vcf.gz myfile.processed.tsv

# Generate multi samples summary plots in a file named myfileFig.png
# Rscript myfile.processed.tsv myfileFig

generate_stats <- function(ficstats){
   if (is.null(ficstats)) return()

    #determine the number of samples from the header
    df = read.table(ficstats, nrows=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    ncols  = dim(df)[2]
    
    nsamples = (dim(df)[2] - 6)/5
    startSample =  1 
    startLocusInfo = nsamples * 5 + 1

    #Read Pi, Maf and Missing values
    coltypes = rep('_', ncols)
    coltypes[c( ncols-2, ncols-1, ncols)]=c('n','n','i')
    df <- readr::read_tsv(ficstats, col_types = paste(coltypes, collapse=""),na = c("", "NA","."))
    nsites = dim(df)[1]
    
    SitePi = hist(df$Pi, plot=FALSE)
    maf = hist(as.numeric(df$MAF), breaks=20, plot=FALSE) 
    freq=df$MAF
    SiteMissingness = hist(df$Miss/nsamples, plot=FALSE)
    rm(df)
    gc()

    # Read All DP values
    coltypes = rep('_', ncols)
    DPCols = seq(2,startLocusInfo -4, 5)

    coltypes[DPCols] = 'i'
    df <- readr::read_tsv(ficstats, col_types = paste(coltypes, collapse=""),na = c("", "NA","."))
    colnames(df)[1:nsamples] = stringr::str_sub(colnames(df)[1:nsamples], 1, -4) # remove trailing .DP in sample names
    if (all(is.na(df))) depth=NULL else    depth = boxplot(df[1:min(1000000, nsites),], plot=FALSE) # limit to the first 1 million snp
    rm(df)
    gc()
  
    # Read GT values
    GTCols = seq(1,startLocusInfo - 4, 5)
    coltypes = rep('_', ncols)
    coltypes[GTCols] = 'i'
    df <- readr::read_tsv(ficstats, col_types = paste(coltypes, collapse=""))

    colnames(df)=stringr::str_sub(colnames(df), 1, -4)
    samples = colnames(df)
    N_missing_sites = apply(df, 2, function(x) sum(is.na(x) ))
    Missingness =  N_missing_sites / nsites
    rm(N_missing_sites)
    #  the inbreeding coefficient F from --het i.e. ([observed hom. count] - [expected hom. count]) / ([total observations] - [expected hom. count]))
    # From vcftools  variant_file_output.cpp#L165
    #// P(Homo) = F + (1-F)P(Homo by chance)
    # // P(Homo by chance) = p^2+q^2 for a biallelic locus.
    # // For an individual with N genotyped loci, we
    # //   1. count the total observed number of loci which are homozygous (O),
    # //   2. calculate the total expected number of loci homozygous by chance (E)
    # // Then, using the method of moments, we have
    # //    O = NF + (1-F)E
    # // Which rearranges to give
    # //    F = (O-E)/(N-E)


    biallelic_sites = apply(df, 1, function(x) {nballeles = length(unique(unlist(strsplit(names(table(x)),"")))); if (nballeles == 2) return(TRUE) else return(FALSE) } )
    df =  df[biallelic_sites,] #only bi-allelic
    freq = freq[biallelic_sites] #this is from maf 
    rm(biallelic_sites)
    gc()

    N_non_missing_site = apply(df, 1, function(x) sum( ! is.na(x)) )
    site_expected_hom  = 1.0 - (2.0 * freq * (1.0 - freq) * (2*N_non_missing_site / (2*N_non_missing_site - 1.0)));
    rm(freq,N_non_missing_site)
    gc()
    N_obs_hom          = apply(df, 2, function(x) sum(x %in% c(0,11))) # sum(x == "0/0" | x=="0|0" | x == "1/1" | x == "1|1"))
    N_expected_hom     = apply(df, 2, function(x, expected_hom) { nonMissing = ! is.na(x); return(sum(expected_hom[nonMissing], na.rm=T)) }, expected_hom = site_expected_hom )
    N_non_missing_ind  = apply(df, 2, function(x) sum( ! is.na(x) ) )
    rm(df); gc()
    F = (N_obs_hom - N_expected_hom) / (N_non_missing_ind - N_expected_hom);

    return(list(samples = samples,  depth=depth,missingness=Missingness, maf=maf, sitemissingness = SiteMissingness, F=F, SitePi=SitePi ) )
}

args = commandArgs(trailingOnly=TRUE)

ficstats = args[1] #a file resulting from parseMultisampleVCF.R
basefigres=args[2] #base name for figures (will not call plotVCF.R if its empty ) 

vcftools_summary = generate_stats(ficstats) 
samples = 1:length(vcftools_summary$samples)
names(samples) = vcftools_summary$samples

depth=vcftools_summary$depth
missingness = vcftools_summary$missingness
maf = vcftools_summary$maf
sitemissingness = vcftools_summary$sitemissingness
SitePi = vcftools_summary$SitePi
F = vcftools_summary$F
if(is.null(depth)) { 
  depth=boxplot(c(0,0,0,0),plot=F)
  titre="No sample DP available! "
} else titre=paste0("Per sample SNP Depth (DP) distrib. for the first 1-e6 snp")

png(file=paste0(basefigres,"_vcftools_plots_mqc.png"), width = 1200, height = 1000)
par(mfrow=c(3,2))
bxp(depth, outline=FALSE, ylab=titre,  boxfill=2:8, las=3 )
plot(maf, main="Histogram of minor allele frequency across all SNP", xlab="maf frequency", col= rgb(1,0,0,1/4))
barplot(missingness, ylab=paste0("Per-sample missingness (%)"),  col=2:8, las=3 )

plot(sitemissingness, main="Per-site missingness", xlab="% missing",col=rgb(0,0,1,1/4))
barplot(as.numeric(F), ylab=paste0("Per-sample F (inbreeding Coef using a method of moments) "),  col=2:8, las=3 )

plot(SitePi, main="Per-site Pi", xlab="Pi",col=rgb(0,0,1,1/4))
mtext(ficstats, side = 3, line = -2, outer = TRUE)
par()
dev.off()