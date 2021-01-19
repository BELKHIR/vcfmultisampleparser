#This version will use calls to vcftools to get  statistics ovec multi samples vcf file

generate_stats <- function(fic){
    vcf.fn = fic
    
    cmd = paste0("vcftools --gzvcf ",vcf.fn," --geno-depth --stdout ")
    ttt = system(cmd, intern=T)
    dp = data.frame(matrix(unlist(strsplit(ttt,'\t')), nrow=length(ttt), byrow=T),  stringsAsFactors=F)
    colnames(dp) = dp[1,]
    dp=dp[-1,-c(1,2)]
    dp[dp == -1] <- NA
    depth = boxplot(data.matrix(dp), plot=F)
    
    cmd = paste0("vcftools --gzvcf ",vcf.fn," --missing-indv  --stdout | cut -f5")
    ttt = system(cmd, intern=T)
    #Missingness <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) })
    #Missingness <- Missingness/nrow(dp)
    Missingness = as.numeric(ttt[-1])

   #  the inbreeding coefficient F from --het i.e. ([observed hom. count] - [expected hom. count]) / ([total observations] - [expected hom. count]))
    cmd = paste0("vcftools --gzvcf ",vcf.fn," --het  --stdout ")
    ttt = system(cmd, intern=T)
    F = data.frame(matrix(unlist(strsplit(ttt,'\t')), nrow=length(ttt), byrow=T),  stringsAsFactors=F)
    colnames(F) = F[1,]
    F=F[-1,5]
    samples = length(F)

    cmd = paste0("vcftools --gzvcf ",vcf.fn," --missing-site  --stdout | cut -f6")
    ttt = system(cmd, intern=T)
    snps = length(ttt) - 1

    SiteMissingness = hist(as.numeric(ttt[-1]) )

    cmd = paste0("vcftools --gzvcf ",vcf.fn," --site-pi  --stdout | cut -f3")
    ttt = system(cmd, intern=T)
    SitePi = hist(as.numeric(ttt[-1]) )


    cmd = paste0("vcftools --gzvcf ",vcf.fn," --freq2 --stdout | awk '{if(NR>1) {if($5<$6) {print $5} else {print $6} } }'")
    maf = system(cmd, intern=T)

    return(list(samples = samples, snps = snps, depth=depth,missingness=Missingness, maf=as.numeric(maf), sitemissingness = SiteMissingness, F=as.numeric(F), SitePi=SitePi ) )
}

args = commandArgs(trailingOnly=TRUE)

vcf.fn = args[1] #a vcf file optionnaly gziped
basefigres=args[2] #base name for figures (will not call plotVCF.R if if empty ) 

resu = generate_stats(vcf.fn) 
depth=resu$depth
missingness = resu$missingness
maf = resu$maf
sitemissingness = resu$sitemissingness
SitePi = resu$SitePi
F = resu$F
samples = resu$samples
snps = resu$snps

minmafF = 0
maxmissingF = 0
png(file=paste0(basefigres,"_vcftools_plots_mqc.png"), width = 1200, height = 1000)
par(mfrow=c(3,2)) 
bxp(depth, outline=FALSE, main=paste0("Per sample SNP Depth (DP) distrib. minMAF:", minmafF, " maxMissing:", maxmissingF),  boxfill=2:8, las=3 )
hist(maf, breaks=20, xlim=c(min(maf),0.5), main="Histogram of minor allele frequency across all SNP", xlab="maf frequency", col= rgb(1,0,0,1/4)) 
barplot(missingness, main=paste0("Per-sample missingness (%)"),  col=2:8, las=3 )
plot(sitemissingness, main="Per-site missingness", xlab="% missing",col=rgb(0,0,1,1/4))
barplot(F, main=paste0("Per-sample F (inbreeding Coef using a method of moments) "),  col=2:8, las=3 )
plot(SitePi, main="Per-site Pi", xlab="Pi",col=rgb(0,0,1,1/4))
par()
dev.off()