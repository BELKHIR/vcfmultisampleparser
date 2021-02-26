library(ggplot2)
library(tidyverse)
library("patchwork")

#vérifier l'usage de ALT dans un contexte multi-samples !
Draw <- function(inFile, sample, Z=NULL)
{   
            if (is.null(Z)) {
            print(paste0("reading data: ", inFile) )
            #determine the number of samples from the header
            df = read.table(inFile, nrows=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
            nsamples = (dim(df)[2] - 6)/5
            #filter data for the given sample
            startSample = (sample -1) * 5 + 1 
            startLocusInfo = nsamples * 5 + 1
            needCols = c(startSample:(startSample + 4), startLocusInfo:(startLocusInfo +2) )
            #AllData=read.table(inFile, head=T, sep="\t", stringsAsFactors=F)
            #Z = fread(inFile, header=T, select = needCols, sep="\t", stringsAsFactors=F)

            coltypes = rep('_', dim(df)[2])
            coltypes[needCols] = c('i','n','c','n','c','c','c','c')
            #Z = vroom::vroom(inFile, col_types = paste(coltypes, collapse=""),  delim="\t")
            Z <- readr::read_tsv(inFile, col_types = paste(coltypes, collapse=""),na = c("", "NA","."))

            }
            
            #get sample name
            sampleName = stringr::str_sub(colnames(Z)[1], 1, -4) # remove trailing .GT in sample name

            colnames(Z) = c('GT','DP','AD','GQ','TYPE','QUAL','REF','ALT')
            if (all(is.na(Z$GQ))) {
            p1 = ggplot() + ggtitle("Genotype quality not available") + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
            } else {
            GQ99= quantile(Z$GQ, 0.99, na.rm=T)
            binWidth = ifelse(GQ99>30,round(GQ99/30), 1)
            GQvals = select(Z, GQ) %>% filter( ! is.na(GQ)) %>% mutate(GQ.class=cut(GQ,breaks=c(seq(0,GQ99-1,binWidth), max(GQ,na.rm=T)+1), include.lowest = T)) %>% group_by(GQ.class) %>% summarise(count=n())


            labs = levels(GQvals$GQ.class)
            GQmean = mean(Z$GQ,na.rm=T)
            x=as.data.frame (cbind(rowmean=1:length(labs), lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ), upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))) %>% filter(lower < GQmean & upper >= GQmean)

            p1 = ggplot(GQvals,aes(GQ.class,count))+ geom_bar(stat="identity",width=1, fill = "#0072B2") + xlab("GQ") + ggtitle("Genotype quality") + 
                theme(plot.title = element_text(hjust = 0.5),  axis.text.x  = element_text(angle=90, vjust=0.5), panel.border = element_rect(colour = "black", fill=NA, size=2)) + 
                geom_vline(xintercept = x$rowmean, color="red" ) + 
                geom_text(aes(x=x$rowmean, label="\nmean", y=round(max(count)/2)), colour="red", angle=90) 
            rm(GQvals)    
            }

            #Z$gt <-  mutate(Z$gt, x_new = ifelse(gt_DP > 100, 100, gt_DP))

            if (all(is.na(Z$DP))){
            p2 = ggplot() + ggtitle("Depth not available") + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
            } else {
            DP99 = quantile(as.numeric(Z$DP), 0.99, na.rm=T)
            binWidth = ifelse(DP99 > 30, round(DP99/30), 1)
            DPvals <- select(Z, DP) %>%  mutate(DP=as.numeric(DP) ) %>% filter( ! is.na(DP))  %>% mutate(DP.class=cut(DP,breaks=c(seq(0,DP99-1,binWidth), max(DP)+1), include.lowest = T)) %>% group_by(DP.class) %>% summarise(count=n())

            labs = levels(DPvals$DP.class)
            DPmean = mean(as.numeric(Z$DP),na.rm=T)
            x=as.data.frame (cbind(rowmean=1:length(labs), lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ), upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))) %>% filter(lower < DPmean & upper >= DPmean)

            p2 = ggplot(DPvals,aes(DP.class,count))+geom_bar(stat="identity",width=1, fill = "#0072B2")  + xlab("DP") + ggtitle("Depth") + 
                theme(plot.title = element_text(hjust = 0.5),  axis.text.x  = element_text(angle=90, vjust=0.5), panel.border = element_rect(colour = "black", fill=NA, size=2))  +
                geom_vline(xintercept = x$rowmean, color="red" ) + 
                geom_text(aes(x=x$rowmean, label="\nmean", y=round(max(count)/2)), colour="red", angle=90) 
            rm(DPvals)  
            #mean_class = 2 trouver le moyen repérer la classe qui contient la moyenne
            #p2 <- p + geom_vline(aes(xintercept=mean_class), color="blue", linetype="dashed", size=1)

            }
            
            if (all(is.na(Z$QUAL)) | all(Z$QUAL==".")){
            p3 = ggplot() + ggtitle("QUAL not available") + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
            }else {
            QUAL99 = quantile(as.numeric(Z$QUAL), 0.99, na.rm=T)
            binWidth = ifelse(QUAL99 > 30, round(QUAL99/30), 1)
            QUALvals <-  select(Z, QUAL)  %>% mutate(QUAL=as.numeric(QUAL) ) %>% filter( ! is.na(QUAL)) %>% mutate(qual.class = cut(QUAL,breaks=c(seq(0,QUAL99-1,binWidth),max(QUAL)+1), include.lowest = T )) %>% group_by(qual.class) %>% summarise(count=n())

            labs = levels(QUALvals$qual.class)
            QUALmean = mean(as.numeric(Z$QUAL),na.rm=T)
            x=as.data.frame (cbind(rowmean=1:length(labs), lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ), upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))) %>% filter(lower < QUALmean & upper >= QUALmean)

            p3 = ggplot(QUALvals, aes(qual.class, count)) + geom_bar(stat="identity",width=1,fill = "#0072B2")  + xlab("QUAL") + ggtitle("Quality score") + 
                theme(plot.title = element_text(hjust = 0.5),  axis.text.x  = element_text(angle=90, vjust=0.5),panel.border = element_rect(colour = "black", fill=NA, size=2)) +
                geom_vline(xintercept = x$rowmean, color="red" ) + 
                geom_text(aes(x=x$rowmean, label="\nmean", y=round(max(count)/2)), colour="red", angle=90) 
            rm(QUALvals)
            }
            p4 = Z %>%  group_by(TYPE) %>% count() %>% arrange(n) %>% ggplot(., aes(x=TYPE, y=n, label=n, fill=TYPE)) + geom_bar(stat = "identity")+ geom_text(size = 3, position = position_stack(vjust = 0.5),col="black") + ggtitle(paste0(sampleName,": Variant types")) + theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=2)) + scale_fill_discrete(name="Variant types",  breaks=c("bins", "bdel", "snp",  "snpmulti", "mnp","multiIcmpl", "ref"),  labels=c("Biallelic_insertion", "Biallelic_deletion", "Biallelic snp", "Multi snp", "mnp", "MultiIndel","RefCall"))
            
 
gc()

if (all(is.na(Z$GT)) | all(is.na(Z$AD))){
  p5 = ggplot() + ggtitle("Ref (0/0)") + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
  p6 = ggplot() + ggtitle("Het (0/x)") + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
  p7 = ggplot() + ggtitle("Hom (x/x)") + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
  p8 = ggplot() + ggtitle("Het two variants (x/y)") + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
} else {
  # The histograms show the variant allele frequency (VAF) distributions for different genotypes. Black guiding lines are shown to indicate the theoretical VAF for the main genotypes. For example heterozygous variants should have about as many variant-supporting as reference-supporting reads, for a VAF of 0.5. The reference calls will not usually show a VAF as low as 0 because otherwise they wouldn't have been flagged as candidates in the first place. The genotypes are based on the GT sub-column and consolidated. For example, 0/1 and 0/2 both become Het (0/x). 1/1 and 3/3 are Hom (x/x). Het - both variants (x/y) includes all calls with two different alternate alleles, such as 1/2 or 3/5.

  # attention DP peut être > à la somme des AD (AD ne contient que les comptes des allèles sûr (filtrés) alors que DP compte tous les allèles)
  # https://gatk.broadinstitute.org/hc/en-us/articles/360035532252-Allele-Depth-AD-is-lower-than-expected
  vaf <- function(ad, gt) {
    totalAlleles = sum(as.numeric(unlist(strsplit(ad,",")) ))
    z = 0
    if (gt==0){ #(gt == "0/0") {
        col = 2;
        x = unlist(strsplit(ad,","))[col]
        #z <- as.numeric(x)/as.numeric(dp)
        z <-  as.numeric(x)/ totalAlleles
    }
    else{
      #als = as.numeric(unique(unlist(strsplit(gt,"/"))  ))
      als = as.numeric(unique(unlist(strsplit(as.character(gt),""))  ))
      als = als[als != 0]
      x = unlist(strsplit(ad,","))[als+1]
      # si deux allèles on somme les freq. !!!
      #z <- sum(as.numeric(x))/as.numeric(dp)
      z <-  sum(as.numeric(x)) / totalAlleles
    }
    #if (z > 1) cat (ad, dp, gt,"\n", sep="\t")
    return(z)
  }

  # genotype are coded as integers : 0/0 -> 00, 0/1 or 1/0 -> 01, 1/1 -> 11, 1/2 -> 12 
  # no values like 10 or 20 ... as alleles are sorted
  # the phase is lost !
  # if readed as integers : 0 REF REF, 1-9 REF ALT, 11 22 homz. ALT , 12 13 23 hetero 2 for ALT 
  # homzygote 0/0 Have we to include ./. ??  GATK ne sort pas les 0/0 par défaut !!
  df = Z %>%   filter(GT == 0  & TYPE=="ref") %>% select(GT, AD)
  if(nrow(df) > 0){
  vaf00 = mapply(vaf, df$AD, df$GT) %>% tibble::enframe(name=NULL, value="vaf00") %>% filter( ! is.nan(vaf00))
  p5 = ggplot(vaf00, aes(vaf00)) +  geom_histogram(binwidth=0.02,fill = "#0072B2") + ggtitle("Ref (0/0)") + 
       theme(plot.title = element_text(hjust = 0.5),panel.border = element_rect(colour = "black", fill=NA, size=2)) + xlab("VAF") +
       geom_vline(xintercept = 0, color="red" )  
  rm(vaf00)     
  }   else{
      p5 = ggplot() + ggtitle("No Ref (0/0) found ") + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
  }
     

  # hetero 0/X (x = 1 or 2 or 3)
  #df = Z %>%   filter((GT == "0/1" | GT == "0/2" | GT == "0/3" | GT == "0/4") & TYPE=="snp")  %>% select(GT, AD, DP)
  df = Z %>%   filter((GT >0 & GT <= 9) & TYPE=="snp")  %>% select(GT, AD)
  if(nrow(df) > 0){
  vaf0X = mapply(vaf,  df$AD, df$GT) %>% tibble::enframe(name=NULL, value="vaf0X") %>% filter( ! is.nan(vaf0X))
  p6 = ggplot(vaf0X, aes(vaf0X)) +  geom_histogram(binwidth=0.02,fill = "#0072B2") + ggtitle("Het (0/x)") +
       theme(plot.title = element_text(hjust = 0.5),panel.border = element_rect(colour = "black", fill=NA, size=2)) + xlab("VAF") +
       geom_vline(xintercept = 0.5, color="red" ) 
  rm(vaf0X)      
  }  else{
      p6 = ggplot() + ggtitle("No Het (0/x) found ") + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
  } 

  # homozygotes X/X  (x = 1 or 2 or 3 ...)
  homoALT = seq(11, 99, 11) # 11, 22, 33 .. 99
  df = Z %>%   filter((GT %in% homoALT)  & TYPE=="snp" ) %>% select(GT, AD)
  if(nrow(df) > 0){
  vafXX = mapply(vaf,  df$AD, df$GT) %>% tibble::enframe(name=NULL, value="vafXX") %>% filter( ! is.nan(vafXX))
  p7 = ggplot(vafXX, aes(vafXX)) +  geom_histogram(binwidth=0.02,fill = "#0072B2") + ggtitle("Hom (x/x)") + 
       theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=2)) + xlab("VAF") +
       geom_vline(xintercept = 1, color="red" )
  rm(vafXX)  
  }  else{
        p7 = ggplot() + ggtitle("No Hom (x/x) found ") + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))
  }

  # hetero X/Y (x and y = 1 or 2 or 3 ...)
  # GT = 12, 13, 14, ... 23, 24 ... 45, ... 89
  df = Z %>%   filter( (GT > 9 &  !(GT %in% homoALT)   )  & TYPE=="snp" ) %>% select(GT, AD)
  if(nrow(df) > 0){
  vafXY = mapply(vaf,  df$AD, df$GT) %>% tibble::enframe(name=NULL, value="vafXY") %>% filter( ! is.nan(vafXY))

  p8 = ggplot(vafXY, aes(vafXY)) +  geom_histogram(binwidth=0.02,fill = "#0072B2") + ggtitle("Het two variants (x/y)") + theme(plot.title = element_text(hjust = 0.5),panel.border = element_rect(colour = "black", fill=NA, size=2)) + xlab("VAF")
  rm(vafXY)
  }  else{
  p8 = ggplot() + ggtitle("No Het two variants (x/y) found ") + theme(panel.border = element_rect(colour = "black", fill=NA, size=2))

  }
  # png(file=paste0(figs,"1_mqc.png"), width = 1000, height = 600)
  # (p4 | p1 | p2 | p3) / (p5 | p6 | p7 | p8)
  # dev.off()
}

plot_list = list()
allBases = data.frame(ALT= c("A", "C", "G", "T"))
for (base in c("A", "C", "G", "T")) {
  p = Z %>%   filter(TYPE == "snp" & REF==base) %>% select(GT,REF, ALT) %>% group_by(ALT) %>% summarise(count=n()) %>% dplyr::right_join(allBases) %>% replace_na (list(count=0)) %>% ggplot(., aes(x=ALT, y=count, label=count, color = ALT, fill = ALT)) + geom_bar(stat = "identity")+ geom_text(size = 3, position = position_stack(vjust = 0.5),col="white") + ggtitle(paste0("from ",base, " To ALT" ) ) + theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(colour = "black", fill=NA, size=2))
  if (base != "T") {
      p = p  + theme(legend.position = "none")
  } #enlever la légende 

 plot_list[[base]] = p
}



purines    = c("A", "G")
pyrimidines= c("C", "T")

Ti = Z %>%   filter(TYPE == "snp" & ( (REF %in% purines & ALT %in% purines) | (REF %in% pyrimidines & ALT %in% pyrimidines) ) ) %>% summarise(count=n())
Tv = Z %>%   filter(TYPE == "snp" & ( (REF %in% purines & ALT %in% pyrimidines) | (REF %in% pyrimidines & ALT %in% purines) ) ) %>% summarise(count=n())
ratio = Ti/Tv
p13=rbind(Transition=Ti,Transversion=Tv) %>% ggplot(., aes(x=rownames(.), y=count, label=count)) + geom_bar(stat = "identity")+ geom_text(size = 3, position = position_stack(vjust = 0.5),col="white") + ggtitle(paste("Bialllelic Ti/Tv ratio: ",ratio) ) + theme(plot.title = element_text(hjust = 0.5),panel.border = element_rect(colour = "black", fill=NA, size=2)) + labs(x="Biallelic change")


#indel ldistrib
bindel = Z %>%   filter(TYPE == "bdel" | TYPE == "bins") %>% mutate(size = nchar(as.character(ALT)) - nchar(as.character(REF) )) %>% select (TYPE, size)

p14 = ggplot(bindel, aes(size, fill = TYPE)) + geom_histogram(alpha = 0.5, binwidth=1) + labs(color="TYPE") + scale_fill_discrete(name="Indel types", breaks=c("bins", "bdel"), labels=c("Biallelic insertion", "Biallelic deletion")) + ggtitle("Bialllelic Indel size distribution " ) + theme(plot.title = element_text(hjust = 0.5),panel.border = element_rect(colour = "black", fill=NA, size=2))

p15 = ggplot(bindel, aes(size, fill = TYPE)) + geom_histogram(alpha = 0.5, binwidth=1) + labs(color="TYPE") + scale_fill_discrete(name="Indel types", breaks=c("bins", "bdel"), labels=c("Biallelic insertion", "Biallelic deletion")) + ggtitle("Bialllelic Indel size logscaled distribution " ) + theme(plot.title = element_text(hjust = 0.5),panel.border = element_rect(colour = "black", fill=NA, size=2)) + scale_y_log10()

# png(file=paste0(figs,"2_mqc.png"), width = 1000, height = 600)
# (plot_list[["A"]] | plot_list[["C"]] | plot_list[["G"]] | plot_list[["T"]])/(p13|p14|p15)
# dev.off()

rm(Z); rm(df);  gc()
(p4 | p1 | p2 | p3) / (p5 | p6 | p7 | p8) /
(plot_list[["A"]] | plot_list[["C"]] | plot_list[["G"]] | plot_list[["T"]])/(p13|p14|p15)


}
