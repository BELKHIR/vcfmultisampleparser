library("patchwork")
library(ggplot2)
library(tidyverse)

args = commandArgs(trailingOnly=FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", args[grep(file.arg.name, args)])
script.dir <- dirname(script.name)

args = commandArgs(trailingOnly=TRUE)
data = args[1]
figs = args[2]
sample = as.integer(args[3])

source(paste0(script.dir,"/Draw_fonction_MultiSamples.R"))

png(file=paste0("sample_",sample,"_",figs,"_vcf_plots_mqc.png"), width = 1200, height = 1000)
Draw(data, sample)
dev.off()
