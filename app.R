#Generated application by SAGE MBB plateform

#@author khalid.belkhir@univ-montp2.fr
#vcfparser is GPLv3 software, authored and maintained by Khalid Belkhir
#Rscript -e 'shiny::runApp(".", host="127.0.0.1", port=4123)' &
    
library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)
library(ggplot2)
library(tidyverse)
library("patchwork")
library(shinycssloaders)
library(shinyFiles)

library(data.table)
options(encoding = 'UTF-8', shiny.maxRequestSize=500*1024^2) #130MB


source ("Draw_fonction_MultiSamples.R")

# Define UI for random distribution app ----
UI <- fluidPage(

  # App title ----
  titlePanel("ShinyVCFparser"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
        #fileInput("vcf_file", label = "Your VCF file :", multiple = TRUE, accept = c()),
        shinyFilesButton("servervcffile" ,label="Select a VCF in the server", title="", multiple=FALSE),
        numericInput("sample", h3("Draw details for sample :"), value = 1) ,
        width = 2 
    ),
   
    # Main panel for displaying outputs ----
    mainPanel(
      width = 10,
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("MultisampleSummary", plotOutput("Plots2", height="900px") %>% withSpinner(color="#0dc5c1")),
                  tabPanel("SampleDetails", plotOutput("Plots1", height="700px") %>% withSpinner(color="#0dc5c1"),
                           DT::dataTableOutput("IntermediateState") %>% withSpinner(color="#0dc5c1")
                          )
                 )
      )
    )
)

summarizeVCF <- function(vcfFile)
{
  ficout = tempfile()  
  x <- system(paste0('bash ./parseMultiSamplesVCF.sh ', vcfFile,'  ',ficout), intern=TRUE)
  return( ficout)
}    

generate_stats <- function(fic, ficstats){
    vcf.fn = fic
    
    if (is.null(ficstats)) return()

    #determine the number of samples from the header
    df = read.table(ficstats, nrows=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    ncols  = dim(df)[2]
    
    nsamples = (dim(df)[2] - 6)/5
    startSample =  1 
    startLocusInfo = nsamples * 5 + 1
    DPCols = seq(2,startLocusInfo -4, 5)
    neededCols = c(DPCols, ncols-2, ncols-1, ncols)
    coltypes = rep('_', dim(df)[2])
    coltypes[neededCols] = 'n'
    df <- vroom::vroom(ficstats, col_types = paste(coltypes, collapse=""),  delim="\t")
    nsites = dim(df)[1]
    colnames(df)[1:nsamples] = str_sub(colnames(df)[1:nsamples], 1, -4) # remove trailing .DP in sample names
    colmissing = dim(df)[2]
    colmaf = colmissing -1
    colpi = colmaf -1

    if (all(is.na(df[,c(-colpi, -colmissing, -colmaf)]))) depth=NULL else    depth = boxplot(df[,c(-colpi, -colmissing, -colmaf)], plot=F)

    SitePi = hist(df$Pi)
    maf = as.numeric(df$MAF)
    SiteMissingness = hist(df$Miss)
   
    GTCols = seq(1,startLocusInfo - 4, 5)
    coltypes = rep('_', ncols)
    coltypes[GTCols] = 'c'
    df <- vroom::vroom(ficstats, col_types = paste(coltypes, collapse=""),  delim="\t")
    colnames(df)=str_sub(colnames(df), 1, -4)

    N_non_missing_sites = apply(df, 2, function(x) sum(x =="./." | x ==".|.") )
    Missingness = N_non_missing_sites / nsites

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

    #Frequency of non-reference allele for bi-allelic sites
    nb_all_sites = apply(df, 1, function(x) { all = which(unique(unlist(strsplit(x,"[/|]"))) != "."); return(length(all))} )
    df =  df[nb_all_sites==2,] #only bi-allelic
    freq = apply(df, 1, function(x) {
        f=table(x); 
        totInd = sum(f) #sum(f[c("0/0","0/1","1/0","1/1", "0|0","0|1","1|0","1|1")], na.rm=T);  
        althom = sum(f[c("1/1", "1|1")], na.rm=T)
        althet = sum(f[c("0/1","1/0","0|1","1|0")] , na.rm=T)
        return( ((2* althom) + althet)/(2*totInd) )
    })

    N_non_missing_site = apply(df, 1, function(x) sum(x !="./." & x !=".|.") )
    site_expected_hom = 1.0 - (2.0 * freq * (1.0 - freq) * (N_non_missing_site / (N_non_missing_site - 1.0)));
  
    N_obs_hom      = apply(df, 2, function(x) sum(x == "0/0" | x=="0|0" | x == "1/1" | x == "1|1"))
    N_expected_hom = apply(df, 2, function(x, expected_hom) { nonMissing = (x !="./." & x !=".|."); return(sum(expected_hom[nonMissing], na.rm=T)) }, expected_hom = site_expected_hom )
    N_non_missing_ind = apply(df, 2, function(x) sum(x !="./." & x !=".|.") )

    F = (N_obs_hom - N_expected_hom) / (N_non_missing_ind - N_expected_hom);


    return(list(samples = nsamples,  depth=depth,missingness=Missingness, maf=as.numeric(maf), sitemissingness = SiteMissingness, F=as.numeric(F), SitePi=SitePi ) )
}

SERVER <- function( input, output, session) {
            data <- NULL
            vcf.fn <<- ""
            VCFsummary=NULL
            shinyFileChoose(input, "servervcffile", root=c(Data="/home/",temp="/tmp/"),filetypes=c('vcf', 'gz'), session = session)
 

            VCFsummary <- reactive({ 
              if (is.null(input$vcf_file$datapath)){
                fics = parseFilePaths(c(Data="/home/",temp="/tmp/"),input$servervcffile)
                if (nrow(fics)>0) {
                  vcf.fn <<- fics$datapath[1]
                  ficstats = summarizeVCF(vcf.fn)
                  return(ficstats)
                }
                else return(NULL)
              }   

            } )

            output$IntermediateState= DT::renderDataTable({
                inFile <- VCFsummary()
                
                if (is.null(inFile)) return()
                #determine the number of samples from the header
                 df = read.table(inFile, nrows=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
                 nsamples = (dim(df)[2] - 6)/5
                 startSample = (input$sample -1) * 5 + 1 
                 startLocusInfo = nsamples * 5 + 1
                 needCols = c(startSample:(startSample + 4), startLocusInfo:(startLocusInfo +2) )
                
                coltypes = rep('_', dim(df)[2])
                coltypes[needCols] = c('c','n','c','n','c','c','c','c')
                data <- vroom::vroom(inFile, col_types = paste(coltypes, collapse=""),  delim="\t")
                output$Plots1 <- renderPlot({
                  inFile <- VCFsummary()
                  if (is.null(inFile)) return() 
                  Draw(inFile, isolate(input$sample), Z=data)
                })
                data 
            })


            
            output$Plots2 <- renderPlot({
                if (! is.null(VCFsummary()) )#(! is.null(input$vcf_file)) 
                {
                    vcftools_summary <- generate_stats(vcf.fn, VCFsummary() )#(input$vcf_file$datapath) 
                    
                    depth=vcftools_summary$depth
                    missingness = vcftools_summary$missingness
                    maf = vcftools_summary$maf
                    sitemissingness = vcftools_summary$sitemissingness
                    SitePi = vcftools_summary$SitePi
                    F = vcftools_summary$F
                    samples = vcftools_summary$samples
                    

                    minmafF = 0
                    maxmissingF = 0
                    if(is.null(depth)) { 
                      depth=boxplot(c(0,0,0,0),plot=F)
                      titre="No sample DP available! "
                    } else titre=paste0("Per sample SNP Depth (DP) distrib.")

                    par(mfrow=c(3,2))
                    bxp(depth, outline=FALSE, ylab=titre,  boxfill=2:8, las=3 )
                    hist(maf, breaks=20, xlim=c(min(maf, na.rm=T),0.5), main="Histogram of minor allele frequency across all SNP", xlab="maf frequency", col= rgb(1,0,0,1/4)) 
                    barplot(missingness, ylab=paste0("Per-sample missingness (%)"),  col=2:8, las=3 )
                    plot(sitemissingness, main="Per-site missingness", xlab="% missing",col=rgb(0,0,1,1/4))
                    barplot(as.numeric(F), ylab=paste0("Per-sample F (inbreeding Coef using a method of moments) "),  col=2:8, las=3 )
                    plot(SitePi, main="Per-site Pi", xlab="Pi",col=rgb(0,0,1,1/4))
                    mtext(vcf.fn, side = 3, line = -2, outer = TRUE)
                    par()
                }
            })
}

shinyApp(ui = UI, server = SERVER)
