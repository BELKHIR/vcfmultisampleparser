#Generated application by SAGE MBB plateform

#@author khalid.belkhir@umontpellier.fr
#vcfmultisampleparser is GPLv3 software, authored and maintained by Khalid Belkhir
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

source ("Draw_fonction_MultiSamples.R")

# Define UI
UI <- fluidPage(

  # App title ----
  titlePanel("ShinyVCFmultiSampleparser"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
        shinyFilesButton("servervcffile" ,label="Select a VCF in the server", title="", multiple=FALSE),
        #numericInput("sample", h3("Draw details for sample :"), value = 1) ,
        selectInput("sample",  h4("Draw details for sample :"), choices =character(0)),
        width = 2 
    ),
   
    # Main panel for displaying outputs ----
    mainPanel(
      width = 10,
      
      tabsetPanel(id = "inTabset", type = "tabs",
                  tabPanel(value = "Multi", "Multisample Summary", plotOutput("Plots2", height="900px") %>% withSpinner(color="#0dc5c1")),
                  tabPanel(value = "Mono", "Sample Details", plotOutput("Plots1", height="900px") %>% withSpinner(color="#0dc5c1")
                    #, DT::dataTableOutput("IntermediateState") %>% withSpinner(color="#0dc5c1")
                          )
                 )
      )
    )
)

summarizeVCF <- function(vcfFile)
{
  ficout = tempfile()  
  x <- system(paste0('bash ./parseMultiSamplesVCF.sh ', vcfFile,' ',ficout), intern=TRUE)
  return( ficout)
}    

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
    colnames(df)[1:nsamples] = str_sub(colnames(df)[1:nsamples], 1, -4) # remove trailing .DP in sample names
    if (all(is.na(df))) depth=NULL else    depth = boxplot(df[1:min(1000000, nsites),], plot=FALSE) # limit to the first 1 million snp
    rm(df)
    gc()
  
    # Read GT values
    GTCols = seq(1,startLocusInfo - 4, 5)
    coltypes = rep('_', ncols)
    coltypes[GTCols] = 'i'
    df <- readr::read_tsv(ficstats, col_types = paste(coltypes, collapse=""))

    colnames(df)=str_sub(colnames(df), 1, -4)
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

    #Frequency of non-reference allele for bi-allelic sites
    # nb_all_sites = apply(df, 1, function(x) { all = which(unique(unlist(strsplit(x,"[/|]"))) != "."); return(length(all))} )
    # df =  df[nb_all_sites==2,] #only bi-allelic
    # freq = freq[nb_all_sites==2] #this is from maf 
    # rm(nb_all_sites)
    # gc()

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

SERVER <- function( input, output, session) {
            data <- NULL
            vcf.fn <<- ""
            VCFsummary=NULL
            shinyFileChoose(input, "servervcffile", root=c(home="/home/",temp="/tmp/"),filetypes=c('vcf', 'gz'), session = session)

            VCFsummary <- reactive({ 
              updateTabsetPanel(session, "inTabset", selected = "Multi"  )

              if (is.null(input$vcf_file$datapath)){
                fics = parseFilePaths(c(home="/home/",temp="/tmp/"),input$servervcffile)
                if (nrow(fics)>0) {
                  vcf.fn <<- fics$datapath[1]
                  ficstats = summarizeVCF(vcf.fn)
                  return(ficstats)
                }
                else return(NULL)
              }   

            } )

            output$Plots1 <- renderPlot({
                inFile <- VCFsummary()
                if (is.null(inFile)) return()

                #determine the number of samples from the header
                # df = read.table(inFile, nrows=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
                # nsamples = (dim(df)[2] - 6)/5
                # print(input$sample)
                # startSample = (as.integer(input$sample) -1) * 5 + 1 
                # startLocusInfo = nsamples * 5 + 1
                # needCols = c(startSample:(startSample + 4), startLocusInfo:(startLocusInfo +2) )
                
                # coltypes = rep('_', dim(df)[2])
                # coltypes[needCols] = c('c','n','c','n','c','c','c','c')
                # #data <- vroom::vroom(inFile, col_types = paste(coltypes, collapse=""),  delim="\t") 
                # data <- readr::read_tsv(inFile, col_types = paste(coltypes, collapse=""))
                data = NULL
                p = Draw(inFile, as.integer(input$sample), Z=data)
                if (! is.null(data)) rm(data); gc()  
                p
            })
            
           
            output$Plots2 <- renderPlot({
                if (! is.null(VCFsummary()) )
                {
                    vcftools_summary <- generate_stats(VCFsummary() )
                    gc()
                    samples = 1:length(vcftools_summary$samples)
                    names(samples) = vcftools_summary$samples
                    updateSelectInput(session, "sample", label = NULL, choices =  samples ,selected = 1)

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

                    par(mfrow=c(3,2))
                    bxp(depth, outline=FALSE, ylab=titre,  boxfill=2:8, las=3 )
                    #hist(maf, breaks=20, xlim=c(min(maf, na.rm=T),0.5), main="Histogram of minor allele frequency across all SNP", xlab="maf frequency", col= rgb(1,0,0,1/4)) 
                    plot(maf, main="Histogram of minor allele frequency across all SNP", xlab="maf frequency", col= rgb(1,0,0,1/4))
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
