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

generate_stats <- function(fic){
    vcf.fn = fic
    
    cmd = paste0("vcftools --gzvcf ",vcf.fn," --geno-depth --stdout ")
    ttt = system(cmd, intern=T)
    dp = data.frame(matrix(unlist(strsplit(ttt,'\t')), nrow=length(ttt), byrow=T),  stringsAsFactors=F)
    colnames(dp) = dp[1,]
    dp=dp[-1,-c(1,2)]
    dp[dp == -1] <- NA
    dp = apply(dp, 2, function(x) as.numeric(x) )
    if (all(is.na(dp))) depth=NULL else    depth = boxplot(dp, plot=F)
    
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

SERVER <- function( input, output, session) {
            data <- NULL
            vcf.fn <<- ""
            VCFsummary=NULL
            shinyFileChoose(input, "servervcffile", root=c(Data="/home/",Results="/home/"),filetypes=c('vcf', 'gz'), session = session)
 

            VCFsummary <- reactive({ 
              if (is.null(input$vcf_file$datapath)){
                fics = parseFilePaths(c(Data="/home/",Results="/home/"),input$servervcffile)
                print(fics)
                if (nrow(fics)>0) {
                  vcf.fn <<- fics$datapath[1]
                  summarizeVCF(vcf.fn)
                }
                else return(NULL)
              }   
            #if (is.null(input$vcf_file$datapath)) return(NULL);  
            #summarizeVCF(input$vcf_file$datapath) 
            } )

            output$IntermediateState= DT::renderDataTable({
                inFile <- VCFsummary()
                
                if (is.null(inFile)) return()
                #determine the number of samples from the header
                 df = read.table(inFile, nrows=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
                 nsamples = (dim(df)[2] - 3)/5
                 startSample = (input$sample -1) * 5 + 1 
                 startLocusInfo = nsamples * 5 + 1
                 needCols = c(startSample:(startSample + 4), startLocusInfo:(startLocusInfo +2) )
                
                #data = read.table(inFile, head=T, colClasses = colClasses,  sep="\t", stringsAsFactors=F)
                #data = fread(inFile, header=T, sep="\t", select=needCols, stringsAsFactors=F)
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
                    vcftools_summary <- generate_stats(vcf.fn )#(input$vcf_file$datapath) 
                    
                    depth=vcftools_summary$depth
                    missingness = vcftools_summary$missingness
                    maf = vcftools_summary$maf
                    sitemissingness = vcftools_summary$sitemissingness
                    SitePi = vcftools_summary$SitePi
                    F = vcftools_summary$F
                    samples = vcftools_summary$samples
                    snps = vcftools_summary$snps

                    minmafF = 0
                    maxmissingF = 0
                    if(is.null(depth)) { 
                      depth=boxplot(c(0,0,0,0),plot=F)
                      titre="No sample DP available! "
                    } else titre=paste0("Per sample SNP Depth (DP) distrib. minMAF:", minmafF, " maxMissing:", maxmissingF)

                    par(mfrow=c(3,2))
                    bxp(depth, outline=FALSE, main=titre,  boxfill=2:8, las=3 )
                    hist(maf, breaks=20, xlim=c(min(maf, na.rm=T),0.5), main="Histogram of minor allele frequency across all SNP", xlab="maf frequency", col= rgb(1,0,0,1/4)) 
                    barplot(missingness, main=paste0("Per-sample missingness (%)"),  col=2:8, las=3 )
                    plot(sitemissingness, main="Per-site missingness", xlab="% missing",col=rgb(0,0,1,1/4))
                    barplot(as.numeric(F), main=paste0("Per-sample F (inbreeding Coef using a method of moments) "),  col=2:8, las=3 )
                    plot(SitePi, main="Per-site Pi", xlab="Pi",col=rgb(0,0,1,1/4))
                    mtext(vcf.fn, side = 3, line = -2, outer = TRUE)
                    par()
                }
            })
}

shinyApp(ui = UI, server = SERVER)
