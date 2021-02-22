scriptDir=`dirname "$0"`

vcf=$1 #a vcf file optionnaly gziped
outfile=$2 #intermediate outfile
basefigres=$3 #base name for figures (will not call plotVCF.R if if empty ) 

gzip -t $vcf 2>/dev/null

[[ $? -eq 0 ]]&&   cmd="zcat" || cmd="cat"

$cmd $vcf | \
awk 'BEGIN{nb=0; OFS="\t"}
/^#CHROM/{ 
    nbsamples = NF - 9;
    Header = ""
    for (i = 10; i<= NF; i++){
        samplesNames[i-9] = $i
        Header = Header $i".GT\t"$i".DP\t"$i".AD\t"$i".GQ\t"$i".TYPE\t"
        }
    print Header"QUAL\tREF\tALT"    
}
/^[^#]/{
    nb++
    if (nb == 1) #we suppose that FORMAT does not change from snp pos
    {
         
        #DP AD GQ and GT cols in FORMAT
        
        DP=0
        AD=0
        GQ=0
        GT=0
        #freebayes output The AO and RO fields reflect the number of read observations for each alternate and the reference allele respectively.
        #this is a replacment of AD. AO+RO <= DP in general. This is because DP counts all reads covering the variant site, but AO and RO only include reads deemed good enough for an allele call.
        AO=0
        RO=0
        ll = split($9,a,":")
        for (i=1; i <= ll; i++)
        {
            if (a[i] == "DP") DP=i;
            if (a[i] == "AD") AD=i;
            if (a[i] == "GQ") GQ=i;
            if (a[i] == "GT") GT=i;
            
            if (a[i] == "AO") AO=i;
            if (a[i] == "RO") RO=i;
        }
    }

    sortie = ""

    for (i = 10; i<= NF; i++)
    {

    split($i, FORMAT, ":")
    x = split($5,a,",");

    if ($4 == $5 || $5 == ".") {refCall++} # we consider that a "." in the ALT field is a ref !
    else{

    if ( FORMAT[GT] == "./.") {TYPE="noCall"}
    else {
         if (length($5)==1) {
            if(FORMAT[GT] == "0/0" )    {TYPE="ref"} # we consider that a "." in GT field is a like ref !
            else {
                if (length($4) == 1) {TYPE="snp"}
                else TYPE="bdel"; # ATTTAC   -> A
            }
         }
         else {

             if (x >= 2 && length($5) == x+(x-1) )
             {
                if(FORMAT[GT] == "0/0" || FORMAT[GT] == "./.") { TYPE="ref"}
                else {
                    if (length($4) == 1) TYPE="snpmulti";
                    else TYPE="multiIcmpl";
                }
             }

             if (x >= 2 && length($5) != x+(x-1)) {  TYPE="multiIcmpl";  }

             if (x == 1) {

               if ( length($4) == length($5) ) if(FORMAT[GT] == "0/0" || FORMAT[GT] == "./.")   {TYPE="ref";} else  {TYPE="mnp";}
               else {
                    if(FORMAT[GT] == "0/0" || FORMAT[GT] == "./.")  {TYPE="ref"; }
                    else {
                    if (length($4) < length($5))  {TYPE="bins";} else  {TYPE="bdel";} }
                   }
             }
         }
        }
    }    

	# Missing values
	if (DP == 0 || (FORMAT[DP]=="") ) { 
        DPVAL="NA" 
        #if we have AD or (AO and RO) we can sum the alleles count to get a pseudo DP
        if (AD != 0 && (FORMAT[AD] !="")) {split(FORMAT[AD],counts,","); for(al in counts) DPVAL = DPVAL+counts[al]}
        if (AO!=0 || RO!=0) {DPVAL=0;ADVAL=FORMAT[RO]","FORMAT[AO]; split(ADVAL,counts,","); for(al in counts) DPVAL = DPVAL+counts[al] }
    }
    else {  DPVAL=FORMAT[DP] }

	if (AD == 0 || (FORMAT[AD] =="")) { 
        #test if we have AO and RO
        if(AO!=0 || RO!=0) {ADVAL=FORMAT[RO]","FORMAT[AO]}
        else {ADVAL="NA" }
    }
    else { ADVAL=FORMAT[AD] }

    if (GQ == 0 || (FORMAT[GQ] == "")) { GQVAL="NA" }
        else { GQVAL=FORMAT[GQ] }

	if (GT == 0 || (FORMAT[GT]=="")) { GTVAL="NA" }
        else { GTVAL=FORMAT[GT] }
    
    sortie = sortie GTVAL"\t"DPVAL"\t"ADVAL"\t"GQVAL"\t"TYPE"\t" 

    }# for each sample
    print sortie $6"\t"$4"\t"$5
     
    }
    END    {
      #  print bidel, biins, mnp, multindelcomplex, refCall, snp, snpmulti
    }  '       > $outfile


if [ "$basefigres" != '' ] 
then
  #Rscript ${scriptDir}/plotMultiSamplesVCF.R $outfile $basefigres 1 #plot for the first sample

  Rscript ${scriptDir}/vcftools_plots.R $vcf $basefigres  #multisample vcftools  plot
fi



 
