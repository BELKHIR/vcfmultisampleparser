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
    print Header"QUAL\tREF\tALT\tPi\tMAF\tMiss"    
}
/^[^#]/{
    nb++
    if (nb == 1) #we suppose that FORMAT does not change from snp to snp !
    {
                 
        #check DP AD GQ and GT cols in FORMAT field
        DP=0
        AD=0
        GQ=0
        GT=0
        #freebayes output The AO and RO fields reflect the number of read observations for each alternate and the reference allele respectively.
        #this is a replacment of AD. AO+RO <= DP in general. This is because DP counts all reads covering the variant site, but AO and RO only include reads deemed good enough for an allele call.
        AO=0
        RO=0
        ll = split($9,a,":") #FORMAT field
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

    x = split($5,a,",");
    split("", alcounts) #clear the array
    missing_site = 0
    totGamets = 0

    for (i = 10; i<= NF; i++)
    {
    
    split($i, FORMAT, ":")
   
    if ($4 == $5 ) {refCall++}
    else{

    if ( index(FORMAT[GT],".") != 0) {TYPE="noCall"; missing_site += 1}
    else {
         if (length($5)==1) {
            if(FORMAT[GT] == "0/0" || FORMAT[GT] == "0|0" )    {TYPE="ref"} 
            else {
                if (length($4) == 1) {TYPE="snp"}
                else TYPE="bdel"; # ATTTAC   -> A
            }
         }
         else {

             if (x >= 2 && length($5) == x+(x-1) )
             {
                if(FORMAT[GT] == "0/0" || FORMAT[GT] == "0|0") { TYPE="ref"}
                else {
                    if (length($4) == 1) TYPE="snpmulti";
                    else TYPE="multiIcmpl";
                }
             }

             if (x >= 2 && length($5) != x+(x-1)) {  TYPE="multiIcmpl";  }

             if (x == 1) {

               if ( length($4) == length($5) ) if(FORMAT[GT] == "0/0" || FORMAT[GT] == "0|0")   {TYPE="ref";} else  {TYPE="mnp";}
               else {
                    if(FORMAT[GT] == "0/0" || FORMAT[GT] == "0|0")  {TYPE="ref"; }
                    else {
                    if (length($4) < length($5))  {TYPE="bins";} else  {TYPE="bdel";} }
                   }
             }
         }
        }
    }    


	
    if ( index(FORMAT[GT],".") == 0) {
        #calc all. freq
        ploidy = split(FORMAT[GT],allels, "[/|]")
        if (ploidy == 2)
        {
            totGamets = totGamets + 2
            if (allels[2] <  allels[1]) {al1 = allels[2]; al2= allels[1]} else {al1 = allels[1]; al2= allels[2]}
            alcounts[al1]++
            alcounts[al2]++  
            #code genotype as integer : 0/0 -> 00, 0/1 or 1/0 -> 01, 1/1 -> 11, 1/2 -> 12 
            # no values like 10 or 20 ... as alleles are sorted
            # the phase is lost !
            #if readed as integers : 0 REF REF, 1-9 REF ALT, 11 22 homz. ALT , 12 13 23 hetero ALT 
            FORMAT[GT] = al1 al2           
        }
        
    } else FORMAT[GT] = "NA"


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
    else { GTVAL=FORMAT[GT]     }
    
    sortie = sortie GTVAL"\t"DPVAL"\t"ADVAL"\t"GQVAL"\t"TYPE"\t" 

    }# for each sample

    # MAF
    total_alleles=0
    if (totGamets > 0)
    {
        MAF=999999
        for (al in alcounts) {
            if (alcounts[al] <MAF) MAF = alcounts[al]
            total_alleles += alcounts[al]
        }
        #MAF= MAF/totGamets
        MAF= MAF/total_alleles
        if (MAF == 1) MAF=0
    }
    else MAF="NA"

    # Pi
    mismatches = 0
    for (al in alcounts) {
        other_alleles_count = (total_alleles - alcounts[al]);
		mismatches += (alcounts[al] * other_alleles_count);
    }
    pairs = (total_alleles * (total_alleles - 1));
	pi = mismatches/pairs;

    print sortie $6"\t"$4"\t"$5"\t"pi"\t"MAF"\t"missing_site
     
    }
    END    {
      #  print bidel, biins, mnp, multindelcomplex, refCall, snp, snpmulti
    }  '     > $outfile


if [ "$basefigres" != '' ] 
then
  #Rscript ${scriptDir}/plotMultiSamplesVCF.R $outfile $basefigres 1 #plot for the first sample

  Rscript ${scriptDir}/vcftools_plots.R $vcf $basefigres  #multisample vcftools  plot
fi



 
