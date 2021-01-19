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
    for (i = 10; i<= NF; i++){samplesNames[i-9] = $i}
}
/^[^#]/{
    nb++
    if (nb == 1) #we suppose that FORMAT does not change from snp pos
    {
        print "GT\tDP\tAD\tQUAL\tGQ\tTYPE\tREF\tALT"
        #DP AD GQ and GT cols in FORMAT
        
        DP=0
        AD=0
        GQ=0
        GT=0

        ll = split($9,a,":")
        for (i=1; i <= ll; i++)
        {
            if (a[i] == "DP") DP=i;
            if (a[i] == "AD") AD=i;
            if (a[i] == "GQ") GQ=i;
            if (a[i] == "GT") GT=i;
        }
    }

    for (i = 10; i<= NF; i++)
    {
    split($i, FORMAT, ":")
    x = split($5,a,",");

    if ($4 == $5 || $5 == ".") {refCall++}
    else {
         if (length($5)==1) {
            if(FORMAT[GT] == "0/0" || FORMAT[GT] == "./.")    {TYPE="ref"}
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

	# Missing values
	if (DP == 0) { DPVAL="NA" }
        else { DPVAL=FORMAT[DP] }

	if (AD == 0) { ADVAL="NA" }
        else { ADVAL=FORMAT[AD] }

        if (GQ == 0) { GQVAL="NA" }
        else { GQVAL=FORMAT[GQ] }

	if (GT == 0) { GTVAL="NA" }
        else { GTVAL=FORMAT[GT] }
    }# for each sample

    print GTVAL,DPVAL, ADVAL,$6,GQVAL , TYPE , $4, $5
    }
    END    {
      #  print bidel, biins, mnp, multindelcomplex, refCall, snp, snpmulti
    }  '       > $outfile


if [ "$basefigres" != '' ] 
then
  Rscript ${scriptDir}/plotVCF.R $outfile $basefigres 
fi
