
This application summarize VCF files visualy. 

First, commands are used to extract multisamples infos and generate their plots.
![Multi samples plot](multisamples.png)

If desired one can get more details about selected sample in a separate tab.

Detailed statistics are mostly inspired by deepvariant VCF stats report : https://github.com/google/deepvariant/blob/r0.10/docs/deepvariant-vcf-stats-report.md.

![Detailed sample plot](SampleDetails.png)

First the VCF file is parsed to extract various informations (GT, DP, GQ, ...) that are then ploted.

Freebayes users have to use --genotype-qualities option to get the GQ values.

Bcftools mpileup must be run  at least with the option -a AD,DP

This can be achieved either via a command line with the 2 scripts or via a provided shiny application.

vcfmultisampleparser is GPLv3 software, authored and maintained by Khalid Belkhir