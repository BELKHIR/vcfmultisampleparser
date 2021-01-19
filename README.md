
This application summarize VCF files visualy. 

First vcftools commands are used to extract multisamples infos and generate their plots.

If desired one can get more details about selected sample in a separate tab.

Statistc mostly inspired by deepvariant VCF stats report : https://github.com/google/deepvariant/blob/r0.10/docs/deepvariant-vcf-stats-report.md are generated.

First the VCF file is parsed to extract various information (GT, DP, GQ, ...) that are then ploted.

This can be achieved either via a command line with the 2 scripts or via a shiny application.

vcfparser is GPLv3 software, authored and maintained by Khalid Belkhir