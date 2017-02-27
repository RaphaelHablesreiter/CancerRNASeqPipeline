# CancerRNASeq - A pipeline for transcriptomics datasets

**CancerRNAseq** is implemented in Perl and uses R to generate various figures to visulaize the results.
**CancerRNAseq** consists of three modules: (i) **preprocessing**, (ii) **mapping** and (iii) **analysis and visualization**.

The proprocessing consists of an adapter trimming, a quality trimming and a sequencing error correction step with three successively executed tools (CutAdapt, Trimmomatic and rCorrector). In the mapping step one of three alignment tools (TopHat2, HISAT2 or STAR) can be chosen to map the preprocessed reads to the reference genome. The calculation of differentially expressed transcripts with the Cufflinks suite (Cufflinks, Cuffmerge and Cuffdiff) as well as to generate graphics of the resulting gene lists is performed in the analysis and visualization step. 

**CancerRNAseq** is a command line tool that can be executed by calling the `CancerRNASeq.pl`-file with the ability to add parameters directly to the execution command (*e.g.,* `$ perl CancerRNASeq.pl -config /path/to/file/userconfig.ini -silent -overwrite -threads 20`).

Paramaters that can be added to the execution command of **CancerRNAseq** are:

* `config <User Configuration File>`
* `help`
* `overwrite (default = FALSE)`
* `silent (default = FALSE)`
* `threads <Number of threads> (default = 1)`

The only necessary parameter to execute the pipeline successfully is the path to the user configuration file in INI-format (*e.g.,* `$ perl CancerRNASeq.pl -config /path/to/file/userconfig.ini`). This file contains information such as input file paths and parameters (*e.g.,* user name, project name, alignment tool).
