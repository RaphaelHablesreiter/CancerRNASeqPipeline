# CancerRNASeq

##What is **CancerRNAseq**?
**CancerRNAseq** is a pipeline for transcriptomics datasets.
**CancerRNAseq** is implemented in Perl and uses R to generate various figures to visulaize the results and consists of three modules: (i) **preprocessing**, (ii) **mapping** and (iii) **analysis and visualization**.

The proprocessing consists of an adapter trimming, a quality trimming and a sequencing error correction step with three successively executed tools ([CutAdapt](http://cutadapt.readthedocs.io/en/stable/index.html), [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) and [rCorrector](https://github.com/mourisl/Rcorrector)). In the mapping step one of three alignment tools ([TopHat2](http://www.ccb.jhu.edu/software/tophat/index.shtml), [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) or [STAR](https://github.com/alexdobin/STAR)) can be chosen to map the preprocessed reads to the reference genome. The calculation of differentially expressed transcripts with the Cufflinks suite ([Cufflinks, Cuffmerge and Cuffdiff](http://cole-trapnell-lab.github.io/cufflinks/)) as well as to generate graphics of the resulting gene lists is performed in the analysis and visualization step. 

##Install
1. Clone the [GitHub repo](https://github.com/RaphaelHablesreiter/CancerRNASeqPipeline), *e.g.* `git clone https://github.com/RaphaelHablesreiter/CancerRNASeqPipeline`
2. To run **CancerRNAseq** call the `CancerRNASeq.pl`-file, *e.g.* `$ perl CancerRNASeq.pl -config /path/to/file/userconfig.ini`)

##Usage

**CancerRNAseq** is a command line tool that can be executed by calling the `CancerRNASeq.pl`-file with the ability to add parameters directly to the execution command (*e.g.,* `$ perl CancerRNASeq.pl -config /path/to/file/userconfig.ini -silent -overwrite -threads 20`).

Paramaters that can be added to the execution command of **CancerRNAseq** are:

* `config <User Configuration File>`
* `help`
* `overwrite (default = FALSE)`
* `silent (default = FALSE)`
* `threads <Number of threads> (default = 1)`

The only necessary parameter to execute the pipeline successfully is the path to the user configuration file in INI-format (*e.g.,* `$ perl CancerRNASeq.pl -config /path/to/file/userconfig.ini`). This file contains information such as input file paths and user parameters (example user configuration file is attached).

##Output

**CancerRNAseq** creates a folder with the user name, containing subfolders with all the intermediate files and results for the pipeline run.

##Terms of use

Copyright (C) 2017  Raphael Hablesreiter

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by  the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

