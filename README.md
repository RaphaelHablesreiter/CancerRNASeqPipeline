# CancerRNASeq - A pipeline for transcriptomics datasets

CancerRNAseq is implemented in Perl and uses R to generate various figures to visulaize the results.
CancerRNAseq consists of three modules: (i) **preprocessing**, (ii) **mapping** and (iii) **analysis and visualization**.
The proprocessing consists of an adapter trimming, a quality trimming and a correction of sequencing error with three successively executed tools (CutAdapt, Trimmomatic and rCorrector).  
