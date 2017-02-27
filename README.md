# CancerRNASeq - A pipeline for transcriptomics datasets

CancerRNAseq is implemented in Perl and uses R to generate various figures to visulaize the results.
CancerRNAseq consists of three modules: (i) **preprocessing**, (ii) **mapping** and (iii) **analysis and visualization**.
The proprocessing consists of an adapter trimming, a quality trimming and a sequencing error correction step with three successively executed tools (CutAdapt, Trimmomatic and rCorrector).  
In the mapping step one of three alignment tools (TopHat2, HISAT2 or STAR) can be chosen to map the preprocessed reads to the reference genome. The calculation of differentially expressed transcripts with the Cufflinks suite (Cufflinks, Cuffmerge and Cuffdiff) as well as to generate graphics of the resulting gene lists is performed in the analysis and visualization step. 
