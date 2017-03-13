#------------------------------------------------------------------------------
# CancerRNASeq - a pipeline for transcriptomics datasets
# Author: Raphael Hablesreiter
#
# ANALYZING.pm
#          
# Description: Analyzing Module
#              In the analyzing module the calculation of DE expressed genes
#              between the two datasets is performed.
#              The calculation of differential expression is performed with 
#              the Cufflinks suite (i.e., Cufflinks, Cuffmerge and Cuffdiff).
#              The output of Cuffdiff is processed with R and the package
#              CummeRbund to generate graphics and tables of the results.
#
#-------------------------------------------------------------------------------
#

package ANALYZING;

use 5.010;
use strict;
use warnings;

use Backticks;

use PIPELINETOOLS;

#-------------------------------------------------------------------------------
# Main: Starts one of the 3 parts of the analyzing module.
#       Input: 1,2 or 3 (is provided by the Driver module)
#       Options: 1 -> Cufflinks
#                2 -> Cuffmerge
#                3 -> Cuffdiff
#-------------------------------------------------------------------------------
sub Main{
    
    my $i = $_[0];
    
    given ($i)
    {
        when (1) {ANALYZING::CufflinksRun();}
        when (2) {ANALYZING::CuffmergeRun();}
        when (3) {ANALYZING::CuffdiffRun();}
        when (4) {ANALYZING::CummerbundRun();}
        default {die "ERROR:Cant find out where to start the analyziation!\n";}
    }
    
    return;
}

#-------------------------------------------------------------------------------
# CutadaptRun: The function creates the command for a proper call of
#              Cufflinks and then starts Cufflinks.
#-------------------------------------------------------------------------------
sub CufflinksRun{
    PIPELINETOOLS::WriteLogChapter("Cufflinks");
    PIPELINETOOLS::PrintTask("Started Cufflinks");
    
    my @aligned_reads;
    
    if (($::global_cfg -> param("required.MappingType")) eq "Tophat")
    {
        @aligned_reads = PIPELINETOOLS::ReadDirectory(
            ($::global_cfg -> param("results.alignment")),"accepted_hits.bam");
    }
    elsif (($::global_cfg -> param("required.MappingType")) eq "STAR")
    {
        @aligned_reads = PIPELINETOOLS::ReadDirectory(
            ($::global_cfg -> param("results.alignment")),"aligned.bam");        
    } 
    else
    {
        @aligned_reads = PIPELINETOOLS::ReadDirectory(
            ($::global_cfg -> param("results.alignment")),".bam");        
    }
    
    my $alignedreads_length = scalar(@aligned_reads);
    my @aligned_seq1 = ();
    my @aligned_seq2 = ();
    
    
    if (@aligned_reads)
    {
        my $check_succ = 0;
        while (@aligned_reads)
        {   
            my $tmp = shift(@aligned_reads);
            
            if ($tmp =~ ($::global_cfg -> param("required.sequence_1")) &&
                $tmp !~ ".bai")
            {
                push(@aligned_seq1,
                    ($::global_cfg -> param("results.alignment"))."/".$tmp);
                $check_succ = $check_succ + 1;
            }
            elsif ($tmp =~ ($::global_cfg -> param("required.sequence_2")) &&
                $tmp !~ ".bai")
            {
                push(@aligned_seq2,
                    ($::global_cfg -> param("results.alignment"))."/".$tmp);
                $check_succ = $check_succ + 1;
            }  
        }
        
        @aligned_seq1 = sort(@aligned_seq1);
        @aligned_seq2 = sort(@aligned_seq2);
    }
    else
    {
        die "ERROR: Can't find the alignment files!\n";
    }
    
    if ($alignedreads_length != 0)
    {     
        my $assemblies_txt = ($::global_cfg -> param("results.tmp"))."/".
            "assemblies.txt";
            
        open(ASSEMBLIES, ">", $assemblies_txt)
            or die "ERROR: Can't create assemblies.txt!\n";
        
        my $succ_count = 0;
        
        my $i = 1;
        my $j = 1;
        my $length_seq1 = scalar(@aligned_seq1);
        while ($i <= $alignedreads_length)
        {
            my $tmp_alignment;
            my $alignment_f;
            
            
            if (@aligned_seq1)
            {
                $tmp_alignment = shift(@aligned_seq1);
                $alignment_f =
                    ($::global_cfg -> param("required.sequence_1")).
                    "_".$j."_Cufflinks";
            }
            elsif (@aligned_seq2)
            {
                $tmp_alignment = shift(@aligned_seq2);
                $alignment_f =
                    ($::global_cfg -> param("required.sequence_2")).
                    "_".$j."_Cufflinks";                
            }
            
            
            if ($::global_cfg -> param("optional.delete_mt") &&
                ($::global_cfg -> param("optional.delete_mt")) eq "yes")
            {
                PIPELINETOOLS::PrintTask("Started Deleting MT");
                my $index_cmd = "samtools index ".$tmp_alignment;
                
                PIPELINETOOLS::PrintStd($index_cmd);
                
                my $cmd = `$index_cmd`;
                if ($cmd -> success() != 1)
                {
                    PIPELINETOOLS::WriteLogOut($cmd -> merged);
                    PIPELINETOOLS::PrintStd($cmd -> merged);
                    die "ERROR: Can't create index for $tmp_alignment !\n";
                }                
                    
                my $deletemt_cmd = "samtools idxstats ".
                    $tmp_alignment.
                    " | cut -f 1 | grep -v MT | xargs samtools view -b ".
                    $tmp_alignment.
                    " > ".
                    ($::global_cfg -> param("results.alignment"))."/".
                    "out_wo_MT.bam";
                    
                PIPELINETOOLS::PrintStd($deletemt_cmd);
                
                $cmd = `$deletemt_cmd`;
                if ($cmd -> success() != 1)
                {
                    PIPELINETOOLS::WriteLogOut($cmd -> merged);
                    PIPELINETOOLS::PrintStd($cmd -> merged);
                    die "ERROR: Can't create .bam without ".
                        "ChrMT for $tmp_alignment!\n";
                }

                PIPELINETOOLS::WriteLogOut($cmd -> stdout);
                PIPELINETOOLS::PrintStd($cmd -> stdout);
                
                unlink ($tmp_alignment.".bai");
                    
                unlink $tmp_alignment;
                
                rename (($::global_cfg -> param("results.alignment"))."/".
                    "out_wo_MT.bam",
                    $tmp_alignment);               
            }
            
            my $mkdir_cmd = system("mkdir -p ".
                ($::global_cfg -> param("results.tmp"))."/".
                $alignment_f)
                and die "ERROR: mkdir\n";
                
            my $cufflinks_cmd;
            
            if ($::global_cfg -> param("required.MappingType") eq "STAR")
            {
                $cufflinks_cmd =
                    ($::global_cfg -> param("sprograms.Cufflinks")).
                    " -q -p ".($::global_cfg -> param("cmdline.threads")).
                    #" --library-type fr-firststrand".                    
                    " -o ".($::global_cfg -> param("results.tmp"))."/".
                    $alignment_f.
                    " ".$tmp_alignment;
            }
            else
            {
                $cufflinks_cmd =
                    ($::global_cfg -> param("sprograms.Cufflinks")).
                    " -q -p ".($::global_cfg -> param("cmdline.threads")).
                    " -o ".($::global_cfg -> param("results.tmp"))."/".
                    $alignment_f.
                    " ".$tmp_alignment;                
            }
            
                
            PIPELINETOOLS::PrintStd($cufflinks_cmd);
            
            my $cmd1 = `$cufflinks_cmd`;
            if ($cmd1 -> success != 1)
            {
                PIPELINETOOLS::WriteLogOut($cmd1 -> merged);
                PIPELINETOOLS::PrintStd($cmd1 -> merged);
                die "ERROR: Can't run cufflinks on $tmp_alignment!\n";
            }
            
            PIPELINETOOLS::WriteLogOut($cmd1 -> stdout);
            PIPELINETOOLS::PrintStd($cmd1 -> stdout);
            
            print ASSEMBLIES ($::global_cfg -> param("results.tmp"))."/".
                $alignment_f."/transcripts.gtf\n";     
            
            $succ_count = $succ_count + 1;
            
            $i++;
            
            if ($length_seq1 == $j)
            {
                $j = 1;
            }
            else
            {
                $j++;
            }
            
            
        }
        
        if ($succ_count != $alignedreads_length)
        {
            die "ERROR: Can't find all alignments!\n";
        }
        
    }
    else
    {
        die "ERROR: Can't find aligned read files!\n";
    }
    
    close ASSEMBLIES;
    PIPELINETOOLS::CreateSuccess("Cufflinks"); 
    PIPELINETOOLS::PrintTask("Finished Cufflinks");
    
    return;
}

#-------------------------------------------------------------------------------
# CutadaptRun: The function creates the command for a proper call of
#              Cuffmerge and then starts Cuffmerge.
#-------------------------------------------------------------------------------
sub CuffmergeRun{
    PIPELINETOOLS::WriteLogChapter("Cuffmerge");
    PIPELINETOOLS::PrintTask("Started Cuffmerge");
    
    my $assemblies_txt = ($::global_cfg -> param("results.tmp"))."/".
        "assemblies.txt";
    my $ref_genome = $::global_cfg -> param("required.GenomeFASTA");    
    my $ref_gtf = $::global_cfg -> param("required.GenomeGTF");
    
        
    if (-e $assemblies_txt and -e $ref_gtf and -e $ref_genome)
    {
        my $cuffmerge_cmd = ($::global_cfg -> param("sprograms.Cuffmerge")).
            " -g ".$ref_gtf.
            " -s ".$ref_genome.
            " -p ".($::global_cfg -> param("cmdline.threads")).
            " -o ".($::global_cfg -> param("results.tmp"))."/merged_asm".
            " ".$assemblies_txt;
        
        PIPELINETOOLS::PrintStd($cuffmerge_cmd);
        
        my $cmd = `$cuffmerge_cmd`;
        if ($cmd -> success() != 1)
        {
            PIPELINETOOLS::PrintStd($cmd -> merged);
            PIPELINETOOLS::WriteLogOut(($cmd -> merged));
            die "ERROR: Can't run cuffmerge!\n";
        }
        
        PIPELINETOOLS::PrintStd($cmd -> stdout);
        PIPELINETOOLS::WriteLogOut($cmd -> stdout);
    }
    elsif (!-e $assemblies_txt)
    {
        die "ERROR: Can't find assemblies *.txt - file";
    }
    elsif (!-e $ref_genome)
    {
        die "ERROR: Can't find reference annotations *.gtf - file";    
    }
    else
    {
        die "ERROR: Can't find reference genome *.fasta - file";
    }
    
    PIPELINETOOLS::CreateSuccess("Cuffmerge");
    PIPELINETOOLS::PrintTask("Finished Cuffmerge");
    
    return;
}

#-------------------------------------------------------------------------------
# CutadaptRun: The function creates the command for a proper call of
#              Cuffdiff and then starts Cuffdiff.
#-------------------------------------------------------------------------------
sub CuffdiffRun{
    PIPELINETOOLS::WriteLogChapter("Cuffdiff");
    PIPELINETOOLS::PrintTask("Started Cuffdiff");
    
    my $ref_genome = $::global_cfg -> param("required.GenomeFASTA");
    my $merged_gtf = ($::global_cfg -> param("results.tmp")).
        "/merged_asm/merged.gtf";
    
    my @aligned_reads;
    
    if (($::global_cfg -> param("required.MappingType")) eq "Tophat")
    {
        @aligned_reads = PIPELINETOOLS::ReadDirectory(
            ($::global_cfg -> param("results.alignment")),"accepted_hits.bam");
    }
    else
    {
        @aligned_reads = PIPELINETOOLS::ReadDirectory(
            ($::global_cfg -> param("results.alignment")),".bam");        
    }
    
    my @aligned_seq1 = ();
    my @aligned_seq2 = ();
    
    if (@aligned_reads)
    {
        my $check_succ = 0;
        while (@aligned_reads)
        {   
            my $tmp = shift(@aligned_reads);
            
            if ($tmp =~ ($::global_cfg -> param("required.sequence_1")) &&
                $tmp !~ ".bai")
            {
                push(@aligned_seq1,
                    ($::global_cfg -> param("results.alignment"))."/".$tmp);
                $check_succ = $check_succ + 1;
            }
            elsif ($tmp =~ ($::global_cfg -> param("required.sequence_2")) &&
                $tmp !~ ".bai")
            {
                push(@aligned_seq2,
                    ($::global_cfg -> param("results.alignment"))."/".$tmp);
                $check_succ = $check_succ + 1;
            }  
        }
        
        @aligned_seq1 = sort(@aligned_seq1);
        @aligned_seq2 = sort(@aligned_seq2);
    }
    else
    {
        die "ERROR: Can't find the *accepted_hits.bam-files!\n";
    }
    
  
        
    if (-e $ref_genome and
        -e $merged_gtf)
    {
        if (($::global_cfg -> param("required.analysistype")) eq
            "pooled")
        {
            my $aligned_seqs1;
            my $aligned_seqs2;
    
            while (@aligned_seq1)
            {                
                $aligned_seqs1 = $aligned_seqs1.shift(@aligned_seq1);
                
                if (@aligned_seq1)
                {
                    $aligned_seqs1 = $aligned_seqs1.",";
                }
            }
            
            while (@aligned_seq2)
            {                
                $aligned_seqs2 = $aligned_seqs2.shift(@aligned_seq2);
                
                if (@aligned_seq2)
                {
                    $aligned_seqs2 = $aligned_seqs2.",";
                }
            }
            
            my $cuffdiff_cmd = ($::global_cfg -> param("sprograms.Cuffdiff")).
                " -q -o ".($::global_cfg -> param("results.tmp"))."/cuffdiff_output".
                " -b ".$ref_genome.
                " -p ".($::global_cfg -> param("cmdline.threads")).
                " -L ".($::global_cfg -> param("required.sequence_1")).",".
                ($::global_cfg -> param("required.sequence_2")).
                " -u ".$merged_gtf.
                " ".$aligned_seqs1." ".$aligned_seqs2;
            
            PIPELINETOOLS::PrintStd($cuffdiff_cmd);
            
            my $cmd = `$cuffdiff_cmd`;
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: Can't run cuffdiff!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);
        }
        elsif (($::global_cfg -> param("required.analysistype")) eq
            "pairedindividual")
        {
            my $i = 1;
            while (@aligned_seq1)
            {
                 my $cuffdiff_cmd =
                    ($::global_cfg -> param("sprograms.Cuffdiff")).
                    " -q -o ".($::global_cfg -> param("results.tmp"))
                        ."/cuffdiff_output_$i".
                    " -b ".$ref_genome.
                    " -p ".($::global_cfg -> param("cmdline.threads")).
                    " -L ".($::global_cfg -> param("required.sequence_1")).",".
                    ($::global_cfg -> param("required.sequence_2")).
                    " -u ".$merged_gtf.
                    " ".shift(@aligned_seq1)." ".shift(@aligned_seq2);
                
                PIPELINETOOLS::PrintStd($cuffdiff_cmd);
                
                my $cmd = `$cuffdiff_cmd`;
                if ($cmd -> success() != 1)
                {
                    PIPELINETOOLS::PrintStd($cmd -> merged);
                    PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                    die "ERROR: Can't run cuffdiff!\n";
                }
                
                PIPELINETOOLS::PrintStd($cmd -> stdout);
                PIPELINETOOLS::WriteLogOut($cmd -> stdout);
                
                $i++;
            }
        }
        else
        {
            die "ERROR: Paired Cuffdiff not available yet!\n";
        }        
    }
    elsif (!-e $ref_genome)
    {
        die "ERROR: Can't find reference genome *.fasta - file";
    }
    elsif (!-e $merged_gtf)
    {
        die "ERROR: Can't find reference merged *.gtf- file";
    }
    
    if (-e (($::global_cfg -> param("results.tmp"))."/merged_asm/merged.gtf"))
    {
        my $mv_cmd = "mv ".
            ($::global_cfg -> param("results.tmp"))."/merged_asm/merged.gtf".
            " ".($::global_cfg -> param("results.results"));
        my $cmd = `$mv_cmd`;
        if ($cmd -> success() != 1)
        {
            PIPELINETOOLS::PrintStd($cmd -> merged);
            PIPELINETOOLS::WriteLogOut(($cmd -> merged));
            print "Warning: Not able to move merged.gtf-file to
                Results folder!\n";
        }
    }
    

    
    
    PIPELINETOOLS::CreateSuccess("Cuffdiff");
    PIPELINETOOLS::PrintTask("Finished Cuffdiff");
    
    return;
} 

#-------------------------------------------------------------------------------
# CutadaptRun: The function creates the command for a proper call of
#              Cufflinks and then starts Cufflinks.
#-------------------------------------------------------------------------------
sub CummerbundRun{
    PIPELINETOOLS::WriteLogChapter("Cummerbund");
    PIPELINETOOLS::PrintTask("Started CummeRbund");
    
    my $cummerbund_cmd = ($::global_cfg -> param("sprograms.Rscript")).
        " ".($::global_cfg -> param("required.sequence_1")).
        " ".($::global_cfg -> param("required.sequence_2")).
        " ".($::global_cfg -> param("results.tmp"))."/cuffdiff_output".
        " ".($::global_cfg -> param("results.tmp")).
        " ".($::global_cfg -> param("results.graphics")).
        " ".($::global_cfg -> param("results.results"));
    
    PIPELINETOOLS::PrintStd($cummerbund_cmd);    
    my $cmd = `$cummerbund_cmd`;
    if ($cmd -> success() != 1)
    {
        PIPELINETOOLS::PrintStd($cmd -> merged);
        PIPELINETOOLS::WriteLogOut(($cmd -> merged));
        die "ERROR: Can't run CummeRbund!\n";
    }
    
    PIPELINETOOLS::PrintStd($cmd -> stdout);
    PIPELINETOOLS::WriteLogOut($cmd -> stdout);
    

    PIPELINETOOLS::CreateSuccess("CummeRbund"); 
    PIPELINETOOLS::PrintTask("Finished CummeRbund");
    
    return;
}


1;