#------------------------------------------------------------------------------
# CancerRNASeq - a pipeline for transcriptomics datasets
# Author: Raphael Hablesreiter
#
# DRIVER.pm
#          
# Description: Driver Module
#              This module checks, if the pipeline has to start at the begining
#              or if the pipeline has to start at a certain module.
#              There are 8 different options available:
#              (1. Preprocessing)
#                 1.1. Cutadapt
#                 1.2. Trimmomatic
#                 1.3. rCorrector
#              2. Mapping
#              (3. Analyzing)
#                 3.1. Cufflinks
#                 3.2. Cuffmerge
#                 3.3. Cuffdiff
#                 3.4. CummeRbund 
#
#-------------------------------------------------------------------------------
#

package DRIVER;

use 5.010;
use strict;
use warnings;
use Config::Simple;


use PIPELINETOOLS;
use PREPROCESSING;
use MAPPING;
use ANALYZING;

#-------------------------------------------------------------------------------
# Main: Checks the progress of the Pipeline and after that
#       starts all the different modules.
#-------------------------------------------------------------------------------
sub Main{
    
    my $end_condition = 0;
    my $succ = DRIVER::CheckSuccess();

    if ($::global_cfg -> param("optional.end") &&
        $::global_cfg -> param("optional.end") eq $succ)
    {
        $end_condition = 1;
    }
    elsif ($succ eq "Cummerbund")
    {
        $end_condition = 1;
    }
        
    #Pipeline Run
    while ($end_condition == 0)
    {        
        given ($succ)
        {
            when ("Start") {PREPROCESSING::Main(1); $succ = "Cutadapt";}
            when ("Cutadapt") {PREPROCESSING::Main(2); $succ = "Trimmomatic";}
            when ("Trimmomatic") {PREPROCESSING::Main(3); $succ = "rCorrector";}
            when ("rCorrector") {MAPPING::Main(1); $succ = "Mapping";}
            when ("Mapping") {MAPPING::Main(2); $succ = "QualiMap";}
            when ("QualiMap") {ANALYZING::Main(1); $succ = "Cufflinks";}
            when ("Cufflinks") {ANALYZING::Main(2); $succ = "Cuffmerge";}
            when ("Cuffmerge") {ANALYZING::Main(3); $succ = "Cuffdiff";}
            when ("Cuffdiff") {ANALYZING::Main(4); $succ = "CummeRbund";}
            when ("CummeRbund") {}
            when ("ERROR") {die "ERROR: Can't read success files!\n";}
            default {die "ERROR: In SUCCES Routine!\n";}
        
        }
        
        # $::global_cfg -> param("optional.end") can also have
        # the same values as the return values of CheckSuccess are
        # if value of $::global_cfg -> param("optional.end") and
        # pipeline progress is equal the pipeline stops
        if ($::global_cfg -> param("optional.end") &&
            $::global_cfg -> param("optional.end") eq $succ)
        {
            $end_condition = 1;
        }
        elsif ($succ eq "Cummerbund")
        {
            $end_condition = 1;
        }
        
    }
    
    PIPELINETOOLS::PrintStd("Pipeline Finished");
    PIPELINETOOLS::WriteLogOut("Pipeline Finished");
        
    return;
}

#-------------------------------------------------------------------------------
# CheckSuccess: Checks the progress of the Pipeline and returns last finished
#               pipeline module.
#               return values: "Start" -> starts with Cutadapt
#                              "Cutadapt" -> starts with Trimmomatic
#                              "Trimmomatic" -> starts with rCorrector
#                              "rCorrector" -> starts with Mapping
#                              "Mapping" -> starts with Cufflinks
#                              "Cufflinks" -> starts with Cuffmerge
#                              "Cuffmerge" -> starts with Cuffdiff
#                              "Cuffdiff" -> starts with CummeRbund
# 
#-------------------------------------------------------------------------------
sub CheckSuccess{
    my $succ_dir = $::global_cfg -> param("results.tmp");
	my @succ_files = PIPELINETOOLS::ReadDirectory($succ_dir,".succ");
	
    if (@succ_files)
    {
        @succ_files = sort(@succ_files);
        my $last_succ_name = pop(@succ_files);
        
        my $last_succ = new Config::Simple($succ_dir."/".$last_succ_name);
        
        if ($last_succ)
        {
            my $last_cmd = $last_succ -> param("success.cmd");            
            
            if ($last_cmd eq "Cutadapt" or
                $last_cmd eq "Trimmomatic" or
                $last_cmd eq "rCorrector")
            {
                @::sequences_1 = ();
                @::sequences_2 = ();
                
                my @tmp_files = PIPELINETOOLS::ReadDirectory(
                    $succ_dir,".fastq.gz");
                my $namead;
                
                given($last_cmd)
                {
                    when("Cutadapt"){$namead = "_CA.fastq.gz";}
                    when("Trimmomatic"){$namead = "_CA_Trim_P.fastq.gz";}
                    when("rCorrector"){$namead = "_CA_Trim_P.cor.fq.gz";}
                    default{die "ERROR! In reading successfiles!\n";}
                }
                
                my @sequences = PIPELINETOOLS::ReadDirectory($succ_dir,$namead);
                
                while (@sequences)
                {
                    my $tmp_seq = shift(@sequences);
                    
                    if ($tmp_seq =~ $::global_cfg -> param("required.sequence_1"))
                    {
                        push(@::sequences_1, $succ_dir."/".$tmp_seq);
                    }
                    elsif ($tmp_seq =~ $::global_cfg -> param("required.sequence_2"))
                    {
                        push(@::sequences_2, $succ_dir."/".$tmp_seq);
                    }
                }
                
                @::sequences_1 = sort(@::sequences_1);
                @::sequences_2 = sort(@::sequences_2);
                
                $::global_cfg -> param("required.seqlength1", scalar(@::sequences_1));
                $::global_cfg -> param("required.seqlength2", scalar(@::sequences_2));
                
                
               
            }
            
            return $last_cmd;
        }
        else
        {
            die "ERROR: Wrong success-file in the tmp-folder!\n";
        }
    }
    else
    {
        return "Start";
    }    
}


1;