#------------------------------------------------------------------------------
# CancerRNASeq - a pipeline for transcriptomics datasets
# Author: Raphael Hablesreiter
#
# Preprocessing.pm
#          
# Description: Preprocessing Module
#              In this module the raw reads are preprocessed, with three 
#              consecuitively executed tools (i.e., CutAdapt, Trimmomatic,
#              rCorrector). CutAdapt removes adapter sequences from the raw
#              reads. Trimmomatic is used to remove bases with a quality under
#              a certain threshold and too short sequences. The last tool,
#              rCorrector, is used to repair sequencing errors.
#
#-------------------------------------------------------------------------------
#

package PREPROCESSING;

use 5.010;
use strict;
use warnings;
	
use Bio::SeqIO;
use File::chdir;
use Backticks;

use PIPELINETOOLS;

#-------------------------------------------------------------------------------
# Main: Starts one of the 3 parts of the preprocessing module.
#       Input is provided by the Driver module.
#       Options: 1 -> Cutadapt
#                2 -> Trimmomatic
#                3 -> rCorrector
#-------------------------------------------------------------------------------
sub Main{
    
    my $i = $_[0];
    
    given ($i)
    {
        when (1) {PREPROCESSING::CutadaptRun();}
        when (2) {PREPROCESSING::TrimmomaticRun();}
        when (3) {PREPROCESSING::RCorrectorRun();}
        default {die "ERROR:Cant find out where to start the preprocessing!\n";}
    }
    
    return;
}

#-------------------------------------------------------------------------------
# CutadaptRun: The function creates the command for a proper call of
#              Cutadapt and then starts Cutadapt.
#-------------------------------------------------------------------------------
sub CutadaptRun{
    
    PIPELINETOOLS::WriteLogChapter("Cutadapt");
    PIPELINETOOLS::PrintTask("Started Cutadapt");
    
    my $s_adapterscmd = "";
    my $s_qualitycmd = " -m 20";
    
    my @in;
    $in[0] = $::global_cfg -> param("results.tmp");
    $in[1] = $::global_cfg -> param("sfiles.AdapterFasta");
    $in[2] = $::global_cfg -> param("optional.cutadapt_cmds");
    $in[3] = $::global_cfg -> param("required.seqtype");
    $in[4] = $::global_cfg -> param("required.sequence_1");
    $in[5] = $::global_cfg -> param("required.sequence_2");          
    
    if ($in[2])
    {
        if ($in[2] =~ "-m")
        {
            $s_qualitycmd = " ".$in[2];
        }
        else
        {
            $s_qualitycmd = $s_qualitycmd." ".$in[2];
        }
    } 
    
    
    my @a_adaptersequences =
        PREPROCESSING::ReadAdapterSequences("$in[1]");
    
    #http://stackoverflow.com/questions/16626891/determine-if-a-string-starts-with-an-upper-case-character-in-perl
    while (@a_adaptersequences)
    {
        my $a_type = shift(@a_adaptersequences);
        
        if ($in[3] eq "SE" && $a_type =~ /^[[:upper:]]/ )
        {
            shift(@a_adaptersequences);        
        }
        else
        {
            $s_adapterscmd = $s_adapterscmd."-".$a_type.
                " ".shift(@a_adaptersequences)." ";
        }
    }

    PIPELINETOOLS::ChangeWorkingDir($::global_cfg -> param("results.tmp"));
    
    if ($in[3] eq "PE")
    {
        my $m = 0;
        my $sample_nr = 1;
        while ($m < ($::global_cfg -> param("required.seqlength1")))
        {
            my $seq_1 = $::sequences_1[$m];
            my $seq_2 = $::sequences_1[$m + 1];          
        
            PREPROCESSING::CreateFastqc($seq_1);
            PREPROCESSING::CreateFastqc($seq_2);
                
            if (int($sample_nr) < 10)
            {
                $sample_nr = "0".int($sample_nr);
            }                
                
            my $s_cutadaptcmd = "cutadapt ".$s_adapterscmd."$s_qualitycmd"." -o ".
                $in[4]."_".$sample_nr."_1_CA.fastq.gz -p ".$in[4]."_".$sample_nr.
                "_2_CA.fastq.gz ".$seq_1." ".$seq_2;
            
            my $cmd = `$s_cutadaptcmd`; 
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: Problem with 1. pair of Sequences in CUTADAPT!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);
            
            $::sequences_1[$m] = $in[4]."_".$sample_nr."_1_CA.fastq.gz";
            $::sequences_1[$m + 1] = $in[4]."_".$sample_nr."_2_CA.fastq.gz ";
            
            PREPROCESSING::CreateFastqc($::sequences_1[$m]);
            PREPROCESSING::CreateFastqc($::sequences_1[$m + 1]);
        
            $m = $m + 2;
            $sample_nr++;
        }
        
        $m = 0;
        $sample_nr  = 1;
        while ($m < ($::global_cfg -> param("required.seqlength2")))
        {
            my $seq_1 = $::sequences_2[$m];
            my $seq_2 = $::sequences_2[$m + 1];           
            
            PREPROCESSING::CreateFastqc($seq_1);
            PREPROCESSING::CreateFastqc($seq_2);
            
            if (int($sample_nr) < 10)
            {
                $sample_nr = "0".int($sample_nr);
            }
            
            my $s_cutadaptcmd = "cutadapt ".$s_adapterscmd."$s_qualitycmd"." -o ".
                $in[5]."_".$sample_nr."_1_CA.fastq.gz -p ".$in[5]."_".$sample_nr.
                "_2_CA.fastq.gz ".$seq_1." ".$seq_2;
            
            my $cmd = `$s_cutadaptcmd`; 
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: Problem with 1. pair of Sequences in CUTADAPT!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);
            
            $::sequences_2[$m] = $in[5]."_".$sample_nr."_1_CA.fastq.gz";
            $::sequences_2[$m + 1] = $in[5]."_".$sample_nr."_2_CA.fastq.gz ";
            
            PREPROCESSING::CreateFastqc($::sequences_2[$m]);
            PREPROCESSING::CreateFastqc($::sequences_2[$m + 1]);
        
            $m = $m + 2;
            $sample_nr++;
        }
        
    }
    else
    {
        my $m = 0;
        my $sample_nr = 1;
        while ($m < ($::global_cfg -> param("required.seqlength1")))
        {
            my $seq_1 = $::sequences_1[$m];       
            
            PREPROCESSING::CreateFastqc($seq_1);
            
            if (int($sample_nr) < 10)
            {
                $sample_nr = "0".int($sample_nr);
            }
            
            my $s_cutadaptcmd = "cutadapt ".$s_adapterscmd."$s_qualitycmd"." -o ".
                $in[4]."_".$sample_nr."_CA.fastq.gz ".$seq_1;
                
            my $cmd = `$s_cutadaptcmd`; 
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: Problem with ".$seq_1." in CUTADAPT!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);
            
            $::sequences_1[$m] = $in[4]."_".$sample_nr."_CA.fastq.gz";
            
            PREPROCESSING::CreateFastqc($::sequences_1[$m]);

            $m++;
            $sample_nr++;
        }
        
        $m = 0;
        $sample_nr = 1;
        while ($m < ($::global_cfg -> param("required.seqlength2")))
        {
            my $seq_2 = $::sequences_2[$m];       
            
            PREPROCESSING::CreateFastqc($seq_2);
            
            if (int($sample_nr) < 10)
            {
                $sample_nr = "0".int($sample_nr);
            }
            
            my $s_cutadaptcmd = "cutadapt ".$s_adapterscmd."$s_qualitycmd"." -o ".
                $in[5]."_".$sample_nr."_CA.fastq.gz ".$seq_2;
                
            my $cmd = `$s_cutadaptcmd`; 
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: Problem with ".$seq_2." in CUTADAPT!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);
            
            $::sequences_2[$m] = $in[5]."_".$sample_nr."_CA.fastq.gz";
            
            PREPROCESSING::CreateFastqc($::sequences_2[$m]);

            $m++;    
        }       

    }
    
    PIPELINETOOLS::CreateSuccess("Cutadapt");
    PIPELINETOOLS::PrintTask("Finished Cutadapt");    
    
    return 0;
}

#-------------------------------------------------------------------------------
# TrimmomaticRun: The function creates the command for a proper call of
#                 Trimmomatic and then starts Trimmomatic.
#-------------------------------------------------------------------------------
sub TrimmomaticRun{
    PIPELINETOOLS::WriteLogChapter("Trimmomatic");
    PIPELINETOOLS::PrintTask("Started Trimmomatic");
    
    my $qualityscore = $::global_cfg -> param("required.qscore");
    if ($qualityscore ne "phred33" or $qualityscore ne "phred64")
    {
        $::global_cfg -> param("required.qscore","phred33");    
    }
    
    my @in;
    $in[0] = $::global_cfg -> param("sprograms.Trimmomatic");
    $in[1] = $::global_cfg -> param("scmd.Trimmomatic");
    $in[2] = $::global_cfg -> param("required.qscore");
    $in[3] = $::global_cfg -> param("cmdline.threads");
    $in[4] = ($::global_cfg -> param("required.sequence_1"));
    $in[5] = ($::global_cfg -> param("required.sequence_2"));
  
    
    PIPELINETOOLS::ChangeWorkingDir($::global_cfg -> param("results.tmp"));
    
    if ($::global_cfg -> param("required.seqtype") eq "PE")
    {
        my $m = 0;
        my $sample_nr = 1;
        while ($m < ($::global_cfg -> param("required.seqlength1")))
        {
            if (int($sample_nr) < 10)
            {
                $sample_nr = "0".int($sample_nr);
            }
            
            my $trimmomatic_cmd = "java -jar ".$in[0]." PE".
                " -threads ".$in[3].
                " -".$in[2].
                " ".$::sequences_1[$m]." ".$::sequences_1[$m + 1].
                " ".$in[4]."_".$sample_nr."_1_CA_Trim_P.fastq.gz".
                " ".$in[4]."_".$sample_nr."_1_CA_Trim_U.fastq.gz".
                " ".$in[4]."_".$sample_nr."_2_CA_Trim_P.fastq.gz".
                " ".$in[4]."_".$sample_nr."_2_CA_Trim_U.fastq.gz".
                " ".$in[1];
            
            PIPELINETOOLS::PrintStd($trimmomatic_cmd);
            
            my $cmd = `$trimmomatic_cmd`;
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: Can't execude Trimmomatic for ".
                    $::sequences_1[$m]." & ".$::sequences_1[$m]."!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);

            PIPELINETOOLS::MoveFile($::sequences_1[$m],
                ($::global_cfg -> param("results.preprocess"))."/"); 
            PIPELINETOOLS::MoveFile($::sequences_1[$m + 1],
                ($::global_cfg -> param("results.preprocess"))."/");
            
            $::sequences_1[$m] = $in[4]."_".$sample_nr."_1_CA_Trim_P.fastq.gz";
            $::sequences_1[$m + 1] = $in[4]."_".$sample_nr."_2_CA_Trim_P.fastq.gz";
            
            PREPROCESSING::CreateFastqc($::sequences_1[$m]);
            PREPROCESSING::CreateFastqc($::sequences_1[$m + 1]);
            
            $m = $m + 2;
            $sample_nr++;
        }
        
        $m = 0;
        $sample_nr = 1;
        while ($m < ($::global_cfg -> param("required.seqlength2")))
        {
            if (int($sample_nr) < 10)
            {
                $sample_nr = "0".int($sample_nr);
            }
            
            my $trimmomatic_cmd = "java -jar ".$in[0]." PE".
                " -threads ".$in[3].
                " -".$in[2].
                " ".$::sequences_2[$m]." ".$::sequences_2[$m + 1].
                " ".$in[5]."_".$sample_nr."_1_CA_Trim_P.fastq.gz".
                " ".$in[5]."_".$sample_nr."_1_CA_Trim_U.fastq.gz".
                " ".$in[5]."_".$sample_nr."_2_CA_Trim_P.fastq.gz".
                " ".$in[5]."_".$sample_nr."_2_CA_Trim_U.fastq.gz".
                " ".$in[1];
            
            PIPELINETOOLS::PrintStd($trimmomatic_cmd);
            
            my $cmd = `$trimmomatic_cmd`;
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: Can't execude Trimmomatic for ".
                    $::sequences_2[$m]." & ".$::sequences_2[$m]."!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);

            PIPELINETOOLS::MoveFile($::sequences_2[$m],
                ($::global_cfg -> param("results.preprocess"))."/"); 
            PIPELINETOOLS::MoveFile($::sequences_2[$m + 1],
                ($::global_cfg -> param("results.preprocess"))."/");
            
            $::sequences_2[$m] = $in[5]."_".$sample_nr."_1_CA_Trim_P.fastq.gz";
            $::sequences_2[$m + 1] = $in[5]."_".$sample_nr."_2_CA_Trim_P.fastq.gz";
            
            PREPROCESSING::CreateFastqc($::sequences_2[$m]);
            PREPROCESSING::CreateFastqc($::sequences_2[$m + 1]);
            
            $m = $m + 2;
            $sample_nr++;
        }   
    }
    else
    {
        my $m = 0;
        my $sample_nr = 1;
        
        while ($m < ($::global_cfg -> param("required.seqlength1")))
        {
            if (int($sample_nr) < 10)
            {
                $sample_nr = "0".int($sample_nr);
            }
            
            my $trimmomatic_cmd = "java -jar ".$in[0]." SE".
                " -threads ".$in[3].
                " - ".$in[2].
                " ".$::sequences_1[$m].
                " ".$in[4]."_".$sample_nr."_CA_Trim_P.fastq.gz".
                " ".$in[1];

            PIPELINETOOLS::PrintStd($trimmomatic_cmd);
            
            my $cmd = `$trimmomatic_cmd`;
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: Can't execude Trimmomatic for".
                    $::sequences_1[$m]."!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);
              
            PIPELINETOOLS::MoveFile($::sequences_1[$m],
                ($::global_cfg -> param("results.preprocess"))."/");
             
            $::sequences_1[$m] = $in[4]."_".$sample_nr."_CA_Trim_P.fastq.gz";
            
            PREPROCESSING::CreateFastqc($::sequences_1[$m]);
            
            $m++;
            $sample_nr++;
        }
        
        $m = 0;
        $sample_nr = 1;
        while ($m < ($::global_cfg -> param("required.seqlength2")))
        {
            if (int($sample_nr) < 10)
            {
                $sample_nr = "0".int($sample_nr);
            }
            
            my $trimmomatic_cmd = "java -jar ".$in[0]." SE".
                " -threads ".$in[3].
                " - ".$in[2].
                " ".$::sequences_2[$m].
                " ".$in[5]."_".$sample_nr."_CA_Trim_P.fastq.gz".
                " ".$in[1];

            PIPELINETOOLS::PrintStd($trimmomatic_cmd);
            
            my $cmd = `$trimmomatic_cmd`;
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: Can't execude Trimmomatic for".
                    $::sequences_2[$m]."!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);

            PIPELINETOOLS::MoveFile($::sequences_2[$m],
                ($::global_cfg -> param("results.preprocess"))."/");
            
            $::sequences_2[$m] = $in[5]."_".$sample_nr."_CA_Trim_P.fastq.gz";
            
            PREPROCESSING::CreateFastqc($::sequences_2[$m]);
            
            $m++;
        }
    }
    
    PIPELINETOOLS::CreateSuccess("Trimmomatic");
    PIPELINETOOLS::PrintTask("Finished Trimmomatic");
    
    return;
}
#-------------------------------------------------------------------------------
# RCorrectorRun: The function creates the command for a proper call of
#                rCorrector and then starts rCorrector.
#-------------------------------------------------------------------------------
sub RCorrectorRun
{
    PIPELINETOOLS::WriteLogChapter("rCorrector");
    PIPELINETOOLS::PrintTask("Started rCorrector");
    
    my @in;
    $in[0] = $::global_cfg -> param("sprograms.rCorrector");
    $in[1] = $::global_cfg -> param("cmdline.threads");
    $in[2] = ($::global_cfg -> param("required.sequence_1"));
    $in[3] = ($::global_cfg -> param("required.sequence_2"));

    PIPELINETOOLS::ChangeWorkingDir($::global_cfg -> param("results.tmp"));
     
    if ($::global_cfg -> param("required.seqtype") eq "PE")
    {
        my $m = 0;
        my $sample_nr = 1;
        
        while ($m < ($::global_cfg -> param("required.seqlength1")))
        {
            if (int($sample_nr) < 10)
            {
                $sample_nr = "0".int($sample_nr);
            }
            
            my $rcorrector_cmd =
                "perl ".$in[0].
                " -1 ".$::sequences_1[$m].
                " -2 ".$::sequences_1[$m + 1].
                " -t ".$in[1];        
            
            PIPELINETOOLS::PrintStd($rcorrector_cmd);
            
            my $cmd = `$rcorrector_cmd`;
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: rCorrector problem with Sequences".
                    $::sequences_1[$m]." & ".$::sequences_1[$m + 1]."!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);

            PIPELINETOOLS::MoveFile($::sequences_1[$m],
                ($::global_cfg -> param("results.preprocess"))."/"); 
            PIPELINETOOLS::MoveFile($::sequences_1[$m + 1],
                ($::global_cfg -> param("results.preprocess"))."/");
            
            $::sequences_1[$m] =
                $in[2]."_".$sample_nr."_1_CA_Trim_P.cor.fq.gz";
            $::sequences_1[$m + 1] =
                $in[2]."_".$sample_nr."_2_CA_Trim_P.cor.fq.gz";
            
            PREPROCESSING::CreateFastqc($::sequences_1[$m]);
            PREPROCESSING::CreateFastqc($::sequences_1[$m + 1]);
            
            $m = $m + 2;
            $sample_nr++;
        }
        
        $m = 0;
        $sample_nr = 1;
        while ($m < ($::global_cfg -> param("required.seqlength2")))
        {
            if (int($sample_nr) < 10)
            {
                $sample_nr = "0".int($sample_nr);
            }
            
            my $rcorrector_cmd =
                "perl ".$in[0].
                " -1 ".$::sequences_2[$m].
                " -2 ".$::sequences_2[$m + 1].
                " -t ".$in[1];        
            
            PIPELINETOOLS::PrintStd($rcorrector_cmd);
            
            my $cmd = `$rcorrector_cmd`;
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: rCorrector problem with Sequences".
                    $::sequences_2[$m]." & ".$::sequences_2[$m + 1]."!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);

            PIPELINETOOLS::MoveFile($::sequences_2[$m],
                ($::global_cfg -> param("results.preprocess"))."/"); 
            PIPELINETOOLS::MoveFile($::sequences_2[$m + 1],
                ($::global_cfg -> param("results.preprocess"))."/");
            
            $::sequences_2[$m] = $in[3]."_".$sample_nr."_1_CA_Trim_P.cor.fq.gz";
            $::sequences_2[$m + 1] = $in[3]."_".$sample_nr."_2_CA_Trim_P.cor.fq.gz";
            
            PREPROCESSING::CreateFastqc($::sequences_2[$m]);
            PREPROCESSING::CreateFastqc($::sequences_2[$m + 1]);
            
            $m = $m + 2;
            $sample_nr++;
        }       
    }
    else
    {
        
        my $m = 0;
        my $sample_nr = 1;
        while ($m < ($::global_cfg -> param("required.seqlength1")))
        {
            if (int($sample_nr) < 10)
            {
                $sample_nr = "0".int($sample_nr);
            }
            
            my $rcorrector_cmd =
                "perl ".$in[0].
                " -s ".$::sequences_1[$m].
                " -t ".$in[1];       
            
            PIPELINETOOLS::PrintStd($rcorrector_cmd);
            
            my $cmd = `$rcorrector_cmd`;
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: rCorrector problem with Sequence".
                    $::sequences_1[$m]."!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);

            PIPELINETOOLS::MoveFile($::sequences_1[$m],
                ($::global_cfg -> param("results.preprocess"))."/"); 
            
            $::sequences_1[$m] =
                $in[2]."_".$sample_nr."_CA_Trim_P.cor.fq.gz";
            
            PREPROCESSING::CreateFastqc($::sequences_1[$m]);

            $m++;
            $sample_nr++;
        }
        
        $m = 0;
        $sample_nr = 1;
        while ($m < ($::global_cfg -> param("required.seqlength2")))
        {
            if (int($sample_nr) < 10)
            {
                $sample_nr = "0".int($sample_nr);
            }
            
            my $rcorrector_cmd =
                "perl ".$in[0].
                " -s ".$::sequences_2[$m].
                " -t ".$in[1];       
            
            PIPELINETOOLS::PrintStd($rcorrector_cmd);
            
            my $cmd = `$rcorrector_cmd`;
            if ($cmd -> success() != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                die "ERROR: rCorrector problem with Sequence".
                    $::sequences_2[$m]."!\n";
            }
            
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout);

            PIPELINETOOLS::MoveFile($::sequences_2[$m],
                ($::global_cfg -> param("results.preprocess"))."/"); 
            
            $::sequences_2[$m] =
                $in[3]."_".$sample_nr."_CA_Trim_P.cor.fq.gz";
            
            PREPROCESSING::CreateFastqc($::sequences_2[$m]);

            $m++;
            $sample_nr++;
        }        
    }
    PIPELINETOOLS::CreateSuccess("rCorrector");    
    PIPELINETOOLS::PrintTask("Finished rCorrector");

    return;
}

#-------------------------------------------------------------------------------
# CreateFastqc: The function creates the command for a proper call of
#               FASTQC and then starts FASTQC.
#               Input: *.fastq(.gz)-File
#-------------------------------------------------------------------------------
sub CreateFastqc{
    
    my $file = $_[0];
    my $dir = $::global_cfg -> param("results.qreports");
    
    my $cmd_fastqc =  "fastqc -o \"$dir\" $file";
    
    my $system_call = `$cmd_fastqc`;
    if ($system_call -> success() != 1)
    {
        PIPELINETOOLS::PrintStd($system_call -> merged);
        PIPELINETOOLS::WriteLogOut($system_call -> merged);
    }
    PIPELINETOOLS::WriteLogOut($system_call -> stdout);
    
    my @zipfiles = PIPELINETOOLS::ReadDirectory($dir,".zip");
    
    while (@zipfiles)
    {
        unlink $dir."/".pop(@zipfiles);
    }
    
    return;
}

#-------------------------------------------------------------------------------
# ReadAdapterSequences: This function reads a fasta file, which provides the
#                       adapter sequences needed for cutadapt. 
#                       Input: *.fasta-File
#                       Return: Array, with the command parameter and sequences
#                               alternating  
#-------------------------------------------------------------------------------
sub ReadAdapterSequences{
    
    my $file = $_[0];
    my $a_sequencefile;
    my @a_adaptersequences;
    
    if(-e $file)
    {
        $a_sequencefile = Bio::SeqIO -> new (-file => $file,
            -format => 'Fasta');        
    }
    else
    {
        die "ERROR: default Fasta File not found!";
    }
    
    while (my $s_seq = $a_sequencefile -> next_seq)
    {    
        push(@a_adaptersequences, $s_seq -> id());
        push(@a_adaptersequences, $s_seq -> seq());    
    }
    
    return @a_adaptersequences;
}

1;