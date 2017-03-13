#------------------------------------------------------------------------------
# CancerRNASeq - a pipeline for transcriptomics datasets
# Author: Raphael Hablesreiter
#
# CONFIGURATION.pm
#          
# Description: Configuration Module
#              This module handles the configuration of the pipeline, which
#              means that all the necassery parameters for the pipeline are set
#              in this module.
#
#-------------------------------------------------------------------------------
#

package CONFIGURE;

use 5.010;
use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use Config::Simple;
use Backticks;
use PIPELINETOOLS;
 
sub Main{
    my $check_value = 0;
    
    #store values of the commandline input into the global variable
    $check_value = CONFIGURE::GetCmdOptions();
    
    if ($check_value == 0)
    {
        PIPELINETOOLS::PrintTask("Started To Check *.config-Files");
        #store values of the system configuration file into the global variable
        CONFIGURE::ReadConfigFile("/export/home/rhablesreiter_914/Pipeline_Version1/Pipeline/Program/sys.ini");
        #store values of the user configuration file into the global variable
        CONFIGURE::ReadConfigFile($::global_cfg -> param("cmdline.usrconfig"));
        CONFIGURE::CheckUsrConfig();
        PIPELINETOOLS::PrintTask("Finished To Check *.config-Files");
        return 0;
    }
    else
    {
        return -1;    
    }
}

#-------------------------------------------------------------------------------
# GetCmdOptions: This function reads all the commandline parameters, with the
#                the help of the package Getopt::Long, and saves the
#                parameters into the global variable
#-------------------------------------------------------------------------------
sub GetCmdOptions{
    #From
	#https://metacpan.org/pod/Getopt::Long
	
	#GetOptions() parses the command line from @ARGV,
	#recognizing and removing specified options and their possible values.
	
	#pass_through (default: disabled)
	#With pass_through anything that is unknown, ambiguous
	#or supplied with an invalid option will not be flagged as an error.
	#Instead the unknown option(s) will be passed to the catchall <> if present,
	#otherwise through to @ARGV.
	Getopt::Long::Configure("pass_through");
    
    
    my %input_options;
	my $values_chk = 0;
	my $key;
	my $value;
	
    
	GetOptions(
		'config:s' => \$input_options{usrconfig},
		'help' => \$input_options{help},
		'silent' => \$input_options{silent},
        'threads:i' => \$input_options{threads},
        'overwrite'=> \$input_options{overwrite}
	);
	
	
	#From
	#http://www.perlhowto.com/iterate_through_a_hash
	while (($key, $value) = each %input_options)
	{
		if ($key eq "help")
		{
            if ($value)
            {
                $values_chk += 1;
            }
        }
        elsif($key eq "silent")
        {
            if ($value)
            {
                $::global_cfg -> param("cmdline.$key","1");
            }        
        }
        elsif($key eq "overwrite")
        {
            if ($value)
            {
                $::global_cfg -> param("cmdline.$key","1");
            }        
        }
        
        elsif($key eq "threads")
        {
            if ($value)
            {
                if ($value > 0)
                {
                    $::global_cfg -> param("cmdline.$key",$value);
                }
                else
                {
                    $::global_cfg -> param("cmdline.$key","1");
                }
            }
            else
            {
                $::global_cfg -> param("cmdline.$key","1");    
            }
        }
        elsif($value eq "")
        {
            $values_chk += 1;
        }
        else
        {
            $::global_cfg -> param("cmdline.$key","$value");
        }
	}
    
    if ($values_chk != 0 && $input_options{help})
    {
        #help has to be printed
        CONFIGURE::PrintHelp();
        return -1;    
    }
    elsif ($values_chk != 0)
    {
        CONFIGURE::PrintUsage();
        return -1;
    }
    else
    {   
        return 0;
    }
    
}

#-------------------------------------------------------------------------------
# GetCmdOptions: This function takes as input the name of a configuration file.
#                The Configuration file has to be a INI-file, which is saved
#                into the global variable.
#                See Cinfig::Simple
#-------------------------------------------------------------------------------
sub ReadConfigFile{
    #See
    #http://search.cpan.org/~sherzodr/Config-Simple-4.59/Simple.pm
    my $cfg_name = $_[0];
    
	if (-e $cfg_name)
    {
        
        my $system_cfg = new Config::Simple($cfg_name);
        my @cfg_layout = $system_cfg -> get_block();
        
        
        while (@cfg_layout)
        {
            my $current_block = shift(@cfg_layout);
            $::global_cfg -> param(-block => $current_block,
                -values => ($system_cfg -> param(-block => $current_block)));  
        }    
    }
    else
    {
        die "Config-File not found.";
    }
    
    return;
}

#-------------------------------------------------------------------------------
# CheckUsrConfig: This function checks the necessary parameters for a smooth
#                 run of the pipeline.
#-------------------------------------------------------------------------------
sub CheckUsrConfig{
    my $fail = "";
    my $dir_raw = $::global_cfg -> param("sdirectories.RawData");

    
    my @tc_info = ("User", "ProjectName", "sequence_1",
        "sequence_2", "ReferenceGenome");
    my @tc_files = ("GenomeFASTA", "GenomeGTF");
    
    while (@tc_info)
    {
        my $tmp_name = shift(@tc_info);
        
        if (! $::global_cfg -> param("required.$tmp_name"))
        {
            $fail = $fail."$tmp_name is not provided!\n";
        }
        else
        {
            my $tmp = $::global_cfg -> param("required.$tmp_name");
            $tmp =~ s/ /_/g;
            $::global_cfg -> param("required.$tmp_name",$tmp);
        }
    }
    
    if (! $::global_cfg -> param("required.seqtype"))
    {
        $fail = $fail."Sequences type is not provided!\n";
    }
    elsif($::global_cfg -> param("required.seqtype") ne "PE" and
          $::global_cfg -> param("required.seqtype") ne "SE")
    {
        $fail = $fail."Sequences type is not known!\n";
    }   
    
    if (! $::global_cfg -> param("required.MappingType"))
    {
        $fail = $fail."Mapping type is not provided!\n";
    }
    elsif($::global_cfg -> param("required.MappingType") ne "Tophat" and
          $::global_cfg -> param("required.MappingType") ne "STAR" and
           $::global_cfg -> param("required.MappingType") ne "HISAT2")
    {
        $fail = $fail."Mapping type is not known!\n";
    }       

    if (! $::global_cfg -> param("required.savemode"))
    {
        $fail = $fail."Mapping type is not provided!\n";
    }
    elsif($::global_cfg -> param("required.savemode") == 1 and
          $::global_cfg -> param("required.savemode") == 2)
    {
        $fail = $fail."Savemode is not known!\n";
    }

    if (! $::global_cfg -> param("required.qscore"))
    {
        $fail = $fail."Mapping type is not provided!\n";
    }
    elsif($::global_cfg -> param("required.qscore") ne "phred33" and
          $::global_cfg -> param("required.qscore") ne "phred64")
    {
        $fail = $fail."Quality scoring type is not known!\n";
    }
    
    

    for (my $j = 1; $j <= 2; $j++)
    {
        my $end_seq = 0;
        my $seq_nr = 1;
        while ($end_seq == 0)
        {
            my $seqname = "sequence".$j."_".$seq_nr."_";
            
            if ($::global_cfg -> param("required.seqtype") eq "PE")
            {
                if (($::global_cfg -> param("sequences.".$seqname."1")) &&
                    ($::global_cfg -> param("sequences.".$seqname."2")))
                {
                    if ($j == 1)
                    {
                        push(@::sequences_1,
                             ($::global_cfg -> param("sequences.".$seqname."1")));
                        push(@::sequences_1,
                             ($::global_cfg -> param("sequences.".$seqname."2")));
                    }
                    else
                    {
                        push(@::sequences_2,
                             ($::global_cfg -> param("sequences.".$seqname."1")));
                        push(@::sequences_2,
                             ($::global_cfg -> param("sequences.".$seqname."2")));                       
                    }

                }
                else
                {
                    if ($j == 1)
                    {
                        $::global_cfg -> param("required.seqlength1",
                            scalar(@::sequences_1));
                    }
                    else
                    {
                        $::global_cfg -> param("required.seqlength2",
                            scalar(@::sequences_1));
                    }
                    
                    $end_seq++; 
                }
            }
            else
            {

                if (($::global_cfg -> param("sequences.".$seqname."1")))
                {
                    if ($j == 1)
                    {
                        push(@::sequences_1,
                             ($::global_cfg -> param("sequences.".$seqname."1")));
                    }
                    else
                    {
                        push(@::sequences_2,
                             ($::global_cfg -> param("sequences.".$seqname."1")));
                    }
                }
                else
                {
                    if ($j == 1)
                    {
                        $::global_cfg -> param("required.seqlength1",
                            scalar(@::sequences_1));
                    }
                    else
                    {
                        $::global_cfg -> param("required.seqlength2",
                            scalar(@::sequences_2));
                    }
                    $end_seq++;
                }
            }
            
            $seq_nr++;
        }

        if ($j == 1)
        {
            for (my $z = 0; $z < scalar(@::sequences_1); $z++)
            {
                my $file = $dir_raw."/Sequences/".($::sequences_1[$z]);
                
                if (!-e $file)
                {
                    $fail = $fail."Can't find ".$file."!\n";
                }
                
            }
        }
        else
        {
            for (my $z = 0; $z < scalar(@::sequences_2); $z++)
            {
                my $file = $dir_raw."/Sequences/".($::sequences_2[$z]);
                
                if (!-e $file)
                {
                    $fail = $fail."Can't find ".$file."!\n";
                }
            }            
        }
    }
    
    
    if ($fail ne "")
    {
        print $fail;
        die "ERROR: INFORMATION is required!\n";
    }
    
    CONFIGURE::CreateProjectFolder();
    
    my $dir_tmp = $::global_cfg -> param("results.tmp");
    
    my $files_check = 0;
    
    while (@tc_files)
    {
        my $tmp_file = shift(@tc_files);
        
        if ($::global_cfg -> param("required.$tmp_file"))
        {
            my $tmp = $::global_cfg -> param("required.$tmp_file");
            my $counf = 0;
            
            given ($tmp_file)
            {
                when ("GenomeFASTA")
                {
                    if (!-e $dir_raw."/FASTA/".$tmp) {$counf += 1;}
                }
                when ("GenomeGTF")
                {
                    if (!-e $dir_raw."/GTF/".$tmp) {$counf += 1;}
                }
                default
                {
                    if (!-e $dir_raw."/Sequences/".$tmp) {$counf += 1;}
                }
            }
            
            if ($counf != 0)
            {
                $fail = $fail."Cant' find  $tmp!\n";
                $files_check += 1;
            }    
        }        
    }
    
    if ($files_check == 0)
    {
        @tc_files = ("GenomeFASTA", "GenomeGTF");
        
        while (@tc_files)
        {
            my $tmp_file = shift(@tc_files);
            my $tmp = $::global_cfg -> param("required.$tmp_file");
            
            if (($tmp_file eq "GenomeGTF" or $tmp_file eq "GenomeFASTA") and
                   -e ($::global_cfg -> param("results.tmp"))."/".(substr($tmp,0,-3)))
            {
                 $::global_cfg -> param("required.$tmp_file",
                    $dir_tmp."/".(substr($tmp,0,-3)));                
            }
            elsif (($tmp_file eq "GenomeGTF" or $tmp_file eq "GenomeFASTA") and
                   -e ($::global_cfg -> param("results.tmp"))."/".$tmp)
            {
                 $::global_cfg -> param("required.$tmp_file",
                    $dir_tmp."/".$tmp);                
            }
            elsif ($tmp_file eq "GenomeFASTA" and
                $tmp =~ ".gz" and
                !-e $dir_tmp."/".(substr($tmp,0,-3)))
            {
                my $cp_cmd = "cp ".$dir_raw."/FASTA/".$tmp." ".$dir_tmp;

                my $cmd = `$cp_cmd`;
                if ($cmd -> success != 1)
                {
                    PIPELINETOOLS::PrintStd($cmd -> stdout);
                    PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                    die "ERROR: Can't copy $tmp to tmp-folder!\n";
                }
                
                my $gz_cmd = "gunzip ".$dir_tmp."/".$tmp;
                $cmd = `$gz_cmd`;
                if ($cmd -> success != 1)
                {
                    PIPELINETOOLS::PrintStd($cmd -> stdout);
                    PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                    die "ERROR: Can't gunzip $tmp in tmp-folder!\n";
                }
                
                $::global_cfg -> param("required.$tmp_file",
                    $dir_tmp."/".(substr($tmp,0,-3)));
                
            }
            elsif ($tmp_file eq "GenomeFASTA" and
                   !-e ($::global_cfg -> param("results.tmp"))."/".$tmp)
            {
                my $cp_cmd = "cp ".$dir_raw."/FASTA/".$tmp." ".$dir_tmp;
                my $cmd = `$cp_cmd`;
                if ($cmd -> success != 1)
                {
                    PIPELINETOOLS::PrintStd($cmd -> stdout);
                    PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                    die "ERROR: Can't copy $tmp to tmp-folder!\n";
                }
                
                $::global_cfg -> param("required.$tmp_file",
                    $dir_tmp."/".$tmp);                
            }
            elsif ($tmp_file eq "GenomeGTF" and
                $tmp =~ ".gz" and
                !-e $dir_tmp."/".(substr($tmp,0,-3)))
            {
                my $cp_cmd = "cp ".$dir_raw."/GTF/".$tmp." ".$dir_tmp;
                my $cmd = `$cp_cmd`;
                if ($cmd -> success != 1)
                {
                    PIPELINETOOLS::PrintStd($cmd -> stdout);
                    PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                    die "ERROR: Can't copy $tmp to tmp-folder!\n";
                }
                
                my $gz_cmd = "gunzip ".$dir_tmp."/".$tmp;
                $cmd = `$gz_cmd`;
                if ($cmd -> success != 1)
                {
                    PIPELINETOOLS::PrintStd($cmd -> stdout);
                    PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                    die "ERROR: Can't gunzip $tmp in tmp-folder!\n";
                }
                
                $::global_cfg -> param("required.$tmp_file",
                    $dir_tmp."/".(substr($tmp,0,-3)));
            }
            elsif ($tmp_file eq "GenomeGTF" and
                   !-e ($::global_cfg -> param("results.tmp"))."/".$tmp)
            {
                my $cp_cmd = "cp ".$dir_raw."/GTF/".$tmp." ".$dir_tmp;
                my $cmd = `$cp_cmd`;
                if ($cmd -> success != 1)
                {
                    PIPELINETOOLS::PrintStd($cmd -> stdout);
                    PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                    die "ERROR: Can't copy $tmp to tmp-folder!\n";
                }
                
                $::global_cfg -> param("required.$tmp_file",
                    $dir_tmp."/".$tmp);                
            }
            elsif ($tmp_file ne "GenomeFASTA" and $tmp_file ne "GenomeGTF" and
                   !-e ($::global_cfg -> param("results.tmp"))."/".$tmp)
            {
                my $cp_cmd = "cp ".$dir_raw."/Sequences/".$tmp." ".$dir_tmp;
                my $cmd = `$cp_cmd`;
                if ($cmd -> success != 1)
                {
                    PIPELINETOOLS::PrintStd($cmd -> stdout);
                    PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                    die "ERROR: Can't copy $tmp to tmp-folder!\n";
                }
                
                $::global_cfg -> param("required.$tmp_file",
                    $dir_tmp."/".$tmp);                 
            }
            
        }
        
        my @init_succfiles = PIPELINETOOLS::ReadDirectory($dir_tmp,".succ");
        
        if (!@init_succfiles)
        {
            my $m = 0;
            while ($m < scalar(@::sequences_1))
            {
                if (!-e ($dir_tmp."/".($::sequences_1[$m])))
                {
                    my $cp_cmd = "cp ".$dir_raw."/Sequences/".($::sequences_1[$m]).
                        " ".$dir_tmp;
                    my $cmd = `$cp_cmd`;
                    if ($cmd -> success != 1)
                    {
                        PIPELINETOOLS::PrintStd($cmd -> stdout);
                        PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                        die "ERROR: Can't copy ".$::sequences_1[$m]." to tmp-folder!\n";
                    }
                }
    
                $::sequences_1[$m] = $dir_tmp."/".$::sequences_1[$m];
                
                $m++;
            }
            
            my $n = 0;
            while ($n < scalar(@::sequences_2))
            {
                #my @test = ();
                my $tmp = $dir_raw."/Sequences/".($::sequences_2[$n])."/*";
                #if (#@test == <$tmp> &&
                if(!-e ($dir_tmp."/".($::sequences_2[$n])))
                {
                    my $cp_cmd = "cp ".$dir_raw."/Sequences/".($::sequences_2[$n]).
                        " ".$dir_tmp;
                    my $cmd = `$cp_cmd`;
                    if ($cmd -> success != 1)
                    {
                        PIPELINETOOLS::PrintStd($cmd -> stdout);
                        PIPELINETOOLS::WriteLogOut(($cmd -> merged));
                        die "ERROR: Can't copy ".$::sequences_2[$n]."to tmp-folder!\n";
                    }
                }
    
                $::sequences_2[$n] = $dir_tmp."/".$::sequences_2[$n];
                
                $n++;
            }            
        }
    }
    
    return;
}

#-------------------------------------------------------------------------------
# CheckUsrConfig: This function creates the necassary folders for the project.
#-------------------------------------------------------------------------------
sub CreateProjectFolder{
    #Creates the folder architecture for the project
    my $project_name = $::global_cfg -> param("required.ProjectName");
    my $user_name = $::global_cfg -> param("required.User");
    my $dir_results = $::global_cfg -> param("sdirectories.Results");
    
    #create user folder
    $dir_results = $dir_results."/".$user_name;
    CONFIGURE::CreateFolder($dir_results);
    
    #create project folder
    my $dir_prj = $dir_results."/".$project_name;
    
    if ($::global_cfg -> param("cmdline.overwrite"))
    {
        my $uinput = "";
        print "Do you like to overwrite folder: $dir_prj !\n
            yes(y) or no(n): ";
        
        
        
        while ($uinput ne "yes" && $uinput ne "no" &&
               $uinput ne "y" && $uinput ne "n" &&
               $uinput ne "Yes" && $uinput ne "No" &&
               $uinput ne "Y" && $uinput ne "N" &&
               $uinput ne "YES" && $uinput ne "NO")
        {
            chomp($uinput = <STDIN>);
        }
        
        if ($uinput eq "yes" || $uinput eq "YES" ||
            $uinput eq "Y" || $uinput eq "Yes" ||
            $uinput eq "y")
        {
            my $delete_cmd = "rm -R -f ".$dir_prj;
            
            PIPELINETOOLS::PrintStd($delete_cmd);
            
            my $cmd = `$delete_cmd`;
            if ($cmd -> success != 1)
            {
                PIPELINETOOLS::PrintStd($cmd -> merged);
                PIPELINETOOLS::WriteLogOut($cmd -> merged);
                die "ERROR: Can't overwrite $dir_results!\n";
            }
            PIPELINETOOLS::PrintStd($cmd -> stdout);
            PIPELINETOOLS::WriteLogOut($cmd -> stdout)
        }     
    }
    
    CONFIGURE::CreateFolder($dir_prj);
    
    #create subfolders for the projet
    CONFIGURE::CreateFolder($dir_prj."/tmp");
    CONFIGURE::CreateFolder($dir_prj."/Results");
    CONFIGURE::CreateFolder($dir_prj."/QualityReports");
    CONFIGURE::CreateFolder($dir_prj."/Preprocessing");
    CONFIGURE::CreateFolder($dir_prj."/Graphics");
    CONFIGURE::CreateFolder($dir_prj."/Alignment");
    
    #adding directories to the global variable
    $::global_cfg -> param("results.std", $dir_prj);
    $::global_cfg -> param("results.tmp", $dir_prj."/tmp"); 
    $::global_cfg -> param("results.results", $dir_prj."/Results");
    $::global_cfg -> param("results.qreports", $dir_prj."/QualityReports");
    $::global_cfg -> param("results.preprocess", $dir_prj."/Preprocessing");
    $::global_cfg -> param("results.graphics", $dir_prj."/Graphics");
    $::global_cfg -> param("results.alignment", $dir_prj."/Alignment");    
    
    return;
}

#-------------------------------------------------------------------------------
# CreateFolder: This function creates a folder.
#               Input: folder path
#-------------------------------------------------------------------------------
sub CreateFolder{
	use File::Path qw(make_path);

	my $folder_name = $_[0];
	
	if (!-d $folder_name)
	{
        make_path($folder_name) or
            die "ERROR! Can't create folders for Results!\n";	
	}
	
    return;
}

#-------------------------------------------------------------------------------
# CheckUsrConfig: This function prints the help message.
#-------------------------------------------------------------------------------
sub PrintHelp{
	#subroutine for printing
	my $line = ("#" x 80)."\n";
    print "Help: \n";
    print $line;
    CONFIGURE::PrintUsage();
    print $line."\n";
    print "--threads <n>: This can be used to set a higher number of threads\n";
    print "         used by the pipeline.";
    print "--overwrite: If overwrite is set to TRUE, it  will start the,\n";
    print "           pipeline always at the beginning.\n";
    print "           That means if there is already a folder with the same\n";
    print "           name it will be deleted and after that the pipeline\n";
    print "           starts.\n";
    print "--silent: This option supresses the STDOUT, and only prints the\n";
    print "        crucial output.\n";
    print "--help: Prints this message!\n\n";
    print "Created by Raphael Hablesreiter\n";
    
    return;
}

#-------------------------------------------------------------------------------
# CheckUsrConfig: This function prints the usage message.
#-------------------------------------------------------------------------------
sub PrintUsage{
    
    print "Usage: perl Main.pl --config <user configuration file> [options]\n";
    print "    Options:\n";
    print "        --threads n (n = number of threads, default: n = 1)\n";
    print "        --overwrite (default: FALSE)\n";
    print "        --silent\n";
    print "        --help\n";
    
    return;
}

1;