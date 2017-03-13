#-------------------------------------------------------------------------------
# PIPELINETOOLS.pm
#
# Author: Raphael Hablesreiter (1130649)
#          
# Description: 
# The PIPELINETOOLS - module consists of all the essential function, which 
# are used more often in the entire pipeline.
#
#------------------------------------------------------------------------------
#

package PIPELINETOOLS; 

use strict;
use warnings;

use Backticks;
use File::chdir;
use File::Copy;

#-------------------------------------------------------------------------------
# WriteToLogFile: This function takes the input as the text, that has to be
#                 written to the log.file. If the log.file doesn't exist it
#                 creates the log.file in the project folder.
#                 Input: text for log.file [string]
#-------------------------------------------------------------------------------

sub WriteToLogFile {
	
	my $log_name = ($::global_cfg -> param("results.std"))."/log.file";
	my $log_text = $_[0];
	
    #creates the log.file, attaches the header and adds the first log entry
	if ($log_text ne "" && !(-e $log_name))
	{
	    open(LOGFILE, ">", $log_name)
			or die "ERROR: Can't create log-file!\n";	
		print LOGFILE CreateLogHeader();
		print LOGFILE $log_text;
		close LOGFILE;
	}
    #opens the log-file for append and attaches the log entry
	elsif ($log_text ne "" && -e $log_name)
	{
        my $log_number = $::global_cfg -> param("tmp.lognumber");
		
		open(LOGFILE, ">>", $log_name)
			or die "ERROR: Can't open log-file!\n";	
		print LOGFILE $log_text;		
		close LOGFILE;	
	}	
	else
	{
		die "ERROR: Problem with the log-file!\n";
	}
	
	return;
}


#-------------------------------------------------------------------------------
# WriteToLogChapter: This function takes the input as the text,
#                    that has to be written to the log.file as a chapter.
#             Input: text for log.file chapter [string]
#-------------------------------------------------------------------------------
sub WriteLogChapter{
	my $chapter = $_[0];
	my $space1 = int((78 - length($chapter))/2);
	my $time = "Started at ".(PIPELINETOOLS::GetTime(1));
	my $space2 = int((78 - length($time))/2);
	
	
	my $chapter_txt = ("#" x 80)."\n";
	$chapter_txt = $chapter_txt."#".(" " x 78)."#\n";
	
	$chapter_txt = $chapter_txt."#".(" " x $space1).$chapter;
	$chapter_txt = $chapter_txt.(" " x (78 - length($chapter) - $space1))."#\n";
	
	$chapter_txt = $chapter_txt."#".(" " x $space2).$time;
	$chapter_txt = $chapter_txt.(" " x (78 - length($time) - $space2))."#\n";	
	$chapter_txt = $chapter_txt."#".(" " x 78)."#\n";
	$chapter_txt = $chapter_txt.("#" x 80)."\n";
	
	PIPELINETOOLS::WriteToLogFile($chapter_txt);
	
	return;
}

#-------------------------------------------------------------------------------
# WriteLogOut: This function writes the std::Out of a function into the log.file
#              and finish this with line that provides the finishing time.
#       Input: std::Out of a function for log.file [string]
#-------------------------------------------------------------------------------
sub WriteLogOut{
	
	my $p_out = $_[0]."\n";
	my $ad = "*" x 27;
	$p_out = $p_out.$ad.
		"Task finished at ".(PIPELINETOOLS::GetTime(2)).$ad."*\n\n\n";
	
	PIPELINETOOLS::WriteToLogFile($p_out);
	
	return;
}

#-------------------------------------------------------------------------------
# CreateLogHeader: This function writes the header into the log.file.
#-------------------------------------------------------------------------------
sub CreateLogHeader {
    #This subprogramm is creating a header for the log-file
	
	
	my $author = "Raphael Hablesreiter";
	my $version = "V0.1";
	my $institute =
		"Graz University of Technology - Institute of Molecular Biotechnology";
	
	my $seqtype;
	if (($::global_cfg -> param("required.sequence_1")) eq "PE")
	{
        $seqtype = "Paired End Reads";
    }
	else
	{
		$seqtype = "Single End Reads";
	}
    
	
	my $header = ("#" x 80)."\n";
	$header = $header."Author: $author\n"
		."Institute: $institute\n"
		."Version: $version\n"
		."Created: SS16\n"
		.("#" x 80)."\n"
		."Run Information:\n\n"
		."    User: ".($::global_cfg -> param("required.User"))."\n"
		."    Projectname: ".($::global_cfg -> param("required.ProjectName"))."\n\n"
		."    1. Sequence Name: ".($::global_cfg -> param("required.sequence_1"))."\n"
		."    2. Sequence Name: ".($::global_cfg -> param("required.sequence_2"))."\n"
		."    Sequences Type: ".$seqtype."\n\n"
		."    Mapping Type: ".($::global_cfg -> param("required.MappingType"))."\n"
		.("#" x 80)."\n\n";
	
	return $header;
}

#-------------------------------------------------------------------------------
# GetTime: This function provides two types of time and date.
#          Input = 1 -> Return: $mday.$mon.$year $hour:$min:$sec
#          Input = 2 -> Return: $hour:$min:$sec
#-------------------------------------------------------------------------------
sub GetTime{
    
    #FROM: Perl Cookbook S.71
    #time - returns the number of seconds that have psased sunce the last Epoch
    #localtime - converts the output of time into distinct values
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	$year += 1900;
    #the range for months and days are starting from 0
	$mon += 1;
    my $time;
    
    $sec = sprintf('%02s',$sec);
    $hour = sprintf('%02s',$hour);
    $min = sprintf('%02s',$min);
    $mday = sprintf('%02s',$mday);
	$mon = sprintf('%02s',$mon);
	
    if ($_[0] == 1)
    {
        $time = "$mday.$mon.$year $hour:$min:$sec";    
    }
    elsif ($_[0] == 2)
    {
        $time = "$hour:$min:$sec";  
    }
    else
    {
        $time = "$hour:$min:$sec";
    }
    
    return $time;
    
}

#-------------------------------------------------------------------------------
# ReadDirectory: This function provides an array with all the files in a folder.
#                Input: $_[0] = directory
#                       $_[1] = string, that should be in the filename,
#                               if $_[1] = "", all files of a folder
#                               are returned.
#                Return: Array with all the files in a folder that matches a
#                        certain criteria.
#-------------------------------------------------------------------------------
sub ReadDirectory{
	
    my $directory = $_[0]; #directory
	my $file_name;
    my @output = ();   
    
    my $name_containing = "";
    if (scalar(@_) > 1)
    {
        $name_containing = $_[1]; #name, that has to be in the files
    }
    

    opendir(my $dir, $directory) or die "Cannot open $directory: $!";
    while (defined ($file_name = readdir($dir)))
    {
        next if $file_name =~ /^\./; #skips files with .. and .
        next if $file_name !~ /\./;
        
        if ($name_containing ne "")
        {
            next if $file_name !~ /$name_containing/;
        }
        
        push @output, $file_name;
        
    }
    closedir $dir;
    
    if (@output)
	{
        return @output;
    }
	else
	{
		@output = ();
		return @output;
	}
    
	
    
    
}
#-------------------------------------------------------------------------------
# ChangeWorkingDir: This function changes the working directory.
#                   Input: $_[0] = directory
#-------------------------------------------------------------------------------
sub ChangeWorkingDir{
	
	my $path = $_[0];
	$CWD = $path;
	
	return;	
}

#-------------------------------------------------------------------------------
# MoveFile: This function changes the working directory.
#           Input: $_[0] = sourcefile
#                  $_[1] = destinationfolder/-file
#-------------------------------------------------------------------------------
sub MoveFile{
	my $sourcefile = $_[0];
	my $destinationfile = $_[1];
	
	my $mv_cmd = "mv ".$sourcefile." ".$destinationfile;
	my $cmd = `$mv_cmd`;
	if ($cmd -> success != 1)
	{
		print "ERROR: Moving $sourcefile to $destinationfile failed!\n";
	}	
    
	return;
}

#-------------------------------------------------------------------------------
# CreateSuccess: This function creates a successfile (part###.succ) and this
#                successfile is in the style of an INI-file.
#                see Config::Simple
#
#                Input: $_[0] = name of part of pipeline finished
#                               "Cutadapt" -> for Cutadapt has finished
#                               "Trimmomatic" -> for Trimmomatic has finished
#                               "rCorrector" -> for rCorrector has finished
#                               "Mapping" -> for Mapping has finished
#                               "Cufflinks" -> for Cufflinks has finished
#                               "Cuffmerge" -> for Cuffmerge has finished
#                               "Cuffdiff" -> for Cuffdiff has finished
#                               "CummeRbund" -> for CummeRbund has finished
#
#   Layout sample:
#   ;Config::Simple 4.58
#   ;Wed Apr 20 09:28:07 2016
#
#   [success]
#   done=okay
#   cmd=Cutadapt
#   time=20.04.2016 09:28:07
#-------------------------------------------------------------------------------
sub CreateSuccess{
    #call with the wanted command-parameter

    my $prgm_cmd = $_[0];	
	
	my $succ_name;
	my $i = 0;
    my $time = PIPELINETOOLS::GetTime(1);
	my $succ_dir = $::global_cfg -> param("results.tmp");
    
	my $succ = new Config::Simple(syntax=>'ini');
	$succ->param("success.cmd","$prgm_cmd");
	$succ->param("success.time","$time");
	$succ->param("success.done","okay");
				 
	my @succ_files = PIPELINETOOLS::ReadDirectory($succ_dir,".succ");
	
	if (@succ_files)
	{
        #getting the number(i) of the last successfile and creating new
        #successfilename(i+1)
		@succ_files = sort @succ_files;
        my $last_succ = pop(@succ_files);
		$last_succ = substr($last_succ,4,-5);
		print "Last Success = $last_succ\n";
		$last_succ = int($last_succ) + 1;
		
		$succ_name = $succ_dir."/part".sprintf("%03s",$last_succ).".succ";
    }
    else
	{
		$succ_name = "part001.succ";	
	}
	
	$succ -> write($succ_name) or die "ERROR";
	
	
	return 0;
}

#-------------------------------------------------------------------------------
# MoveFile: This prints to Std::out, if pipeline parameter silent is not TRUE.
#           Input: $_[0] = std::out of a cmd program.
#-------------------------------------------------------------------------------
sub PrintStd{
	my $stdoutput = $_[0];
	$stdoutput = $stdoutput."\n";
	
	if (!($::global_cfg -> param("cmdline.silent")))
	{
        print $stdoutput;
    }
	
	return;
}
#-------------------------------------------------------------------------------
# MoveFile: This prints to Std::out, if pipeline parameter silent is not TRUE.
#           Input: $_[0] = pipeline part
#-------------------------------------------------------------------------------
sub PrintTask{
	my $stdoutput = $_[0];
	$stdoutput = $stdoutput." @ ".(PIPELINETOOLS::GetTime(1))."\n";
	
	print $stdoutput;
	
	return;
}

1;