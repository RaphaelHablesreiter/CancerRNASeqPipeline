#------------------------------------------------------------------------------
# CancerRNASeq - a pipeline for transcriptomics datasets
# Author: Raphael Hablesreiter
#
# CancerRNASeq.pl
#          
# Description: Use this to call CancerRNASeq.
#
#-----------------------------------------------------------------------------


use File::Find;
use strict;
use warnings;
use Cwd;
use File::Path qw(make_path);
use PREPROCESSING;
use CONFIGURE;
use DRIVER;
use Config::Simple;

our $global_cfg = new Config::Simple(syntax => 'ini');
our @sequences_1 = (); 
our @sequences_2 = ();

PIPELINETOOLS::PrintTask("Starting Pipeline");
    
if(CONFIGURE::Main() == 0)
{
    
    DRIVER::Main();
    PIPELINETOOLS::PrintTask("Finished Pipeline");
}

exit; 