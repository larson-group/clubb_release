#!/usr/bin/perl

###########################################################################
# Plotgen v3.0
#
#  
###########################################################################

package Plotgen;

use strict;
use CaseReader;
use Cwd;

my $plotgenDir = cwd;

# Loop through each .case file so the case can be plotted
my @cases = <cases/*>;
foreach my $file (@cases) 
{
    # Read the case file. If there is an error, exit.
    if (my $err = CaseReader::ReadCase($file))
    {
        print(STDERR $err, "\n");
        exit(1);
    }

    print("Successfully opened: " . $file . "\n");
    print("Plotting case: " . $CASE::CASE{'name'} . "\n");
    print("nightly_output: " . $CASE::CASE{'test'}[2] . "\n");
}
