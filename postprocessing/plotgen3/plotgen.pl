#!/usr/bin/perl

###########################################################################
# Plotgen v3.0
#
#  
###########################################################################

package Plotgen;

use strict;
use CaseReader;
use OutputWriter;
use Cwd;
use Switch;

# Argument list

# If running in nightly mode, this value should be set to 1
my $nightly = 0;

# Specifies to overwrite a directory (Default: false)
my $overwrite = 0;
my $dir1 = "";
my $dir2 = "";
my $output = "";
my $plotLes = 0;
my $plotGolaz = 0;
my $plotHOC = 0;

my $plotgenDir = cwd;

readArgs();

my $outputIndex = $output . "/index.html";

OutputWriter->writeHeader($outputIndex);
runCases();
OutputWriter->writeFooter($outputIndex);

###############################################################################
# Runs all of the .case file in the cases folder.
###############################################################################
sub runCases()
{
	print("Dir 1: " . $dir1 . "\n");
	print("Dir 2: " . $dir2 . "\n");
	print("Output: " . $output . "\n");

	# Loop through each .case file so the case can be plotted
	my @cases = <cases/*>;
	foreach my $file (@cases) 
	{
		print("File: " . $file . "\n");
		# Read the case file. If there is an error, exit.
		if (my $err = CaseReader->readCase($file))
		{
	    		print(STDERR $err, "\n");
	    		exit(1);
	    	}
		
		# Print the case title to the HTML page
		OutputWriter->writeCaseTitle($outputIndex, $CASE::CASE{'headerText'});
	}
}

###############################################################################
# Reads any arguments passed in.
###############################################################################
sub readArgs()
{
	my $numArgs = $#ARGV + 1;

	# This variable is used for parsing out the folder path variables
	my $pos = 0;

	if($numArgs == 0)
	{
		printHelp();
	}

	foreach my $argnum (0 .. $#ARGV) 
	{
		switch(@ARGV[$argnum])
		{
			case "--help"
			{
				printHelp();
			}
			case "--nightly"
			{
				$nightly = 1;
			}
			case "-r"
			{
				$overwrite = 1;
			}
			case "--replace"
			{
				$overwrite = 1;
			}
			else
			{
				$pos++;

				if($pos == 1)
				{
					$dir1 = @ARGV[$argnum];
				}
				elsif($pos == 2)
				{
					$output = @ARGV[$argnum];
				}
				elsif($pos == 3)
				{
					$dir2 = $output;
					$output = @ARGV[$argnum];
				}
			}
		}
	}
	
	if($pos < 2)
	{
		printHelp();
	}
}

###############################################################################
# Prints the help message
###############################################################################
sub printHelp()
{
	print("Plotgen" . "\n");
	print("Usage: plotgen [ options ... ]" . "\n");
	exit(0);
}
