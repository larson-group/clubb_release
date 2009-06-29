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
use Cwd 'abs_path';
use Switch;
use Getopt::Std;

# Argument list

# If running in nightly mode, this value should be set to 1
my $nightly = 0;

# Specifies to overwrite a directory (Default: false)
my $overwrite = 0;
my @inputDirs;
my $output = "";
my $plotLes = 0;
my $plotBest = 0;
my $plotDec = 0;

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
	print("Input Dirs: @inputDirs\n");
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
	my $argsSet = 0;

	if($numArgs == 0)
	{
		printHelp();
	}

	my %option = ();
	getopts("rlbdan?", \%option);

	if ($option{r})
	{
		$overwrite = 1;
		$argsSet ++;
	}

	if ($option{l})
	{
		$plotLes = 1;
		$argsSet ++;
	}

	if ($option{b})
	{
		$plotBest = 1;
		$argsSet ++;
	}

	if ($option{d})
	{
		$plotDec = 1;
		$argsSet ++;
	}

	if ($option{a})
	{
		$plotLes = 1;
		$plotBest = 1;
		$plotDec = 1;
		$argsSet ++;
	}

	if ($option{n})
	{
		$nightly = 1;
		$argsSet ++;
	}

	if ($option{h})
	{
		printHelp();
	}

	my $currentCount = 0;

	# Reset the number of arguments
	$numArgs = $#ARGV + 1;

	# Parse any additional arguments
	foreach my $argnum (0 .. $numArgs) 
	{
		print(@ARGV[$argnum] . "\n");
		# If the argument does not start with '-' and if $output was not set
		if(!(@ARGV[$argnum] =~ m/^-/) && !$output)
		{
			if((($numArgs - 1) - $currentCount) > 0)
			{
				my $inputDirToAdd = abs_path(@ARGV[$argnum]);

				if(-d $inputDirToAdd)
				{
					push(@inputDirs, $inputDirToAdd);
				}
				else
				{
					print("The input folder: $inputDirToAdd does not exist.\n");
					exit(1);
				}
			}
			else
			{
				$output = abs_path(@ARGV[$argnum]);
			}
			
			$currentCount++;
		}
	}

	# Finally, check to see if the output folder exists. If it does, and
	# '-r' was not passed in, exit. Otherwise, create it.
	if(-d $output && $overwrite == 0)
	{
		print("Output folder already exists. To overwrite, use the -r option.\n");
		exit(1);
	}
	else
	{
		mkdir $output unless -d "$output";
	}
}

###############################################################################
# Prints the help message
###############################################################################
sub printHelp()
{
	print("Usage: plotgen [OPTION]... INPUT... OUTPUT\n");
	print("  -r\t\tIf the output folder already exists, replace the contents\n");	
	print("  -l\t\tPlot LES data for comparison.\n");	
	print("  -b\t\tPlot Best Ever data for comparison.\n");	
	print("  -d\t\tPlot December data for comparison.\n");	
	print("  -a\t\tSame as -lbd. Plots LES, Best Ever, and December data for comparison.\n");
	print("  -n\t\tRuns in nightly mode.\n");
	print("  -h\t\tPrints this help message.\n");
	exit(0);
}
