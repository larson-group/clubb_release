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
	my $argCount = 0;

	# Parse any additional arguments
	foreach my $argnum (0 .. $numArgs) 
	{
		# If the argument does not start with '-' and if $output was not set
		if(!(@ARGV[$argnum] =~ m/^-/) && !$output)
		{
			my $count = $currentCount + $argCount;
			print("Count: " . $count . "\n");
			if(($numArgs - 1) - ($currentCount + $argCount) > 1)
			{
				push(@inputDirs, @ARGV[$argnum]);
			}
			else
			{
				$output = @ARGV[$argnum];
			}
			
			$currentCount++;
		}
		else
		{
			$argCount += length(@ARGV[$argnum]) - 1;
		}
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
	print("  -a\t\t\n");	
	print("  -h\t\t\n");
	exit(0);
}
