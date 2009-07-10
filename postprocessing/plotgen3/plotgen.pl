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
use File::Basename;

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
	print("Input Folders: @inputDirs\n");
	print("Output Folder: $output \n");

	my $plotNumber = 0;

	# Loop through each .case file so the case can be plotted
	my @cases = <cases/*.cas*>;
	foreach my $file (@cases) 
	{
		# Read the case file. If there is an error, exit.
		if (my $err = CaseReader->readCase($file))
		{
	    		print(STDERR $err, "\n");
	    		exit(1);
	    	}
	
		if(dataExists($CASE::CASE{'name'}))
		{	
			# Print the case title to the HTML page
			OutputWriter->writeCaseTitle($outputIndex, $CASE::CASE{'headerText'});
	
			# Print any additional text/html specified
			if($nightly == 1) # If in nightly mode
			{
				my $nightlySubText = $CASE::CASE{'nightlyOutput'}{'subText'};
				my $nightlySubHtml = $CASE::CASE{'nightlyOutput'}{'subHtml'};
				
				# Check to see if there was any additional text specified. If there was,
				# write it to the HTML file.
	
				if($nightlySubText)
				{
					OutputWriter->writeSubHeader($outputIndex, $nightlySubText);
				}
				
				if($nightlySubHtml)
				{
					OutputWriter->writeSubHtml($outputIndex, $nightlySubHtml);
				}
			}
			else # If not in nightly mode
			{
				my $subText = $CASE::CASE{'additionalOutput'}{'subText'};
				my $subHtml = $CASE::CASE{'additionalOutput'}{'subHtml'};
				
				if($subText)
				{
					OutputWriter->writeSubHeader($outputIndex, $subText);
				}
				
				if($subHtml)
				{
					OutputWriter->writeSubHtml($outputIndex, $subHtml);
				}
	
			}
	
			# Parse Case information
			my $caseName =  $CASE::CASE{'name'};
			my $plotTitle =  $CASE::CASE{'plotTitle'};
			my $startTime =  $CASE::CASE{'startTime'};
			my $endTime =  $CASE::CASE{'endTime'};
			my $startHeight =  $CASE::CASE{'startHeight'};
			my $endHeight =  $CASE::CASE{'endHeight'};
			my $units =  $CASE::CASE{'units'};
			my $type = $CASE::CASE{'type'};
			my $tickCount = 12344324234532;

			my $matlabArgs = $caseName . " " . $plotTitle . " " . $plotNumber . " " . $type . " " . $startTime . " " . $endTime . " " . $startHeight . " " . $endHeight . " " . $units . " " . $tickCount . " ";
		
			# Parse variables from .case file
			for(my $count = 0; $count < $CASE::CASE{'numVars'}; $count++)
			{
				my $file = $CASE::CASE{'variables'}[$count]{'file'};
				my $name = $CASE::CASE{'variables'}[$count]{'name'};
				my $expression = $CASE::CASE{'variables'}[$count]{'expression'};
				my $title = $CASE::CASE{'variables'}[$count]{'title'};
				my $unit = $CASE::CASE{'variables'}[$count]{'units'};
				my $lineWidth = $CASE::CASE{'variables'}[$count]{'lineWidth'};
				my $lineType = $CASE::CASE{'variables'}[$count]{'lineType'};
				my $lineColor = $CASE::CASE{'variables'}[$count]{'lineColor'};

				$matlabArgs = "$matlabArgs $file $name $expression $title $unit $lineWidth $lineType $lineColor";
			}

			$plotNumber++;

			print("Matlab args: $matlabArgs \n");
		}
		else
		{
			print("Not plotting $CASE::CASE{'name'}\n");
		}
	
		# Convert the eps files to jpq
		convertEps();

		# Add image file to HTML page
		OutputWriter->placeImage($outputIndex, "jpg/cloud_feedback_s6_page1.eps.jpg");
	}
}

sub convertEps()
{
	mkdir "$output/jpg";
	my @epsFiles = <$output/eps/*eps>;
	foreach my $eps (@epsFiles)
	{
		my $filename = basename($eps);
		system("convert -density 90 $eps $output/jpg/$filename.jpg");
	}
}

###############################################################################
# Checks all input directories to see if one of them contains the current case.
# Will return true if at least on input folder contains data, otherwise false.
# #############################################################################
sub dataExists()
{
	my $dataFile = shift(@_);
	my $retValue = 0;
	
	foreach my $currInput (@inputDirs)
	{
		my @files = <$currInput/$dataFile*>;
		if(@files)
		{
			$retValue = 1;
		}
	}

	return $retValue;
}

###############################################################################
# Reads any arguments passed in.
###############################################################################
sub readArgs()
{
	my $numArgs = $#ARGV + 1;

	if($numArgs == 0)
	{
		printHelp();
	}

	my %option = ();
	getopts("rlbdan?", \%option);

	if ($option{r})
	{
		$overwrite = 1;
	}

	if ($option{l})
	{
		$plotLes = 1;
	}

	if ($option{b})
	{
		$plotBest = 1;
	}

	if ($option{d})
	{
		$plotDec = 1;
	}

	if ($option{a})
	{
		$plotLes = 1;
		$plotBest = 1;
		$plotDec = 1;
	}

	if ($option{n})
	{
		$nightly = 1;
	}

	if ($option{h})
	{
		printHelp();
	}

	my $currentCount = 0;

	# Reset the number of arguments
	$numArgs = $#ARGV + 1;

	# Parse any additional arguments
	foreach (@ARGV)
	{
		# If the argument does not start with '-' and if $output was not set
		if(!$output)
		{
			my $currentDir = abs_path($_);

			if((($numArgs - 1) - $currentCount) > 0)
			{
				if(-d $currentDir)
				{
					push(@inputDirs, $currentDir);
				}
				else
				{
					print("The input folder: $currentDir does not exist.\n");
					exit(1);
				}
			}
			else
			{
				if($_ =~ m/\/$/)
				{
					$currentDir = abs_path(substr($_, 0, length($_) - 1));
				}

				$output = $currentDir;
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
	print("  -r\tIf the output folder already exists, replace the contents\n");	
	print("  -l\tPlot LES data for comparison.\n");	
	print("  -b\tPlot Best Ever data for comparison.\n");	
	print("  -d\tPlot December data for comparison.\n");	
	print("  -a\tSame as -lbd. Plots LES, Best Ever, and December data for comparison.\n");
	print("  -n\tRuns in nightly mode.\n");
	print("  -h\tPrints this help message.\n");
	exit(0);
}
