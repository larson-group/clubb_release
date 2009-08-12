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
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Path;

# The location of MATLAB. External users probably will not need the 
# sudo -u matlabuser part.
my $MATLAB = "sudo -u matlabuser /usr/local/bin/matlab -nodisplay -nodesktop";

# Plotgen Version Number
my $VERSION = 3.0;

# Used to create a "random" output directory so multiple runs
# don't overwrite each other.
my $randInt = int(rand(999999999999999)) + 1;

# Image Conversion Settings
my $DPI = 300;
my $QUALITY = 100;

# Argument list

# If running in nightly mode, this value should be set to 1
my $nightly = 0;

# Specifies to overwrite a directory (Default: false)
my $overwrite = 0;
my @inputDirs;
my $output = "";
my $outputTemp = "";
my $plotLes = 0;
my $plotBest = 0;
my $plotDec = 0;

# Custom Color Definitions for "CLUBB_current" and "CLUBB_previous"
my $lt_blue = "[ 0.00, 0.63, 1.00 ]";
my $orange = "[ 0.94, 0.50, 0.16 ]";
# Arrays to cycle through when auto is set for lines
my @lineStyles = ("--", "-", "-.");
my @lineColors = ($orange, $lt_blue, "green", "red", "blue", "cyan", "yellow", "magenta");
my @lineWidths = (4, 3, 2.5, 2, 1.5, 1);

# Counters for automatic lines
my $lineStyleCounter = 0;
my $lineColorCounter = 0;
my $lineWidthCounter = 0;

my $outputIndex = "";

my $plotCount = 0;

# Do not keep permissions on file copies
$File::Copy::Recursive::KeepMode = 0;

# This is needed so we can restore the display after plotgen runs.
# We need to clear this environment variable before running, otherwise
# matlab will fail.
my $sessionType = $ENV{'DISPLAY'}; 
$SIG{INT} = "cleanup";
$SIG{CHLD} = "IGNORE";

main();

###############################################################################
# Plotgen 3 Main subroutine
###############################################################################
sub main()
{
	readArgs();

	# Ensure that Matlab can write to the temp output folder
	my $mode = 0777; chmod $mode, $outputTemp;

	print("Input Folders: @inputDirs\n");
	print("Output Folder: $output \n");

	$outputIndex = $outputTemp . "/index.html";

	# Fork to make MATLAB faster
	my $pid = fork();

	if($pid == 0) # Child
	{
		$ENV{'DISPLAY'} = '';

		system("mkfifo matlab_pipe");
		system("$MATLAB <> matlab_pipe");

		# Convert the eps files to jpq
		convertEps();

		OutputWriter->writeFooter($outputIndex);
	
		cleanup();
	
		print("Done! To display the plots, open: \n$output/index.html \nin your web browser.\n\n");
		print("Hit [Enter] to continue.");

		exit(0);
	}
	else # Parent
	{
		OutputWriter->writeHeader($outputIndex);
		runCases();

		# Quit matlab
		system("echo quit > matlab_pipe");
		system("rm matlab_pipe");

		print("\n");
		exit(0);
	}
}

###############################################################################
# Runs all of the .case file in the cases folder.
###############################################################################
sub runCases()
{
	# Counter used to place images
	my $count = 0;

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

		if(dataExists($CASE::CASE{'name'}) && ($CASE::CASE{'enabled'} ne 'false'))
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
			
			# Check to see if this is a budget plot or standard plot
			if($CASE::CASE{'type'} eq "budget")
			{
				buildMatlabStringBudget($CASE::CASE, $count);
				
				# Convert the eps files to jpq
				convertEps($CASE::CASE{'name'} . "_" . $count . "_budget");

				# Add image file to HTML page
				placeImages($CASE::CASE{'name'} . "_" . $count . "_budget");
			}
			else
			{
				buildMatlabStringStd($CASE::CASE);
				
				# Add image file to HTML page
				placeImages($CASE::CASE{'name'}, $plotCount);
			}

			$count++;
		}
		else
		{
			#print("Not plotting $CASE::CASE{'headerText'}\n");
		}
	}
}

###############################################################################
# Converts all EPS files to JPG file
###############################################################################
sub convertEps()
{
	print("\nConverting images...\n");
	mkdir "$outputTemp/jpg" unless -d "$outputTemp/jpg";
	my @epsFiles = <$outputTemp/*eps>;
	foreach my $eps (@epsFiles)
	{
		my $filename = basename($eps);
		system("convert -density $DPI -quality $QUALITY -colorspace RGB $eps $outputTemp/jpg/$filename.jpg");

		unlink($eps);
	}
}

sub placeImages()
{
	my $caseName = shift(@_);
	my $numImages = shift(@_);

	OutputWriter->printDivCenter($outputIndex);

	for(my $x = 0; $x < $numImages; $x++)
	{
		OutputWriter->placeImage($outputIndex, "jpg/$caseName" . "_" . "$x.eps.jpg");
	}

	OutputWriter->printCloseDivCenter($outputIndex);
}

###############################################################################
# Generates argument lists needed for the Matlab part of Plotgen. This is for
# budget plots that will plot each input directory on a seperate plot.
###############################################################################
sub buildMatlabStringBudget()
{
	my $CASE = shift(@_);
	my $budgetNum = shift(@_);
	
	my $totPlotNum = 0;
	
	# Get Common Case information
	my $caseName =  $CASE::CASE{'name'} . "_" . $budgetNum . "_budget";
	my $startTime =  $CASE::CASE{'startTime'};
	my $endTime =  $CASE::CASE{'endTime'};
	my $startHeight =  $CASE::CASE{'startHeight'};
	my $endHeight =  $CASE::CASE{'endHeight'};
	
	# Get array of plots from case file
	my @plots;
	push(@plots, @{$CASE::CASE{'plots'}});

	my $units = $plots[0]{'axisLabel'};
	my $type = $plots[0]{'type'};
	
	# Loop through each input folder and create a plot for each folder
	for(my $plotCount = 0; $plotCount < @inputDirs; $plotCount++)
	{
		my $plotTitle = basename($inputDirs[$plotCount]);
		
		# Loop through each plot in the Case file
		for(my $plotNum = 0; $plotNum < @plots; $plotNum++)
		{
			# Counters for automatic lines
			$lineStyleCounter = 0;
			$lineColorCounter = 0;
			$lineWidthCounter = 0;
			my $matlabArgs = "\'$caseName\', \'$plotTitle\', $totPlotNum, \'$type\', $startTime, $endTime, $startHeight, $endHeight, \'$units\', $randInt";
			my $tempMatlabArgs = $matlabArgs;
		
			my @lines;
			push(@lines, @{$plots[$plotNum]{'lines'}});
			
			for(my $lineNum = 0; $lineNum < @lines; $lineNum++)
			{
				my $file = "$inputDirs[$plotCount]/$lines[$lineNum]{'filename'}";
				
				my $name = $lines[$lineNum]{'name'};
				my $expression = $lines[$lineNum]{'expression'};
				my $type = $lines[$lineNum]{'type'};

				if(-e $file)
				{
					my $lineWidth = $lineWidths[$lineWidthCounter];
					my $lineStyle = $lineStyles[$lineStyleCounter];
					my $lineColor = $lineColors[$lineColorCounter];
							
					$matlabArgs = "$matlabArgs, \'$file\', \'$expression\', \'$name\', $lineWidth, \'$lineStyle\', \'$lineColor\'";		

					incrementLineTypes();
				}
			}
			
			# Check to see if there are lines to be plotted:
			if($matlabArgs eq $tempMatlabArgs)
			{
				print(STDERR "No valid data available to plot.\n");
			}
			else
			{
				executeMatlab($matlabArgs);
				$totPlotNum ++;
			}
		}
	}
}

###############################################################################
# Generates the argument list needed for the Matlab part of Plotgen. This is
# for standard runs that plots variables from each input folder on the same
# plot.
###############################################################################
sub buildMatlabStringStd()
{
	my $CASE = shift(@_);

	$plotCount = 0;
	
	# Get Common Case information
	my $caseName =  $CASE::CASE{'name'};
	my $startTime =  $CASE::CASE{'startTime'};
	my $endTime =  $CASE::CASE{'endTime'};
	my $startHeight =  $CASE::CASE{'startHeight'};
	my $endHeight =  $CASE::CASE{'endHeight'};
	
	# Get array of plots from case file
	my @plots;
	push(@plots, @{$CASE::CASE{'plots'}});

	# Get plots from .case file
	for(my $count = 0; $count < @plots; $count++)
	{
		# Counters for automatic lines
		$lineStyleCounter = 0;
		$lineColorCounter = 0;
		$lineWidthCounter = 0;

		my $plotTitle = $plots[$count]{'plotTitle'};
		my $units = $plots[$count]{'axisLabel'};
		my $type = $plots[$count]{'type'};

		my $matlabArgs = "\'$caseName\', \'$plotTitle\', $count, \'$type\', $startTime, $endTime, $startHeight, $endHeight, \'$units\', $randInt";
		my $tempMatlabArgs = $matlabArgs;

		my @lines;
		push(@lines, @{$plots[$count]{'lines'}});

		for(my $lineNum = 0; $lineNum < @lines; $lineNum++)
		{
			my $name = $lines[$lineNum]{'name'};
			my $expression = $lines[$lineNum]{'expression'};

			my $type = $lines[$lineNum]{'type'};

			if(($type eq "les" && $plotLes == 1) || ($type eq "dec17" && $plotDec) || ($type eq "bestever" && $plotBest))
			{
				my $file = "$lines[$lineNum]{'filename'}";
				if(-e $file)
				{
					my $title;

					if($type eq "les")
					{
						$title = "LES";
					}
					elsif($type eq "dec17")
					{
						$title = "HOC 12/17/2005";
					}
					elsif($type eq "bestever")
					{
						$title = "HOC \"best-ever\"";
					}

					my $lineWidth = $lines[$lineNum]{'lineWidth'};
					my $lineStyle = $lines[$lineNum]{'lineType'};
					my $lineColor = $lines[$lineNum]{'lineColor'};
	
					$matlabArgs = "$matlabArgs, \'$file\', \'$expression\', \'$title\', $lineWidth, \'$lineStyle\', \'$lineColor\'";
				}
			}
			elsif($type eq "auto")
			{
				foreach (@inputDirs)
				{
					my $file = "$_/$lines[$lineNum]{'filename'}";

					if(-e $file)
					{
						my $title = basename($_);
						
						my $lineWidth = $lineWidths[$lineWidthCounter];
						my $lineStyle = $lineStyles[$lineStyleCounter];
						my $lineColor = $lineColors[$lineColorCounter];
						
						$matlabArgs = "$matlabArgs, \'$file\', \'$expression\', \'$title\', $lineWidth, \'$lineStyle\', \'$lineColor\'";

						incrementLineTypes();
					}
				}
			}
		}

		# Check to see if there are lines to be plotted:
		if($matlabArgs eq $tempMatlabArgs)
		{
			#print(STDERR "No valid data available to plot.\n");
		}
		else
		{
			$plotCount++;
			executeMatlab($matlabArgs);	
		}
	}
}

sub executeMatlab()
{
	my $matlabArgs = shift(@_);

	#my $args = "echo \"quit\" | $MATLAB -nodisplay -nodesktop -r PlotCreator\"($matlabArgs)\"";
	my $args = "PlotCreator\"($matlabArgs)\"";

	#print $args;

	system("echo $args > matlab_pipe");
}

###############################################################################
# Increments the counters for the automatic line types
###############################################################################
sub incrementLineTypes()
{
	# Increment counters
	if($lineColorCounter + 1 >= @lineColors)
	{
		$lineColorCounter = 0;

		if($lineStyleCounter + 1 >= @lineStyles)
		{
			$lineStyleCounter = 0;

			if($lineWidthCounter + 1 >= @lineWidths)
			{
				$lineWidthCounter = 0;
			}
			else
			{
				$lineWidthCounter ++;
			}
		}
		else
		{
			$lineStyleCounter ++;

			if($lineWidthCounter + 1 >= @lineWidths)
			{
				$lineWidthCounter = 0;
			}
			else
			{
				$lineWidthCounter ++;
			}
		}
	}
	else
	{
		$lineColorCounter++;

		if($lineStyleCounter + 1 >= @lineStyles)
		{
			$lineStyleCounter = 0;
	
			if($lineWidthCounter + 1 >= @lineWidths)
			{
				$lineWidthCounter = 0;
			}
			else
			{
				$lineWidthCounter ++;
			}
		}
		else
		{
			$lineStyleCounter ++;

			if($lineWidthCounter + 1 >= @lineWidths)
			{
				$lineWidthCounter = 0;
			}
			else
			{
				$lineWidthCounter ++;
			}
		}
	}
}

sub cleanup()
{
	# Copy temp. output folder to actual output location and remove the temp. folder
	dircopy($outputTemp, $output);
	rmtree($outputTemp);

	$ENV{'DISPLAY'} = $sessionType;
}

###############################################################################
# Checks all input directories to see if one of them contains the current case.
# Will return true if at least on input folder contains data, otherwise false.
# #############################################################################
sub dataExists()
{
	my $dataFile = shift(@_);
	my $retValue = 0;
	
	foreach (@inputDirs)
	{
		my @files = <$_/$dataFile*>;
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
		main::HELP_MESSAGE();
	}

	my %option = ();
	getopts("rlbdanh?", \%option);

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
		main::HELP_MESSAGE();
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
				$outputTemp = "/tmp/output" . "_" . "$randInt";
			}
			
			$currentCount++;
		}
	}

	if(@inputDirs == 0)
	{
		main::HELP_MESSAGE();
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
		rmtree($output);
		mkdir $output;
		mkdir $outputTemp;
	}
}

###############################################################################
# Prints the help message
###############################################################################
sub main::HELP_MESSAGE()
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

###############################################################################
# Prints the version number of Plotgen.
###############################################################################
sub main::VERSION_MESSAGE()
{
	print("Plotgen version $VERSION, Copyright (c) 2009 Larson Group.\n");
}
