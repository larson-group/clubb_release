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

# Used to create a "random" output directory so multiple runs
# don't overwrite each other.
my $randInt = rand(999999999999999);

my $DPI = 75;

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

# Arrays to cycle through when auto is set for lines
my @lineStyles = ("-", "--", "-.");
my @lineColors = ("blue", "green", "red", "cyan", "yellow", "black", "magenta");
my @lineWidths = (5.5, 5, 4.5, 4, 3.5, 3, 2.5, 2, 1.5, 1);

my $outputIndex = "";

# Do not keep permissions on file copies
$File::Copy::Recursive::KeepMode = 0;

# This is needed so we can restore the display after plotgen runs.
# We need to clear this environment variable before running, otherwise
# matlab will fail.
my $sessionType = $ENV{'DISPLAY'}; 
$SIG{INT} = "cleanup";

# Call the entry point to Plotgen
main();

###############################################################################
# Plotgen 3 Main subroutine
###############################################################################
sub main()
{
	$ENV{'DISPLAY'} = '';

	readArgs();

	# Ensure that Matlab can write to the temp output folder
	chmod(0777, $outputTemp);

	print("Input Folders: @inputDirs\n");
	print("Output Folder: $output \n");

	$outputIndex = $outputTemp . "/index.html";

	OutputWriter->writeHeader($outputIndex);
	runCases();
	OutputWriter->writeFooter($outputIndex);

	cleanup();

	print("\n");

	exit(0);
}

###############################################################################
# Runs all of the .case file in the cases folder.
###############################################################################
sub runCases()
{
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
			
			callMatlab($CASE::CASE);

			# Convert the eps files to jpq
			convertEps($CASE::CASE{'name'});

			# Add image file to HTML page
			placeImages($CASE::CASE{'name'});
		}
		else
		{
			print("Not plotting $CASE::CASE{'name'}\n");
		}
	}
}

sub placeImages()
{
	my $caseName = shift(@_);

	print("Case Name: $caseName\n");

	OutputWriter->printDivCenter($outputIndex);

	my @imgFiles = <$outputTemp/jpg/$caseName*.jpg>;

	for(my $x = 0; $x < @imgFiles; $x++)
	{
		OutputWriter->placeImage($outputIndex, "jpg/$caseName" . "_" . "$x.eps.jpg");
	}

	OutputWriter->printCloseDivCenter($outputIndex);
}

###############################################################################
# Generates the argument list needed for the Matlab part of Plotgen. Also,
# Matlab is executed using sudo
# #############################################################################
sub callMatlab()
{
	my $CASE = shift(@_);

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
		my $lineStyleCounter = 0;
		my $lineColorCounter = 0;
		my $lineWidthCounter = 0;

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

			if(($type eq "les" && $plotLes == 1) || ($type eq "dec17" && $plotDec) || ($type eq "bestEver" && $plotBest))
			{
				my $file = "$lines[$lineNum]{'filename'}";
				if(-e $file)
				{
					my $title = $type;
					my $lineWidth = $lines[$lineNum]{'lineWidth'};
					my $lineStyle = $lines[$lineNum]{'lineType'};
					my $lineColor = $lines[$lineNum]{'lineColor'};
	
					$matlabArgs = "$matlabArgs, \'$file\', \'$name\', \'$expression\', \'$title\', $lineWidth, \'$lineStyle\', \'$lineColor\'";
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
	
						$matlabArgs = "$matlabArgs, \'$file\', \'$name\', \'$expression\', \'$title\', $lineWidth, \'$lineStyle\', \'$lineColor\'";

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
				}
			}
		}

		# Check to see if there are lines to be plotted:
		if($matlabArgs eq $tempMatlabArgs)
		{
			print(STDERR "No valid data available to plot.\n");
#			cleanup();
#			exit(1);
		}
		else
		{
			my $args = "echo \"quit\" | sudo -u matlabuser /usr/local/bin/matlab -nodisplay -nodesktop -r PlotCreator\"($matlabArgs)\"";
			print("Call: $args\n\n");
			
			system($args);
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
# Converts all EPS files to JPG file
###############################################################################
sub convertEps()
{
	mkdir "$outputTemp/jpg" unless -d "$outputTemp/jpg";
	my @epsFiles = <$outputTemp/*eps>;
	foreach my $eps (@epsFiles)
	{
		my $filename = basename($eps);
		system("convert -density $DPI $eps $outputTemp/jpg/$filename.jpg");

		unlink($eps);
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
				$outputTemp = "$output" . "_" . "$randInt";
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
		mkdir $outputTemp;
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
