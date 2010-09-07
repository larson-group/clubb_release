#!/usr/bin/perl
###########################################################################
# Plotgen v3.41
#
# Documentation is available here:
# http://larson-group.com/twiki/bin/view.pl/Documentation/CarsonDoc/Plotgen3
###########################################################################

package Plotgen;

my $dirPrefix = "";
# This is only needed when installing globally to ensure additional perl 
# modules can be found.
BEGIN
{
    # Only do this if plotgen was executed from a path other than the plotgen directory
    if($0 ne "./plotgen.pl")
    {
        push @INC,"/home/matlabuser/plotgen";
    }
}

# Add the path to the plotgen directory for executing external scripts
if($0 ne "./plotgen.pl")
{
    $dirPrefix = "/home/matlabuser/plotgen/";
}

use strict;
use CaseReader;
use OutputWriter;
use Cwd 'abs_path';
use Cwd;
use Switch;
use Getopt::Std;
use File::Basename;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Path;

# The location of MATLAB. External users probably will not need the 
# sudo -u matlabuser part.
my $MATLAB = "sudo -u matlabuser /usr/local/bin/matlab -nodisplay -nodesktop";

# The pipe name used to communicate with MATLAB
my $matlabPipe;

# The lock file used for image conversion
my $imageConversionLock;

# Plotgen Version Number
my $VERSION = 3.41;

# Used to create a "random" output directory so multiple runs
# don't overwrite each other.
my $randInt = int(rand(999999999999999)) + 1;

my $DPI = 300;
my $QUALITY = 100;

# Field to hold the total number of cases plotted. This will
# be used when automatically specifying image quality.
my $caseCount = 0;

# If running in nightly mode, this value should be set to 1
my $nightly = 0;

# Which mode to run in. Default to plotgen (CLUBB). This will specify 
# what directory to look for .case files in. 
# Valid modes: 
#   plotgen
#   splotgen
#   wrfgen 
my $plotgenMode = "plotgen";

# The type of data file to use.
# Valid types:
#   grads
#   netcdf
my $dataFileType = "grads";

# Specifies to overwrite a directory (Default: false)
my $overwrite = 0;
my @inputDirs;
my $output = "";
my $outputTemp = "";
my $plotLes = 0;
my $plotBest = 0;
my $plotDec = 0;

# Specifies high quality images. If this is true, keep
# image quality at default.
my $highQuality = 0;

# Whether or not to delete the EPS images. Default (0) is to delete them.
my $keepEps = 0;

# Whether or not to compress output into a maff file. Default (0) is no.
my $outputAsMaff = 0;

# Custom Color Definitions for "CLUBB_current" and "CLUBB_previous"
my $lt_blue = "[ 0.00, 0.63, 1.00 ]";
my $orange = "[ 0.94, 0.50, 0.16 ]";
# Custom Colors (taken from splotgen)
my $goldenRod = "[ 0.80, 0.80, 0.10 ]";
my $grey = "[ 0.40, 0.40, 0.40 ]";
my $purple = "[ 0.50, 0., 0.50 ]";
my $peach = "[ 1.00, 0.90, 0.40 ]";
my $darkGreen = "[ 0.00, 0.40, 0.00 ]";
# Arrays to cycle through when auto is set for lines
my @lineStyles = ("--", "-", "-.", "-");
my @lineColors = ($orange, $lt_blue, $purple, "blue", $peach, $grey, $goldenRod, $darkGreen, "cyan", "yellow", "magenta", "green", "red");
my @lineWidths = (6, 4, 3, 2);

# Counters for automatic lines
my $lineStyleCounter = 0;
my $lineColorCounter = 0;
my $lineWidthCounter = 0;

my $outputIndex = "";
my $navigationPage = "navigation.html";
my $indexPage = "index.html";

my $plotCount = 0;

my $consoleOutput = "$dirPrefix./console_output.pl";

# Do not keep permissions on file copies
$File::Copy::Recursive::KeepMode = 0;

# This is needed so we can restore the display after plotgen runs.
# We need to clear this environment variable before running, otherwise
# matlab will fail.
my $sessionType = $ENV{'DISPLAY'}; 
$SIG{INT} = "cleanup";

# This might be needed to allow plotgen to run child processes. So far, it works without it.
#sub REAPER {
#    my $waitedpid = wait;
#    # loathe SysV: it makes us not only reinstate
#    # the handler, but place it after the wait
#    $SIG{CHLD} = \&REAPER;
#}
#$SIG{CHLD} = \&REAPER;

my @casesExecuted;

main();

###############################################################################
# Plotgen 3 Main subroutine
# Arguments:
#   None.
###############################################################################
sub main()
{
    readArgs();

    $matlabPipe = "matlab_pipe_$randInt";

    $imageConversionLock = "$outputTemp/.lock";

    # Ensure that Matlab can write to the temp output folder
    my $mode = 0777; chmod $mode, $outputTemp;

    print("Input Folders: @inputDirs\n");
    print("Output Folder: $output \n");

    $outputIndex = "$outputTemp/plots.html";

    OutputWriter->writeNavPageStart("$outputTemp/$navigationPage");
    OutputWriter->writeIndex("$outputTemp/index.html");

    # Fork to make MATLAB faster
    my $pid = fork();

    if($pid == 0) # Child
    {
        $ENV{'DISPLAY'} = '';

        system("mkfifo $matlabPipe");
        system("$MATLAB <> $matlabPipe");

        system("rm $imageConversionLock");

        print("\nPlease wait while the remaining images are converted...\n");

        # Wait until all images are converted
        my @epsFiles = <$outputTemp/*eps>;
        my @jpgFiles = <$outputTemp/jpg/*jpg>;
        my $arraySizeEps = @epsFiles;
        my $arraySizeJpg = @jpgFiles;

        if($keepEps == 0)
        {
            # Block this thread from moving on until there are no more eps files
            while($arraySizeEps != 0)
            {
                @epsFiles = <$outputTemp/*eps>;
                $arraySizeEps = @epsFiles;
                sleep(1);
            }
        }
        else
        {
            # Block this thread from moving on until the number of JPGs equal the
            # numper of EPS files.
            while($arraySizeEps != $arraySizeJpg)
            {
                @epsFiles = <$outputTemp/*eps>;
                @jpgFiles = <$outputTemp/jpg/*jpg>;
                $arraySizeEps = @epsFiles;
                $arraySizeJpg = @jpgFiles;

                sleep(1);
            }
        }

        OutputWriter->writeFooter($outputIndex);
    
        cleanup();
    
        if($outputAsMaff == 0)
        {
            print("Done! To display the plots, open: \n$output/index.html \nin your web browser.\n");
        }
        else
        {
            print("Done! To display the plots, open: \n" . $output . ".maff\nin your maff compatible web browser.\n");
        }

        exit(0);
    }
    else # Parent
    {
        # Now fork to create images in the background. This should hopefully
        # speed things up a little
        my $convertPid = fork();

        if($convertPid == 0) # Image Conversion Child
        {
            system("touch $imageConversionLock");
            convertEps();
            exit(0);
        }
        else # Main thread
        {
            if($plotgenMode eq "splotgen")
            {
                OutputWriter->writeHeader($outputIndex, "Splotgen");
                
                # Only print the SAM_CLUBB Equiv. table if in nightly mode.
                if($nightly == 1)
                {
                    OutputWriter->writeSamSubHeader($outputIndex);
                }
            }
            elsif($plotgenMode eq "plotgen")
            {
                OutputWriter->writeHeader($outputIndex, "Plotgen");
            }
	          elsif($plotgenMode eq "wrfgen")
	          {
		            OutputWriter->writeHeader($outputIndex, "WRFGen");
	          }
    
            runCases();

            # Quit MATLAB
            system("echo quit > $matlabPipe");
          
            # Wait for image conversion process to exit
            wait;
        }

        # Wait for the MATLAB child process to exit
        wait;

        exit(0);
    }
}

###############################################################################
# Returns the path to where case files should be looked for in.
# Arguments:
#   None.
###############################################################################
sub getCasePath()
{
    my $casePath = "cases";

    if($plotgenMode eq "plotgen")
    {
        return "$casePath/clubb";
    }
    elsif($plotgenMode eq "splotgen")
    {
        return "$casePath/sam_clubb";
    }
    elsif($plotgenMode eq "wrfgen")
    {
	return"$casePath/wrf";
    }
}

###############################################################################
# Runs all of the .case file in the cases folder.
# Arguments:
#   None.
###############################################################################
sub runCases()
{
    # Counter used to place images
    my $count = 0;

    # Get the path to look for case files in
    my $casePath = getCasePath();

    my @cases = <$casePath/*.cas*>;

    # Loop through each .case file so the case can be plotted    
    foreach my $file (@cases) 
    {
        # The first part of this if statement will run all case files
        # if the nightly flag was passed in. The second part will only
        # run the case file if the nightly flag was not passed in and
        # the filename does not contain _nightly. 
        # We need to do this because we have special condition case files.
        # For example, gabls2 plots special variables at night. Therefore,
        # there is a special gabls2_nightly.case file that only includes 
        # those plots. 
        if($nightly == 1 || ($nightly == 0 && $file !~ m/_nightly/))
        {
            #print("Case file: $file\n");
            my $runCase = 'true';
            # Read the case file. If there is an error, display it and continue without plotting it.
            if (my $err = CaseReader->readCase($file))
            {
                system("$consoleOutput -w \"$err\"");
                $runCase = 'false';
            }
            else
            {
                # Check to see if case was plotted already. This fixes the infinite loop problem
                # when converting images. If the case was already plotted, do not do it again.
                foreach my $chkCase (@casesExecuted)
                {
                    if($chkCase eq  $CASE::CASE{'name'})
                    {
                        $runCase = 'false';
                    }
                }
            }
    
            if($runCase eq 'true' && dataExists($CASE::CASE) && ($CASE::CASE{'enabled'} ne 'false'))
            {
                push(@casesExecuted, $CASE::CASE{'name'});

                # Print the case title to the HTML page
                OutputWriter->writeCaseTitle($outputIndex, $CASE::CASE{'headerText'});
                OutputWriter->writeNavPageCase("$outputTemp/$navigationPage", $CASE::CASE{'name'}, $CASE::CASE{'headerText'});
        
                # Print any additional text/html specified
                if($nightly == 1) # If in nightly mode
                {
                    my $nightlySubText = $CASE::CASE{'nightlyOutput'}{'subText'};
                    my $nightlySubHtml = $CASE::CASE{'nightlyOutput'}{'subHtml'};
                    
                    # Check to see if there was any additional text specified. If there was,
                    # write it to the HTML file.
                    if($nightlySubText ne "")
                    {
                        OutputWriter->writeSubHeader($outputIndex, $nightlySubText);
                    }
                    
                    if($nightlySubHtml ne "")
                    {
                        OutputWriter->writeSubHtml($outputIndex, $nightlySubHtml);
                    }
                }
                else # If not in nightly mode
                {
                    my $subText = $CASE::CASE{'additionalOutput'}{'subText'};
                    my $subHtml = $CASE::CASE{'additionalOutput'}{'subHtml'};
                    
                    if($subText ne "")
                    {
                        OutputWriter->writeSubHeader($outputIndex, $subText);
                    }
                    
                    if($subHtml ne "")
                    {
                        OutputWriter->writeSubHtml($outputIndex, $subHtml);
                    }
        
                }
                
                # Check to see if this is a budget plot or standard plot
                if($CASE::CASE{'type'} eq "budget")
                {
                    buildMatlabStringBudget($CASE::CASE, $count);
                    
                    # Add image file to HTML page
                    placeImages($CASE::CASE{'name'} . "_" . $count . "_budget", $plotCount, 0);
                }
                else
                {
                    my $nightlyCase = 0;

                    if($nightly == 1 && $file =~ m/_nightly/)
                    {
                        $nightlyCase = 1;
                    }
                   
                    buildMatlabStringStd($CASE::CASE, $nightlyCase);
                    
                    # Add image file to HTML page
                    placeImages($CASE::CASE{'name'}, $plotCount, $nightlyCase);
                }
    
                $count++;
            }
            else
            {
                #print("Not plotting $CASE::CASE{'headerText'}\n");
            }
        }
    }

    OutputWriter->writeNavPageClose("$outputTemp/navigation.html");
}

###############################################################################
# Converts all EPS files to JPG file
# Arguments:
#   None.
###############################################################################
sub convertEps()
{
    mkdir "$outputTemp/jpg" unless -d "$outputTemp/jpg";
    
    my @epsFiles = <$outputTemp/*eps>;
    my @jpgFiles = <$outputTemp/jpg/*jpg>;
    my $arraySizeEps = @epsFiles;
    my $arraySizeJpg = @jpgFiles;

    # This will keep calling convertEps() if the image conversion lock exists OR:
    #   If -e was not passed in:
    #       Until there are no more eps files left
    #   If -e was passed in:
    #       Until the number of eps files equals the number of jpg files       
    if(($keepEps == 0 && $arraySizeEps != 0) || ($keepEps == 1 && ($arraySizeEps != $arraySizeJpg)) || -e "$imageConversionLock")
    {
        # Set the image scale if -q was not passed in
        if($highQuality == 0)
        {
            $DPI = 120;
            $QUALITY = 80;
        }    

        foreach my $eps (@epsFiles)
        {
            my $filename = basename($eps);

            if(! -e "$outputTemp/jpg/$filename.jpg")
            {
                # First convert the image to png and trim all white space
                system("convert -density $DPI -colorspace RGB -trim $eps $outputTemp/jpg/$filename.png");

                # Then convert (and scale if not in high quality mode) to jpg
                if($highQuality == 0)
                {
                    system("convert -geometry 324x312\\! -quality $QUALITY $outputTemp/jpg/$filename.png $outputTemp/jpg/$filename.jpg");
                }
                else
                {
                    system("convert -quality $QUALITY $outputTemp/jpg/$filename.png $outputTemp/jpg/$filename.jpg");
                }
        
                unlink("$outputTemp/jpg/$filename.png");
                if($keepEps == 0)
                {
                    unlink($eps);
                }
            }
        }
       
        sleep(1);

        # Call convert again. If the doesn't exist anymore, and there are no more .eps files,
        # this will be the last time the subroutine is called.
        convertEps();
    }
}

###############################################################################
# Places the jpg images on the HTML page
# Arguments:
#   placeImages(String caseName, int numImages, boolean nightlyCaseFile)
#     - caseName: The case name to plot images for
#     - numImages: The number of images for the case
#     - nightlyCaseFile: If the case file ends in _nightly.case, this should be
#                        1, otherwise 0.
###############################################################################
sub placeImages()
{
    my $caseName = shift(@_);
    my $numImages = shift(@_);
    my $nightlyCaseFile = shift(@_);

    OutputWriter->printDivCenter($outputIndex);

    for(my $x = 0; $x < $numImages; $x++)
    {
        if($nightlyCaseFile == 1)
        {
            OutputWriter->placeImage($outputIndex, "jpg/$caseName" . "_nightly_" . "$x.eps.jpg");
        }
        else
        {
            OutputWriter->placeImage($outputIndex, "jpg/$caseName" . "_" . "$x.eps.jpg");
        }
    }

    OutputWriter->printCloseDivCenter($outputIndex);
}

###############################################################################
# Generates argument lists needed for the Matlab part of Plotgen. This is for
# budget plots that will plot each input directory on a seperate plot.
# Arguments:
#   buildMatlabStringBudget(Hash CASE, int budgetNumber)
#     - CASE: The case file
#     - budgetNumber: The current budget number. This is incremented after each
#                     budget is plotted.
###############################################################################
sub buildMatlabStringBudget()
{
    $plotCount = 0;
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
    
    # Loop through each input folder and create a plot for each folder
    for(my $bgtPlotCount = 0; $bgtPlotCount < @inputDirs; $bgtPlotCount++)
    {
        my $plotTitle = basename($inputDirs[$bgtPlotCount]);
        
        # Loop through each plot in the Case file
        for(my $plotNum = 0; $plotNum < @plots; $plotNum++)
        {

	    my $units = $plots[$plotNum]{'axisLabel'};
	    my $type = $plots[$plotNum]{'type'};

            # Counters for automatic lines
            $lineStyleCounter = 0;
            $lineColorCounter = 0;
            $lineWidthCounter = 0;

            # Check to see if start and end time was specified for the specific plot. This is usually done for timeseries plots. If it isn't
            # use case defined times.
            my $startTimeOverride = $plots[$plotNum]{'startTime'};
            my $endTimeOverride = $plots[$plotNum]{'endTime'};

            if($startTimeOverride)
            {
                $startTime = $startTimeOverride;
            }

            if($endTimeOverride)
            {
                $endTime = $endTimeOverride;
            }


            my $matlabArgs = "\'$caseName\', \'$plotTitle\', $totPlotNum, \'$type\', $startTime, $endTime, $startHeight, $endHeight, \'$units\', $randInt";
            my $tempMatlabArgs = $matlabArgs;
        
            my @lines;
            push(@lines, @{$plots[$plotNum]{'lines'}});
            
            for(my $lineNum = 0; $lineNum < @lines; $lineNum++)
            {
                my $file = "$inputDirs[$bgtPlotCount]/$lines[$lineNum]{'filename'}";
                
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
                #print(STDERR "No valid data available to plot.\n");
            }
            else
            {
                executeMatlab($matlabArgs);
                $totPlotNum ++;
                $plotCount++;
            }
        }
    }
}

###############################################################################
# Generates the argument list needed for the Matlab part of Plotgen. This is
# for standard runs that plots variables from each input folder on the same
# plot.
#
# Arguments:
#   buildMatlabStringStd(Hash case, boolean nightlyCaseFile)
#     - case: The case file hash
#     - nightlyCaseFile: If the case file ends in _nightly.case, this should be
#                        1, otherwise 0.
###############################################################################
sub buildMatlabStringStd()
{
    my $CASE = shift(@_);
    my $nightlyCaseFile = shift(@_);

    $plotCount = 0;
    
    # Get Common Case information
    my $caseName;

    # Modify the case name if this is a nightly case
    if($nightlyCaseFile == 1)
    {
        $caseName =  $CASE::CASE{'name'} . "_nightly";
    }
    else
    {
        $caseName = $CASE::CASE{'name'};
    }
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

        # Check to see if start and end time was specified for the specific plot. This is usually done for timeseries plots. If it isn't
        # use case defined times.
        my $startTimeOverride = $plots[$count]{'startTime'};
        my $endTimeOverride = $plots[$count]{'endTime'};

        if($startTimeOverride ne "")
        {
            $startTime = $startTimeOverride;
        }

        if($endTimeOverride ne "")
        {
            $endTime = $endTimeOverride;
        }

        my $matlabArgs = "\'$caseName\', \'$plotTitle\', $count, \'$type\', $startTime, $endTime, $startHeight, $endHeight, \'$units\', $randInt";
        my $tempMatlabArgs = $matlabArgs;

        my @lines;
        push(@lines, @{$plots[$count]{'lines'}});

        for(my $lineNum = 0; $lineNum < @lines; $lineNum++)
        {
            my $name = $lines[$lineNum]{'name'};
            my $expression = $lines[$lineNum]{'expression'};
            my $type = $lines[$lineNum]{'type'};

            if($type eq "auto")
            {
                foreach (@inputDirs)
                {
                    my $file = "$_/$lines[$lineNum]{'filename'}";

                    if($dataFileType eq "netcdf")
                    {
                        # Replace all .ctl file extensions with .nc
                        # By default, the case files have .ctl extensions
                        $file =~ s/.ctl/.nc/g;
                    }

                    if(-e $file)
                    {
                        my $title;
                        my $folderName = basename($_);
                        $folderName =~ s/_/\\_/g; # Replace all '_' with '\_' for MATLAB

                        if($name eq "auto")
                        {
                            $title = $folderName;
                        }
                        else
                        {
                            $title = $name;

                            # Replace any '{0}' with the folder name
                            $title =~ s/\{0\}/$folderName/;
                        }
                        
                        my $lineWidth = $lineWidths[$lineWidthCounter];
                        my $lineStyle = $lineStyles[$lineStyleCounter];
                        my $lineColor = $lineColors[$lineColorCounter];
                        
                        $matlabArgs = "$matlabArgs, \'$file\', \'$expression\', \'$title\', $lineWidth, \'$lineStyle\', \'$lineColor\'";
                    }

                    incrementLineTypes();
                }
            }
            elsif(($type eq "les" && $plotLes == 1) || ($type eq "dec17" && $plotDec) || ($type eq "bestever" && $plotBest))
            {
                my $file = "$lines[$lineNum]{'filename'}";
                if(-e $file)
                {
                    my $title = $name;

                    my $lineWidth = $lines[$lineNum]{'lineWidth'};
                    my $lineStyle = $lines[$lineNum]{'lineType'};
                    my $lineColor = $lines[$lineNum]{'lineColor'};
    
                    $matlabArgs = "$matlabArgs, \'$file\', \'$expression\', \'$title\', $lineWidth, \'$lineStyle\', \'$lineColor\'";
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

###############################################################################
# Executes the PlotCreator.
# Arguments:
#   executeMatlab(String matlabArgs)
#     - matlabArgs: The arguments to be passed to the plot creator
###############################################################################
sub executeMatlab()
{
    my $matlabArgs = shift(@_);

    my $args = "PlotCreator\"($matlabArgs)\"";

    #print("\n$args\n\n");

    system("echo $args > $matlabPipe");
}

###############################################################################
# Changes to the next line style, width, and color.
# Arguments:
#   None.
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

###############################################################################
# Does necessary cleanup code.
# Arguments:
#   None.
###############################################################################
sub cleanup()
{
    # Copy temp. output folder to actual output location
    if($outputAsMaff == 1)
    {
        print("\nCompressing output to maff file:\n");
        my $outputName = $outputTemp;
        
        # Shortens the path name to just the lowest directory, if this is not done,
        # the UNIX zip command makes a zip including the absolute file structure.
        while($outputName =~ m/\/[^\/]*\//)
        {
            $outputName =~ s/\/[^\/]*//;
        }
        $outputName = substr $outputName, 1;
        
        # Remove possible .maff from output name to avoid "*.maff.maff" file names
        if($output =~ m/.maff$/)
        {
            $output = substr $output, 0, -5
        }

        system("cd $outputTemp/.. && zip -r " . $output . ".maff $outputName/");
    }
    else
    {
        dircopy($outputTemp, $output);
    }
    
    # Remove the temp. folder
    rmtree($outputTemp);

    # Remove matlab pipe
    system("rm $matlabPipe");

    $ENV{'DISPLAY'} = $sessionType;
}

###############################################################################
# Checks all input directories to see if one of them contains the current case.
# Will return true if at least on input folder contains data, otherwise false.
# Arguments:
#   dataExists(Hash CASE)
#     - CASE: The case file to check
###############################################################################
sub dataExists()
{
    my $CASE = shift(@_);
    
    # We need to figure out if any of the input directories contain
    # data to plot. I am assuming that input directories are of type
    # "auto" for now. If the type is auto, look in each input directory
    # to see if some data exists and return 1. 
    
    # Get array of plots from case file
    my @plots;
    push(@plots, @{$CASE::CASE{'plots'}});

    for(my $count = 0; $count < @plots; $count++)
    {
        my @lines;
        push(@lines, @{$plots[$count]{'lines'}});
        
        for(my $lineNum = 0; $lineNum < @lines; $lineNum++)
        {
            my $lineType = $lines[$lineNum]{'type'};
            my $filename = $lines[$lineNum]{'filename'};
            
            if($dataFileType eq "netcdf")
            {
                # Replace all .ctl file extensions with .nc
                $filename =~ s/.ctl/.nc/g;
            }

            if($lineType eq 'auto')
            {
                foreach(@inputDirs)
                {
                    my @inputFiles = <$_/$filename*>;

                    if(@inputFiles)
                    {
                        return 1;
                    }
                }
            }
        }
    }

    return 0;
}

###############################################################################
# Reads any arguments passed in.
# Arguments:
#   None.
###############################################################################
sub readArgs()
{
    my $numArgs = $#ARGV + 1;

    if($numArgs == 0)
    {
        main::HELP_MESSAGE();
    }

    my %option = ();
    my $result = getopts("rlbdanqemhcswtg?", \%option);

    # A 1 will be returned from getopts if there weren't any
    # invalid options passed in.
    if($result != 1)
    {
        main::HELP_MESSAGE();
    }

    if ($option{r}) # Option to replace data if it already exists
    {
        $overwrite = 1;
    }

    if ($option{l}) # Option to plot LES data
    {
        $plotLes = 1;
    }

    if ($option{b}) # Option to plot Chris Golaz Best Ever Data
    {
        $plotBest = 1;
    }

    if ($option{d}) # Option to plot HOC Dec. 17
    {
        $plotDec = 1;
    }

    if ($option{a}) # Plot LES, CGBE, and HOC Dec. 17
    {
        $plotLes = 1;
        $plotBest = 1;
        $plotDec = 1;
    }

    if ($option{n}) # Run in nightly mode
    {
        $nightly = 1;
    }

    if ($option{q}) # Force high quality images
    {
        $highQuality = 1;
    }

    if ($option{e}) # Keep EPS files
    {
        $keepEps = 1;
    }
    
    if ($option{g})
    {
        $dataFileType = "grads";
    }

    if ($option{t})
    {
        $dataFileType = "netcdf";
    }

    if ($option{m}) # Output as maff file
    {
    
      if ($option{e}) # It's inconvenient for maff if eps files are kept
      {
      print("Argument conflict: Please do not simultaneously choose to 
           save .eps files (-e) and create a .maff file (-m)\n");
      exit(1);
      }
      
       $outputAsMaff = 1;
    }

    if($option{s})
    {
        $plotgenMode = "splotgen";
    }

    if($option{c})
    {
        $plotgenMode = "plotgen";
    }

    if($option{w})
    {
	    $plotgenMode = "wrfgen";
    }

    if ($option{h}) # Print the help message
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
                    print("The input folder $currentDir does not exist.\n");
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

    # If there were no input directories passed in, print the help and exit
    if(@inputDirs == 0)
    {
        main::HELP_MESSAGE();
    }

    # Finally, check to see if the output folder exists. If it does, and
    # '-r' was not passed in, exit. Otherwise, create it.
    
    if($outputAsMaff == 0 && -e $output && $overwrite == 0)
    {
        system("$consoleOutput -s \"Output folder already exists. To overwrite, use the -r option.\"");
        exit(1);
    }
    elsif($outputAsMaff == 1 && -e $output . ".maff" && $overwrite == 0)
    {
        system("$consoleOutput -s \"maff file already exists. To overwrite, use the -r option.\"");
        exit(1);
    }
    else
    {
        if($outputAsMaff == 0)
        {
            rmtree($output);
            mkdir $output;
        }
        mkdir $outputTemp;

        # The following was added to allow us to install plotgen globally.
        # This will follow the symlink to the actual plotgen.pl script.
        my $plotgenDirectory = readlink($0);

        print("$plotgenDirectory\n");

        if($plotgenDirectory ne "")
        {
            $plotgenDirectory = dirname(abs_path($plotgenDirectory));
            
            # In case we aren't already in the directory that plotgen.pl is, goto it.
            chdir($plotgenDirectory);
        }
    }
}

###############################################################################
# Prints the help message
# Arguments:
#   None.
###############################################################################
sub main::HELP_MESSAGE()
{
    print("Usage: plotgen [OPTION]... INPUT... OUTPUT\n");
    print("  -c\tPlot CLUBB cases [DEFAULT] (equiv to plotgen)\n");
    print("  -s\tPlot SAM_CLUBB cases (equiv to splotgen)\n");
    print("  -w\tPlot WRF_CLUBB cases\n");
    print("  -r\tIf the output folder already exists, replace the contents\n");    
    print("  -l\tPlot LES data for comparison.\n");    
    print("  -b\tPlot HOC Best Ever data for comparison.\n");    
    print("  -d\tPlot HOC 12/17/2005 data for comparison.\n");    
    print("  -a\tSame as -lbd. Plots LES, Best Ever, and 12/17/2005 data for comparison.\n");
    print("  -n\tRuns in nightly mode.\n");
    print("  -q\tOutputs high quality images (does not auto scale).\n");
    print("  -e\tDoes not delete EPS images after conversion.\n");
    print("  -m\tOutputs plots compressed inside a .maff directory.\n");
    print("  -g\tUses GrADS data files. [DEFAULT]\n");
    print("  -t\tUses NetCDF data files.\n");
    print("  -h\tPrints this help message.\n");
    exit(0);
}

###############################################################################
# Prints the version number of Plotgen.
# Arguments:
#   None.
###############################################################################
sub main::VERSION_MESSAGE()
{
    print("Plotgen version $VERSION, Copyright (c) 2010 Larson Group.\n");
}
