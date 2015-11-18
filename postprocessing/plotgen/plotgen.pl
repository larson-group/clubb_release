#!/usr/bin/perl
###########################################################################
# Plotgen v3.104
#
# Documentation is available here:
# http://larson-group.com/twiki/bin/view.pl/Documentation/CarsonDoc/Plotgen3
###########################################################################

package Plotgen;

##This was added to allow this script to function from outside its own directory.
##See other:274
my $dirPrefix;
use Cwd;
use Cwd 'abs_path';
use File::Basename;

BEGIN {
    ( undef, $dirPrefix ) = fileparse( abs_path($0) );
##print STDERR "absolute path is " . $dirPrefix;
    push @INC, $dirPrefix;

    #my $dirOrig = getcwd();
##print STDERR "original folder is " . $dirOrig;
}

use strict;
use Switch;
use Getopt::Std;
use File::Path;
use CaseReader;
use OutputWriter;
use File::Copy::Recursive qw(fcopy rcopy dircopy fmove rmove dirmove);
use File::Copy;

# The location of MATLAB. External users probably will not need the
# sudo -u matlabuser part.
my $MATLAB =
  "sudo -H -u matlabuser /usr/local/bin/matlab -nodisplay -nodesktop";

if ( $^O eq "darwin" ) {
    $MATLAB = "/Applications/MATLAB_R2014a.app/bin/matlab -nodisplay -nodesktop"; #For Macs in W434
}

# The pipe name used to communicate with MATLAB
my $matlabPipe;

# The lock file used for image conversion
my $imageConversionLock;

# Plotgen Version Number
my $VERSION = 3.110;

# Used to create a "random" output directory so multiple runs
# don't overwrite each other.
my $randInt = int( rand(999999999999999) ) + 1;

my $DPI     = 300;
my $QUALITY = 100;

# Field to hold the total number of cases plotted. This will
# be used when automatically specifying image quality.
my $caseCount = 0;

# If running in nightly mode, this value should be set to 1
my $nightly = 0;

# If we are doing a "difference" run, this value will become a 1.
my $diffRun = 0;

# Which mode to run in. Default to plotgen (CLUBB). This will specify
# what directory to look for .case files in.
# Valid modes:
#   plotgen
#   splotgen
#   wrfgen
#   camgen
#   gfdlgen
my $plotgenMode = "plotgen";

# The type of data file to use.
# Valid types:
#   grads
#   netcdf
my $dataFileType = "grads";

# Flag for if budget plots should be plotted (Default: false)
# Set by -bu
my $plotBudgets = 0;

# Flag for if Morrison budget plots should be plotted (Default: false)
# Set by -bumorr
my $plotMorrBudgets = 0;

# Specifies to overwrite a directory (Default: false)
my $overwrite = 0;
my @inputDirs;
my $output     = "";
my $outputTemp = "";
my $plotLes    = 0;
my $plotBest   = 0;
my $plotDec    = 0;

# Specifies high quality images. If this is true, keep
# image quality at default.
my $highQuality = 0;

# Whether or not to delete the EPS images. Default (0) is to delete them.
my $keepEps = 0;

# Whether or not to compress output into a maff file. Default (0) is no.
my $outputAsMaff = 0;

# Whether or not to include setup text files.  Default (1) is yes.
my $includeSetup = 1;

# Whether or not to display legens.  Default (1) is yes.
my $displayLegend = 1;

# Whether or not to use smaller lines for ensemble runs.  Default (0) is no.
my $thinLines = 0;

# Is used for plotting emsemble tuner runs, it will show a set of plots as the same color.
# Default (0) is no
my $ensembleTuner = 0;

# Arrays to cycle through when auto is set for lines
# Brian:  for this branch, I have altered the order to match the PDF plots.
my @lineStyles = ( "-", "--", "-.", "-" );

# These colors were stolen from http://www.colorbrewer2.org/
# Brian:  for this branch, I have altered the colors to match the PDF plots.
my @lineColors = (
    "[ 0.000, 0.000, 1.000 ]",    # blue
    "[ 0.000, 1.000, 0.000 ]",    # green
    "[ 1.000, 0.000, 1.000 ]",    # magenta
    "[ 0.471, 0.471, 0.471 ]",    # gray
    "[ 0.969, 0.506, 0.749 ]",    # pink
    "[ 1.000, 0.498, 0.000 ]",    # orange
    "[ 0.651, 0.337, 0.157 ]",    # brown
    "[ 0.894, 0.102, 0.110 ]",    # red
    "[ 1.000, 1.000, 0.200 ]"     # yellow
);

my @lineWidthsNormal = ( 3.75, 3, 2.25, 1.5 );

# We apply a constant budget width for budget cases.
my @lineWidthsBudget = (3);
my @lineWidths;

# Counters for automatic lines
my $lineStyleCounter = 0;
my $lineColorCounter = 0;
my $lineWidthCounter = 0;

my $outputIndex    = "";
my $navigationPage = "navigation.html";
my $indexPage      = "index.html";

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
sub main() {
    readArgs();

    # If thinLines is enabled then make it so we only use thin solid lines
    if ( $thinLines == 1 ) {
        @lineStyles = ("-");
        @lineWidths = (1);
    }

    # If ensembleTuner is enabled then change the colors to match that of
    # the other ensemble plotter
    if ( $ensembleTuner == 1 ) {
        @lineColors = ( "[1,.65,0.]", "[0.,0.,1]", "[0.,1,0.]", "[.63,.13,.94]",
            "[0,0,0]" );
    }

    $matlabPipe = "matlab_pipe_$randInt";

    $imageConversionLock = "$outputTemp/.lock";

    # Ensure that Matlab can write to the temp output folder
    my $mode = 0777;
    chmod $mode, $outputTemp;

    print("Input Folders: @inputDirs\n");
    print("Output Folder: $output \n");

    $outputIndex = "$outputTemp/plots.html";

    OutputWriter->writeNavPageStart("$outputTemp/$navigationPage");
    OutputWriter->writeIndex("$outputTemp/index.html");

    # Fork to make MATLAB faster
    my $pid = fork();

    if ( $pid == 0 )    # Child
    {
        $ENV{'DISPLAY'} = '';

        system("mkfifo $matlabPipe");
        system("$MATLAB <> $matlabPipe");

        system("rm $imageConversionLock");

        print("\nPlease wait while the remaining images are converted...\n");

        # Wait until all images are converted
        my @epsFiles     = <$outputTemp/*eps>;
        my @jpgFiles     = <$outputTemp/jpg/*jpg>;
        my $arraySizeEps = @epsFiles;
        my $arraySizeJpg = @jpgFiles;

        if ( $keepEps == 0 ) {

            # Block this thread from moving on until there are no more eps files
            while ( $arraySizeEps != 0 ) {
                @epsFiles     = <$outputTemp/*eps>;
                $arraySizeEps = @epsFiles;
                sleep(1);
            }
        }
        else {

           # Block this thread from moving on until the number of JPGs equal the
           # numper of EPS files.
            while ( $arraySizeEps != $arraySizeJpg ) {
                @epsFiles     = <$outputTemp/*eps>;
                @jpgFiles     = <$outputTemp/jpg/*jpg>;
                $arraySizeEps = @epsFiles;
                $arraySizeJpg = @jpgFiles;

                sleep(1);
            }
        }

        OutputWriter->writeFooter($outputIndex);

        cleanup();

        if ( $outputAsMaff == 0 ) {
            print(
"Done! To display the plots, open: \n$output/index.html \nin your web browser.\n"
            );
        }
        else {
            print(  "Done! To display the plots, open: \n" 
                  . $output
                  . ".maff\nin your maff compatible web browser.\n" );
        }

        exit(0);
    }
    else    # Parent
    {
        sleep(1);

        # Now fork to create images in the background. This should hopefully
        # speed things up a little
        my $convertPid = fork();

        if ( $convertPid == 0 )    # Image Conversion Child
        {
            system("touch $imageConversionLock");
            convertEps();
            exit(0);
        }
        else                       # Main thread
        {
            if ( $plotgenMode eq "splotgen" ) {
                OutputWriter->writeHeader( $outputIndex, "Splotgen" );

                # Only print the SAM_CLUBB Equiv. table if in nightly mode.
                if ( $nightly == 1 ) {
                    OutputWriter->writeSamSubHeader($outputIndex);
                }
                if ( $plotBudgets == 1 ) {
                    OutputWriter->writeSamBudgetSubHeader($outputIndex);
                }
            }
            elsif ( $plotgenMode eq "plotgen" ) {
                OutputWriter->writeHeader( $outputIndex, "Plotgen" );
            }
            elsif ( $plotgenMode eq "wrfgen" ) {
                OutputWriter->writeHeader( $outputIndex, "WRFGen" );
            }
            elsif ( $plotgenMode eq "camgen" ) {
                OutputWriter->writeHeader( $outputIndex, "CAMGen" );
            }
            elsif ( $plotgenMode eq "gfdlgen" ) {
                OutputWriter->writeHeader( $outputIndex, "GFDLGen" );

                #Only print the GFDL header if it is in nightly mode
                if ( $nightly == 1 ) {
                    OutputWriter->writeGfdlHeader($outputIndex);
                }
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
sub getCasePath() {
    my $casePath = $dirPrefix . "cases";

    if ( $plotgenMode eq "plotgen" ) {
        return "$casePath/clubb";
    }
    elsif ( $plotgenMode eq "splotgen" ) {
        return "$casePath/sam_clubb";
    }
    elsif ( $plotgenMode eq "wrfgen" ) {
        return "$casePath/wrf";
    }
    elsif ( $plotgenMode eq "camgen" ) {
        return "$casePath/cam";
    }
    elsif ( $plotgenMode eq "gfdlgen" ) {
        return "$casePath/gfdl";
    }
}

###############################################################################
# Runs all of the .case file in the cases folder.
# Arguments:
#   None.
###############################################################################
sub runCases() {

    # Counter used to place images
    my $count = 0;

    # Get the path to look for case files in
    my $casePath = getCasePath();

    my @cases = <$casePath/*.cas*>;

    # Loop through each .case file so the case can be plotted
    foreach my $file (@cases) {

        # The first part of this if statement will run all case files
        # if the nightly flag was passed in. The second part will only
        # run the case file if the nightly flag was not passed in and
        # the filename does not contain _nightly.
        # We need to do this because we have special condition case files.
        # For example, gabls2 plots special variables at night. Therefore,
        # there is a special gabls2_nightly.case file that only includes
        # those plots.
        if ( $nightly == 1 || ( $nightly == 0 && $file !~ m/_nightly/ ) ) {

            #print("Case file: $file\n");
            my $runCase = 'true';

# Read the case file. If there is an error, display it and continue without plotting it.
            if ( my $err = CaseReader->readCase($file) ) {
                system("$consoleOutput -w \"$err\"");
                $runCase = 'false';
            }
            else {

# Check to see if case was plotted already. This fixes the infinite loop problem
# when converting images. If the case was already plotted, do not do it again.
                foreach my $chkCase (@casesExecuted) {
                    if ( $chkCase eq $CASE::CASE{'name'} ) {
                        $runCase = 'false';
                    }
                }
            }

            if (   $runCase eq 'true'
                && dataExists($CASE::CASE)
                && ( $CASE::CASE{'enabled'} ne 'false' ) )
            {
                push( @casesExecuted, $CASE::CASE{'name'} );
                if (
                    ( $CASE::CASE{'type'} eq "budget" && $plotBudgets == 1 )
                    || (   $CASE::CASE{'type'} eq "morrbudget"
                        && $plotMorrBudgets == 1 )
                    || ( $CASE::CASE{'type'} eq "standard" )
                  )
                {

                    # Print the case title to the HTML page
                    OutputWriter->writeCaseTitle( $outputIndex,
                        $CASE::CASE{'headerText'} );
                    OutputWriter->writeNavPageCase(
                        "$outputTemp/$navigationPage", $CASE::CASE{'name'},
                        $CASE::CASE{'headerText'} );

                    # Print any additional text/html specified
                    if ( $nightly == 1 )    # If in nightly mode
                    {
                        my $nightlySubText =
                          $CASE::CASE{'nightlyOutput'}{'subText'};
                        my $nightlySubHtml =
                          $CASE::CASE{'nightlyOutput'}{'subHtml'};

        # Check to see if there was any additional text specified. If there was,
        # write it to the HTML file.
                        if ( $nightlySubText ne "" ) {
                            OutputWriter->writeSubHeader( $outputIndex,
                                $nightlySubText );
                        }

                        if ( $nightlySubHtml ne "" ) {
                            OutputWriter->writeSubHtml( $outputIndex,
                                $nightlySubHtml );
                        }
                    }
                    else    # If not in nightly mode
                    {
                        my $subText =
                          $CASE::CASE{'additionalOutput'}{'subText'};
                        my $subHtml =
                          $CASE::CASE{'additionalOutput'}{'subHtml'};

                        if ( $CASE::CASE{'type'} eq "morrbudget" ) {
                            OutputWriter->writeMorrBudgetSubHeader(
                                $outputIndex, $subText );
                        }
                        if ( $subText ne "" ) {
                            OutputWriter->writeSubHeader( $outputIndex,
                                $subText );
                        }

                        if ( $subHtml ne "" ) {
                            OutputWriter->writeSubHtml( $outputIndex,
                                $subHtml );
                        }

                    }
        # Check for includeSetup and print links for each input folder
                    if ( $includeSetup == 1 ) {
                        my $i = 1;
                        foreach my $inDir (@inputDirs) {
                            my $setupFileName = "${inDir}/$CASE::CASE{'name'}_setup.txt";
                            if ( -e $setupFileName ) {
                                OutputWriter->writeSetupLink( $outputIndex,
                                    $i, basename($inDir), $CASE::CASE{'name'} );
                                $i++;
                            }
                        }
                    }
                }

                # Check to see if this is a budget plot or standard plot
                if (   $CASE::CASE{'type'} eq "budget" && $plotBudgets
                    || $CASE::CASE{'type'} eq "morrbudget" && $plotMorrBudgets )
                {

                    #buildMatlabStringBudget($CASE::CASE, $count);

    # Add image file to HTML page
    #placeImages($CASE::CASE{'name'} . "_" . $count . "_budget", $plotCount, 0);

                    my $nightlyCase = 0;

                    if ( $nightly == 1 && $file =~ m/_nightly/ ) {
                        $nightlyCase = 1;
                    }

                    buildMatlabStringStd( $CASE::CASE, $nightlyCase );

                    # Add image file to HTML page
                    placeImages( $CASE::CASE{'name'}, $plotCount,
                        $nightlyCase );
                }
                if ( $CASE::CASE{'type'} eq "standard" ) {
                    my $nightlyCase = 0;

                    if ( $nightly == 1 && $file =~ m/_nightly/ ) {
                        $nightlyCase = 1;
                    }

                    buildMatlabStringStd( $CASE::CASE, $nightlyCase );

                    # Add image file to HTML page
                    placeImages( $CASE::CASE{'name'}, $plotCount,
                        $nightlyCase );
                }

                $count++;
            }
            else {

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
sub convertEps() {
    mkdir "$outputTemp/jpg" unless -d "$outputTemp/jpg";

    my @epsFiles     = <$outputTemp/*eps>;
    my @jpgFiles     = <$outputTemp/jpg/*jpg>;
    my $arraySizeEps = @epsFiles;
    my $arraySizeJpg = @jpgFiles;

   # This will keep calling convertEps() if the image conversion lock exists OR:
   #   If -e was not passed in:
   #       Until there are no more eps files left
   #   If -e was passed in:
   #       Until the number of eps files equals the number of jpg files
    if (   ( $keepEps == 0 && $arraySizeEps != 0 )
        || ( $keepEps == 1 && ( $arraySizeEps != $arraySizeJpg ) )
        || -e "$imageConversionLock" )
    {

        # Set the image scale if -q was not passed in
        if ( $highQuality == 0 ) {
            $DPI     = 120;
            $QUALITY = 80;
        }

        foreach my $eps (@epsFiles) {
            my $filename = basename($eps);

            if ( !-e "$outputTemp/jpg/$filename.jpg" ) {

                # First convert the image to png and trim all white space
                system(
"convert -density $DPI -colorspace RGB -trim $eps $outputTemp/jpg/$filename.png"
                );

                # Then convert (and scale if not in high quality mode) to jpg
                if ( $highQuality == 0 ) {
                    system(
"convert -geometry 324x312\\! -quality $QUALITY $outputTemp/jpg/$filename.png $outputTemp/jpg/$filename.jpg"
                    );
                }
                else {
                    system(
"convert -quality $QUALITY $outputTemp/jpg/$filename.png $outputTemp/jpg/$filename.jpg"
                    );
                }

                unlink("$outputTemp/jpg/$filename.png");
                if ( $keepEps == 0 ) {
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
sub placeImages() {
    my $caseName        = shift(@_);
    my $numImages       = shift(@_);
    my $nightlyCaseFile = shift(@_);

    OutputWriter->printDivCenter($outputIndex);

    for ( my $x = 0 ; $x < $numImages ; $x++ ) {
        if ( $nightlyCaseFile == 1 ) {
            OutputWriter->placeImage( $outputIndex,
                "jpg/$caseName" . "_nightly_" . "$x.eps.jpg" );
        }
        else {
            OutputWriter->placeImage( $outputIndex,
                "jpg/$caseName" . "_" . "$x.eps.jpg" );
        }
    }

    OutputWriter->printCloseDivCenter($outputIndex);
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
sub buildMatlabStringStd() {
    my $CASE            = shift(@_);
    my $nightlyCaseFile = shift(@_);

    $plotCount = 0;
    my $bool = 0;

    # Get Common Case information
    my $caseName;

    # Modify the case name if this is a nightly case
    if ( $nightlyCaseFile == 1 ) {
        $caseName = $CASE::CASE{'name'} . "_nightly";
    }
    else {
        $caseName = $CASE::CASE{'name'};
    }

    my $startHeight = $CASE::CASE{'startHeight'};
    my $endHeight   = $CASE::CASE{'endHeight'};

    # Get array of plots from case file
    my @plots;
    push( @plots, @{ $CASE::CASE{'plots'} } );

    # Get plots from .case file
    for ( my $count = 0 ; $count < @plots ; $count++ ) {

        # Adjust the widths if not using thinLines
        if ( $thinLines == 0 ) {
            if (   $CASE::CASE{'type'} eq "budget"
                || $CASE::CASE{'type'} eq "morrbudget" )
            {
                @lineWidths = @lineWidthsBudget;
            }
            else {
                @lineWidths = @lineWidthsNormal;
            }
        }

        # Counters for automatic lines
        $lineStyleCounter = 0;
        $lineColorCounter = 0;
        $lineWidthCounter = 0;
        my $fileCounter = 0;

        my $plotTitle = $plots[$count]{'plotTitle'};
        my $units     = $plots[$count]{'axisLabel'};
        my $type      = $plots[$count]{'type'};

# Define the startTime and endTime variables as the global variables in the case file.
        my $startTime = $CASE::CASE{'startTime'};
        my $endTime   = $CASE::CASE{'endTime'};

# Check to see if start and end time was specified for the specific plot. This is usually done for timeseries plots. If it isn't
# use case defined times.
        my $startTimeOverride = $plots[$count]{'startTime'};
        my $endTimeOverride   = $plots[$count]{'endTime'};

        if ( $startTimeOverride ne "" ) {
            $startTime = $startTimeOverride;
        }

        if ( $endTimeOverride ne "" ) {
            $endTime = $endTimeOverride;
        }

        my $matlabArgs =
"\'$caseName\', \'$CASE::CASE{'type'}\', \'$plotTitle\', $count, \'$type\', $startTime, $endTime, $startHeight, $endHeight, \'$units\', $randInt, $displayLegend";
        my $tempMatlabArgs = $matlabArgs;

        my @lines;
        push( @lines, @{ $plots[$count]{'lines'} } );
        for ( my $lineNum = 0 ; $lineNum < @lines ; $lineNum++ ) {
            my $name       = $lines[$lineNum]{'name'};
            my $expression = $lines[$lineNum]{'expression'};
            my $type       = $lines[$lineNum]{'type'};

            if ( $type eq "auto" ) {
                if ( $diffRun == 0 ) {
                    foreach (@inputDirs) {
                        my $file = "$_/$lines[$lineNum]{'filename'}";

                        if ( $dataFileType eq "netcdf" ) {

                            # Replace all .ctl file extensions with .nc
                            # By default, the case files have .ctl extensions
                            $file =~ s/.ctl/.nc/g;
                        }

                        if ( -e $file ) {
                            my $title;
                            my $folderName = basename($_);
                            $folderName =~
                              s/_/\\_/g;  # Replace all '_' with '\_' for MATLAB

                            if ( $name eq "auto" ) {
                                $title = $folderName;
                            }
                            else {
                                $title = $name;

                                # Replace any '{0}' with the folder name
                                $title =~ s/\{0\}/$folderName/;
                            }

                            my $lineWidth = $lineWidths[$lineWidthCounter];
                            my $lineStyle = $lineStyles[$lineStyleCounter];
                            my $lineColor = $lineColors[$lineColorCounter];

                            $matlabArgs =
"$matlabArgs, \'$file\', \'$file\', \'$expression\', \'$title\', $lineWidth, \'$lineStyle\', \'$lineColor\' ...\n";

            # Used for ensemble runs to make a group of plots all the same color
            # The five(5) is used to set the number of files per group
                            if ( $ensembleTuner == 1 ) {
                                if (   ( ( ( $fileCounter + 1 ) % 10 ) == 0 )
                                    && ( $fileCounter != 0 ) )
                                {
                                    incrementLineTypes();
                                }
                                $fileCounter++;
                            }
                            else {
                                incrementLineTypes();
                            }
                        }
                        else {
                            if ( $ensembleTuner == 1 ) {
                                if (   ( ( ( $fileCounter + 1 ) % 10 ) == 0 )
                                    && ( $fileCounter != 0 ) )
                                {
                                    incrementLineTypes();
                                }
                                $fileCounter++;
                            }
                        }
                    }
                }
                else    # This is a "difference run"
                {
                    my $file1 = "$inputDirs[0]/$lines[$lineNum]{'filename'}";
                    my $file2 = "$inputDirs[1]/$lines[$lineNum]{'filename'}";

#Ticket 543
#If a CaseName_lh_sfc file is found in the output, the lh_morr_rain_rate should be output
#for the Surface rainfall rate plot
                    my $lhSfcFile1 =
                      $inputDirs[0] . "/" . $caseName . "_lh_sfc.ctl";
                    my $lhSfcFile2 =
                      $inputDirs[1] . "/" . $caseName . "_lh_sfc.ctl";
                    if (   $plotTitle eq "Surface rainfall rate"
                        && -e $lhSfcFile1
                        && -e $lhSfcFile2 )
                    {
                        $file1      = $lhSfcFile1;
                        $file2      = $lhSfcFile2;
                        $expression = "lh_morr_rain_rate";
                    }

                    if ( $dataFileType eq "netcdf" ) {

                        # Replace all .ctl file extensions with .nc
                        # By default, the case files have .ctl extensions
                        $file1 =~ s/.ctl/.nc/g;
                        $file2 =~ s/.ctl/.nc/g;
                    }

                    if ( -e $file1 && -e $file2 ) {
                        my $title;
                        my $folderName1 = basename( $inputDirs[0] );
                        my $folderName2 = basename( $inputDirs[1] );
                        my $folderName  = $folderName2 . "-" . $folderName1;
                        $folderName =~
                          s/_/\\_/g;    # Replace all '_' with '\_' for MATLAB

                        if ( $name eq "auto" ) {
                            $title = $folderName;
                        }
                        else {
                            $title = $name;

                            # Replace any '{0}' with the folder name
                            $title =~ s/\{0\}/$folderName/;
                        }

                        my $lineWidth = $lineWidths[$lineWidthCounter];
                        my $lineStyle = $lineStyles[$lineStyleCounter];
                        my $lineColor = $lineColors[$lineColorCounter];
                        $matlabArgs =
"$matlabArgs, \'$file1\', \'$file2\', \'$expression\', \'$title\', $lineWidth, \'$lineStyle\', \'$lineColor\' ...\n";
                        incrementLineTypes();
                    }
                }
            }
            elsif (( $type eq "les" && $plotLes == 1 )
                || ( $type eq "dec17" && $plotDec )
                || ( $type eq "bestever" && $plotBest ) )
            {
                my $file = "$lines[$lineNum]{'filename'}";
                unless ( -e $file ) {

                    # See if the user specified a relative path.
                    # Try the first input directory.
                    $file = "$inputDirs[0]/$lines[$lineNum]{'filename'}";
                }
                if ( -e $file ) {
                    my $title = $name;

                    my $lineWidth = $lines[$lineNum]{'lineWidth'};
                    my $lineStyle = $lines[$lineNum]{'lineType'};
                    my $lineColor = $lines[$lineNum]{'lineColor'};

                    if ( $lineWidth eq "auto" ) {
                        $lineWidth = $lineWidths[$lineWidthCounter];
                    }
                    if ( $lineStyle eq "auto" ) {
                        $lineStyle = $lineStyles[$lineStyleCounter];
                    }
                    if ( $lineColor eq "auto" ) {
                        $lineColor = $lineColors[$lineColorCounter];
                    }

                    if (   ( "$lines[$lineNum]{'lineWidth'}" eq "auto" )
                        || ( "$lines[$lineNum]{'lineType'}"  eq "auto" )
                        || ( "$lines[$lineNum]{'lineColor'}" eq "auto" ) )
                    {
                        incrementLineTypes();
                    }

                    $matlabArgs =
"$matlabArgs, \'$file\', \'$file\', \'$expression\', \'$title\', $lineWidth, \'$lineStyle\', \'$lineColor\' ...\n";
                }
            }
        }

        # Check to see if there are lines to be plotted:
        if ( $matlabArgs eq $tempMatlabArgs ) {

            #print(STDERR "No valid data available to plot.\n");
        }
        else {
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
sub executeMatlab() {
    my $matlabArgs = shift(@_);

    my $args = "PlotCreator\"($matlabArgs)\"";

    #    print("\n$args\n\n");

    system("echo $args > $matlabPipe");
}

###############################################################################
# Changes to the next line style, width, and color.
# Arguments:
#   None.
###############################################################################
sub incrementLineTypes() {

    # Increment counters
    if ( $lineColorCounter + 1 >= @lineColors ) {
        $lineColorCounter = 0;

        if ( $lineStyleCounter + 1 >= @lineStyles ) {
            $lineStyleCounter = 0;

            if ( $lineWidthCounter + 1 >= @lineWidths ) {
                $lineWidthCounter = 0;
            }
            else {
                $lineWidthCounter++;
            }
        }
        else {
            $lineStyleCounter++;

            if ( $lineWidthCounter + 1 >= @lineWidths ) {
                $lineWidthCounter = 0;
            }
            else {
                $lineWidthCounter++;
            }
        }
    }
    else {
        $lineColorCounter++;

        if ( $lineStyleCounter + 1 >= @lineStyles ) {
            $lineStyleCounter = 0;

            if ( $lineWidthCounter + 1 >= @lineWidths ) {
                $lineWidthCounter = 0;
            }
            else {
                $lineWidthCounter++;
            }
        }
        else {
            $lineStyleCounter++;

            if ( $lineWidthCounter + 1 >= @lineWidths ) {
                $lineWidthCounter = 0;
            }
            else {
                $lineWidthCounter++;
            }
        }
    }
}

###############################################################################
# Does necessary cleanup code.
# Arguments:
#   None.
###############################################################################
sub cleanup() {

    # Copy temp. output folder to actual output location
    if ( $outputAsMaff == 1 ) {
        print("\nCompressing output to maff file:\n");
        my $outputName = $outputTemp;

     # Shortens the path name to just the lowest directory, if this is not done,
     # the UNIX zip command makes a zip including the absolute file structure.
        while ( $outputName =~ m/\/[^\/]*\// ) {
            $outputName =~ s/\/[^\/]*//;
        }
        $outputName = substr $outputName, 1;

      # Remove possible .maff from output name to avoid "*.maff.maff" file names
        if ( $output =~ m/.maff$/ ) {
            $output = substr $output, 0, -5;
        }

        system(
            "cd $outputTemp/.. && zip -r " . $output . ".maff $outputName/" );
    }
    else {
        dircopy( $outputTemp, $output );
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
sub dataExists() {
    my $CASE = shift(@_);

    # We need to figure out if any of the input directories contain
    # data to plot. I am assuming that input directories are of type
    # "auto" for now. If the type is auto, look in each input directory
    # to see if some data exists and return 1.

    # Get array of plots from case file
    my @plots;
    push( @plots, @{ $CASE::CASE{'plots'} } );

    for ( my $count = 0 ; $count < @plots ; $count++ ) {
        my @lines;
        push( @lines, @{ $plots[$count]{'lines'} } );

        for ( my $lineNum = 0 ; $lineNum < @lines ; $lineNum++ ) {
            my $lineType = $lines[$lineNum]{'type'};
            my $filename = $lines[$lineNum]{'filename'};

            if ( $dataFileType eq "netcdf" ) {

                # Replace all .ctl file extensions with .nc
                $filename =~ s/.ctl/.nc/g;
            }

            if ( $lineType eq 'auto' ) {
                foreach (@inputDirs) {
                    my @inputFiles = <$_/$filename*>;

                    if (@inputFiles) {
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
sub readArgs() {
    use Switch;

    my $numArgs = $#ARGV + 1;

    if ( $numArgs == 0 ) {
        main::HELP_MESSAGE();
    }

    my $option     = '';
    my $numOptions = 0;

# Loop through each agrument and see if it has a '-', then see if the remaining letters
# are a valid option
    foreach (@ARGV) {
        if ( ( substr $_, 0, 1 ) eq '-' ) {
            $option = substr $_, 1;
            $numOptions++;

            switch ($option) {
                case 'r' {
                    $overwrite = 1;
                }    # Option to overwrite file or folder if it already exsists
                case 'l' { $plotLes = 1; }    # Option to plot LES data
                case 'b' {
                    $plotBest = 1;
                }    # Option to plot Chris Golaz Best Ever Data
                case 'd' { $plotDec = 1; }    # Option to plot HOC Dec. 17
                case 'a' {
                    $plotLes  = 1;
                    $plotBest = 1;
                    $plotDec  = 1;
                }                             # Plot LES, CGBE, and HOC Dec. 17
                case 'n' { $nightly      = 1; }          # Run in nightly mode
                case 'q' { $highQuality  = 1; }          # Run in nightly mode
                case 'e' { $keepEps      = 1; }          # Keep EPS files
                case 'g' { $dataFileType = "grads"; }
                case 't' { $dataFileType = "netcdf"; }
                case "m" {                               # Output as maff file
                    if ( $option eq
                        "e" ) # It's inconvenient for maff if eps files are kept
                    {
                        print(
"Argument conflict: Please do not simultaneously choose to 
           save .eps files (-e) and create a .maff file (-m)\n"
                        );
                        exit(1);
                    }
                    $outputAsMaff = 1;
                }
                case 'thin'     { $thinLines       = 1; }
                case 'dfrnce'   { $diffRun         = 1; }
                case 'ensemble' { $ensembleTuner   = 1; }
                case 'nolegend' { $displayLegend   = 0; }
                case 'nosetup'  { $includeSetup    = 0; }
                case 's'        { $plotgenMode     = "splotgen"; }
                case 'c'        { $plotgenMode     = "plotgen"; }
                case 'w'        { $plotgenMode     = "wrfgen"; }
                case 'cam'      { $plotgenMode     = "camgen"; }
                case 'gfdl'     { $plotgenMode     = "gfdlgen"; }
                case 'bu'       { $plotBudgets     = 1; }
                case 'bumorr'   { $plotMorrBudgets = 1; }
                case 'h' { main::HELP_MESSAGE(); }    #Prints help message
                else     { main::HELP_MESSAGE(); }
            }
        }
    }

    # Create a second array that only has the folders in it
    my @fileArgs;
    my $fileNum = $numArgs - $numOptions;
    my $i;
    for ( $i = 0 ; $i < $fileNum ; $i++ ) {
        @fileArgs[$i] = @ARGV[$numOptions];
        $numOptions++;
    }

    if ( $ensembleTuner == 1 ) {
        $thinLines     = 1;
        $displayLegend = 0;
    }

    my $currentCount = 0;

    # Reset the number of arguments
    $numArgs = $#ARGV + 1;

    # Parse any additional arguments
    foreach (@fileArgs) {

        # If the argument does not start with '-' and if $output was not set
        if ( !$output ) {
            my $currentDir = abs_path($_);

            if ( ( ( $fileNum - 1 ) - $currentCount ) > 0 ) {
                if ( -d $currentDir ) {
                    push( @inputDirs, $currentDir );
                }
                else {
                    print("The input folder $currentDir does not exist.\n");
                    exit(1);
                }
            }
            else {
                if ( $_ =~ m/\/$/ ) {
                    $currentDir = abs_path( substr( $_, 0, length($_) - 1 ) );
                }

                $output     = $currentDir;
                $outputTemp = "/tmp/output" . "_" . "$randInt";
            }

            $currentCount++;
        }
    }

    # If there were no input directories passed in, print the help and exit
    if ( @inputDirs == 0 ) {
        main::HELP_MESSAGE();
    }

    # Finally, check to see if the output folder exists. If it does, and
    # '-r' was not passed in, exit. Otherwise, create it.

    if ( $outputAsMaff == 0 && -e $output && $overwrite == 0 ) {
        system(
"$consoleOutput -s \"Output folder already exists. To overwrite, use the -r option.\""
        );
        exit(1);
    }
    elsif ( $outputAsMaff == 1 && -e $output . ".maff" && $overwrite == 0 ) {
        system(
"$consoleOutput -s \"maff file already exists. To overwrite, use the -r option.\""
        );
        exit(1);
    }
    else {
        if ( $outputAsMaff == 0 ) {
            rmtree($output);
            mkdir $output;
        }
        mkdir $outputTemp;

        # includeSetup only works in clubb standalone.
        if ( $plotgenMode ne "plotgen" ) {
            $includeSetup = 0;
        }

        # Copies all (case)_setup.txt files to outputTemp folder.
        if ( $includeSetup == 1 ) {
            my $i = 1;
            foreach my $inputD ( @inputDirs ) {
                mkpath("${outputTemp}/setup/${i}");
                my @files = glob("$inputD/*_setup.txt");
                for my $file (@files) {
                    my $fileOut = basename($file);
                    copy("$file", "${outputTemp}/setup/$i/${fileOut}");
                }
            $i++;
            }
        }

        my $plotgenDirectory = $dirPrefix;

       # In case we aren't already in the directory that plotgen.pl is, goto it.
        chdir($plotgenDirectory);
    }
}

###############################################################################
# Prints the help message
# Arguments:
#   None.
###############################################################################
sub main::HELP_MESSAGE() {
    print("Usage: plotgen [OPTION]... INPUT... OUTPUT\n");
    print("  -c\tPlot CLUBB cases [DEFAULT] (equiv to plotgen)\n");
    print("  -s\tPlot SAM_CLUBB cases (equiv to splotgen)\n");
    print("  -w\tPlot WRF_CLUBB cases\n");
    print("  -cam\tPlot CAM cases\n");
    print("  -gfdl\tPlot GFDL cases\n");
    print("  -r\tIf the output folder already exists, replace the contents\n");
    print("  -l\tPlot LES data for comparison.\n");
    print("  -b\tPlot HOC Best Ever data for comparison.\n");
    print("  -d\tPlot HOC 12/17/2005 data for comparison.\n");
    print(
"  -a\tSame as -lbd. Plots LES, Best Ever, and 12/17/2005 data for comparison.\n"
    );
    print("  -n\tRuns in nightly mode.\n");
    print("  -q\tOutputs high quality images (does not auto scale).\n");
    print("  -e\tDoes not delete EPS images after conversion.\n");
    print("  -m\tOutputs plots compressed inside a .maff directory.\n");
    print("  -g\tUses GrADS data files. [DEFAULT]\n");
    print("  -t\tUses NetCDF data files.\n");
    print("  -thin\tUses thin solid lines\n");
    print("  -dfrnce\tPerforms a 'difference' plot of two output folders\n");
    print("  -nolegend\tPlot without legends\n");
    print("  -nosetup\tPlot without including (case)_setup.txt\n");
    print("  -ensemble\tUsed for plotting ensemble tuner runs\n");
    print("  -bu\tUsed to plot standard budget plots\n");
    print("  -bumorr\tUsed to plot Morrison budget plots\n");
    print("  -h\tPrints this help message.\n");
    print("Each option must be seperate, eg -r -a not -ra\n");
    exit(0);
}

###############################################################################
# Prints the version number of Plotgen.
# Arguments:
#   None.
###############################################################################
sub main::VERSION_MESSAGE() {
    print("Plotgen version $VERSION, Copyright (C) 2013 Larson Group.\n");
}
