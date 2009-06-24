#!/usr/bin/perl

###########################################################################
# Plotgen v3.0
#
#  
###########################################################################

package plotgen;

use strict;
use CaseReader;
use Cwd;
use Switch;

my $plotgenVersion = 3.0;
my $plotgenDir = cwd;

readArgs();

# Open a new index.html file.
open(INDEXFILE, ">index.html");

writeHeader();
runCases();
writeFooter();

close(INDEXFILE);

###############################################################################
# Runs all of the .case file in the cases folder.
###############################################################################
sub runCases()
{
	# Loop through each .case file so the case can be plotted
	my @cases = <cases/*>;
	foreach my $file (@cases) 
	{
		# Read the case file. If there is an error, exit.
		if (my $err = CaseReader::readCase($file))
		{
	    		print(STDERR $err, "\n");
	    		exit(1);
	    	}
		
		# Print the case title to the HTML page
		writeCaseTitle($CASE::CASE{'name'});
	    	print("Successfully opened: " . $file . "\n");
	    	print("Plotting case: " . $CASE::CASE{'name'} . "\n");
	    	print("nightly_output: " . $CASE::CASE{'test'}[2] . "\n");
	}
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

	foreach my $argnum (0 .. $#ARGV) 
	{
		switch(@ARGV[$argnum])
		{
			case "--nightly"
			{
				print("Running in nightly mode" . "\n");
			}
			else
			{
				printHelp();
			}
		}
	}
}

###############################################################################
# Prints the help message
###############################################################################
sub printHelp()
{
	print("Plotgen v. " . $plotgenVersion . "\n");
	print("Usage: plotgen [ options ... ]" . "\n");
	exit(0);
}

###############################################################################
# Writes a case title to the HTML file
###############################################################################
sub writeCaseTitle()
{
	my($CASE) = @_;
	print INDEXFILE <<"EOF";
	<div align="CENTER">
		<b><font size="+2">
		<font color="#0000ff"> $CASE </font> </font></b>
	</div>
EOF
}

###############################################################################
# Writes the HTML header information
###############################################################################
sub writeHeader()
{
	print INDEXFILE <<"EOF";
<html>
<head>
	<title>Plotgen $plotgenVersion</title>
	<div align="CENTER">
		<h1>Plotgen</h1>
	</div>

</head>
<body>
EOF
}

###############################################################################
# Writes the HTML footer information and closes the file.
###############################################################################
sub writeFooter()
{
	print INDEXFILE <<"EOF";
	<br /> <br /> <br /> <br />
	<hr noshade size=5 width=70%>
	<div align="CENTER">
		<font size="-2">
		Copyright (c) 2009 Larson Group. All rights reserved. 
		</font>
	</div>
</body>
</htlm>
EOF
}
