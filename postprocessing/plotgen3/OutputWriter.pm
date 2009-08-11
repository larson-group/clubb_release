#!/usr/bin/perl

package OutputWriter;

use strict;

use warnings;

our $VERSION = '1.00';

###############################################################################
# Writes a case title to the HTML file
###############################################################################
sub writeCaseTitle()
{
	shift(@_);
	my $fh = shift(@_);

	open(FILE, ">> $fh");

	my($CASE) = shift;
	print FILE <<"EOF";
	<a name="$CASE"></a>
	<div align="CENTER">
		<font size="+2">
		<font color="#0000ff"> <a href="#$CASE">$CASE</a> </font> </font>
	</div>
EOF
	close(FILE);
}

###############################################################################
# Writes sub text under the header
###############################################################################
sub writeSubHeader()
{
	shift(@_);
	my $fh = shift(@_);
	open(FILE, ">> $fh");

	my($TEXT) = shift;

	print FILE <<"EOF";
	<div align="CENTER">
	        <b><font size="-1" color="#430e9a"> $TEXT </font></b>
	</div>
EOF
	close(FILE);
}

###############################################################################
# Writes HTML under the case header
###############################################################################
sub writeSubHtml()
{
	shift(@_);
	my $fh = shift(@_);
	open(FILE, ">> $fh");

	my($HTML) = shift;

	print FILE <<"EOF";
	<div align="CENTER">
		$HTML
	</div>
EOF
	close(FILE);
}

###############################################################################
# Writes the HTML header information
###############################################################################
sub writeHeader()
{
	
	(my $sec, my $min, my $hour, my $mday, my $mon, my $year, my$wday,
	my $yday, my $isdst)=localtime(time);

	my $currYear = $year + 1900;
	my $currMon = $mon + 1;
	my $currDay = $mday;
	my $currHour = $hour;
	my $currMin = $min;
	my $currSec = $sec;

	shift(@_);
	my $fh = shift(@_);
	open(FILE, "> $fh");

	print FILE <<"EOF";
<html>
<head>
	<title>Plotgen</title>
	<p><div align="CENTER">
		<font size="+3" color="#811212">Plotgen</font>
		<font size="-1"><br/>$currMon/$currDay/$currYear</font></p>
	</div></p>

</head>
<body>
EOF
	close(FILE);
}

###############################################################################
# Writes the HTML footer information and closes the file.
###############################################################################
sub writeFooter()
{
	shift(@_);
	my $fh = shift(@_);
	open(FILE, ">> $fh");

	print FILE <<"EOF";
	<br /> <br /> <br /> <br />
	<hr noshade size=5 width=70%>
	<div align="CENTER">
		<font size="-2">
		Copyright (c) 2009 Larson Group. All rights reserved. 
		</font>
	</div>
</body>
</html>
EOF
	close(FILE);
}

###############################################################################
# Inserts an image
###############################################################################
sub placeImage()
{
	shift(@_);
	my $fh = shift(@_);
	open(FILE, ">> $fh");

	my($img) = shift;

	print FILE <<"EOF";
		<img width="328" height="312" align="BOTTOM" border="0" src="$img" alt="$img" />
EOF
	close(FILE);
}

sub printDivCenter()
{
	shift(@_);
	my $fh = shift(@_);
	open(FILE, ">> $fh");

	print FILE <<"EOF";
	<div align="CENTER">
EOF
	close(FILE);
}

sub printCloseDivCenter()
{
	shift(@_);
	my $fh = shift(@_);
	open(FILE, ">> $fh");

	print FILE <<"EOF";
	</div>
EOF
	close(FILE);
}

1;
