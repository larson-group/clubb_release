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
	<div align="CENTER">
		<b><font size="+2">
		<font color="#0000ff"> $CASE </font> </font></b>
	</div>
EOF
	close(FILE);
}

###############################################################################
# Writes the HTML header information
###############################################################################
sub writeHeader()
{
	shift(@_);
	my $fh = shift(@_);
	open(FILE, "> $fh");

	print FILE <<"EOF";
<html>
<head>
	<title>Plotgen</title>
	<div align="CENTER">
		<h1><font color="#811212">Plotgen</font></h1>
	</div>

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
</htlm>
EOF
	close(FILE);
}

1;
