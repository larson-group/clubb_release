#!/usr/bin/perl

package OutputWriter;

use strict;

use warnings;

our $VERSION = '1.00';

###############################################################################
# Writes the beginning of the navigation page HTML. This should be called 
# before the cases are written to the navigation page.
###############################################################################
sub writeNavPageStart()
{
    shift(@_);
    my $fh = shift(@_);

    open(FILE, "> $fh");
    print FILE <<"EOF";
    <html>
        <h4>Cases:</h4>
EOF
}

###############################################################################
# Writes a case to the navigation page
###############################################################################
sub writeNavPageCase()
{
    shift(@_);
    my $fh = shift(@_);
    my $case = shift(@_);
    my $link = shift(@_);

    open(FILE, ">> $fh");
    print FILE <<"EOF";
        <a href="plots.html#$link" target="plots">$case</a><br/>
EOF
}

###############################################################################
# Writes HTML to the end of the navigation page. This should be called after
# all cases were added.
###############################################################################
sub writeNavPageClose()
{
    shift(@_);
    my $fh = shift(@_);

    open(FILE, ">> $fh");
    print FILE <<"EOF";
    </html>
EOF
}

###############################################################################
# Writes the index page that contains the two frames.
###############################################################################
sub writeIndex()
{
    shift(@_);
    my $fh = shift(@_);

    open(FILE, "> $fh");
    print FILE <<"EOF";
    <html>
        <title>Plotgen</title>
        <frameset cols="180,*">
            <frame src="navigation.html" frameborder="0" name="nav">
            <frame src="plots.html" frameborder="0" name="plots">
        </frameset>
    </html>
EOF
}

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
# Writes the HTML header information for CLUBB
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
    my $mode = shift(@_);

	open(FILE, "> $fh");

	print FILE <<"EOF";
<html>
<head>
	<title>$mode</title>
	<p><div align="CENTER">
		<font size="+3" color="#811212">$mode</font>
		<font size="-1"><br/>$currMon/$currDay/$currYear</font></p>
	</div></p>

</head>
<body>
EOF
	close(FILE);
}

###############################################################################
# Writes the SAM_CLUBB Variable Equivalence Table
###############################################################################
sub writeSamSubHeader()
{
    shift(@_);
    my $fh = shift(@_);

    OutputWriter->writeSubHeader($fh, "2D SAM_CLUBB runs use a 64-km horizontal domain and a 10-s timestep, with CLUBB called every 6th SAM timestep. 3D SAM_CLUBB runs use 4x4 columns with a 16-km horizontal grid spacing and a 10-s timestep with CLUBB called every 6th SAM timestep. All SAM_CLUBB runs except LBA use Morrison microphysics. CLUBB standalone runs use a 10-s timestep and the Morrison microphysics.");
    OutputWriter->writeSubHeader($fh, "When two variables are listed, the first variable is the SAM-CLUBB variable and the second is the SAM-Standalone variable.");
     
    my $text = <<HTML;
    <br />
    <DIV ALIGN="CENTER"><TABLE CELLPADDING=3 BORDER="1">
    <TR>
        <TD ALIGN="CENTER" COLSPAN=11><B>Variable Equivalence Table</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER"><B>CLUBB</B></TD>
        <TD ALIGN="CENTER"><B>SAM CLUBB</B></TD>
        <TD ALIGN="CENTER"><B>SAM Standalone</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">wpthlp</TD>
        <TD ALIGN="CENTER">wpthlp+tlflux</TD>
        <TD ALIGN="CENTER">tlflux</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">wprtp</TD>
        <TD ALIGN="CENTER">wprtp+qtflux</TD>
        <TD ALIGN="CENTER">qtflux</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">thlp2</TD>
        <TD ALIGN="CENTER">thlp2+tl2</TD>
        <TD ALIGN="CENTER">tl2</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">rtp2</TD>
        <TD ALIGN="CENTER">rtp2 + qt2</TD>
        <TD ALIGN="CENTER">qt2</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">upwp</TD>
        <TD ALIGN="CENTER">uw + upwp</TD>
        <TD ALIGN="CENTER">uw</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">vpwp</TD>
        <TD ALIGN="CENTER">vw + vpwp</TD>
        <TD ALIGN="CENTER">vw</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">up2</TD>
        <TD ALIGN="CENTER">up2 + u2</TD>
        <TD ALIGN="CENTER">u2</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">vp2</TD>
        <TD ALIGN="CENTER">vp2 + v2</TD>
        <TD ALIGN="CENTER">v2</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">wp2</TD>
        <TD ALIGN="CENTER">wp2 + w2</TD>
        <TD ALIGN="CENTER">w2</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">wp3</TD>
        <TD ALIGN="CENTER">wp3 + w3</TD>
        <TD ALIGN="CENTER">w3</TD>
    </TR>
</TABLE>
<br />
</DIV>
HTML
    
    OutputWriter->writeSubHtml($fh, $text);
}

###############################################################################
# Writes the HTML header for WRF
###############################################################################
sub writeWrfHeader()
{
     shift(@_);
     my $fh = shift(@_);

    OutputWriter->writeSubHeader($fh, "WRF_CLUBB runs use a 100-km horizontal grid spacing, 3x3 grid columns, and a 60-s timestep.");
    
    my $text = <<HTML;
<br />
<br />
HTML
    
    OutputWriter->writeSubHtml($fh, $text);
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

    my $imgWidth = 324;
    my $imgHeight = 312;

	print FILE <<"EOF";
		<img width="$imgWidth" height="$imgHeight" align="BOTTOM" border="0" style="padding: 5px;" src="$img" alt="$img" />
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
