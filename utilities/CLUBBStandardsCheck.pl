#!/usr/bin/perl

#$Id: CLUBBStandardsCheck.pl,v 1.12 2008-08-01 14:59:29 faschinj Exp $

##################################################################
#	
#	Program: CLUBBStandardsCheck.pl
#	Purpose: This program checks valid Fortran 95 files for:
#
#		(1) Missing implicit none statements.  
#		Warns if "implicit note" is missing from programs, 
#		modules, subroutines, or functions.
#
#		(2) Missing "only" lists in use statements.
#		Displays the line and number that an unrestricted 
#		"use" statement occurs on.  
#		
#		(3) Default Private declarations in modules
#		Warns if "private" is missing from modules.
#
#		(4) $ Id $ comment tags at the top of the file.
#		Warns if the file does not contain one.
#
#		(5) Lines that are longer than a specified size.
#
#		(6) Use of deprecated keywords such as .ge.
#
#               This perl script assumes that 
#               the Fortran code compiles!!
#
##################################################################

# Enforce strict checking
use strict;

# Includes
use Getopt::Long;
use Thread;

# Global Variables

# Print verbose messages?
our $verbose = 0;

# Name of the program
our $programName = "CLUBBStandardsCheck.pl";

# Line Separator
our $horizontal = "--------------------------------------------------------------------------------\n";

# Regular Expressions for Fortran statements
# See http://perldoc.perl.org/perlre.html for more information

our $implicitNoneRegEx = qr/^				# Bind to beginning of line
	      		\s*?				# Zero or more spaces before statement
			\b				# Bind to front of statement
			implicit\snone			# Implicit None
			\b				# Bind to end of statement
			\s*?				# Zero or more spaces after statement
			/ix;				# Whole expression is case insensitive

our $functionRegEx = 	qr/^				# Bind to beginning of line
        		\s*? 				# Zero or more spaces before statement
          	 	(real|double\sprecision|complex|
			logical|character|integer|pure|
			elemental|recursive|\s*?)*?	# Zero or more specifications.
							# Note the expression does not fit	
							# Fortran syntax exactly but we can
							# rely on the compiler to verify 
							# appropriate uses of specifications.
	        	\s?				# Zero or One space before statement
                	\b				# Bind to front of statement
			function			# Function
			\b				# Bind to back of statement
			/ix;				# Whole expression is case insensitive

our $subroutineRegEx = 	qr/^				# Bind to beginning of line
	      		\s*?				# Zero or more spaces before statement
	           	(pure|elemental|recursive|\s*?)*?	# Zero or more Specifications
		        \s?				# Zero or One Space before statement
			\b				# Bind to front of statement
			subroutine			# Subroutine
			\b				# Bind to end of statement
			/ix;				# Whole expression is case insensitive

our $moduleRegEx =	qr/^				# Bind to beginning of line
	      		\s*?				# Zero or more spaces before statement
			\b				# Bind to front of statement
			module				# Module
			\b				# Bind to end of statement
			\s+?				# One or more spaces
			\w+?				# One or more word characters 
			\s*?				# Zero or more spaces after statement
			(!.*?)*?			# After statement, accept comments and nothing else
			$				# Bind to end of line
			/ix;				# Whole expression is case insensitive
		
our $programRegEx =	qr/^				# Bind to beginning of line
	      		\s*?				# Zero or more spaces before statement
			\b				# Bind to beginning of statement
			program				# Program
			\b				# Bind to end of statement
			/ix;				# Whole expression is case insensitve

our $useRegEx =		qr/^				# Bind to beginning of line 
			\s*?				# Zero or more spaces before statement
			use				# Use
			\s+?				# One or more spaces after statement
			\w+?				# One or more word characters
			\s*?				# Zero or more spaces
			(!.*?)*?			# After statement, accept comments and nothing else
			$				# Bind to end of line 
			/ix;				# Whole expression is case insensitive

our $privateRegEx =	qr/^				# Bind to beginning of line
	      		\s*?				# Zero or more spaces before statement
			\bprivate\b			# Private
			\s*? 				# Zero or more spaces after statement
			(!.*?)*?			# Accepts only comments after statement
			$				# Bind to end of line
			/ix;				# Whole expression is case insensitive

our $IdTagRegEx =	qr/^				# Bind to beginning of line
			\s*?				# Zero or More spaces
			!				# Comment Character
			\s*?				# Zero or More spaces
			\$Id				# $Id
			:?                              # Zero or one colon
			\s?                             # Zero or one spaces
			((\w|\.)*)                      # Capture Zero or more alphanumeric characters
			.*				# Zero or More of any character
			\$				# Ending Dollar Sign
			\s*?				# Zero or more spaces
			$				# Bind to end of line
			/ix;				# Whole Expression is case insensitve

our $deprecatedRegEx = qr/
			(
			\.le\.
			|
			\.ge\.
			|
			\.lt\.
			|
			\.gt\.
			|
			\.eq\.
			|
			\.ne\.
			)
			/ix;

# Captures verbose command switch.
GetOptions ('v|verbose' => \$verbose);

# Test to see if at least one argument is present.
if(@ARGV < 1 )
{
	warn "Not enough arguments"	
}
#### BEGIN MAIN PROGRAM ####
else{
	warn "CLUBBStandardsCheck.pl has begun.\n";

	# Declare Local Variables
	our $thr;
	our $file;
	
	# For Every File
	foreach $file (@ARGV)
	{
		$thr = new Thread \&fileThread, $file, $thr;
	}
	
	while ( my(@list)=threads->list()) {
		grep { $_->join } @list;
	};

	warn "CLUBBStandardsCheck.pl has finished.\n";
}
#### END MAIN PROGRAM ####

sub fileThread
{
	my($file) = shift(@_);
		
	my(@input);

	# Open the file
	open FILE, $file or die "Bad Filename";
	
	# Store the lines to an array
	@input = <FILE>;

	# Check for missing implicit nones
	if( ! &implicitCheck( $verbose, @input ) )
	{
		warn "$file\n";	
	        warn $horizontal;

	}

	# Check for use statements without only
	if( ! &useCheck( $verbose, @input ) )
	{
		warn "$file\n";
	        warn $horizontal;
	}

	# Check for default private statements
	if( ! &privateCheck( $verbose, @input ) )
	{
		warn "$file\n";
		warn $horizontal;
	}
		
	# Check for missing $ Id $ tags
	if( ! &idCheck( $verbose, $file, @input ) )
	{
		warn "$file\n";
		warn $horizontal;
	}

	# Check for use of forbidden keywords
	if( ! &deprecateCheck( $verbose, @input ) )
	{
		warn "$file\n";
		warn $horizontal;
	}

	# Check for long lines
	if( ! &lineCheck( $verbose, @input ) )
	{
		warn "$file\n";
		warn $horizontal;
	}

	# Close File		
	close FILE;
	
}

#####################################################################
sub implicitCheck
#
#     &implicitCheck( $verbose, @input ) 
#
#     Description: This subroutine verifies that the number of program 
#     block and implicit none statements in a FORTRAN file match. 
#
#        This subroutine works by testing whether
#           ( number of programs,modules,subroutines, and functions )  
#                     ==  ( number of 'implicit none' statements ).
#     
#     Arguments:
#     	Sverbose - Prints verbose messages when true.
#     	@input   - Lines of a Fortran source file.
#
#     Returns
#     	True if the file has no missing implicit none statements.
#####################################################################
{
	# Grab first argument
	my( $verbose ) = shift( @_ );

	# Grab second argument
	my( @input ) = @_;

	# Declare Local Variables
	my( $line, $implicitCount, $statementCount, $result );

	# Number of instances of implicit none.
	$implicitCount = 0;

	# Number of instances of statements requiring implicit none
	$statementCount = 0;

	if( $verbose  )
	{
		print "----------------Implicit None Test----------------\n";
	}

	# Parse every line for a program block or implicit none statements.	
	foreach $line (@input)
	{
		# Is it a function declaration?
	        if( $line =~ $functionRegEx )		
		{
			# Print if verbose
			if( $verbose )
			{	
				print "$programName comment: Found Function\n$line";
			}
			$statementCount++;
		}
		# Is it a subroutine declaration?	
		elsif($line =~ $subroutineRegEx)
		{	
			# Print if verbose
			if( $verbose )
			{
				print "$programName comment: Found Subroutine\n$line";
			}
			$statementCount++;
			
		}
		# Is it a module declaration?
		elsif($line =~ $moduleRegEx)           	
		{
			# Print if verbose	
			if( $verbose )
			{
				print "$programName comment: Found Module\n$line";
			}
			$statementCount++;
		}
		# Is it a program declaration?
		elsif( $line =~ $programRegEx )
		{	
			if( $verbose )
			{
				print "$programName comment: Found Program\n$line";
			}
			$statementCount++;
		}
		# Is it an implicit none declaration?	
		elsif( $line =~ $implicitNoneRegEx )
		{
			# Print if verbose	
			if($verbose)
			{
				print "$programName comment: Implicit none found\n$line";
			}
			$implicitCount++;
		}
	}

	# Check to see if the statements match up
	if($implicitCount == $statementCount)
	{
		$result = 1; # Passed check	
	}
	else
	{
	
		warn $horizontal;
		warn "$programName WARNING: Missing 'implicit none' statements. 'Implicit None' check FAILED! \n";  
		warn "$programName Add a line containing 'implicit none' to each program, module, subroutine, and function.\n";
		$result = 0; # Failed check
	}
	return $result;
}
####################################################################
sub useCheck
#	
#	&useCheck( $verbose, @input );
#
#     Description: This subroutine verifies that the file does not contain
#     any "use" statements without being followed by an "only" restriction.
#
#     	Basically, this checks for any line of the form
#       	use module_name  ! Comment
#     	and if so, returns an error because there is no "only" statement.
#     
#     Arguments:
#     	$verbose - Prints verbose messages when true.
#     	@input   - Lines of a Fortran source file.
#
#     Returns
#     	True if the file has no unrestricted "use" statements.
####################################################################
{
	# Grab first argument, assign to "verbose".
        my( $verbose ) = shift(@_);

	# Grab second argument, assign to "input".
	my( @input ) = @_;

	# Local variables
	my( $line, $result, $lineNumber,@warnings );

	$result = 1; # Assume file passed use test unless proven wrong below

	# Line number	
	$lineNumber = 1;

	if( $verbose )
	{
		print "---------------- Use Check ----------------\n";
	}

	# Check whether (true if) each line matches the format described in the comments below.
	foreach $line ( @input )
	{
		if( $line =~ $useRegEx )
		{       # If the above is true, execute the statements below (i.e. return an error):
			push(@warnings,"$programName WARNING: 'use' statement w/o 'only' found in the following line:\n");
			push(@warnings, "$lineNumber : $line");
			$result = 0; # Failed check
		}
		$lineNumber++;
	}

	if ( ! $result )
	{
		push (@warnings, "$programName WARNING: Use check FAILED!\n");
		push (@warnings, "$programName WARNING: Check that comma is on same line as 'use', as CLUBB requires.\n");

		warn $horizontal;

		warn @warnings;
	}
	return $result;
}

#####################################################################
sub privateCheck
#
#	&privateCheck( $verbose, @input );
#
#     Description: This subroutine verifies that a file's modules
#     have a corresponding private statement to set the default
#     scope of the module to private.
#
#     This subroutine works by testing whether
#     	( number of module statements )  ==  ( number of private statements ).
#     
#     Arguments:
#     	Sverbose - Prints verbose messages when true.
#     	@input   - Lines of a Fortran source file.
#
#     Returns
#     	True if the module is set to private default scope.  
####################################################################
{
	# Grab first argument
	my( $verbose ) = shift(@_);

	# Grab the second argument
	my( @input ) = @_;

	# Local variables
	my( $line, $result, $moduleCount, $privateCount );

	# Print if verbose
	if( $verbose )
	{
		print "---------------- Private Check  ----------------\n";
	}

	# For every line of the file
	foreach $line (	@input )
	{
		# Does the line contain a module declaration?
		if( $line =~ $moduleRegEx )                           
		{
			# Print if verbose
			if($verbose)
			{
				print "$programName comment: Module found\n $line";
			}
			$moduleCount++;	
		}
		# Does the line contain a private declaration?
		elsif($line =~ $privateRegEx )
		{	
			# Print if verbose
			if($verbose)
			{
				print "$programName comment: Private found\n$line";
			}
			$privateCount++;
		}		
	}

	if($privateCount == $moduleCount )
	{
		$result = 1; # Check passed
	}
	else
	{
		warn $horizontal;
		warn "$programName WARNING: Number of \"private\" statements does not not match number of modules.\n";
		warn "$programName WARNING: Private Test failed!\n";
		$result = 0; # Check failed
	}	

	return $result;
}
####################################################################
sub idCheck
#
#     &idCheck( $verbose, $filepath, @input ) 
#
#     Description: This subroutine verifies that an $ Id $ comment 
#     exists somewhere in the file. 
#
#        This subroutine works by testing each line for the presence
#        of a line containing a $ Id $.
#     
#     Arguments:
#     	Sverbose - Prints verbose messages when true.
#     	$filepath - Used to determine if the Id tag matches the filename.
#     	@input   - Lines of a Fortran source file.
#
#     Returns
#     	True if the file has an $ Id $ comment in it.
#####################################################################
{
	# Grab first argument
	my( $verbose ) = shift(@_);

        # Grab the filepath
	my( $filepath ) = shift(@_);

	# Grab Second argument
	my( @input ) = @_;

	# Local Variables
	my( $line, $result, @split_path, $filename, $id_filename, @warnings);

	# Default initialization to false
	$result = 0;
        
	# Getting the name of the file
	@split_path = split(/\//,$filepath);
	$filename = @split_path[-1];

	# Print if verbose
	if( $verbose )
	{
		print "---------------- \$Id\$ Check  ----------------\n";	
	}

	# For every line of the file
	foreach $line (@input)
	{
		# If it contains an $ Id $ comment
		if( $line =~ $IdTagRegEx )
		{
			# Store the filename from the Id tag
		        $id_filename = $1;
		        if( ! ($filename =~ m/$id_filename/) )
			{
		           push(@warnings,"$programName CAUTION: \"$filename\" does not match \"$id_filename\" in Id tag.\n");
			   
			}	
			if( $verbose )
			{
				print "$programName comment: Id tag found\n$line";
			}
			$result = 1; # Remember an $ Id $ tag was found
		}
	}

	# If no $ Id $ tags were found
	if( ! $result )
	{
		push(@warnings, "$programName WARNING: Missing \$Id\$ Tag. \$Id\$ check FAILED!\n");
		push(@warnings, "Add ! \$Id\$ to top of file.\n");
	}

	if($#warnings > -1)
	{
		warn $horizontal;
		warn @warnings;
		$result = 0;
        }

	return $result;
}
#####################################################################
sub lineCheck
#
#     &lineCheck( $verbose, @input ) 
#
#     Description: This subroutine  verifies that no line in the file
#     is longer than the maxLength. 
#
#        This subroutine works by comparing the length of each line
#        to the maxLength
#     
#     Arguments:
#     	Sverbose - Prints verbose messages when true.
#     	@input   - Lines of a Fortran source file.
#
#     Returns
#     	True if no line in the file is greater than the maxLength.
#####################################################################
{
	
	# Grab first argument
	my( $verbose ) = shift(@_);
	
	# Grab Second argument
	my( @input ) = @_;

	# Local Variables
	my( $line, $result, $lineNumber, $maxLength, @warnings);

	# The maximum character length a line is allowed to be.
	$maxLength = 100;
	
	# Default initialization
	$lineNumber = 0;

	$result = 1;

	# For every line of the file
	for $line (@input)
	{
		# Keep track of the line number
		$lineNumber++;

		# If the number of characters in the line is greater than the max length.
		if( length($line) - 1 > $maxLength )
		{
			# Start building the warning message
			push(@warnings,"$programName WARNING: Line has exceeded $maxLength characters.\n");
			push(@warnings, "$lineNumber : $line");

			$result = 0; # Check Failed!
		}
	}

	# Did the check fail?
	if( ! $result )
	{
		push(@warnings, "$programName WARNING: File has lines that exceed $maxLength characters.\n");
		
		# Show the warning messages	
		warn $horizontal;
		warn @warnings;
	}

	return $result;
}
#####################################################################
sub deprecateCheck
#
#     &deprecateCheck( $verbose, @input ) 
#
#     Description: This subroutine verifies that there are no deprecated
#       tokes used in the Fortran source file.
#
#     Arguments:
#     	Sverbose - Prints verbose messages when true.
#     	@input   - Lines of a Fortran source file.
#
#     Returns
#     	True if no line has deprecated tokens.
#####################################################################
{
	
	# Grab first argument
	my( $verbose ) = shift(@_);
	
	# Grab Second argument
	my( @input ) = @_;

	# Local Variables
	my( $line, $result, $lineNumber, $maxLength, @warnings,@test, $prev, %seen, @unique);

	# The maximum character length a line is allowed to be.
	$maxLength = 100;
	
	# Default initialization
	$lineNumber = 0;

	$result = 1;

	# For every line of the file
	for $line (@input)
	{
		# Keep track of the line number
		$lineNumber++;
                
		# If the number of characters in the line is greater than the max length.
		$line =~ m/^([^!]*)!?/;
		if( @test = $1 =~ m/$deprecatedRegEx/ig )
		{
			# Start building the warning message
			%seen = ();

                        @unique = grep { ! $seen{$_} ++ } @test;
			push(@warnings,"$programName WARNING: Line has deprecated elements:  @unique .\n");
			push(@warnings, "$lineNumber : $line");

			$result = 0; # Check Failed!
		}
	}

	# Did the check fail?
	if( ! $result )
	{
		push(@warnings, "$programName WARNING: File has deprecated elements.\n");
		
		# Show the warning messages	
		warn $horizontal;
		warn @warnings;
	}

	return $result;
}

