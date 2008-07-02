#$Id: CLUBBStandardsCheck.pl,v 1.1 2008-07-02 21:38:16 faschinj Exp $

#!/usr/bin/perl

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
#
#               This perl script assumes that 
#               the Fortran code compiles!!
#
##################################################################

# Includes
use Getopt::Long;

# Global Variables

# Print verbose messages?
$verbose = 0;

$programName = "CLUBBStandardsCheck.pl";

# Captures verbose command switch.
GetOptions ('v|verbose' => \$verbose);

# Test to see if at least one argument is present.
if(@ARGV < 1 )
{
	warn "Not enough arguments"	
}
#### BEGIN MAIN PROGRAM ####
else{
	# For Every File
	foreach $file (@ARGV)
	{
		# Open the file
		open FILE, $file or die "Bad Filename";
	
		# Store the lines to an array
		@input = <FILE>;

		# Check for missing implicit nones
		if( ! &implicitCheck( $verbose, @input ) )
		{
			warn "$file\n";	
		}

		# Check for use statements without only
		if( ! &useCheck( $verbose, @input ) )
		{
			warn "$file\n";
		}

		# Check for default private statements
		if( ! &privateCheck( $verbose, @input ) )
		{
			warn "$file\n";
		}
		
		# Close File
		close FILE;
	}
}
#### END MAIN PROGRAM ####

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
	        if($line =~ /^				# Bind to beginning of line
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
			/ix				# Whole expression is case insensitive
		)		
		{
			# Print if verbose
			if( $verbose )
			{	
				print "$programName comment: Found Function\n$line";
			}
			$statementCount++;
		}
		# Is it a subroutine declaration?	
		elsif($line =~/^			# Bind to beginning of line
	      		\s*?				# Zero or more spaces before statement
	           	(pure|elemental|recursive|\s*?)*?	# Zero or more Specifications
		        \s?				# Zero or One Space before statement
			\b				# Bind to front of statement
			subroutine			# Subroutine
			\b				# Bind to end of statement
			/ix				# Whole expression is case insensitive
		)
		{	
			# Print if verbose
			if( $verbose )
			{
				print "$programName comment: Found Subroutine\n$line";
			}
			$statementCount++;
			
		}
		# Is it a module declaration?
		elsif($line =~/^			# Bind to beginning of line
	      		\s*?				# Zero or more spaces before statement
			\b				# Bind to front of statement
			module				# Module
			\b				# Bind to end of statement
			\s+?				# One or more spaces
			\w+?				# One or more word characters 
			\s*?				# Zero or more spaces after statement
			(!.*?)*?			# After statement, accept comments and nothing else
			$				# Bind to end of line
			/ix				# Whole expression is case insensitive
		)           	
		{
			# Print if verbose	
			if( $verbose )
			{
				print "$programName comment: Found Module\n$line";
			}
			$statementCount++;
		}
		# Is it a program declaration?
		elsif($line =~/^			# Bind to beginning of line
	      		\s*?				# Zero or more spaces before statement
			\b				# Bind to beginning of statement
			program				# Program
			\b				# Bind to end of statement
			/ix				# Whole expression is case insensitve
		)
		{	
			if( $verbose )
			{
				print "$programName comment: Found Program\n$line";
			}
			$statementCount++;
		}
		# Is it an implicit none declaration?	
		elsif($line =~/^			# Bind to beginning of line
	      		\s*?				# Zero or more spaces before statement
			\b				# Bind to front of statement
			implicit\snone			# Implicit None
			\b				# Bind to end of statement
			\s*?				# Zero or more spaces after statement
			/ix				# Whole expression is case insensitive
		)
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
		warn "$programName warning: Missing 'implicit none' statements.\n";  
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
	my( $line, $result, $lineNumber );

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
		if($line =~ /^				# Bind to beginning of line 
			\s*?				# Zero or more spaces before statement
			use				# Use
			\s+?				# One or more spaces after statement
			\w+?				# One or more word characters
			\s*?				# Zero or more spaces
			(!.*?)*?			# After statement, accept comments and nothing else
			$				# Bind to end of line 
			/ix				# Whole expression is case insensitive
		)
		{       # If the above is true, execute the statements below (i.e. return an error):
			warn " $programName warning: 'use' statement w/o 'only' found in the following line:\n";
			warn "$lineNumber : $line";
			$result = 0; # Failed check
		}
		$lineNumber++;
	}
	if ( ! $result )
	{
		warn "$programName warning: CLUBBStandardsCheck.pl error: Use check failed!\n";
		warn "$programName warning: Check that comma is on same line as 'use', as CLUBB requires.\n";
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
		if($line =~/^				# Bind to beginning of line
	      		\s*?				# Zero or more spaces before statement
			\bmodule\b 			# Module
			\s 				# One space
			\w+? 				# One or more word characters
			\s*?				# Zero or more spaces after statement
			(!.*?)*?			# Accepts only comments after statement
			$				# Bind to end of line
			/ix				# Whole expression is case insensitive
		)                           
		{
			# Print if verbose
			if($verbose)
			{
				print "$programName comment: Module found\n $line";
			}
			$moduleCount++;	
		}
		# Does the line contain a private declaration?
		elsif($line =~/^			# Bind to beginning of line
	      		\s*?				# Zero or more spaces before statement
			\bprivate\b			# Private
			\s*? 				# Zero or more spaces after statement
			(!.*?)*?			# Accepts only comments after statement
			$				# Bind to end of line
			/ix				# Whole expression is case insensitive
		)
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
		warn "$programName warning: Number of \"private\" statements does not not match number of modules.\n";
		warn "$programName warning: Private Test failed!\n";
		$result = 0; # Check failed
	}	

	return $result;
}
