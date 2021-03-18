#!/usr/bin/perl
#
#	This program is used for converting sets of files that used the 
#       file format read by file_functions.F90 into files that can be read in
#       by time_dependent_input.F90
#
#       The old format looked something like this:
#       ------ Varies in time down list ---> 
#     | DATAPOINT DATAPOINT DATAPOINT DATAPOINT DATAPOINT
#  V  | DATAPOINT DATAPOINT DATAPOINT DATAPOINT DATAPOINT
#  a  | DATAPOINT DATAPOINT DATAPOINT DATAPOINT DATAPOINT
#  r  | DATAPOINT DATAPOINT DATAPOINT DATAPOINT DATAPOINT
#  i  | DATAPOINT DATAPOINT DATAPOINT 
#  e  | DATAPOINT DATAPOINT DATAPOINT DATAPOINT DATAPOINT
#  s  | DATAPOINT DATAPOINT DATAPOINT DATAPOINT DATAPOINT
#     | DATAPOINT DATAPOINT DATAPOINT DATAPOINT DATAPOINT
#  z  | DATAPOINT DATAPOINT DATAPOINT DATAPOINT DATAPOINT
#     | DATAPOINT DATAPOINT DATAPOINT 
#     V
#
#       The format being converted to looks like this
#	VARTITLE[UNITS] VARTITLE[UNITS] VARTITLE[UNITS] VARTITLE[UNITS] VARTITLE[UNITS]
#       TIME NZ
#       ZDATAPOINT      DATAPOINT       DATAPOINT       DATAPOINT       DATAPOINT
#       ZDATAPOINT      DATAPOINT       DATAPOINT       DATAPOINT       DATAPOINT
#       ZDATAPOINT      DATAPOINT       DATAPOINT       DATAPOINT       DATAPOINT
#       ZDATAPOINT      DATAPOINT       DATAPOINT       DATAPOINT       DATAPOINT
#       TIME NZ
#       ZDATAPOINT      DATAPOINT       DATAPOINT       DATAPOINT       DATAPOINT
#       ZDATAPOINT      DATAPOINT       DATAPOINT       DATAPOINT       DATAPOINT
#       ZDATAPOINT      DATAPOINT       DATAPOINT       DATAPOINT       DATAPOINT
#       ZDATAPOINT      DATAPOINT       DATAPOINT       DATAPOINT       DATAPOINT
#
#
#       Example usage:
#       ./file_reader_to_time_dependent_converter.pl <PATH TO OUTPUT><PREFIX OF OUTPUT> <INPUT FILES>
#
#------------------------------------------------------------------------------

use strict;

use TextFormat qw/align_table/;
use List::Util qw[min max];

# Necessary for strict
my( %variables,  # Data structure that holds data parsed from input.
    $file,       # Current file being read in.
    $response ); # Response by user to questions.

# Store the first argument to be used as the prefix for output files.
my($case_path) = shift(@ARGV);

# For every input file file
foreach $file (@ARGV){
        
	# Necessary for strict
	my( $dimension, # Dimensionality of file
            @filedata,  # Array of lines from the file
            $line );    # Individual line of the file

	# Assume the dimensionality is zero until found otherwise
	$dimension = 0;

	# Open the file
	open FILE, $file or die "Bad Filename: $file";

	 # Store file contents to memory
	@filedata = <FILE>;

	# Close it
	close FILE;
	
	print "Loaded $file\n";

        # If the file only has one datapoint per line
        if( $filedata[0] =~ m/^\s*[^\s]+\s*$/ ){
		my($line,$var_name,$i);
		print "File is valid.\n";


                print "Processing 1 dimensional file as a column";

		print "What variable does this file represent?\n";

		# Wait for user feedback
		chomp( $var_name = <STDIN> );

		$variables{$var_name}{'dim'} = 1;


                # Store each datapoint in an array associated with a hash which is
		# identified with the name provided by the user
		$i = 0;

		foreach $line (@filedata){

			chomp($line);
			$variables{ $var_name }{'data'}[$i++] = $line;

		}
	}
	# If the file has no more than 5 and no fewer than 1 datapoint on a line
	elsif($filedata[0] =~ m/^\s*([^\s]+\s+){0,4}[^\s]+\s*$/)
	{
		print "File is valid.\n";
                
		# Determine if there are multiple paragraphs of datapoints
		foreach $line (@filedata){
			if( $line !~ m/^\s*([^\s]+\s+){4}[^\s]+\s*$/ )
			{
				$dimension++;
				if($dimension > 1){
					last;
				}
			}
		}
	
		# If the file represents 2 dimensional data by having multiple paragraphs
		if( $dimension == 2 ){
			print "Processing 2 dimensional file\n";
			
			# Required by strict
			my( $height_idx, # Vertical index
                            $time_idx,   # Temporal index
                            $var_name,   # Name of the variable
                            @values,     # Array of values represented by a particular line
                            $i );        # Index

			$height_idx = 0;

			$time_idx = 0;

			print "What variable does this file represent?\n";

			# Wait for user input
			chomp( $var_name = <STDIN> );

			# Assign dimensionality
			$variables{$var_name}{'dim'} = $dimension;
			
			# For every line
			foreach $line (@filedata){
				# Split that line into an array

				$line =~ s/^\s+//;           # Remove all leading spaces
				@values = split /\s+/,$line; # Split up line into entries, delimited by spaces


				# For every value in that array			
				for $i (0 .. $#values) {
					# Assign it to the array in the hash mapped with the user supplied var_name
					$variables{$var_name}{'data'}[$time_idx++][$height_idx] = $values[$i];
				}
				# if there were less than 5 entries on that line
				if( $#values < 4 ){
					# It means that the next line is at a higher altitude
					$height_idx++;
					# And starting back at the first time point
					$time_idx = 0;
				}

			}

		} elsif($dimension = 1){
			print "Processing 1 dimensional file\n";

			my( $idx, @values, $var_name, $entry );

			print "What variable does this file represent?\n";
			
			# Wait for user input
			chomp( $var_name = <STDIN> );

			# Assign dimensionality
			$variables{ $var_name }{ 'dim' } = $dimension;

			$idx = 0;
			
			# For every line
			foreach $line (@filedata){

				# Split the line into an array
				$line =~ s/^\s+//;           # Remove all leading spaces
				@values = split /\s+/,$line; # Split up line into entries, delimited by spaces

				# For every entry in that array
				foreach $entry (@values){
					# Assign it to the array in the hash mapped to the user supplied var_name
					$variables{$var_name}{'data'}[$idx++] = $entry;
				}	
			}

		}
	}
}


# Ask the user how they want to process the data
print "What kind of file do you want to process?\n";
print "1 => sfc\n";
print "2 => forcing\n";
print "3 => both\n";

# Wait for their response
chomp( $response = <STDIN> );

if($response == 1){
	&create_sfc_file();
}
elsif($response == 2){
	&create_forcing_file();
}
elsif($response == 3){
	&create_sfc_file();
	&create_forcing_file();
}

sub create_sfc_file{
#
#  This subroutine attempts to create an <case_path>_sounding.in file using 
#  data provided by the input files.
#
#  Arguments:
#    none
#
#  Returns:
#    none
#
#------------------------------------------------------------------------------

        my($i, $lh, $sh, $thlm, $rt, $press);

	print "Creating sfc file";
	
	# Sanity checking
	if( ! exists $variables{'time'} ){
		print "There is no time variable specified. Exiting.\n";
		exit
	}

        open FILE, ">", $case_path . "_sfc.in";

	# Write out header
        print FILE " Time[s]     LH[W\\m^2]     SH[W\\m^2]     thlm[K]     rt[kg\\kg]     Press[Pa]\n";

	# For every time write out the values for each field at that time
	for $i (0 .. $#{ $variables{ 'time' }{ 'data' } } ){
		$lh = &fetch_val( 'lh', 0, $i );
		$sh = &fetch_val( 'sh', 0, $i );
		$thlm = &fetch_val( 'thlm', 0, $i );
		$rt = &fetch_val( 'rt', 0, $i );
		$press = &fetch_val( 'press', 0, $i );

		print FILE $variables{ 'time' }{ 'data' }[$i] ." ". $lh . " " . $sh . " " . $thlm . " " . $rt . " ". $press . "\n";
	}
	
	close FILE;

	TextFormat->align_table( $case_path . "_sfc.in");
}

sub create_forcing_file{
#
#  This subroutine attempts to create an <case_path>_forcings.in file using 
#  data provided by the input files.
#
#  Arguments:
#    none
#
#  Returns:
#    none
#
#------------------------------------------------------------------------------
	my($key, $i, $j, $thlm_f, $rtm_f, $um_ref,$vm_ref,$um_f,$vm_f, $wm, $ug, $vg, $vthlm_f,$vrtm_f);

	print "Creating Forcing file\n";

	# Sanity Checking
        foreach $key (keys %variables){
		print "Keys: $key\n";
	}

	if( ! exists $variables{'time'} ){
		print "There is no time variable specified. Exiting.\n";
		exit
	}

	if( ! exists $variables{'z'} ){
		print "There is no z variable specified. Exiting.\n";
		exit
	}
	        
        open FILE, ">", $case_path . "_forcings.in";

        # Write out the header
	print FILE "z[m]    thlm_f[K\\s]   rtm_f[kg\\kg\\s]   um_ref[m\\s]   vm_ref[m\\s] " .
                        	"um_f[m\\s^2]   vm_f[m\\s^2]   w[m\\s]  ug[m\\s]  vg[m\\s]\n";

	# For every time
	for $i (0 .. $#{ $variables{'time'}{'data'} } ){

		# Write out the time header
		print FILE $variables{'time'}{'data'}[$i] . " " . @{ $variables{'z'}{'data'} } . "\n";

		# For every altitude
		for $j (0 .. $#{ $variables{'z'}{'data'} }){
	 		# Write out the value of each variable at that time and altitude
			$thlm_f = &fetch_val( 'thlm_f', $i, $j );
			$rtm_f = &fetch_val( 'rtm_f', $i, $j );
			$vthlm_f = &fetch_val( 'vthlm_f', $i, $j );
			$vrtm_f = &fetch_val( 'vrtm_f', $i, $j );
			$um_ref = &fetch_val( 'um_ref', $i, $j );
			$vm_ref = &fetch_val( 'vm_ref', $i, $j );
			$um_f = &fetch_val( 'um_f', $i, $j );
			$vm_f = &fetch_val( 'vm_f', $i, $j );
			$wm = &fetch_val( 'wm', $i, $j );
			$ug = &fetch_val( 'ug', $i, $j );
			$vg = &fetch_val( 'vg', $i, $j );
		
			$thlm_f = $thlm_f + max(0.0, $vthlm_f);
			$rtm_f = $rtm_f + max(0.0, $vrtm_f);

			print FILE $variables{'z'}{'data'}[$j]. " " . $thlm_f . " ". $rtm_f . " " . $um_ref . " " . $vm_ref . " " . $um_f .
				" ". $vm_f . " " . $wm ." ". $ug . " " . $vg . "\n";
		}
	}

	close FILE;

	TextFormat->align_table( $case_path . "_forcings.in" );

}
sub fetch_val{
#
#	This subroutine retrieves a value from the variables has if it exists
#       otherwise it returns a missing value (e.g. - 999.9)
#
#	Arguments:
#		var - Key used to identify variable in the hash
#		i   - Time index of variable
#		j   - Vertical index of variable
#
#	Returns: Value of variable at i,j
#
#------------------------------------------------------------------------------
	my($var) = shift(@_);
	my($i) = shift(@_);
	my($j) = shift(@_);

        my($retVar);

	if( exists $variables{$var}){
		if($variables{$var}{'dim'} == 1) {
			
			$retVar = $variables{$var}{'data'}[$j];
		} elsif($variables{$var}{'dim'} == 2) {
			$retVar = $variables{$var}{'data'}[$i][$j];
		} else {
			warn "$var has invalid dimensions\n";
			$retVar = -999.9;
		}
	}else{
		warn "$var does not exist\n";
		$retVar = -999.9;
	}
	return $retVar;
}
