#!/usr/bin/perl

package TextFormat; 

use strict;

sub align_table{


	my($file) = shift(@_);

	my(@filedata, $itr, $tmp, $line, @longest);

	my($col, $rol, @comments, $entry);

	# Open the file
	open FILE, $file or die "Bad Filename";

	 # Store file contents to memory
	@filedata = <FILE>;

	# Close it
	close FILE;


	# Strip out comment lines
	for $itr ( 0 .. $#filedata ){
		$tmp = shift(@filedata);
		if($tmp =~ /^\s*!.*/i){
			push(@comments,$tmp);
		}else{
			push(@filedata,$tmp);
		}
			
	}


	# Delimit table by single spaces (Essentially setting up the 2-D array.
	foreach $line (@filedata){
		$line =~ s/^\s+//;      # Remove all leading spaces
        	@{$line} = split /\s+/,$line; # Split up line into entries, delimited by spaces
	}


	# Array of the longest element in each column
	@longest = 0;
        
	# Iterate through each column determining which element in the column is the longest.
	for $col (0 .. $#{$filedata[0]}){	# For every column
	  	for $row ( 0 .. $#filedata ) {  # For every row
			if( length($filedata[$row][$col]) > $longest[$col] ){ 
				$longest[$col] = length($filedata[$row][$col]);
			}
		}

		# Pad out for approriate spacing between columns
  		$longest[$col]++;
  		$longest[$col]++;
  		$longest[$col]++;
	}	

	# Open the same file for writing
	open FILE, ">", $file;

	# Construct the corrected table in memory
	for $col (0 .. $#{$filedata[0]}){	# For every column
  		for $row ( 0 .. $#filedata ) {  # For every row
        		$filedata[$row][$col] = $filedata[$row][$col] . ' ' x ($longest[$col] - length($filedata[$row][$col]))
  		}
	}	
	
        # Write out the comments to the file
	foreach $line (@comments){
		print FILE $line;
	}

	# Write out the table to the file
	foreach $line (@filedata){
		foreach $entry (@{$line}){
			print FILE $entry;
		}
		print FILE "\n";
	}

	# Close the file
	close FILE;
}

1;
