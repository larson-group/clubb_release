#!/usr/bin/perl


# This script formats a table of text such that columns 
# are properly spaced with all entities aligned in their column.
#
# Originally written to help ease the transition to 'SAM-like' input files.
#
# Written By: Joshua Fasching
# April 2009


# For every file
foreach $file (@ARGV){
        
	# Open the file
	open FILE, $file or die "Bad Filename";

	 # Store file contents to memory
	@filedata = <FILE>;

	# Close it
	close FILE;

	# Delimit table by single spaces (Essentially setting up the 2-D array.
	foreach $s (@filedata){
		$s =~ s/^\s+//;	
        	@{$s} = split /\s+/,$s;
	}
        
	# Array of the longest element in each column
	@longest = 0;
        
	# Iterate through each column determining which element in the column is the longest.
	for $j (0 .. $#{$filedata[0]}){	
	  	for $i ( 0 .. $#filedata ) {
			if(length($filedata[$i][$j]) > $longest[$j]){
				$longest[$j] = length($filedata[$i][$j]);
		}
	}

	# Pad out for approriate spacing between columns
  	$longest[$j]++;
  	$longest[$j]++;
  	$longest[$j]++;
	}	

	# Open the same file for writing
	open FILE, ">", $file;

	# Construct the corrected table in memory
	for $j (0 .. $#{$filedata[0]}){	
  		for $i ( 0 .. $#filedata ) {
        		$filedata[$i][$j] = $filedata[$i][$j] . ' ' x ($longest[$j] - length($filedata[$i][$j]))
  		}
	}	
	
	# Write out the table to the file
	foreach $a (@filedata){
		foreach $b (@{$a}){
			print FILE "$b";
		}
		print FILE "\n";
	}

	# Close the file
	close FILE;
}
