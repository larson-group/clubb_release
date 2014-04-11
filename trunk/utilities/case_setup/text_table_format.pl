#!/usr/bin/perl

use strict;

use TextFormat qw/align_table/;

my( $file );

# For every file
foreach $file (@ARGV){
        
	TextFormat->align_table( $file );
}
