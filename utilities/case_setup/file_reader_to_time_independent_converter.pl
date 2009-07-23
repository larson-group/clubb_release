#!/usr/bin/perl

use strict;

use TextFormat qw/align_table/;
use List::Util qw[min max];

my(%variables, $file, $response );

# For every file

my($case_path) = shift(@ARGV);


foreach $file (@ARGV){
        
	my($dimension, @filedata, $line);

	$dimension = 0;

	# Open the file
	open FILE, $file or die "Bad Filename";

	 # Store file contents to memory
	@filedata = <FILE>;

	# Close it
	close FILE;
	
	print "Loaded $file\n";

        if( $filedata[0] =~ m/^\s*[^\s]+\s*$/ ){
		my($line,$var_name,$i);
		print "File is valid.\n";

		print "What variable does this file represent?\n";


		chomp( $var_name = <STDIN> );


		$variables{$var_name}{'dim'} = 1;

		$i = 0;
		foreach $line (@filedata){
			chomp($line);
			$variables{ $var_name }{'data'}[$i++] = $line;
		}
	}
	elsif($filedata[0] =~ m/^\s*([^\s]+\s+){0,4}[^\s]+\s*$/)
	{
		print "File is valid.\n";

		foreach $line (@filedata){
			if( $line !~ m/^\s*([^\s]+\s+){4}[^\s]+\s*$/ )
			{
				$dimension++;
				if($dimension > 1){
					last;
				}
			}
		}
	
		if($dimension == 2){
			print "Processing 2 dimensional file\n";
			
			my( $height_idx, $time_idx, $var_name, @values, $i );

			$height_idx = 0;

			$time_idx = 0;

			print "What variable does this file represent?\n";

			chomp( $var_name = <STDIN> );

			$variables{$var_name}{'dim'} = $dimension;

			foreach $line (@filedata){
				$line =~ s/^\s+//;           # Remove all leading spaces
				@values = split /\s+/,$line; # Split up line into entries, delimited by spaces
			
				for $i (0 .. $#values) {
					$variables{$var_name}{'data'}[$time_idx++][$height_idx] = $values[$i];
				}
				if( $#values < 4 ){
					$height_idx++;
					$time_idx = 0;
				}

			}

		} elsif($dimension = 1){
			print "Processing 1 dimensional file\n";

			my( $idx, @values, $var_name, $entry );

			print "What variable does this file represent?\n";
			
			chomp( $var_name = <STDIN> );

			$variables{ $var_name }{ 'dim' } = $dimension;

			$idx = 0;

			foreach $line (@filedata){
				$line =~ s/^\s+//;           # Remove all leading spaces
				@values = split /\s+/,$line; # Split up line into entries, delimited by spaces
				foreach $entry (@values){
					$variables{$var_name}{'data'}[$idx++] = $entry;
				}	
			}

		}
	}
}

print "What kind of file do you want to process?\n";
print "1 => surface\n";
print "2 => forcing\n";
print "3 => both\n";

$response = <STDIN>;

if($response == 1){
	&create_surface_file();
}
elsif($response == 2){
	&create_forcing_file();
}
elsif($response == 3){
	&create_surface_file();
	&create_forcing_file();
}

sub create_surface_file{

        my($i, $lh, $sh, $thlm, $rt, $press);

	print "Creating surface file";

	if( ! exists $variables{'time'} ){
		print "There is no time variable specified. Exiting.\n";
		exit
	}

        open FILE, ">", $case_path . "_surface.in";

        print FILE " Time[s]     LH[W\\m^2]     SH[W\\m^2]     thlm[K]     rt[kg\\kg]     Press[Pa]";

	for $i (0 .. $#{ $variables{'time'}{'data'} } ){
		$lh = &fetch_val( 'lh', $i );
		$sh = &fetch_val( 'sh', $i );
		$thlm = &fetch_val( 'thlm', $i );
		$rt = &fetch_val( 'rt', $i );
		$press = &fetch_val( 'press', $i );
	}
	
	close FILE;

	TextFormat->align_table( $case_path . "_surface.in");
}

sub create_forcing_file{

	my($key, $i, $j, $thlm_f, $rtm_f, $um_ref,$vm_ref,$um_f,$vm_f, $wm, $ug, $vg, $vthlm_f,$vrtm_f);

	print "Creating Forcing file\n";

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

	print FILE "z[m]    thlm_f[K\\s]   rtm_f[kg\\kg\\s]   um_ref[m\\s]   vm_ref[m\\s] " .
                        	"um_f[m\\s^2]   vm_f[m\\s^2]   w[m\\s]  ug[m\\s]  vg[m\\s]\n";

	for $i (0 .. $#{ $variables{'time'}{'data'} } ){

		print FILE $variables{'time'}{'data'}[$i] . " " . @{ $variables{'z'}{'data'} } . "\n";

		for $j (0 .. $#{ $variables{'z'}{'data'} }){
	 		
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

	my($var) = shift(@_);
	my($i) = shift(@_);
	my($j) = shift(@_);

        my($retVar);

	if( exists $variables{$var}){
		if($variables{$var}{'dim'} == 1) {
			$retVar = $variables{$var}{'data'}[$i];
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
