#!/usr/bin/perl -w

open (OUT,">mpace_a_times.dat");

$seconds = 7200; # First time is 2:00 AM on 5 Oct 2004
for ($i=1;$i<140;$i++) {
  print OUT $seconds . " ";
  if (($i % 5) == 0) { print OUT "\n"; } # Print five values then do a newline
  $seconds += 10800; # Three hour time steps
}

close OUT;
