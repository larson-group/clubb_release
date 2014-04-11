#!/usr/bin/perl
###########################################################################
# Plotgen Console Output Writer v1.0
#
# Sample MATLAB calls:
#   unix(['./console_output.pl' ' "Case Name: ' caseName '"']);
#   unix(['./console_output.pl ' '-s ' '"This is an error."']);
###########################################################################

package Plotgen;

use strict;
use Term::ANSIColor;
use Getopt::Std;

my $VERSION = 1.0;

main();

sub main()
{
    readArgs();

    exit(0);
}

sub readArgs()
{
    my $numArgs = $#ARGV + 1;

    if($numArgs == 0)
    {
        main::HELP_MESSAGE();
    }

    if($numArgs == 1)
    {
        my %option = ();
        my $result = getopts("h", \%option);
        my $message = $ARGV[0];

        if($option{h})
        {
            main::HELP_MESSAGE();
        }

        # Assume that this is a message
        message($message);
    }
    elsif($numArgs == 2)
    {
        my %option = ();
        my $result = getopts("swmh?", \%option);
        my $message = $ARGV[0];

        if($option{s})
        {
            severe($message);
        }
    
        if($option{w})
        {
            warning($message);
        }

        if($option{m})
        {
            message($message);
        }
    }
    else # Too many arguments
    {
        main::HELP_MESSAGE();
    }
}

sub message()
{
    my $message = shift(@_);
    
    print("$message\n");
}

sub warning()
{
    my $message = shift(@_);

    print color 'yellow';
    print("WARN: $message\n");

    resetColors();
}

sub severe()
{
    my $message = shift(@_);
    
    print color 'bold red';
    print("SEVERE: $message\n");

    resetColors();
}

sub resetColors()
{
    print color 'reset';
}

sub main::HELP_MESSAGE()
{
    print("Usage: console_output.pl [OPTION] MESSAGE\n");
    print("  -s\tPrint a severe message\n");
    print("  -w\tPrint a warning message\n");
    print("  -m\tPrint a standard message\n");
    print("In no option is specified, a standard message will be printed\n");
    exit(0);
}

sub main::VERSION_MESSAGE()
{
    print("Plotgen Console Output Writer version $VERSION, Copyright (c) 2009 Larson Group.\n");
}
