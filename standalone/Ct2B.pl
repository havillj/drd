#!/usr/bin/perl
# Convert from ct format to Vienna format
# 2 record output - sequence, structure + energy if '=' found in
# header and something follows.
# M. Zuker, February 10, 2009.
#
# INPUT: Use standard input or give file name on the command line.
# OUTPUT: standard output

use strict;
use warnings;

my @Rec ;
my ($energy, $Bp, $Seq);
my $n = 1; my $Last_Seq = '';

while (<>) {
    chomp;
    @Rec = split(' ');
    if ($. % ($n + 1) == 1) {
	$energy = $Seq = $Bp = '';
	$n = $Rec[0];
	my @a = split('=');
	defined($a[1]) && (my @e = split(' ', $a[1]));
	defined($e[0]) && ($energy = "\t(" . $e[0] . ')' );
    } else {
	$Seq .= $Rec[1];
	if ($Rec[4] == 0) {
	    $Bp .= '.';
	} elsif ($Rec[0] < $Rec[4]) {
	    $Bp .= '(';
	} else {
	    $Bp .= ')';
	}
	if ($. % ($n + 1) == 0) {
	    if ($Last_Seq ne $Seq ) {
		printf "%s\n", $Seq;
		$Last_Seq = $Seq;
	    }
	    printf "%s%s\n", $Bp,$energy;
	}
    }
}
