#!/usr/bin/perl -w
use strict;

my $spec_gsize = shift(@ARGV);
my $input = shift(@ARGV);
my $tag_size = shift(@ARGV);
my $center = shift(@ARGV) || 0;
my $output = shift(@ARGV) || $input;
$output =~ s/\.bed/.extended.bed/;

open(IN, "<$spec_gsize") || die "failed open\n";
my %size;
while (my $line = <IN>) {
    chomp $line;
    my @a = split(/\t/, $line);
    $a[0] =~ s/chr//;
    my $t = "chr$a[0]";
    if ($a[0] eq 'MT') {$t = "chrM";}
    $size{$t} = $a[1];
    $size{$a[0]} = $a[1];
}
close (IN);

my $nread=0;

open(IN, "<$input") || die "no $input\n";
open(OUT, ">$output");
while (my $line = <IN>) {
    chomp $line;
    my @a = split(/\s+/, $line);
    if ($a[0] =~ /NT/) {next;}
    if ($a[0] =~ /GL/) {next;}
    $a[0] =~ s/MT/M/;
    $a[3] = 'r_'.$nread;
#  if ($a[4] < 1) {next;}
    if ($a[5] eq '+') {$a[2] = $a[1]+$tag_size;}
    if ($a[5] eq '-') {$a[1] = $a[2] - $tag_size;}
    if($center){
        my $diff = int(($a[2] - $a[1] - $center)/2);
        $a[1] = $a[1] + $diff;
        $a[2] = $a[2] - $diff;
    }
    if ($a[1] < 0) {next;}
    if ($a[2] > $size{$a[0]}) {next;}
    if ($line =~ /^chr/) {
        print OUT join("\t", @a), "\n";
    }
    else {
        print OUT "chr", join("\t", @a), "\n";
    }
    $nread = $nread + 1;
}
close (IN);
close (OUT);
