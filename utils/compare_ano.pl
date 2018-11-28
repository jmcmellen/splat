#!/usr/bin/perl

# compare_ano.pl <file1> <file2>
# 
# compares two .ano files emitted from splat (via the -ano flag). The
# contents don't have to be in the same order, but they should be from
# runs with:
#  a) the same transmitter site
#  b) the same receiver site
#  c) the same settings for the splat compilation (8x8, 2x2, etc)
#  d) the same terrain elevation files
#
# This is primarily to test for differences between linear processing,
# multithreaded-CPU processing, and GPU processing.
#

$|++;

my ($file1, $file2) = @ARGV;
if (not defined $file1 || not defined $file2) {
    printf("compare_ano.pl <file1> <file2>\n");
    exit;
}

open(my $fh1, $file1) or die "Cannot open file: $file1";
open(my $fh2, $file2) or die "Cannot open file: $file2";

printf("Loading first file into memory...\n");

my %ano1;

my $lineno = 0;
while (my $line = <$fh1>) {
    chomp $line;
    my ($lat, $long, $azimuth, $elevation, $loss, $blocked) = split(/,?\s+/, $line);
    chop $lat;
    chop $long;

    my $key = $lat."-".$long;

    $ano1{$key} = $loss;

    if ($lineno % 20000 == 0) {
        printf(".", $lineno);
    }
    $lineno++;
}
close($fh1);


printf("\nComparing to second file...\n");
printf("\nLines Errors (>0.02%%)\n");

$lineno = 0;
my $errcount = 0;
while (my $line = <$fh2>) {
    chomp $line;
    my ($lat, $long, $azimuth, $elevation, $loss2, $blocked) = split(/,?\s+/, $line);
    chop $lat;
    chop $long;


    my $key = $lat."-".$long;
    if (!exists $ano1{$key}) {
        printf("%s %s not found\n", $lat, $long);
        next;
    }

    my $loss1 = $ano1{$key};

    my $pct = 0;
    # prevent divide by zero errors
    if ($loss1 == 0) {
        $loss1 = 0.0001;
    }
    if ($loss2 == 0) {
        $loss2 = 0.0001;
    }

    if ($loss1 >= $loss2) {
        $pct = 100 - (($loss1/$loss2) * 100);
    } else {
        $pct = 100 - (($loss2/$loss1) * 100);
    }

    if (abs($loss1 - $loss2) > 0.02) {
        printf("%s %s: %s <-> %s    %0.2f%%\n", $lat, $long, $loss1, $loss2, abs($pct));
        $errcount++;
    }

    if ($lineno % 50000 == 0) {
        printf("%d   %d\n", $lineno, $errcount);
    }
    $lineno++;

}
printf("\n");

close($fh2);

