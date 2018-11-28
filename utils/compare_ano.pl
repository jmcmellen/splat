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

    $ano1{$lat}{$long} = $loss;

    if ($lineno % 20000 == 0) {
        printf(".", $lineno);
    }
    $lineno++;
}
close($fh1);


printf("\nComparing to second file...\n");

$lineno = 0;
$unfound = 0;
my $errcount = 0;
while (my $line = <$fh2>) {
    chomp $line;
    my ($lat, $long, $azimuth, $elevation, $loss2, $blocked) = split(/,?\s+/, $line);
    chop $lat;
    chop $long;

    # we can be off by 0.000001 in either latitude or longitude and it's ok.
    # So we need to look at the available values and see if there's anything near.
    # first, find exact matches, cuz that's fastest
    if (!exists $ano1{$lat}{$long}) {
        my $found = false;

        #printf("%s %s not found\n", $lat, $long);

        # search at this exact latitude
        foreach my $testlong (keys %{ $ano1{$lat} } ) {
            #printf("trying %s %s: ", $lat, $testlong);
            if (abs($testlong - $long) < 0.000002) {
                #printf("yup!\n");
                $long = $testlong;
                $found = true;
                break;
            }
            #printf("nope\n");
        }

        # search slightly above
        if ($found == false) {
            my $testlat = $lat + 0.000001;
            foreach my $testlong (keys %{ $ano1{$lat} } ) {
                #printf("trying %s %s: ", $lat, $testlong);
                if (abs($testlong - $long) < 0.000002) {
                    #printf("yup!\n");
                    $lat = $testlat;
                    $long = $testlong;
                    $found = true;
                    break;
                }
                #printf("nope\n");
            }
        }

        # search slightly below
        if ($found == false) {
            my $testlat = $lat - 0.000001;
            foreach my $testlong (keys %{ $ano1{$lat} } ) {
                #printf("trying %s %s: ", $lat, $testlong);
                if (abs($testlong - $long) < 0.000002) {
                    #printf("yup!\n");
                    $lat = $testlat;
                    $long = $testlong;
                    $found = true;
                    break;
                }
                #printf("nope\n");
            }
        }

        if ($found == false) {
            printf("%s %s not found\n", $lat, $long);
            $unfound++;
            next;
        }

        #printf("using %s %s\n", $lat, $long);
    }

    my $loss1 = $ano1{$lat}{$long};

    my $delta = abs($loss1 - $loss2);

    if ($delta > 0.5) {
        printf("%s %s: %s <-> %s    %0.2f\n", $lat, $long, $loss1, $loss2, $delta);
        $errcount++;
    }

    if ($lineno % 50000 == 0) {
        printf("%d   %d\n", $lineno, $errcount);
    }
    $lineno++;

}
printf("\n");

printf("%s         %s          %s\n", "total", ">0.5 off", "missing");
printf("%d         %d          %d\n", $lineno, $errcount, $unfound);


close($fh2);

