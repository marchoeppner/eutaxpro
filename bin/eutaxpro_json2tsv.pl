#!/usr/bin/env perl

use strict;
use Getopt::Long;
use JSON::Parse 'read_json';

my $usage = qq{
perl eutaxpro_json2tsv.pl
    Getting help:
    [--help]

    Input:
    [--json filename]
        The name of the JSON file to rea
    Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile     = undef;
my $json        = undef;
my $min_cov     = 10;

my $help;

GetOptions(
    "help" => \$help,
    "json=s" => \$json,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

# We read the JSON file
my $data = read_json($json);

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

# We iterate over the list elements, which are individual samples
foreach my $sample ( @{$data}) {
    
    my %matrix;

    # Keep track of how many reads total were in this sample
    my $sum = 0;

    my $s = $sample->{'sample'} ;
    my $hits = $sample->{'hits'};

    # Iterate over all taxonomic hits
    foreach my $hit (@{$hits}) {
        my $taxon = $hit->{'taxon'};
        my $counts = $hit->{'reads'};
        
        $sum += $counts;
        
        # Add redundant entries, if any
        if (defined $matrix{$taxon}) {
            $matrix{$taxon} += $counts;
        } else {
            $matrix{$taxon} = $counts;
        }
    }
    
    printf $s . "\t" . $sum . "\t" ;

    my @all;

    foreach my $key (sort { $matrix{$b} <=> $matrix{$a} } keys %matrix ){
         my $fcount = $matrix{$key};
        my $perc = sprintf( "%.2f", ($fcount/$sum)*100);
        next if ($perc < 1.0);
        push(@all,"${key}:${perc}");
    }
    printf join(", ",@all);
    printf "\n";

}