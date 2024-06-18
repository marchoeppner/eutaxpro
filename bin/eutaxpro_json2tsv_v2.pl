#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Data::Dumper;
use JSON::XS;

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
my $min_perc    = 0.1;
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

my %levels = (
    "s" => "Species",
    "g" => "Genus",
    "f" => "Family",
    "o" => "Order"
);

my @taxlevels = ( "s", "g", "f");

# We read the JSON file
my $json_in;
{
    local $/; #Enable 'slurp' mode
    open my $fh, "<", $json;
    $json_in = <$fh>;
    close $fh;
}

my $data = decode_json($json_in);

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

printf "Sample\tReads\tPercentage\tTaxon\tTaxon_level\tBootstrap\n" ;

# We iterate over the list elements, which are individual samples
foreach my $entry ( @$data) {
    
    my $reads = $entry->{'reads_total'};
    my $sample = $entry->{'sample'};
    my @hits = @{$entry->{'hits'}};

    my @sorted_hits = reverse sort { $a->{'reads'} <=> $b->{'reads'} } @hits ;

    foreach my $hit (@sorted_hits) {
        
        my $taxdata = $hit->{'taxon'};
        my $species = undef;
        my $count = $hit->{'reads'};
        my $taxgroup = undef;
        my $bs = undef;
        my $perc = sprintf( "%.2f", ($count/$reads)*100);
        
        foreach my $taxlevel (@taxlevels) {
    
            if (!defined $species) {
                if (defined $taxdata->{$taxlevel}) {
                    $taxgroup = $levels{$taxlevel};
                    my $taxon = $taxdata->{$taxlevel};
                    $species = $taxon->{'taxon'};
                    $bs = $taxon->{'bootstrap'};
                }
            }
        }

        next if($perc < $min_perc);
        printf $sample . "\t" . $count . "\t" . $perc . "\t" . $species . "\t" . $taxgroup . "\t" . $bs . "\n";
        
    }

}