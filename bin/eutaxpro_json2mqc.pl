#!/usr/bin/env perl

use strict;
use Getopt::Long;
use JSON::XS;

my $usage = qq{
perl eutaxpro_json2mqc.pl
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

my %mqc;

my %pconfig = ( 'id' => 'custom_data_barplot', 
    'tile' => "Vsearch Sintax",
    'ylab' => "Percentage",
    'xDecimals' => "false"
);
$mqc{'id'} = 'taxonomy';
$mqc{'section_name'} = "Vsearch Sintax";
$mqc{'description'} = "Taxonomic assignments with Vsearch Sintax";
$mqc{'plot_type'} = 'bargraph';
$mqc{'pconfig'} = \%pconfig;
$mqc{'data'} = {};

my @taxa;

# Get a list of all taxa
foreach my $sample (@{$data}) {
    my $hits = $sample->{'hits'};
    foreach my $hit (@{$hits}) {
        my $taxon = $hit->{'taxon'};
        if ( !grep(/^$taxon$/, @taxa)) {
            push(@taxa,$taxon);
        }
    }
}
my @sorted_taxa = sort @taxa;

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
    
    my %bucket;

    foreach my $taxon (@sorted_taxa) {
        if (defined $matrix{$taxon}) {
            my $count = $matrix{$taxon};
            my $perc = sprintf( "%.2f", ($count/$sum)*100);
            $bucket{$taxon} = $perc;
        } else {
            $bucket{$taxon} = 0 ;
        }
    }
    $mqc{'data'}{$s} = \%bucket;

}

my $json_out = encode_json(\%mqc);

printf $json_out ; 