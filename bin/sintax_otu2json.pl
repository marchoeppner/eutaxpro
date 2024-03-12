#!/usr/bin/env perl

use strict;
use Getopt::Long;
use JSON;

my $usage = qq{
perl sintax_otu2tab.pl
    Getting help:
    [--help]

    Input:
    [--otu filename]
        The name of the OTU table to read
    [--sintax filename]
		The name of the Sintax result file to read. 
    Ouput:    
    [--outfile filename]
        The name of the output file. By default the output is the
        standard output
};

my $outfile     = undef;
my $sintax      = undef;
my $otu         = undef;
my $min_cov     = 10;

my %otu_translations;
my %matrix;
my $help;

GetOptions(
    "help" => \$help,
    "otu=s" => \$otu,
    "sintax=s" => \$sintax,
    "outfile=s" => \$outfile);

# Print Help and exit
if ($help) {
    print $usage;
    exit(0);
}

# We read the sintax table
open (my $SINTAX, '<', $sintax) or die "FATAL: Can't open file: $sintax for reading.\n$!\n";

chomp(my @lines = <$SINTAX>);

foreach my $line (@lines) {
    my ($otu,$tax_string,$strand,$final_call) = split("\t", $line);
    $otu =~ s/;size=.*$// ;

    if (defined $final_call) {
        my $taxon = decode_taxstring($final_call);
        # Annotation has species level
        if (defined $taxon->{"s"}) {
            $otu_translations{$otu} = $taxon->{"s"};
        } elsif (defined $taxon->{"g"}) {
            $otu_translations{$otu} = $taxon->{"g"};
        } else {
            $otu_translations{$otu} = "No identification at genus level possible";
        }
    } else {
        $otu_translations{$otu} = "OTU unbekannt.";
    }
}

close($SINTAX);

open (my $OTU, '<', $otu) or die "FATAL: Can't open file: $otu for reading.\n$!\n";

chomp(my @lines = <$OTU>);

# First colum is the OTU name label, don't need that
my $hline = shift @lines;
my @header = (split "\t", $hline);
my $oc = shift @header;
my @data = ();
foreach my $h (@header) {
    push(@data, { "sample" => $h, "hits" => [], "reads_total" => 0 } )
}

foreach my $line (@lines) {

    my @elements = (split "\t",$line);
    # Get OTU label and shift the remaining columns to match header
    my $otu = shift @elements;

    my $taxon = $otu_translations{$otu};

    for (0..$#elements) {
        my $count = @elements[$_];
        my $sample = @header[$_];

        # Skip any assignments that are zero. 
        next if ($count == 0);

        my %payload = ( "taxon" => $taxon, "reads" => $count);
        
        push @{ $data[$_]{'hits'} }, \%payload ;

        $data[$_]{'reads_total'} += $count ;
        
    }

}

my $json = encode_json(\@data);

if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

print $json . "\n";

# Turns the tax string into a hash
sub decode_taxstring {

    my $taxstring = shift;
    my %data;

    foreach my $element (split ",", $taxstring) {
        my ($key,$value) = (split ":", $element);
        $value =~ s/_[0-9]*// ; 
        $data{$key} = $value ;
    }
    return \%data;
}
