#!/usr/bin/env perl

use strict;
use Getopt::Long;

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
            $otu_translations{$otu} = "Keine Auflösung zur Gattung möglich.";
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

foreach my $line (@lines) {

    my @elements = (split "\t",$line);
    # Get OTU label and shift the remaining columns to match header
    my $otu = shift @elements;

    my $taxon = $otu_translations{$otu};

    for (0..$#elements) {
        my $count = @elements[$_];
        my $sample = @header[$_];
        my %payload = ( "taxon" => $taxon, "count" => $count);
        
        if (defined $matrix{$sample}) {
            push @{ $matrix{$sample} },\%payload ;
        } else {
            $matrix{$sample} = [ \%payload ] ;
        }
        
    }

}


if ($outfile) {
    open(STDOUT, ">$outfile") or die("Cannot open $outfile");
}

printf "sample\treads\thits\n";

# One sample, one row of results
foreach my $sample (sort keys %matrix) {

    my $hits = $matrix{$sample};
    my $sum = 0.0 ;
    my @row;

    # We iterate once to get the total sum of reads
    foreach my $hit (@{ $hits} ) {
        my $count = $hit->{"count"};
        my $taxon = $hit->{"taxon"};
        
        $sum += $count;

        if ($count >= $min_cov) {
            push @row, $hit;
        }
    
    }  

    # Sometimes we get slightly divergent OTUs for the same species, 
    # lets add it up for a combined result
    my %results ;

    foreach my $hit(@row) {

        my $taxon = $hit->{"taxon"};
        my $count = $hit->{"count"};
        if (defined $results{$taxon}) {
            $results{$taxon} += $count;
        } else {
            $results{$taxon} = $count;
        }

    }

    my @all;
    printf $sample . "\t" . $sum . "\t";

    foreach my $key (sort { $results{$b} <=> $results{$a} } keys %results ){
        my $fcount = $results{$key};
        my $perc = sprintf( "%.2f", ($fcount/$sum)*100);
        next if ($perc < 1.0);
        push(@all,"${key}:${perc}");
    }

    printf join(", ",@all);
    printf "\n";
    
}

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
