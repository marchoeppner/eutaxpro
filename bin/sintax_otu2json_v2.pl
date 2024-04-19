#!/usr/bin/env perl
use strict;
use Getopt::Long;
use JSON::XS;
use Data::Dumper;

my $usage = qq{
perl sintax_otu2tab_v2.pl
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

my %otu_translations;
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

    my $taxdata;
    
    if (defined $tax_string) {
        $taxdata = decode_taxstring($tax_string);
    } 
    
    $otu_translations{$otu} = $taxdata ;
}

close($SINTAX);

open (my $OTU, '<', $otu) or die "FATAL: Can't open file: $otu for reading.\n$!\n";

chomp(my @lines = <$OTU>);

# First colum is the OTU name label, don't need that
my $hline = shift @lines;
my @header = (split "\t", $hline);
my $oc = shift @header;
my @data = ();

# This assumes that the order of samples is identical throughout but
# saves time since we  not need to initialize anything moving forward. 
foreach my $h (@header) {
    push(@data, { "sample" => $h, "hits" => [], "reads_total" => 0 } )
}

foreach my $line (@lines) {

    my @elements = (split "\t",$line);
    # Get OTU label and shift the remaining columns to match header
    my $otu = shift @elements;

    my $taxdata = $otu_translations{$otu};

    for (0..$#elements) {
        my $count = @elements[$_];
        my $sample = @header[$_];

        # Skip any assignments that are zero. 
        next if ($count == 0);

        my %payload = ( "taxon" => $taxdata, "reads" => $count);
        
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
        my $bootstrap = (split /[),(]/, $taxstring)[-1];
        $value =~ s/_[0-9]*.*// ; 
        $data{$key} = { "taxon" => $value, "bootstrap" => $bootstrap } ;
    }
    return \%data;
}
