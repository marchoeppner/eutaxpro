#!/bin/env ruby

require 'optparse'
require 'ostruct'
require 'json'

### Define modules and classes here
def decode_tax_string(tax_string)

    data = {}
    
    tax_string.split(",").each do |e|
        key,value = e.split(":")
        data[key] = value
    end

    return data

end

def filter_hits(list)
    sum = 0
    bucket = []
    list.map {|l| sum += l["count"]}
    list.each do |l|
        if (l["count"].to_f/sum.to_f)*100 > 1.0
            bucket << l
        end
    end

    return bucket
end

### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "Reads Fastq files from a folder and writes a sample sheet to STDOUT"
opts.separator ""
opts.on("-t","--otutable", "=OTUTABLE","Vsearch OTU table with sample counts") {|argument| options.otutable = argument }
opts.on("-s","--sintax", "=SINTAX","Vsearch sintax assignments for OTUs") {|argument| options.sintax = argument }
opts.on("-h","--help","Display the usage information") {
    puts opts
    exit
}

opts.parse! 

# we ignore any results below this coverage
min_coverage = 10

translations = {}

IO.readlines(options.sintax).each do |line|
    
    line.strip!
    otu,tax_string,strand,final_call = line.split("\t")
    otu.gsub!(/;size=.*$/, "")

    if final_call
        taxon = decode_tax_string(final_call)

        if taxon.has_key?("s")
            translations[otu] = taxon["s"]
        elsif taxon.has_key?("g")
            translations[otu] = taxon["g"]
        else
            translations[otu] = "Nicht bis Gatting aufgeloest!"
        end
    else
        translations[otu] = "Kein Treffer in Taxonomie Datenbank"
    end
end

matrix = {}

lines = IO.readlines(options.otutable)

# the first element is the OTU label, don't need that for the header
header = lines.shift.split("\t")[1..-1]

lines.each do |line|
    line.strip!

    elements = line.split("\t")
    
    otu = elements.shift.gsub(/;size=.*$/, "")
    raise otu unless translations.has_key?(otu)

    taxon = translations[otu]
    
    elements.each_with_index do |e,i|
        count = e.to_i
        if count > min_coverage
            sample = header[i]
            payload = { "taxon" => taxon, "count" => count }
            matrix.has_key?(sample) ? matrix[sample] << payload : matrix[sample] = [ payload ]
        end

    end
end

matrix.each do |sample,hits|
    sum = 0.0
    hits.collect {|h| sum  += h["count"]}

    final_hits = filter_hits(hits)
    result = []
    final_hits.map {|f| result << "#{f['taxon']}:#{((f['count'].to_f/sum)*100).round(2)}%" }
    puts "#{sample}\t#{result.join(',')}"
end


