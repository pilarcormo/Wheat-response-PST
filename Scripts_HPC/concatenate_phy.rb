#!/usr/bin/env ruby


### This script concatenates several files in phylip format into one fasta file. 
## Identifiers for genes will be ignored, the sample names must be the same 
### Usage ruby concatenate_phy.rb list.txt 5

require 'pp'
list_files = ARGV[0] ###text file with all the .3rdcodon.phy files 
number_of_samples = ARGV[1].to_i ####number of samples as integer

phy_files = []



list_files = File.open(list_files, "r")
list_files.each_line { |line| 
	phy_file = line.split("\n")
	phy_file.to_s 
	phy_files << phy_file
}
phy_files.flatten!

concatenate = {}
lengths = []

phy_files.each do |file|
	f = File.open(file, "r")
	lib_names, alignments = [], []
	lib_align = {}
	f.each_with_index { |line, index|
		next if index == 0
		libs = line.split(" ")[0]
		lib_names << libs 
		alig = line.split(" ")[1]
		lengths << alig.length
		alignments << alig
		if lengths.uniq.size.to_i == 1
			lib_align.store(libs, alig)
		end 
	}
	lengths = []
	next if lib_align.keys.length != number_of_samples
	lib_align.each do |lib, align|
		if concatenate.has_key?(lib)
			concatenate[lib] << align
		else 
			concatenate.store(lib, align)
		end 
	end 
end 


#pp concatenate
concatenate.each do |lib, concatenation|
	File.open("concatenate_#{lib}.fasta", "w") do |f|
		f.write(">" + lib + "\n" + concatenation + "\n")
	end 
end 