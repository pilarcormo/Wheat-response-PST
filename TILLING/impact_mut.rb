
require 'csv'
require 'pp'
number = ARGV[0] 

##open csv with all the mutations

mutations = File.open("#{number}", "r")

hash = {}
hash_mod = {}
hash_low = {}
high = 0
moderate = 0
low = 0
gene_effect = {}
effect_array = []
genes = []

mutations.each_line { |line|
 	fields = line.split(",")
	gene = fields[12]
	genes << gene
	eff = fields[13]
	if gene_effect.has_key?(gene)
		gene_effect[gene] << eff
	else
	 	gene_effect[gene] = []
	 	gene_effect[gene] << eff 
	end
}





gene_effect.each do |gene, effects|
	high = 0
	moderate = 0
	low = 0
	effects.each do |effect|
		if effect ==  "transcript_ablation"
			high = high + 1 
		elsif effect == "splice_acceptor_variant"
			high = high + 1 
		elsif effect == "splice_donor_variant"
			high = high + 1 
		elsif effect == "stop_gained"
			high = high + 1 
		elsif effect == "frameshift_variant"
			high = high + 1 
		elsif effect == "stop_lost"
			high = high + 1 
		elsif effect == "start_lost"
			high = high + 1 
		elsif effect == "transcript_amplification"
			high = high + 1 
		elsif effect == "inframe_insertion"
			moderate = moderate + 1
		elsif effect == "inframe_deletion"
			moderate = moderate + 1
		elsif effect == "missense_variant"
			moderate = moderate + 1
		elsif effect == "protein_altering_variant"
			moderate = moderate + 1
		end
	end 
	hash[gene] = high
	hash_mod[gene] = moderate
	# hash_low[gene] = low
end 


hash.reject!{|k| k == "gene" } 
hash_mod.reject!{|k| k == "gene" } 

pp hash
pp hash_mod

hash.each do |gene, high|
	CSV.open("#{number}_total_mutations.csv", 'a') do |csv|
		csv << [gene, high, hash_mod[gene]]
	end 
end 





