require 'csv'
require 'pp'

f = File.open("SNPAndContig.csv", "r")
reference = File.open("LIB12451_Reference_greater_2x.tab", "r")
snps = File.open("SNPs/20X/LIB12451_SNP_freq_20x.txt", 'r')


# f = File.open("test_markers", "r")
# reference = File.open("test_tab", "r")
# snps = File.open("test_20x", 'r')

pos_bases = {}
contig_pos_snp = {}
marker_pos = {}
marker_contig_pos = {}
con_pos_snp_poly = {}
pos_alt = {}

####SNPs with markers file 
f.each_line { |line|
  fields = line.split("\t")
  marker = fields[0].strip 
  ref_alt = fields[1].strip 
  contig = fields[2].strip 
  pos = fields[3].strip
  pos_bases.store(pos, ref_alt)
  marker_pos.store(pos, marker)
  contig_pos_snp.store(contig, pos_bases)
  marker_contig_pos.store(contig, marker_pos)

 }


####20X SNPs from SNP calling
snps.each_line { |snp| 
	polymorphims = snp.split("\t")
  freq = polymorphims[-1].strip
  if freq.to_f >= 0.85 #prefiltering step 
    alt_base = polymorphims[3].strip
    contig_poly = polymorphims[0].strip
    pos_poly = polymorphims[2].strip
    pos_alt.store(pos_poly, alt_base)
    con_pos_snp_poly.store(contig_poly, pos_alt)
  # elsif freq.length > 5 
  #   pp freq
  #   i = i +1
    # if freq 
    #  two_bases = polymorphims[5].strip

  
  end 
}



final_contig_pos = {}
final_pos_bases = {}
bases = []

####Reference tab file 
reference.each_line { |ref| 
	columns = ref.split("\t")
	ref_sample = columns[2].strip
	contig_sample = columns[0].strip
	pos_sample = columns[1].strip
	if con_pos_snp_poly.has_key?(contig_sample)
		if con_pos_snp_poly[contig_sample].has_key?(pos_sample)
      base_alt = con_pos_snp_poly[contig_sample][pos_sample] 
      # bases << ref_sample + base_alt 
      final_pos_bases.store(pos_sample, ref_sample + base_alt)
      final_contig_pos.store(contig_sample, final_pos_bases)
      # pp final_contig_pos
 		end
	end 
}





final_contig_pos.each do |contig, pos_base|
  if contig_pos_snp.has_key?(contig)
    pos_base.each do |pos, ref_alt_snp|
      if contig_pos_snp[contig].has_key?(pos)
        base_ref = contig_pos_snp[contig][pos][0]
        base_alt = contig_pos_snp[contig][pos][-1]
        name_marker = marker_contig_pos[contig][pos]
        if base_alt == ref_alt_snp[-1] && base_ref == ref_alt_snp[0]
            CSV.open("snps_good.csv", 'ab') do |csv|
              csv << [name_marker, contig, pos, base_ref, base_alt, 2]
            end
        elsif base_alt == ref_alt_snp[0]
          CSV.open("snps_good.csv", 'ab') do |csv|
            csv << [name_marker, contig, pos, ref_alt_snp[0], base_alt, 0]
          end
        end 
      end
    end 
  end 
end 



