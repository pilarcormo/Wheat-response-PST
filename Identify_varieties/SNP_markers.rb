require 'csv'
require 'pp'

f = File.open("SNPAndContig.csv", "r")
reference = File.open("LIB4468/LIB4468_Reference_greater_2x.tab", "r")
snps = File.open("LIB4468/SNPs/10X/LIB4468_SNP_freq_10x.txt", 'r')


pos_bases = {}
marker_pos = {}
marker_contig_pos = {}

marker_contig_pos = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
contig_pos_snp = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }

####SNPs with markers file 
f.each_line { |line|
  fields = line.split("\t")
  marker = fields[0].strip 
  ref_alt = fields[1].strip 
  contig = fields[2].strip 
  pos = fields[3].strip
  contig_pos_snp[contig][pos] = ref_alt 
  marker_contig_pos[contig][pos] = marker

 }


con_pos_snp_poly = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
pp "Finished with the SNP markers!"


####20X SNPs from SNP calling
snps.each_line { |snp| 
  polymorphims = snp.split("\t")
  freq = polymorphims[-1].strip
  alt_base = polymorphims[5].strip
  contig_poly = polymorphims[0].strip
  pos_poly = polymorphims[2].strip
  con_pos_snp_poly[contig_poly][pos_poly] = alt_base
}


pp "Finished opening the 10x SNP file"


final_contig_pos = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }

##Reference tab file 
reference.each_line { |ref| 
  columns = ref.split("\t")
  ref_sample = columns[2].strip
  contig_sample = columns[0].strip
  pos_sample = columns[1].strip
  if con_pos_snp_poly.has_key?(contig_sample)
    if con_pos_snp_poly[contig_sample].has_key?(pos_sample)
      base_alt = con_pos_snp_poly[contig_sample][pos_sample] 
      final_contig_pos[contig_sample][pos_sample] = base_alt
    end
  end 
}
pp "Finished combining 10x and reference files"




final_contig_pos.each do |contig, pos_base|
  if contig_pos_snp.has_key?(contig)
    pos_base.each do |pos, ref_alt_snp|
      if contig_pos_snp[contig].has_key?(pos)
        base_ref = contig_pos_snp[contig][pos][0]
        base_alt = contig_pos_snp[contig][pos][-1]
        name_marker = marker_contig_pos[contig][pos]
        if ref_alt_snp[-2] == ref_alt_snp[-1]   
          if base_alt == ref_alt_snp[-1]
            CSV.open("snps_good_last.csv", 'ab') do |csv|
              csv << [name_marker, contig, pos, base_ref, base_alt, 2]
            end
          elsif base_alt == ref_alt_snp[0]

            CSV.open("snps_good_last.csv", 'ab') do |csv|
              csv << [name_marker, contig, pos, ref_alt_snp[0], base_alt, 0]
            end
          end
        elsif ref_alt_snp[-2] != ref_alt_snp[-1]   
          if ref_alt_snp[-1] == base_alt || ref_alt_snp[-2] == base_alt
            CSV.open("snps_good_last.csv", 'ab') do |csv|
              csv << [name_marker, contig, pos, base_ref, base_alt, 1]
            end
          end 
        end 
      end 
    end 
  end 
end 


pp "Done"


