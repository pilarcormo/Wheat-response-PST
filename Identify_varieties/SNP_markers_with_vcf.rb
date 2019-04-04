require 'csv'
require 'pp'
require_relative 'vcf'


dataset = ARGV[0] 
path = ARGV[1] #path to the library 


f = File.open("SNPAndContig.csv", "r")
reference = File.open("2016/#{dataset}_Reference_greater_2x.tab", "r")
#snps = File.open("#{path}/#{dataset}/Trimmed/SNPs/10X/#{dataset}_SNP_freq_10x.txt", 'r')

vcf_file = File.open("2016/#{dataset}.vcf") 

con_pos_snp_poly, vcfsinfo, vcfspos, vcfs_chrom = Vcf.open_vcf(vcf_file)

con_pos_type = Vcf.type_per_pos(vcfsinfo, vcfspos, vcfs_chrom, con_pos_snp_poly)

#pp con_pos_type


marker_contig_pos = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
contig_pos_snp = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }

variants_markers = {}

###markers info
decode.each_line { |marker|
  variants = marker.split("\t")
  name = variants[0]
  refe = variants[1].strip
  hom = variants[2].strip
  het = variants[3].strip
  variants_markers[name] = []
  variants_markers[name] << refe
  variants_markers[name] << hom
  variants_markers[name] << het
}

#pp variants_markers

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

pp marker_contig_pos

final_contig_pos = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }

#Reference tab file 
reference.each_line { |ref| 
  columns = ref.split("\t")
  ref_sample = columns[2].strip
  contig_sample = columns[0].strip
  pos_sample = columns[1].strip
  final_contig_pos[contig_sample][pos_sample] = ref_sample
}



final_contig_pos.each do |contig, pos_base|
  pos_base.each do |pos, alt|
    name_marker = marker_contig_pos[contig][pos]
    if variants_markers.has_key?(name_marker)
      ref = variants_markers[name_marker][0]
      pp ref
      if ref[0] == alt
        CSV.open("snp_markers_#{dataset}.csv", 'a') do |csv|
            csv << [name_marker, contig, pos, ref, 0]
        end
      end
    end 
  end 
end 


contig_pos_snp.each do |contig, pos_base|
  pos_base.each do |pos, ref_alt_snp|
    name_marker = marker_contig_pos[contig][pos]
    if con_pos_snp_poly[contig].has_key?(pos)
      expected = variants_markers[name_marker][1] 
      snp = con_pos_snp_poly[contig][pos]
      if expected[0] == snp
        pos = pos.to_i
         if con_pos_type[contig][pos] == ["HET"]
          CSV.open("snp_markers_#{dataset}.csv", 'a') do |csv|
            csv << [name_marker, contig, pos, snp, 1]
          end
        elsif con_pos_type[contig][pos] == ["HOM"]
          CSV.open("snp_markers_#{dataset}.csv", 'a') do |csv|
            csv << [name_marker, contig, pos, snp, 2]
          end
      end 
    end 
  end 
  end
end 








