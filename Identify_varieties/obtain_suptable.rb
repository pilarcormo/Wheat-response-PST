require 'csv'
require 'pp'

###Reads need to be aligned against the Wheat_SNP_markers 
##Usage ruby SNP_markers.rb LIB10250 2015/new_samples

dataset = ARGV[0] #name of the library
path = ARGV[1] #path to the library 

####change the paths is needed

f = File.open("SNPAndContig.csv", "r")
decode = File.open("SNP_A_B_info.txt", 'r')
code = File.open("varieties_numbercode.csv", "r")

variety = []
code.each_line { |line|
  variety = line.strip.split(",")
    break 
 }
variants_markers = {} #open hash

###markers info. The decode file contains the reference bases and alt bases for each marker. 
decode.each_line { |marker|
  variants = marker.split("\t") #split columns 
  name = variants[0] #first column is the name of the markers
  refe = variants[1] #second column is the reference base (AA)
  hom = variants[2] #third column is the homozygous SNP (TT)
  het = variants[3].strip #fourth column is the heterozygous SNP (AT)
  variants_markers[name] = [] #open empty hash. The nam of the marker as key 
  variants_markers[name] << refe #for each name, all the bases as values
  variants_markers[name] << hom
  variants_markers[name] << het
}


marker_contig_pos = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } #open hash with several keys per value
contig_pos_snp = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) } 

####SNPs with markers file. This file contains:  name of marker, the reference and alternative bases, 
#the contig where the marker is and the position
f.each_line { |line|
  fields = line.split("\t")
  marker = fields[0].strip 
  ref_alt = fields[1].strip 
  contig = fields[2].strip 
  pos = fields[3].strip
  contig_pos_snp[contig][pos] = ref_alt #hash for contig name, position and the ref base and alternative
  marker_contig_pos[contig][pos] = marker #has for contig name, position and the name of the marker

 }


code.each_line { |line|
  variety = line.strip.split(",")
    break 
 }

james_marker_score = {}
scores = []
code.each_with_index { |line, index|
  next if index == 0 
  fields = line.split(",")
  marker = fields[1].strip 
  Array(2..fields.length-1).each do |i|
    scores << fields[i].strip
  end 
  james_marker_score.store(marker, scores)
  scores = []
 }

array = []
output_marker_dic = {}
name_base = {}

contig_pos_snp.each do |contig, pos_base|
  pos_base.each do |pos, ref_alt_snp|
  name_marker = marker_contig_pos[contig][pos]
    if variants_markers.has_key?(name_marker)
      ref = variants_markers[name_marker][0]
      hom = variants_markers[name_marker][1]
      het = variants_markers[name_marker][2]
      array << ref
      array << hom
      name_base.store(name_marker, array)
      array = []
    end 
  end 
end 


base = []



name_base.each do |name_marker, ref_alt|
  hom = name_base[name_marker][1]
  ref = name_base[name_marker][0]
  if james_marker_score.has_key?(name_marker) 
    james_marker_score[name_marker].each do |value|
      val = value.to_i
      if val == 0
        base << ref
      elsif val == 2
        base << hom
      end
      # pp base
      
      
    end 
    output_marker_dic.store(name_marker, base)
    base = []
  end 
end 

#pp output_marker_dic

#         end 
#         output_marker_dic.store(name_marker, base)  
#       end 
#     end
#   end 
# end 
var = variety[1,variety.length]
# ###
File.open("final_scores2.csv", "w+") do |f|
  var.each do |name|
    f << name + ","
  end 
  output_marker_dic.each do |marker, score|
    f << "\n" + marker  
    
    score.each do |value|
      f << "," + value.to_s

    end
  end 

end


# output_marker_dic = {}
# array = []















