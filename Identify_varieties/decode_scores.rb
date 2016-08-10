require 'csv'
require 'pp'

dataset = ARGV[0] 
#path = ARGV[1]
f = File.open("snp_markers_#{dataset}.csv", "r")
code = File.open("varieties_numbercode_example.csv", "r")

sample_marker_score = {}
james_marker_score = {}

variety_marker_score = {}

markers = []
variety = []

code.each_line { |line|
	variety = line.strip.split(",")
  	break 
 }
 
var = variety[1,variety.length]

 f.each_line { |line|
 	fields = line.split(",")
	marker = fields[0]
	score = fields[-1].strip.to_i
	sample_marker_score.store(marker, score)
 }


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

output_marker_dic = {}
array = []

sample_marker_score.each do |marker, score|
	if james_marker_score.has_key?(marker) 
		james_marker_score[marker].each do |value|
			val = value.to_i
			if val == 2
				if score == 2
					array << 1
				elsif score ==1
					array << 0.5
				elsif score == 0
					array << 0
				end
			elsif val == 0  
				if score == 2
					array << 0
				elsif score ==1
					array << 0
				elsif score == 0
					array << 1
				end
			end 
		end 
		output_marker_dic.store(marker, array)	
	end 	
	array = []	

end 

File.open("final_scores_#{dataset}_last.csv", "w+") do |f|
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





