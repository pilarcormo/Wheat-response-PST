require 'csv'
require 'pp'


f = File.open("snps1.csv", "r")
code = File.open("varieties_numbercode.csv", "r")

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
	pp marker
	pp score 
	if james_marker_score.has_key?(marker) 
		james_marker_score[marker].each do |value|
			val = value.to_i
			if val == score
				array << 1
			elsif val != score
				array << 0
			end 
		end 
		output_marker_dic.store(marker, array)	
	end 	
	array = []	

end 



File.open("final_scores.csv", "w+") do |f|
	var.each do |name|
		f << name + ","
	end 
	output_marker_dic.each do |marker, score|
		f << "\n" + marker  
		
		score.each do |value|
			f << "," + value.to_s

		end
		# f << "\n"
	end 

end


# CSV.open("final_scores.csv", 'ab') do |csv|
# 	var = variety[7, -2]
#   csv << [var]
# end

# CSV.open("final_scores.csv", 'ab') do |csv|
#   csv << [i, marker, numeritos 