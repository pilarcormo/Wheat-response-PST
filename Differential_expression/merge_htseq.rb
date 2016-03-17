
require "pp"


ID = ["LIB18260", "LIB19021", "LIB18683", "LIB19042", "LIB19332", "LIB19340", "LIB19530"]

ID.each do |lib|
	f = File.open("I#{lib}_2.htseq.txt", "r")