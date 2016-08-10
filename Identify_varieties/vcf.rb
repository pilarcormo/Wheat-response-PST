class Vcf 
	require 'bio'
	require 'bio-samtools'

	 def self.open_vcf(vcf_file)
	    vcf_ref, vcf_alt, vcfs_chrom, vcfspos, vcfsinfo = [], [], [], [], []
	    con_pos_snp_poly = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
	    con_pos_type = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }

	    File.open(vcf_file, 'r').each do |line|
	    	next if line =~ /^#/
   			polymorphims = line.split("\t")
   			#freq = polymorphims[-1].strip
   			alt_base = polymorphims[4].strip
   			contig_poly = polymorphims[0].strip
   			pos_poly = polymorphims[1].strip
   			con_pos_snp_poly[contig_poly][pos_poly] = alt_base
	  		
	    
	        v = Bio::DB::Vcf.new(line)
	        # vcf_ref << v.ref
	        # vcf_alt << v.alt
	        vcfs_chrom << v.chrom
	        vcfspos << v.pos
	        vcfsinfo << v.info 

	  #       vcfsinfo.each do |hash|
	  #   	 	hash.each do |type, number|
	  #   			if number == "1" 
	  #   			# snps.store(pos_poly, type) 		
	  #   				con_pos_type[contig_poly][pos_poly] = [type]
	  #   			end
	  #   		end
			# end 
	    	#a = line.split("\t")
	        
	        # if v.chrom == "#{chromosome}"
	        # 	new_vcf << line 
	        # end
  			#con_pos_snp_poly[vcfs_chrom][vcfspos] = vcf_alt
	    end 

		return con_pos_snp_poly, vcfsinfo, vcfspos, vcfs_chrom
	end 

	def self.type_per_pos (vcfs_info, vcfs_pos, vcfs_chrom, con_pos_snp_poly)
		con_pos_type = Hash.new{ |h,k| h[k] = Hash.new(&h.default_proc) }
	    #snps = {}

	  	x=0

	    vcfs_info.each do |hash|
	    	hash.each do |type, number|
	    		if number == "1" 
	    			con_pos_type[vcfs_chrom[x]][vcfs_pos[x]] = [type]
	    			#snps.store(vcfs_pos[x], type) 		
	   				x=x+1
	    		end
	    	end

	    end

	    # snps.each do |pos, type|
	    # 	if type == "HET"
	    # 		ht << pos
	    # 	elsif type == "HOM"
	    # 		hm << pos
	    # 	end 
	    # end 
	    return con_pos_type
	end  

	def self.filtering(vcfs_pos_c, snps_p, snps_c, child_chr_vcf)
		short_vcfs_pos_c = vcfs_pos_c
		short_vcfs_pos_c.flatten!
		snps_p.each do |pos, type|
		    if snps_c.has_key?(pos)
		        snps_c.delete(pos) 
		        short_vcfs_pos_c.delete(pos)
		    end 
		end 
		short_child_chr_vcf = []
		child_chr_vcf.each do |line|
		    position = line.split("\t")[1].to_i
		    if short_vcfs_pos_c.include?(position) 
		        short_child_chr_vcf << line
		    end 
		end 
		return short_child_chr_vcf
	end 
end 