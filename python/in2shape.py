import sys
import time
import general
import bed
import metrn
import cetrn
import modencode
import optparse
import numpy
import copy
import ucsc
import os
import pdb

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())

""" Internal operation functions """
def mapplus(x, foo):
	nlist = []
	for bar in foo:
		nlist.append(bar+x)
	return nlist

def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "path from script to files")
	parser.add_option("--infile", action = "store", type = "string", dest = "infile", help = "input file to convert to bed format")
	parser.add_option("--folder", action = "store", type = "string", dest = "folder", help = "folder location of input files")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "mode of operation")
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "output file name")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--nuclear", action = "store", type = "string", dest = "nuclear", help = "Peaks are only nuclear?", default="ON")
	parser.add_option("--parameters", action = "store", type = "string", dest = "parameters", help = "modification parameters")
	parser.add_option("--header", action = "store", type = "string", dest = "header", help = "Is there a header?", default="OFF")
	parser.add_option("--target", action = "store", type = "string", dest = "target", help = "Target identifier", default="feature")
	parser.add_option("--cutChr", action = "store", type = "string", dest = "cutChr", help = "Should first 3 letters (chr) be removed?", default="OFF")
	parser.add_option("--server", action = "store", type = "string", dest = "server", help = "Are we on the server?", default="OFF")
	(option, args) = parser.parse_args()
	
	# import paths:
	if option.server == "OFF":
		path_dict = modencode.configBuild(option.path + "/input/" + "configure_path.txt")
	elif option.server == "ON":
		path_dict = modencode.configBuild(option.path + "/input/" + "configure_server.txt")
	
	# specify input and output paths:
	inpath = path_dict["input"]
	extraspath = path_dict["extras"]
	pythonpath = path_dict["python"]
	scriptspath = path_dict["scripts"]
	downloadpath = path_dict["download"]
	fastqpath = path_dict["fastq"]
	bowtiepath = path_dict["bowtie"]
	bwapath = path_dict["bwa"]
	macspath = path_dict["macs"]
	memepath = path_dict["meme"]
	idrpath = path_dict["idr"]
	igvpath = path_dict["igv"]
	testpath = path_dict["test"]
	processingpath = path_dict["processing"]
	annotationspath = path_dict["annotations"]
	peakspath = path_dict["peaks"]
	gopath = path_dict["go"]
	hotpath = path_dict["hot"]
	qsubpath = path_dict["qsub"]
	
	# standardize paths for analysis:
	alignerpath = bwapath
	indexpath = alignerpath + "index/"
	alignmentpath = alignerpath + "alignment/"
	qcfilterpath = alignerpath + "qcfilter/"
	qcmergepath = alignerpath + "qcmerge/"
	
	# import configuration dictionaries:
	source_dict = modencode.configBuild(inpath + "configure_source.txt")
	method_dict = modencode.configBuild(inpath + "configure_method.txt")
	context_dict = modencode.configBuild(inpath + "configure_context.txt")
	
	# define organism parameters:
	if option.organism == "hs" or option.organism == "h.sapiens":
		organismTag = "hs"
	elif option.organism == "mm" or option.organism == "m.musculus":
		organismTag = "mm"
	elif option.organism == "ce" or option.organism == "c.elegans":
		organismTag = "ce"
	elif option.organism == "dm" or option.organism == "d.melanogaster":
		organismTag = "dm"
	
	# specify genome size file:
	if option.organism != "OFF":
		if option.nuclear == "ON":
			chromosomes = metrn.chromosomes[organismTag]["nuclear"]
			genome_size_file = option.path + "/input/" + metrn.reference[organismTag]["nuclear_sizes"]
			genome_size_dict = general.build_config(genome_size_file, mode="single", separator="\t", spaceReplace=True)
		else:
			chromosomes = metrn.chromosomes[organismTag]["complete"]
			genome_size_file = option.path + "/input/" + metrn.reference[organismTag]["complete_sizes"]
			genome_size_dict = general.build_config(genome_size_file, mode="single", separator="\t", spaceReplace=True)
	
	# determine where the files are stored:
	path = path_dict[option.folder]
	
	# obtain the input file for conversion and mode information:
	infile = path + option.infile
	outfile = annotationspath + "in2shape_" + option.name
	if "windows" == option.mode:
		outfile = annotationspath + "in2shape_" + option.infile.replace(".nh.bed","").replace(".nh.gff","").replace(".bed","").replace(".gff","").replace("in2shape_", "") + "_windows_" + option.name
	
	# read input lines into memory:
	if not "spell:" in option.mode and not "adjust" == option.mode and not "wbGene:protein" == option.mode and not "slopbed" == option.mode and not "waterston:tss" == option.mode:
		indata = open(infile, "U")
	
	# set the output files and adjustment dictionary if necessary:
	if not "adjust" == option.mode and not "slopbed" == option.mode and not "akundaje:report" == option.mode:
		b_outfile = outfile + ".bed"
		g_outfile = outfile + ".gff"
		bh_outfile = outfile + ".nh.bed"
		gh_outfile = outfile + ".nh.gff"
		
		b_output = open(b_outfile, "w")
		g_output = open(g_outfile, "w")
		bh_output = open(bh_outfile, "w")
		gh_output = open(gh_outfile, "w")
	
	elif "adjust" == option.mode:
		adjust_dict = dict()
		adjust_handle = "_adjust"
		for parameter in option.parameters.split(","):
			method, adjustment = parameter.split(":")
			if method == "up":
				adjust_handle += "_up" + adjustment
				adjust_dict[method] = int(adjustment)
			if method == "dn":
				adjust_handle += "_dn" + adjustment
				adjust_dict[method] = int(adjustment)
			if method == "window":
				adjust_handle += "_wd" + adjustment
				adjust_dict["up"] = int(adjustment)
				adjust_dict["dn"] = int(adjustment)
			if method == "center":
				adjust_handle += "_cn" + adjustment
				adjust_dict["cn"] = int(adjustment)
		
		b_output = open(outfile + adjust_handle + ".bed", "w")
		g_output = open(outfile + adjust_handle + ".gff", "w")
		bh_output = open(outfile + adjust_handle + ".nh.bed", "w")
		gh_output = open(outfile + adjust_handle + ".nh.gff", "w")
	
	
	elif "slopbed" == option.mode:
		adjust_handle = "_slopbed_" + option.parameters.replace("-b ", "wd").replace("-l ", "up").replace("-r ", "dn").replace(" ", "_").replace("-s", "").rstrip("_")
		b_outfile = outfile + adjust_handle + ".bed"
		g_outfile = outfile + adjust_handle + ".gff"
		bh_outfile = outfile + adjust_handle + ".nh.bed"
		gh_outfile = outfile + adjust_handle + ".nh.gff"
		
	
	# convert input ensGene annotation data into bed format:
	if "ensGene" in option.mode and not "fly-ensGene" in option.mode:
	
		# set mode and target features:
		mode, target = option.mode.split(":")
		
		# set empty parameters:
		score = "0"
		frame = "."
		group = "touch"
		
		# print header:
		print >>b_output, "\t".join(["#chrm","start","end","feature","score","strand","exons","exon_starts","exon_ends"])
		print >>g_output, "\t".join(["#chrm","source","feature","start","end","score","strand","frame","group","exons","exon_starts","exon_ends"])
		
		# process lines:
		isoforms, genes, features = list(), list(), list()
		inline = indata.readline().replace("\r","\n").strip("\n")
		while inline:
			if not inline == "":
				items = inline.rstrip("\n").split("\t")
				num0, isoform, chrm, strand, start, end, num1, num2, exons, exon_starts, exon_ends, num3, gene, class1, class2, num4 = items
				
				# collect isoforms and genes:
				isoforms.append(isoform)
				genes.append(gene)
				
				# define type of features to store:
				if target == "gene":
					feature = gene
				elif target == "isoform":
					feature = isoform
				
				# rename chromosome to ws220 names:
				chrm = cetrn.chrm2code_dict[chrm]
				
				# process new features not added:
				if not feature in features:
					features.append(feature)
						
					# correct start coordinates from the input ensGene file:
					exon_starts = map(int, exon_starts.rstrip(",").split(","))
					exon_starts = ",".join(map(str, mapplus(1, exon_starts))) + ","
			
					# export data:
					print >>b_output, "\t".join([chrm, str(int(start)+1), end, feature, score, strand, exons, exon_starts, exon_ends])
					print >>g_output, "\t".join([chrm, "in2shape", feature, str(int(start)+1), end, score, strand, frame, group, exons, exon_starts, exon_ends])
					print >>bh_output, "\t".join([chrm, str(int(start)+1), end, feature, score, strand, exons, exon_starts, exon_ends])
					print >>gh_output, "\t".join([chrm, "in2shape", feature, str(int(start)+1), end, score, strand, frame, group, exons, exon_starts, exon_ends])
			
			# reload line
			inline = indata.readline().replace("\r","\n").strip("\n")
			
		print "Genes:", len(set(genes))
		print "Isoforms:", len(set(isoforms))
	
	# extract WormBase annotations from GFF:
	elif "wormbase:" in option.mode:
		
		# preset counters and feature accumulator:
		a, b, c, processed = 0, 0, 0, list()
		
		# print headers:
		print >>b_output, "\t".join(["#chrm","start","end","feature","score","strand","name","gene","feature.type", "feature.key", "feature.source"])
		print >>g_output, "\t".join(["#chrm","source","feature","start","end","score","strand","frame","name","gene","feature.type", "feature.key", "feature.source"])
		
		# define target features of interest:
		if option.mode == "wormbase:complete":
			targets = ['CDS', 'protein_coding_primary_transcript', 'ncRNA_primary_transcript', 'snoRNA', 'Pseudogene', 'tRNA', 'miRNA_primary_transcript', 'snRNA', 'rRNA_primary_transcript', 'intron', 'exon', 'TSS', 'SL1_acceptor_site', 'SL2_acceptor_site', 'transcription_end_site', 'polyA_site', 'polyA_signal_sequence', 'five_prime_UTR', 'three_prime_UTR', 'snlRNA', 'DNAse_I_hypersensitivity']
		elif option.mode == "wormbase:gene":
			targets = ['CDS']
		elif option.mode == "wormbase:transcript" or option.mode == "wormbase:tss" or option.mode == "wormbase:tes":
			targets = ['protein_coding_primary_transcript']
		
		# define target feature count dict:
		feature_count_dict = dict()
		
		# process lines:
		feature_types = list()
		inline = indata.readline().replace("\r","\n").strip("\n")
		while inline:
			items = inline.rstrip("\n").split("\t")
			if not inline == "" and len(items) > 1:
				
				chrm, feature_source, feature_type, start, end, score, strand, frame = items[:8]
				chrm = chrm.replace("CHROMOSOME_","")
				
				if option.cutChr == "ON":
					chrm = chrm.lstrip("chr")
				
				# standardize chromosome names:
				if chrm in cetrn.chrm2code_dict:
					chrm = cetrn.chrm2code_dict[chrm]
				
				# check chromosome names:
				if not chrm in metrn.chromosomes[organismTag]["complete"]:
					print "Error with chromosome name:", chrm
					pdb.set_trace()
						
				# preset processing booleans:
				curated, target, process = False, False, False
				
				# check type of feature coming in:
				if feature_source == "curated":
					a += 1
					curated = True
				if feature_type in targets:
					b += 1
					target = True
				if target:
					c += 1
					process = True
				
				# determine if notes exist
				if len(items) > 8:
					group = items[8]
				else:
					group == ""
					
				# if this is targeted feature, process data	
				if process:
				
					# further filter annotations based on source and characteristics:
					proceed = False
					
					# gather group data:				
					info = dict()
					for data in group.split(" ; "):
						data = [element.strip(" ") for element in data.strip('"').split('"')]
						if not data == [""] and len(data) == 2:
							key, entry = data
							info[key] = entry
					
					# load feature types:
					if not feature_type in feature_count_dict:
						feature_count_dict[feature_type] = 0
						#print feature_type, info
						#pdb.set_trace()
					
					feature = "NA"
					gene = "NA"
					name = "NA"
					
					if feature_type == "CDS" and curated:
						feature_key = feature_type
						feature = info[feature_key]
						gene = info["Gene"]
						if "Locus" in info:
							name = info["Locus"].upper()
						else:
							name = info["CDS"]
						proceed = True
					
					elif feature_type in ["protein_coding_primary_transcript", "ncRNA_primary_transcript", "miRNA_primary_transcript", "rRNA_primary_transcript", "snoRNA", "tRNA", "snRNA", "snlRNA"]:
					    feature = info["Transcript"]
					    feature_key = "transcript"
					    if "Gene" in info:
					    	gene = info["Gene"]
					    proceed = True
					
					elif feature_type in ["intron", "exon","five_prime_UTR", "three_prime_UTR"] and feature_source == "Coding_transcript":
					    feature = info["Transcript"]
					    feature_key = feature_type
					    proceed = True
					    
					elif feature_type in ["TSS", "transcription_end_site"] and feature_source == "RNASeq":
					    feature = info["Feature"]
					    feature_key = feature_type.replace("transcription_end_site", "TES")
					    proceed = True
					
					elif feature_type in ["polyA_site", "polyA_signal_sequence"] and feature_source in ["polyA_site", "polyA_signal_sequence"]: 
					    feature = info["Feature"]
					    feature_key = feature_type
					    proceed = True
					
					elif feature_type in ["SL1_acceptor_site", "SL2_acceptor_site"]:
					    feature = info["Feature"]
					    feature_key = "acceptor_site"
					    proceed = True
					
					elif feature_type in ["Pseudogene"] and feature_source == "Pseudogene":
					    feature = info["Pseudogene"]
					    feature_key = feature_type
					    proceed = True
					
					elif feature_type in ["DNAse_I_hypersensitivity"]:
					    feature = "DHS." + str(feature_count_dict[feature_type])
					    feature_key = "DHS"
					    proceed = True
						
					if proceed:
						
						# update feature count dictionary:
						feature_count_dict[feature_type] += 1
						
						# check to see if we got the feature:
						if not feature == "NA":
							export = True
						
						# check uniqueness
						#if feature in processed:
						#	export = False
						#	print feature, "<-- repeated!"
						#	pdb.set_trace
						#else:
						#	processed.append(feature)
						#	export = True
						
						# make corrections for TSS and TES modes:
						if option.mode == "wormbase:tss":
							if strand == "+":
								start = start
								end = str(int(start) + 1)
							else:
								start = str(int(end) - 1)
								end = str(int(start) + 1)
							
						elif option.mode == "wormbase:tes":
							if strand == "+":
								start = str(int(end) - 1)
								end = str(int(start) + 1)
							else:
								start = start
								end = str(int(start) + 1)
							
						
						# export data:
						if export:
							print >>b_output, "\t".join([chrm, start, end, feature, score, strand, name, gene, feature_type.lower(), feature_key.lower(), feature_source.lower()]).lstrip("#")
							print >>g_output, "\t".join([chrm, feature_source.lower(), feature, start, end, score, strand, frame, name, gene, feature_type.lower(), feature_key])
							print >>bh_output, "\t".join([chrm, start, end, feature, score, strand, name, gene, feature_type.lower(), feature_key.lower(), feature_source.lower()]).lstrip("#")
							print >>gh_output, "\t".join([chrm, feature_source.lower(), feature, start, end, score, strand, frame, name, gene, feature_type.lower(), feature_key])
				
			# reload line
			inline = indata.readline().replace("\r","\n").strip("\n")
			
		print feature_count_dict
		
	
	# parse ensGene mode:
	elif "fly-ensGene" in option.mode:
		
		# process lines:
		inline = indata.readline().replace("\r","\n").strip("\n")
		while inline:
			if not inline == "":
				items = inline.rstrip("\n").split("\t")
				bin, feature, chrm, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name, cdsStartStat, cdsEndStat, exonFrames = items
				if option.cutChr == "ON":
					chrm = chrm.lstrip("chr")
				
				if option.mode == "fly-ensGene:TSS":
					start = txStart
					end = str(int(start) + 1)
					strand = "+"
				
				elif option.mode == "fly-ensGene:TSS":
					end = txEnd
					start = str(int(end) - 1)
					strand = "+"
				
				elif option.mode == "fly-ensGene:Transcript":
					start = txStart
					end = txEnd
				
				elif option.mode == "fly-ensGene:Gene":
					start = cdsStart
					end = cdsEnd
					feature = copy.deepcopy(name)
					name = copy.deepcopy(feature)
				
				# export:
				print >>b_output, "\t".join([chrm, start, end, feature, score, strand, name]).lstrip("#")
				print >>g_output, "\t".join([chrm, "in2shape", feature, start, end, score, strand, frame, group, name])
				print >>bh_output, "\t".join([chrm, start, end, feature, score, strand, name]).lstrip("#")
				print >>gh_output, "\t".join([chrm, "in2shape", feature, start, end, score, strand, frame, group, name])
		
			# reload line
			inline = indata.readline().replace("\r","\n").strip("\n")
	
	
	# load protein sequences:
	elif option.mode == "wbGene:protein":
		
		bed_lines = open(path_dict[option.folder] + "in2shape_wbGene.bed").readlines()
		gff_lines = open(path_dict[option.folder] + "in2shape_wbGene.gff").readlines()
		nh_bed_lines = open(path_dict[option.folder] + "in2shape_wbGene.nh.bed").readlines()
		nh_gff_lines = open(path_dict[option.folder] + "in2shape_wbGene.nh.gff").readlines()
		
		print >>b_output, bed_lines.pop(0).strip() + "\tprotein"
		print >>g_output, gff_lines.pop(0).strip() + "\tprotein"
		
		protein_dict = dict()
		inlines = open(inpath + option.infile, "U").readlines()
		for inline in inlines:
			feature, sequence = inline.strip().split("\t")
			protein_dict[feature] = sequence
		
		for index in range(0,len(bed_lines)):
			chrm, start, end, feature = bed_lines[index].strip().split("\t")[:4]
			if option.cutChr == "ON":
				chrm = chrm.lstrip("chr")
				
			if feature in protein_dict:
				print >>b_output, bed_lines[index].strip() + "\t" + protein_dict[feature]
				print >>g_output, gff_lines[index].strip() + "\t" + protein_dict[feature]
				print >>bh_output, nh_bed_lines[index].strip() + "\t" + protein_dict[feature]
				print >>gh_output, nh_gff_lines[index].strip() + "\t" + protein_dict[feature]
	
	
	# convert integrated transcripts (Waterston Lab) GFF to BED format:
	elif option.mode == "waterston:expression":
		
		# process lines:
		k, processed = 1, list()
		inlines = indata.readlines()
		for inline in inlines:
			if not inline.strip() == "" and not "#" == inline[0]:
				items = inline.rstrip("\n").split("\t")
				chrm, program, feature, start, end, score, strand, frame, group = items[:9]
				if option.cutChr == "ON":
					chrm = chrm.lstrip("chr")
					
				if feature == "transcript":
					
					# recover gene assignment and status information:
					info = dict()
					for element in group.split(";"):
						key, data = element.split("=")
						info[key] = data
					output = list()
					for key in ["ID", "transcribed", "normscore", "overlapping_wormbase_transcript", "prediction_status"]:
						if key in info:
							output.append(info[key])
						else:
							output.append("NONE")
							
					# set feature name:
					parts  = info["ID"].split("_")
					feature = parts[len(parts)-1]
					if not feature in processed:
						processed.append(feature)
					else:
						print feature, info["ID"]
						pdb.set_trace()
					
					# export data:
					print >>b_output, "\t".join([chrm, start, end, feature, score, strand] + output).lstrip("#")
					print >>g_output, "\t".join([chrm, "in2shape", feature, start, end, score, strand, frame, group])
					print >>bh_output, "\t".join([chrm, start, end, feature, score, strand] + output).lstrip("#")
					print >>gh_output, "\t".join([chrm, "in2shape", feature, start, end, score, strand, frame, group])
		
			# reload line
			#inline = indata.readline().replace("\r","\n").strip("\n")
	
	
	# convert integrated TSSs (Waterston Lab) GFF to BED format:
	elif option.mode == "waterston:tss":
		
		# process lines:
		i, j, processed = 0, 0, list()
		inlines = open(infile).readlines()
		#inline = indata.readline().replace("\r","\n").strip("\n")
		for inline in inlines:
			if not inline == "" and not "#" == inline[0]:
				initems = inline.strip().split("\t")
				if not initems == [""]:
					chrm, program, feature, start, end, score, strand, frame, group = initems[:9]
					if option.cutChr == "ON":
						chrm = chrm.lstrip("chr")
				
					# find TSS lines:
					if "TSS" == feature:
						i += 1
						
						# recover gene assignment and status information:
						info = dict()
						for element in group.split(";"):
							key, data = element.split("=")
							info[key] = data
						output = list()
						for key in ["ID", "Parent", "prediction_status"]:
							if key in info:
								output.append(info[key])
							else:
								output.append("NONE")
								
						# set feature name:
						feature = info["ID"]
						
						# filter confirmed TSS:
						if info["prediction_status"] == "Confirmed":
							j += 1
						
							# export data:
							print >>b_output, "\t".join([chrm, start, end, feature, score, strand] + output).lstrip("#")
							print >>g_output, "\t".join([chrm, "in2shape", feature, start, end, score, strand, frame, group])
							print >>bh_output, "\t".join([chrm, start, end, feature, score, strand] + output).lstrip("#")
							print >>gh_output, "\t".join([chrm, "in2shape", feature, start, end, score, strand, frame, group])
						
		print i, j			
	
	# import Gu et al. 2012 TSSs to BED format (all TSSs):
	elif option.mode == "gu2012:all":
		
		# process lines:
		k = 1
		inline = indata.readline().replace("\r","\n").strip("\n")
		inline = indata.readline().replace("\r","\n").strip("\n")
		while inline:
			if not inline == "":
				items = inline.rstrip("\n").split("\t")
				#chromosome	strand	start	reads	distance_to_transcript	transcript	covered_by_CapSeq	transcript type
				chrm, strand, kstart, reads, distance, transcript, capseq, ftype = items[:9]
				if option.cutChr == "ON":
					chrm = chrm.lstrip("chr")
				
				# orient coordinates:
				if strand == "+":
					start = str(kstart)
					end = str(int(kstart) + 1)
				else:
					start = str(int(kstart) - 1)
					end = str(int(kstart))
					
				# rename transcript (if unknown):
				if "NA" == transcript:
					transcript = "gu2012." + str(k)
					k += 1
				geneID = transcript
				score = "0"
				frame = capseq
				group = ftype
				
				# export:
				print >>b_output, "\t".join([chrm, start, end, transcript, score, strand, geneID, frame, group, distance]).lstrip("#")
				print >>g_output, "\t".join([chrm, "in2shape", transcript, start, end, score, strand, frame, group, distance])
				print >>bh_output, "\t".join([chrm, start, end, transcript, score, strand, geneID, frame, group, distance]).lstrip("#")
				print >>gh_output, "\t".join([chrm, "in2shape", transcript, start, end, score, strand, frame, group, distance])
		
			# reload line
			inline = indata.readline().replace("\r","\n").strip("\n")
	
	# import Gu et al. 2012 TSSs to BED format (coding TSSs only):
	elif option.mode == "gu2012:pro":
		
		# process lines:
		k = 1
		inline = indata.readline().replace("\r","\n").strip("\n")
		inline = indata.readline().replace("\r","\n").strip("\n")
		while inline:
			if not inline == "":
				items = inline.rstrip("\n").split("\t")
				#chromosome	strand	start	reads	distance_to_transcript	transcript	covered_by_CapSeq	transcript type
				chrm, strand, kstart, reads, distance, transcript, capseq, ftype = items[:9]
				if option.cutChr == "ON":
					chrm = chrm.lstrip("chr")
				
				# orient coordinates:
				if strand == "+":
					start = str(kstart)
					end = str(int(kstart) + 1)
				else:
					start = str(int(kstart) - 1)
					end = str(int(kstart))
					
				# rename transcript (if unknown):
				if "NA" == transcript:
					transcript = "gu2012." + str(k)
					k += 1
				geneID = transcript
				score = "0"
				frame = capseq
				group = ftype
				
				# export:
				if ftype.lower() == "coding":
					print >>b_output, "\t".join([chrm, start, end, transcript, score, strand, geneID, frame, group, distance]).lstrip("#")
					print >>g_output, "\t".join([chrm, "in2shape", transcript, start, end, score, strand, frame, group, distance])
					print >>bh_output, "\t".join([chrm, start, end, transcript, score, strand, geneID, frame, group, distance]).lstrip("#")
					print >>gh_output, "\t".join([chrm, "in2shape", transcript, start, end, score, strand, frame, group, distance])
		
			# reload line
			inline = indata.readline().replace("\r","\n").strip("\n")
	
	# convert GFF to BED format:
	elif option.mode == "gff2bed":
		
		# process lines:
		inline = indata.readline().replace("\r","\n").strip("\n")
		while inline:
			if not inline == "":
				items = inline.rstrip("\n").split("\t")
				chrm, program, feature, start, end, score, strand, frame, group = items[:9]
				if option.cutChr == "ON":
					chrm = chrm.lstrip("chr")
					
				# export:
				print >>b_output, "\t".join([chrm, start, end, feature, score, strand]).lstrip("#")
				print >>g_output, "\t".join([chrm, "in2shape", feature, start, end, score, strand, frame, group])
				print >>bh_output, "\t".join([chrm, start, end, feature, score, strand]).lstrip("#")
				print >>gh_output, "\t".join([chrm, "in2shape", feature, start, end, score, strand, frame, group])
		
			# reload line
			inline = indata.readline().replace("\r","\n").strip("\n")
	
	# import GENCODE v10 to BED format:
	elif option.mode == "gencode:tss":
		
		# process lines:
		rnaIDs = list()
		inline = indata.readline().replace("\r","\n").strip("\n")
		while inline:
			if not inline == "":
				items = inline.rstrip("\n").split("\t")
				chrm, program, feature, start, end, score, strand, frame, group = items[:9]
				if option.cutChr == "ON":
					chrm = chrm.lstrip("chr")
				
				if "trlist" in group:
					info = group.split(",")[0]
					genekey, geneID, rnakey, rnaID = info.split(" ")
					
					if rnaID in rnaIDs:
						print geneID, rnaID
						pdb.set_trace()
					else:
						rnaIDs.append(rnaID)
					
					# export:
					print >>b_output, "\t".join([chrm, start, end, rnaID, score, strand, geneID]).lstrip("#")
					print >>g_output, "\t".join([chrm, "in2shape", rnaID, start, end, score, strand, frame, group])
					print >>bh_output, "\t".join([chrm, start, end, rnaID, score, strand, geneID]).lstrip("#")
					print >>gh_output, "\t".join([chrm, "in2shape", rnaID, start, end, score, strand, frame, group])
		
			# reload line
			inline = indata.readline().replace("\r","\n").strip("\n")
	
	# import Celniker Lab TSS to BED format:
	elif option.mode == "celniker:tss":
		
		# process lines:
		gene2rna_dict = dict()
		inline = indata.readline().replace("\r","\n").strip("\n")
		while inline:
			if not inline == "":
				items = inline.rstrip("\n").split("\t")
				chrm, program, feature, start, end, score, strand, frame, group = items[:9]
				if option.cutChr == "ON":
					chrm = chrm.lstrip("chr")
				
				if "gene_id" in group:
					geneID = group.split(";")[0].replace('gene_id "', '').replace('"', '')
					if not geneID in gene2rna_dict:
						gene2rna_dict[geneID] = list()
					rnaID = geneID + "." + str(len(gene2rna_dict[geneID]) + 1)
					gene2rna_dict[geneID].append(rnaID)
					
					# export:
					print >>b_output, "\t".join([chrm, start, end, rnaID, score, strand, geneID]).lstrip("#")
					print >>g_output, "\t".join([chrm, "in2shape", rnaID, start, end, score, strand, frame, group])
					print >>bh_output, "\t".join([chrm, start, end, rnaID, score, strand, geneID]).lstrip("#")
					print >>gh_output, "\t".join([chrm, "in2shape", rnaID, start, end, score, strand, frame, group])
		
			# reload line
			inline = indata.readline().replace("\r","\n").strip("\n")
	
	# as-is import BED file mode:
	elif option.mode == "import":
	
		if ".nh.bed" in infile:
			outfile = bh_outfile
		elif ".bed" in infile:
			outfile = b_outfile
		elif ".nh.gff" in infile:
			outfile = gh_outfile
		elif ".gff" in infile:
			outfile = g_outfile
		
		if option.cutChr == "OFF":
			command = "cp " + infile + " " + outfile
			os.system(command)
		
		elif option.cutChr == "ON":
			f_output = open(outfile, "w")
			for inline in open(infile).readlines():
				if "chr" == inline[:3]:
					print >>f_output, inline[3:].strip()
				else:
					print >>f_output, inline.strip()
			f_output.close
		
	
	# import enhancer file mode:
	elif option.mode == "enhancers":
	
		if ".nh.bed" in infile:
			outfile = bh_outfile
		elif ".bed" in infile:
			outfile = b_outfile
		elif ".nh.gff" in infile:
			outfile = gh_outfile
		elif ".gff" in infile:
			outfile = g_outfile
		
		# load input lines:
		inlines = open(infile).readlines()
		if option.header == "ON":
			inlines.pop(0)
		
		# convert to bed file:
		f_output = open(b_outfile, "w")
		index = 1
		for inline in inlines:
			inline = inline.replace(" ", "\t")
			if option.cutChr == "ON" and "chr" == inline[:3]:
				chrm, coord = inline[3:].strip().split("\t")
			else:
				chrm, coord = inline.strip().split("\t")
			print >>f_output, "\t".join(map(str, [chrm, int(coord), int(coord) + 1, "enhancer." + str(index), "0", "+"]))
			index += 1
		f_output.close
		
	
	# mix enhancer and promoters mode:
	elif option.mode == "regulatory.mix":
	
		if ".nh.bed" in infile:
			outfile = bh_outfile
		elif ".bed" in infile:
			outfile = b_outfile
		elif ".nh.gff" in infile:
			outfile = gh_outfile
		elif ".gff" in infile:
			outfile = g_outfile
		
		#
		inheader = general.build_header_dict(infile)
		
		# load input lines:
		inlines = open(infile).readlines()
		if option.header == "ON":
			inlines.pop(0)
		
		# generate to mix file:
		f_output = open(b_outfile, "w")
		index = 1
		for inline in inlines:
			initems = inline.strip().split("\t")
			chrm, start, stop, feature, score, strand = initems[:6]
			target = initems[inheader[option.target]]
			print >>f_output, "\t".join(map(str, [chrm, start, stop, feature, score, strand, target]))
			index += 1
		
		
		# load enhancer lines:
		inlines = open(annotationspath + option.parameters).readlines()
		for inline in inlines:
			initems = inline.strip().split("\t")
			chrm, start, stop, feature, score, strand = initems[:6]
			print >>f_output, "\t".join(map(str, [chrm, start, stop, feature, score, strand, feature.split(".")[0]]))
			index += 1
		
		f_output.close
		
		
	
	# slopBed coordinates mode:
	elif option.mode == "slopbed":
		
		# update input path:
		inpath = path_dict[option.folder]
		
		# process files:
		print
		for infile in sorted(os.listdir(inpath)):
			if option.infile in infile and not "_adjust" in infile and not "_slopbed" in infile and not "_windows_" in infile:
				print infile
				if ".nh.bed" in infile:
					outfile = bh_outfile
				elif ".bed" in infile:
					outfile = b_outfile
				elif ".nh.gff" in infile:
					outfile = gh_outfile
				elif ".gff" in infile:
					outfile = g_outfile
				command = "slopBed -i " + inpath + infile + " -g " + genome_size_file + " " + option.parameters + " > " + outfile
				os.system(command)
					
	# window generator mode:
	elif option.mode == "windows":
		
		# process header info:
		if option.header == "ON" or option.header == "bed": 
		
			# load header dictionary:
			if option.header == "ON":
				hd = general.build_header_dict(infile)
			elif option.header == "bed":
				hd = metrn.bedHeader
				
			# export headers:
			bedHeader = ["chrm","start","end","feature","score","strand"]
			gffHeader = ["chrm","source","feature","start","end","score","strand"]
			for column in general.valuesort(hd):
				if not column in bedHeader:
					bedHeader.append(column)
				if not column in gffHeader:
					gffHeader.append(column)
			bedHeader += ["window", "original.start", "original.end"]
			gffHeader += ["window", "original.start", "original.end"]
			print >>b_output, "\t".join(bedHeader)
			print >>g_output, "\t".join(gffHeader)
			
		# configure processing for filetype:
		indata = open(infile, "U")
		if option.header == "ON":
			inline = indata.readline().replace("\r","\n").strip("\n")
				
		# process input file:
		inline = indata.readline().replace("\r","\n").strip("\n")
		while inline:
			initems = inline.split("\t")
			if not initems == [""]:
				chrm, start, end, strand, feature = initems[hd["chrm"]], int(initems[hd["start"]]), int(initems[hd["end"]]), initems[hd["strand"]], initems[hd[option.target]]
				
				for window in option.parameters.split(","):
					minDist, maxDist = map(int, window.split(":"))
					upWindow, dnWindow = [start-maxDist, start-minDist], [end+minDist, end+maxDist]
					if upWindow[1] + 1 == dnWindow[0]:
						newWindows = [[start-maxDist, end+maxDist]]
					else:
						newWindows = [upWindow, dnWindow]
					
					for newWindow in newWindows:
						winStart, winEnd = newWindow
						if winStart < 1:
							winStart = 1
						if winEnd < 1:
							winEnd = 1
						
						bedOutput, gffOutput = list(), list()
						for column in bedHeader:
							if column == "window":
								bedOutput.append(window)
							elif column == "original.start":
								bedOutput.append(start)
							elif column == "original.end":
								bedOutput.append(end)
							elif column == "feature":
								bedOutput.append(initems[hd[option.target]])
							elif column == "start":
								bedOutput.append(winStart)
							elif column == "end":
								bedOutput.append(winEnd)
							elif column in hd:
								bedOutput.append(initems[hd[column]])
							else:
								bedOutput.append("na")
						
						for column in gffHeader:
							if column == "window":
								gffOutput.append(window)
							elif column == "original.start":
								gffOutput.append(start)
							elif column == "original.end":
								gffOutput.append(end)
							elif column == "feature":
								gffOutput.append(initems[hd[option.target]])
							elif column == "start":
								gffOutput.append(winStart)
							elif column == "end":
								gffOutput.append(winEnd)
							elif column in hd:
								gffOutput.append(initems[hd[column]])
							elif column == "source":
								gffOutput.append("in2shape")
							else:
								gffOutput.append("na")
							
						print >>b_output, "\t".join(map(str, bedOutput))
						print >>bh_output, "\t".join(map(str, bedOutput))
						print >>g_output, "\t".join(map(str, gffOutput))
						print >>gh_output, "\t".join(map(str, gffOutput))
						#print "\t".join(map(str, bedOutput))
						#pdb.set_trace()
			
			# load next line...
			inline = indata.readline().replace("\r","\n").strip("\n")
	
	# adjust coordinates mode:
	elif option.mode == "adjust":
		
		def extendStart(start, end, amount, strand):
			if strand == "+":
				start = start-amount
			elif strand == "-":
				end = end+amount
			if start < 1: 
				start = 1
			if end < 1: 
				end = 1
			return start, end
			
		def extendEnd(start, end, amount, strand):
			if strand == "+":
				end = end+amount
			elif strand == "-":
				start = start-amount
			if start < 1: 
				start = 1
			if end < 1: 
				end = 1
			return start, end
				
		
		# determine where the files are stored:
		inpath = path_dict[option.folder]
	
		# define output file dictionary:
		outfile_dict = dict()
		outfile_dict[".bed"] = [b_output, "bed", True]
		outfile_dict[".gff"] = [g_output, "gff", True] 
		outfile_dict[".nh.bed"] = [bh_output, "bed", False] 
		outfile_dict[".nh.gff"] = [gh_output, "gff", False] 
		
		# prepare a header dict for each file type:
		header_type_dict = dict()
				
		# process files:
		print
		for infile in sorted(os.listdir(inpath)):
			if option.name in infile and not "_adjust" in infile:
				
				# detect filetype:
				key = ""
				for extension in [".bed", ".gff", ".nh.bed", ".nh.gff"]:
					if extension in infile:
						key = extension
				
				# configure processing for filetype:
				indata = open(inpath + infile, "U")
				f_output, filetype, hasheader = outfile_dict[key]
				if hasheader:
					header_type_dict[filetype] = general.build_header_dict(inpath + infile)
					print >>f_output, indata.readline().replace("\r","\n").strip("\n").strip()
				outheader = general.valuesort(header_type_dict[filetype])
				
				print "Processing:", infile
				#print key
				#print filetype
				#print outheader
				#print header_type_dict[filetype]
				#print
				#pdb.set_trace()
				
				# process lines:
				inline = indata.readline().replace("\r","\n").strip("\n")
				while inline:
					if not inline == "":
						initems = inline.rstrip("\n").split("\t")
						
						start = int(initems[header_type_dict[filetype]["start"]])
						end = int(initems[header_type_dict[filetype]["end"]])
						strand = initems[header_type_dict[filetype]["strand"]]
						
						if "up" in adjust_dict:
							start, end = extendStart(start, end, adjust_dict["up"], strand)
							initems[header_type_dict[filetype]["start"]] = start
							initems[header_type_dict[filetype]["end"]] = end
						if "dn" in adjust_dict:
							start, end = extendEnd(start, end, adjust_dict["dn"], strand)
							initems[header_type_dict[filetype]["start"]] = start
							initems[header_type_dict[filetype]["end"]] = end
								
						output = list()
						for header in outheader:
							if len(initems) > header_type_dict[filetype][header]:
								value = initems[header_type_dict[filetype][header]]
							else:
								value = ""
							output.append(value)
							#print header, value
						print >>f_output, "\t".join(map(str, output))
						
					inline = indata.readline().replace("\r","\n").strip("\n")
				
		print
					
	# as-is import report file mode:
	elif option.mode == "akundaje:report":
	
		# load header dictionary:
		hd = general.build_header_dict(infile)
		
		f_output = open(outfile + ".txt", "w")
		print >>f_output, "\t".join(["organism", "strain", "factor", "context", "institute", "method", "filename"])
		inlines = open(infile).readlines()
		inlines.pop(0)
		for inline in inlines:
			initems = inline.strip().split("\t")
			if option.organism == "hs":
				strain = "na"
				factor = initems[hd["HGNC TARGET NAME"]]
				context = initems[hd["CELLTYPE"]].replace("Gm12878", "GM12878")
				institute = initems[hd["LAB"]].lower()
				method = initems[hd["PROTOCOL"]].lower()
				filename = initems[hd["FILENAME"]].replace(".bam", "").replace("AlnRep0", "Aln").replace("AlnRep1", "Aln")
			output = [option.organism, strain, factor, context, institute, method, filename]
			print >>f_output, "\t".join(output)
		f_output.close
	
	
	else:
		print "Error: Unrecognized option.mode parameter!"
	
	
	if not "slopbed" == option.mode and not "akundaje:report" == option.mode:
	
		# close output files:
		b_output.close()
		g_output.close()
		bh_output.close()
		gh_output.close()
		
	

"""
Here is a brief description of the GFF fields:

seqname - The name of the sequence. Must be a chromosome or scaffold.
source - The program that generated this feature.
feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
start - The starting position of the feature in the sequence. The first base is numbered 1.
end - The ending position of the feature (inclusive).
score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ".".
strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
group - All lines with the same group are linked together into a single item.
"""

if __name__ == "__main__":
	main()
	print "Completed:", time.asctime(time.localtime())

#python in2shape.py --path ~/meTRN --infile Gencodev10_TSS_May2012.gff --folder extras --mode gff2bed --name hs_TSS
#python in2shape.py --path ~/meTRN --infile in2shape_hs_TSS --folder annotations --mode slopbed --name hs_TSS --organism hs --parameters "-l 5000 -r 500 -s"
