#!/usr/bin/env python
# Make profile plots from pairs of bed files; based on Doug profilePlot.py

import sys
import time
import optparse
import general
import numpy
import pickle
import pdb
import metrn
import itertools
import modencode
import os
import genomica
import gs

#import HTSeq
#import itertools
#from matplotlib import pyplot

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())

def indexColumns(infile, column, base="", header=True, extension=".idx", append=False, start=1, sep="\t"):
	index = start
	inlines = open(infile).readlines()
	outfile = open(infile + extension, "w")
	if header:
		header_line = inlines.pop(0)
		header_dict = general.build_header_dict(header_line, mode="line")
		column = header_dict[column]
		print >>outfile, header_line.strip()
	for inline in inlines:
		initems = inline.strip().split(sep)
		if append:
			outitems = initems[:column-1] + [base + str(index) + "." + initems[column]] + initems[column:]
		else:
			outitems = initems[:column-1] + [base + str(index)] + initems[column:]
		print >>outfile, sep.join(outitems)
		index += 1
	outfile.close()
				
def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "path from script to files")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "peaks to be used for analysis", default="OFF")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "operations to be performed")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--source", action = "store", type = "string", dest = "source", help = "Path of peaks to be filtered or other files...")
	parser.add_option("--infile", action = "store", type = "string", dest = "infile", help = "Input bed regions to explore")
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "Output file name")
	parser.add_option("--group", action = "store", type = "string", dest = "group", help = "Data-type or modENCODE group", default="regulation")
	parser.add_option("--header", action = "store", type = "string", dest = "header", help = "What files have headers? A, B, AB", default="OFF")
	parser.add_option("--target", action = "store", type = "string", dest = "target", help = "Type of label to generate", default="factor.context")
	parser.add_option("--overwrite", action = "store", type = "string", dest = "overwrite", help = "Overwrite files?", default="OFF")
	parser.add_option("--threshold", action = "store", type = "string", dest = "threshold", help = "Threshold for CNV-seq analysis...", default="OFF")
	parser.add_option("--rename", action = "store", type = "string", dest = "rename", help = "Targets to rename. Comma-separated list of 'target:replacement' pairs to search and replace.", default="OFF")
	parser.add_option("--include", action = "store", type = "string", dest = "include", help = "Targets to include", default="OFF")
	parser.add_option("--exclude", action = "store", type = "string", dest = "exclude", help = "Targets to exclude", default="OFF")
	parser.add_option("--window", action = "store", type = "string", dest = "window", help = "Window size surrounding input bed regions", default="OFF")
	parser.add_option("--up", action = "store", type = "string", dest = "up", help = "Upstream bases", default="OFF")
	parser.add_option("--dn", action = "store", type = "string", dest = "dn", help = "Downstream bases", default="OFF")
	parser.add_option("--spacer", action = "store", type = "int", dest = "spacer", help = "Spacer, window from which to exclude upstream and downstream signals...", default=0)
	parser.add_option("--strand", action = "store", type = "string", dest = "strand", help = "Enforce strandedness", default="ON")
	parser.add_option("--track", action = "store", type = "string", dest = "track", help = "Type of signal to use: 'fold' or 'logL'", default="fold")
	parser.add_option("--metric", action = "store", type = "string", dest = "metric", help = "Signal metric to use: 'signal.rng', 'signal.avg', or other...", default="signal.rng")
	parser.add_option("--multiarray", action = "store", type = "string", dest = "multiarray", help = "Store all values per position to allow calculating per position signal stats (mean, median, stdev)?", default="ON")
	parser.add_option("--format", action = "store", type = "string", dest = "format", help = "File format to use: 'bedGraph' or 'bigWig'", default="bedGraph")
	parser.add_option("--width", action = "store", type = "string", dest = "width", help = "Width of the bins to analyze")
	parser.add_option("--precision", action = "store", type = "int", dest = "precision", help = "Precision in rounding; decimal places...", default=3)
	parser.add_option("--cutChr", action = "store", type = "string", dest = "cutChr", help = "Should first 3 letters (chr) be removed?", default="OFF")
	parser.add_option("--restructure", action = "store", type = "string", dest = "restructure", help = "Restructure signal files to BED format?", default="OFF")
	parser.add_option("--nuclear", action = "store", type = "string", dest = "nuclear", help = "Peaks are only nuclear?", default="ON")
	parser.add_option("--threads", action = "store", type = "string", dest = "threads", help = "multiprocessing threads", default="1")
	parser.add_option("--chunks", action = "store", type = "int", dest = "chunks", help = "", default=1)
	parser.add_option("--qsub", action = "store", type = "string", dest = "qsub", help = "qsub configuration header", default="OFF")
	parser.add_option("--server", action = "store", type = "string", dest = "server", help = "are we on the server?", default="OFF")
	parser.add_option("--job", action = "store", type = "string", dest = "job", help = "job name for cluster", default="OFF")
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
	orthologspath = path_dict["orthologs"]
	peakspath = path_dict["peaks"]
	gopath = path_dict["go"]
	hotpath = path_dict["hot"]
	qsubpath = path_dict["qsub"]
	signalpath = path_dict["signal"]
	profilepath = path_dict["profile"]
	
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
		#organismIGV = "ce6"
	elif option.organism == "mm" or option.organism == "m.musculus":
		organismTag = "mm"
		#organismIGV = "ce6"
	elif option.organism == "ce" or option.organism == "c.elegans":
		organismTag = "ce"
		#organismIGV = "ce6"
	elif option.organism == "dm" or option.organism == "d.melanogaster":
		organismTag = "dm"
		#organismIGV = "dm5"
		
	# update signal path:
	signalpath = signalpath + option.organism + "/" + option.group + "/" + option.track + "/" + option.format + "/"
	
	# specify genome size file:
	if option.nuclear == "ON":
		chromosomes = metrn.chromosomes[organismTag]["nuclear"]
		genome_file = option.path + "/input/" + metrn.reference[organismTag]["genome"]
		genome_size_file = option.path + "/input/" + metrn.reference[organismTag]["nuclear_sizes"]
		genome_size_dict = general.build_config(genome_size_file, mode="single", separator="\t", spaceReplace=True, format="int")
	else:
		chromosomes = metrn.chromosomes[organismTag]["complete"]
		genome_file = option.path + "/input/" + metrn.reference[organismTag]["genome"]
		genome_size_file = option.path + "/input/" + metrn.reference[organismTag]["complete_sizes"]
		genome_size_dict = general.build_config(genome_size_file, mode="single", separator="\t", spaceReplace=True, format="int")
	
	# import bigWig data into bedGraph:
	if option.mode == "import":
		
		# create output folder:
		general.pathGenerator(signalpath)
		
		# load bigWig files:
		inpath = option.path + "/" + option.source + "/"
		infiles = os.listdir(inpath)
		
		print 
		print "Processing bigWig files:", len(infiles)
		for infile in infiles:
			outfile = "na_" + str(infile).replace("FE.bw", option.group + "_stn.fc.signal.bedgraph")
			if option.rename != "OFF":
				for scheme in option.rename.split(","):
					target, replacement = scheme.split(":")
					outfile = outfile.replace(target, replacement)
			command = "bigWigToBedGraph " + inpath + infile + " " + signalpath + outfile
			print command
			os.system(command)
		print
		
	
	# read parsing mode:
	if option.mode == "parser":
	
		# define output path:
		profilepath = profilepath + option.name + "/"
		depthspath = profilepath + "depths/"
		cnvseqpath = profilepath + "cnvseq/"
		parserpath = profilepath + "parser/"
		general.pathGenerator(depthspath)
		general.pathGenerator(cnvseqpath)
		general.pathGenerator(parserpath)
		
		# load (target) SAM files and FASTQ files:
		fastqpath = fastqpath.replace("meTRN", "ceTRN")
		samsamfiles = os.listdir(depthspath)
		fastqsfiles = os.listdir(fastqpath)
		
		print
		print "Scanning SAM and FASTQ files..."
		scanningDict = dict()
		for samsamfile in samsamfiles:
			if ".sam" in samsamfile:
				
				# define output file:
				parserfile = samsamfile.replace("mapprofile_depths_","mapprofile_parser_").replace(".sam",".txt")
			
				# check inclusion tags:
				include = False
				if option.include == "OFF":
					include = True
				else:
					for inclusion in option.include.split(","):
						if inclusion in samsamfile:
							include = True
				
				# check exclusion tags:
				exclude = True
				if option.exclude == "OFF":
					exclude = False
				else:
					exclude = False
					for exclusion in option.exclude.split(","):
						if exclusion in samsamfile:
							exclude = True
				
				# check for re-processing:
				process = True
				if option.overwrite == "ON" or not parserfile in os.listdir(parserpath):
					process = True
				
				# should we process this file?
				if (include and not exclude) and process:
					samsamitems = samsamfile.replace("_q30.sam","").split("_")
					fastqitems = samsamitems[10:]
					fastqitems[3] = "L" + fastqitems[3]
					fastqtarget = "_".join(fastqitems)
					for fastqsfile in fastqsfiles:
						if fastqtarget in fastqsfile:
							scanningDict[samsamfile] = fastqsfile
							#print "Match:", samsamfile, "(" + fastqsfile + ")"
							#break
		
		
		print "Loading target DNA sequence..."
		import fasta
		fastaDict = fasta.buildfile(path_dict[option.source] + option.infile)
		if len(fastaDict) == 1:
			fastaTarget = fastaDict[fastaDict.keys()[0]]
		genomeDict = fasta.buildfile(genome_file)
		
		matched, trackTarget = False, str(option.track)
		while not matched and len(trackTarget) > 8:
			for chrm in genomeDict:
				if trackTarget in genomeDict[chrm]:
					matched = True
					break
			trackTarget = trackTarget[1:]
		print "Search sequence:", trackTarget, "(" + option.track + ")"
		
		fastasfiles = list()
		for samsamfile in sorted(scanningDict.keys()):
			print
			print "Processing:", samsamfile, "(" + scanningDict[samsamfile] + ")"
			
			# define input and output files:
			parserfile = samsamfile.replace("mapprofile_depths_","mapprofile_parser_").replace(".sam",".txt")
			nonmapfile = samsamfile.replace("mapprofile_depths_","mapprofile_parser_").replace(".sam",".non")
			fastasfile = samsamfile.replace("mapprofile_depths_","mapprofile_parser_").replace(".sam",".non.fasta")
			fastqsfile = scanningDict[samsamfile]
			
			# collect FASTA files for de novo assembly:
			fastasfiles.append(parserpath + fastasfile)
			
			# extract unmapped reads:
			if option.overwrite == "ON" or not nonmapfile in os.listdir(parserpath):
				
				mappedIDs = list()
				indata = open(depthspath + samsamfile)
				inline = indata.readline()
				while inline:
					mappedIDs.append(inline.strip().split("\t")[0])
					inline = indata.readline()
				print "Mapped reads:", len(mappedIDs)
				
				inputIDs = dict()
				indata = open(fastqpath + fastqsfile)
				inline = indata.readline()
				while inline:
					if "@" == inline[0]:
						readID = inline.strip().lstrip("@").rstrip("/1")
						sequence = indata.readline().strip()
						inputIDs[readID] = sequence
					inline = indata.readline()
				nonmapIDs = set(inputIDs.keys()).difference(set(mappedIDs))
				print "Unmapped reads:", len(nonmapIDs)
				
				f_output = open(parserpath + nonmapfile, "w")
				for nonmapID in nonmapIDs:
					print >>f_output, "\t".join([nonmapID, inputIDs[nonmapID]])
				f_output.close()
			
			# generate FASTA file reads:
			if option.overwrite == "ON" or not fastasfile in os.listdir(parserpath):
				
				f_output = open(parserpath + fastasfile, "w")
				indata = open(parserpath + nonmapfile)
				inline = indata.readline()
				while inline:
					readID, sequence = inline.strip().split("\t")
					fasta.printout(sequence, key=readID, output_file=f_output, charged="ON")
					inline = indata.readline()
				f_output.close()
			
			# examine unmapped reads:
			if option.overwrite == "OFF" or not parserfile in os.listdir(parserpath):
				
				k = 0
				indata = open(parserpath + nonmapfile)
				inline = indata.readline()
				while inline:
					match = False
					readID, sequence = inline.strip().split("\t")
					if trackTarget in sequence:
						left, right = sequence.split(str(trackTarget))
						sequence = sequence.replace(trackTarget, trackTarget.lower())
						match = True
					elif gs.revcomp(trackTarget) in sequence:
						left, right = sequence.split(gs.revcomp(trackTarget))
						sequence = sequence.replace(gs.revcomp(trackTarget), gs.revcomp(trackTarget).lower())
						match = True
					if match:
						if len(left) > len(right):
							match = left
						else:
							match = right
						if len(match) >= 15:
							for chrm in genomeDict:
								if match in genomeDict[chrm]:
									print chrm, ">", readID, ">", sequence, "::", match, genomeDict[chrm].index(match), "(direct)"
								if gs.revcomp(match) in genomeDict[chrm]:
									print chrm, ">", readID, ">", sequence, "::", match, genomeDict[chrm].index(gs.revcomp(match)), "(revComp)"
								
						k += 1
					inline = indata.readline()
				print "Matched reads:", k	
		
		# define de-novo overlap and assembly files
		ovrlapfile = "mapprofile_parser_" + option.include + ".ovl"
		edenasfile = "mapprofile_parser_" + option.include + ".txt"
		
		print
		print "Launching de novo overlap analysis..."
		command = "edena -p PREFIX -r ".replace("PREFIX", parserpath + "mapprofile_parser_" + option.include) + " ".join(fastasfiles)
		#os.system(command)
		
		print "Launching de novo assembly analysis..."
		command = "edena -p PREFIX -e ".replace("PREFIX", parserpath + "mapprofile_parser_" + option.include) + parserpath + ovrlapfile
		#os.system(command)
		
		print
		
	
	# copy-number analysis mode:
	if option.mode == "cnvseq":
		
		# define output path:
		profilepath = profilepath + option.name + "/"
		depthspath = profilepath + "depths/"
		cnvseqpath = profilepath + "cnvseq/"
		general.pathGenerator(depthspath)
		general.pathGenerator(cnvseqpath)
		
		# define output file:
		cnvseqfile = cnvseqpath + "mapprofile_cnvseq_" + option.target.replace("maprofile_depths_","").replace(".hit",".txt")
		
		print
		print "Executing copy-number analysis..."
		print "Target:", option.target
		print "Reference:", option.source
		print
		os.chdir(cnvseqpath)
		command = "cnv-seq.pl --test TARGETFILE --ref REFERENCEFILE --log2-threshold THRESHOLD --genome-size GENOMESIZE".replace("TARGETFILE", depthspath + option.target).replace("REFERENCEFILE", depthspath + option.source).replace("THRESHOLD", option.threshold).replace("GENOMESIZE", str(sum(genome_size_dict.values())))
		os.system(command)
		print
	
	
	# generate read-depth stuff:
	if option.mode == "depths":
		
		# define output path:
		profilepath = profilepath + option.name + "/"
		depthspath = profilepath + "depths/"
		general.pathGenerator(depthspath)
		
		# tweak read alignment paths:
		indexpath = indexpath.replace("meTRN","ceTRN")
		qcfilterpath = qcfilterpath.replace("meTRN","ceTRN")
		qcfilterfiles = os.listdir(qcfilterpath)
		
		print
		print "Calculating read-depth per base..."
		for qcfilterfile in qcfilterfiles:
		
			# define input/output files:
			filtername = qcfilterfile.replace(".bedGraph", ".txt").replace(".bedgraph", ".txt").replace(".bed", ".txt").replace(".sam", "")
			samsamfile = depthspath + "mapprofile_depths_" + filtername + ".sam"
			indexsfile = depthspath + "mapprofile_depths_" + filtername + ".sai"
			bambamfile = depthspath + "mapprofile_depths_" + filtername + ".bam"
			sortedfile = depthspath + "mapprofile_depths_" + filtername + ".sort"
			depthsfile = depthspath + "mapprofile_depths_" + filtername + ".bedGraph"
			inreadfile = depthspath + "mapprofile_depths_" + filtername + ".bed"
			windowfile = depthspath + "mapprofile_depths_" + filtername + ".win.bed"
			outhitfile = depthspath + "mapprofile_depths_" + filtername + ".hit"
			
			# rename datasets if indicated:
			if option.rename != "OFF":
				for scheme in option.rename.split(","):
					target, replace = scheme.split(":")
					filtername = filtername.replace(target, replace)
					samsamfile = samsamfile.replace(target, replace)
					indexsfile = indexsfile.replace(target, replace)
					bambamfile = bambamfile.replace(target, replace)
					sortedfile = sortedfile.replace(target, replace)
					depthsfile = depthsfile.replace(target, replace)
					inreadfile = inreadfile.replace(target, replace)
					windowfile = windowfile.replace(target, replace)
					outhitfile = outhitfile.replace(target, replace)
			
			# check inclusion tags:
			include = False
			if option.include == "OFF":
				include = True
			else:
				for inclusion in option.include.split(","):
					if inclusion in samsamfile:
						include = True
			
			# check exclusion tags:
			exclude = True
			if option.exclude == "OFF":
				exclude = False
			else:
				exclude = False
				for exclusion in option.exclude.split(","):
					if exclusion in samsamfile:
						exclude = True
			
			# check for re-processing:
			process = True
			if option.overwrite == "ON" or not "mapprofile_depths_" + filtername + ".bedGraph" in os.listdir(depthspath):
				process = True
			
			# should we process this file?
			if (include and not exclude) and process:
				
				print "Processing:", filtername
				
				# copy the SAM file:
				command = "cp " + qcfilterpath + qcfilterfile + " " + samsamfile
				os.system(command)
				
				command = "cp " + indexpath + qcfilterfile.replace("_q30.sam",".sai") + " " + indexsfile
				os.system(command)
				
				# convert to bam format...
				command = "samtools view -bS -T GENOMEFILE -t CHROMINFOFILE INPUTSAMFILE > OUTPUTBAMFILE".replace("GENOMEFILE", genome_file).replace("CHROMINFOFILE", genome_size_file).replace("INPUTSAMFILE", samsamfile).replace("OUTPUTBAMFILE", bambamfile)
				os.system(command)
				
				# sort the bam file...
				command = "samtools sort OUTPUTBAMFILE SORTEDBAMFILE".replace("OUTPUTBAMFILE", bambamfile).replace("SORTEDBAMFILE", sortedfile)
				os.system(command)
				
				# calculate coverage...
				command = "genomeCoverageBed -bg -ibam SORTEDBAMFILE -g CHROMINFOFILE > OUTPUTBEDGRAPH".replace("CHROMINFOFILE", genome_size_file).replace("SORTEDBAMFILE", sortedfile + ".bam").replace("OUTPUTBEDGRAPH", depthsfile)
				os.system(command)
				
				# generate reads BED file...
				command = "bamToBed -i SORTEDBAMFILE > OUTBEDFILE".replace("SORTEDBAMFILE", sortedfile + ".bam").replace("OUTBEDFILE", inreadfile)
				os.system(command)
				
				# generate coverage in windows file...
				command = "coverageBed -a OUTBEDFILE -b INTERVALFILE > WINDOWSFILE".replace("OUTBEDFILE", inreadfile).replace("INTERVALFILE", annotationspath + option.window).replace("WINDOWSFILE", windowfile)
				os.system(command)
				
				# generate best-hits file (for CNV-seq analysis):
				command = str("samtools view -F 4 SORTEDBAMFILE | perl -lane 'print " + '"$F[2]\t$F[3]"' + "' > OUTHITFILE").replace("SORTEDBAMFILE", sortedfile + ".bam").replace("OUTHITFILE", outhitfile)
				os.system(command)
		
		print
	
	
	# scan peak density as a function of distance from features:
	if option.mode == "signal":
	
		print
		print "Importing genomica module..."
		import genomica
		
		# read signal files:
		signalfiles = os.listdir(signalpath)
		
		# define output path:
		profilepath = profilepath + option.name + "/"
		inputspath = profilepath + "inputs/"
		matrixpath = profilepath + "matrix/"
		vectorpath = profilepath + "vector/"
		viewerpath = profilepath + "viewer/"
		general.pathGenerator(inputspath)
		general.pathGenerator(matrixpath)
		general.pathGenerator(vectorpath)
		general.pathGenerator(viewerpath)
		
		# define domain path: 
		domainpath = path_dict[option.source]
		
		# define strand policy:
		if option.strand == "ON":
			strandFlag = True
		else:
			strandFlag = False
		
		# Report per-base statistics?
		# Note: This is where we are at! 
		# If we remove the stats incorporation, 
		# we won't be consuming so much memory!
		# See "multiArray" flag in 'queryArray'
		# function (genomica module).
		# 5/6/2013
		
		# process signal files:
		print "Processing input signal files..."
		print
		for signalfile in signalfiles:
			
			# define input/output files:
			signalname = signalfile.replace(".bedGraph", ".txt").replace(".bedgraph", ".txt").replace(".bed", ".txt")
			vectorname = "mapprofile_signal_" + signalname
			vectorfile = vectorpath + vectorname
			viewerfile = viewerpath + vectorname
			outputfile = matrixpath + "mapprofile_signal_" + signalfile
			domainfile = domainpath + option.infile
			
			# rename datasets if indicated:
			if option.rename != "OFF":
				for scheme in option.rename.split(","):
					target, replace = scheme.split(":")
					signalname = signalname.replace(target, replace)
					vectorname = vectorname.replace(target, replace)
					vectorfile = vectorfile.replace(target, replace)
					outputfile = outputfile.replace(target, replace)
			
			# define end-tag:
			if option.track == "fold":
				endtag = ".fc.signal.txt"
			
			# should we process this file?
			process = False
			if not vectorname in os.listdir(vectorpath) or option.overwrite == "ON":
				process = True
				if option.peaks != "OFF":
					process = False
					peakfiles = os.listdir(peakspath + option.peaks)
					signaltag = signalname.replace(endtag,"")
					for peakfile in peakfiles:
						if signaltag in peakfile:
							process = True
			
			# check inclusion tags:
			include = False
			if option.include == "OFF":
				include = True
			else:
				for inclusion in option.include.split(","):
					if inclusion in signalname:
						include = True
			
			# check exclusion tags:
			exclude = True
			if option.exclude == "OFF":
				exclude = False
			else:
				for exclusion in option.exclude.split(","):
					if exclusion in signalname:
						exclude = True
			
			# re-calculate processing:
			if process:
				if not include or exclude:
					process = False
			
			# should we restructure signal files or remove 'chr' from names? 
			if process and (option.restructure == "ON" or option.cutChr == "ON"):
				sourcefile = inputspath + signalfile
				if not signalfile in os.listdir(inputspath):
					print "Rebuilding signal file:", signalfile
					i_output = open(inputspath + signalfile, "w")
					indata = open(signalpath + signalfile)
					inline = indata.readline()
					feature = 1
					while inline:
					   if option.cutChr == "ON":
					   	output = inline.lstrip("chr").strip()
					   else:
					   	output = inline.strip()
					   if option.restructure == "OFF":
					   	print >>i_output, output
					   else:
					   	output = output.split("\t")
					   	score = output.pop(len(output)-1)
					   	output.append("signal-" + str(feature))
					   	output.append(score)
					   	print >>i_output, "\t".join(output)
					   inline = indata.readline()
					   feature += 1
					i_output.close()
			else:
				sourcefile = signalpath + signalfile
			
			# should we calculate profiles with HTSeq?
			if process and (not vectorname in os.listdir(vectorpath) or option.overwrite == "ON"):
				
				print "Profiling:", vectorname
				print "Started:", time.asctime(time.localtime())
				
				# load signal data:
				print "Loading signal data..."
				genome = genomica.loadArray(infile=sourcefile, sizes=genome_size_dict, index=3, exclude=["MtDNA"], stranded=False, typecode="d")
				
				# extract domain signal:
				print "Extracting domain signals..."
				profile, storage, regions, counter = genomica.queryArray(infile=domainfile, genome=genome, sizes=genome_size_dict, flanks=int(option.window), exclude=["MtDNA"], header=True, multiArray=option.multiarray, multiExport=viewerpath + vectorname, rounding=option.precision)
				
				# setup output file:
				print "Exporting average domain signals..."
				f_output = open(vectorfile, "w")
				print >>f_output, "\t".join(["index", "position", "signal.sum", "signal.max", "signal.avg", "signal.med", "signal.std", "signal.count"])
				k = 0
				profile = profile.tolist()
				for index in range(-int(option.window), int(option.window) + 1):
					output = [k + 1, index]
					
					# process storing per-position stats?
					if option.multiarray == "ON":
						# Note: Here we will use numpy array slicing as follows:
						# x = numpy.array([[1,2,3],[4,5,6],[7,8,9]])
						# x[:,1]
						# array([2, 5, 8])
						values = list(storage[:,k])
						output.append(sum(values))
						output.append(max(values))
						output.append(numpy.mean(values))
						output.append(numpy.median(values))
						output.append(numpy.std(values))
						output.append(len(values))
					
					# process storing per-position sums and mean only?
					else:
						output.append(profile[k])
						output.append(max(profile))
						output.append(float(profile[k])/counter)
						output.append(0)
						output.append(0)
						output.append(len(profile))
						
					print >>f_output, "\t".join(map(str, output))
					k += 1
					
				f_output.close()
				print "Ended:", time.asctime(time.localtime())
				print
				
				# clear items from memory:
				del genome
				del profile
				del storage
				del regions
	
		
	# group peak density as a function of distance from different datasets:
	if option.mode == "matrix":
	
		# define output path:
		profilepath = profilepath + option.name + "/"
		matrixpath = profilepath + "matrix/"
		vectorpath = profilepath + "vector/"
		inputspath = profilepath + "inputs/"
		general.pathGenerator(matrixpath)
		general.pathGenerator(vectorpath)
		general.pathGenerator(inputspath)
		
		# prepare output file:
		f_outfile = open(matrixpath + "mapprofile_matrix_" + option.name +  "_" + option.peaks + "_matrix.txt", "w")
		s_outfile = open(matrixpath + "mapprofile_matrix_" + option.name +  "_" + option.peaks + "_sorted.txt", "w")
		
		# process signal profiles:
		print
		print "Compiling signal per dataset..."
		k, m, outputDict, ratiosDict = 0, 0, dict(), dict()
		for profile in os.listdir(vectorpath):
			if "mapprofile_signal_" in profile:
				dataset = organismTag + "_" + profile.replace("mapprofile_signal_", "").replace(".fc.signal.txt", "")
				
				# rename datasets if indicated:
				if option.rename != "OFF":
					for scheme in option.rename.split(","):
						target, replace = scheme.split(":")
						dataset = dataset.replace(target, replace)
				
				# generate output name:
				datasetID = metrn.labelGenerator(target=option.target, mode="label", dataset=dataset)
				
				# check whether dataset is present in target peaks:
				if option.peaks != "OFF":
					process = False
					peaksfile = dataset + "_peaks.bed"
					if peaksfile in os.listdir(peakspath + option.peaks):
						process = True
						m += 1
				else:
					process = True
				
				# export per-position signal to output matrix:
				if process:
				
					# load lines:
					inlines = open(vectorpath + profile).readlines()
					headerDict = general.build_header_dict(vectorpath + profile)
					header = inlines.pop(0)
					
					# extract max (reference) signal:
					signalDict = general.build2(vectorpath + profile, i="position", x="signal.avg", mode="values")
					refSignal = max(map(float, signalDict.values()))
					minSignal = min(map(float, signalDict.values()))
					
					# initiate output dictionary:
					outputDict[datasetID] = dict()
					
					# export lines:
					N = 0
					Nmax = len(inlines)
					Nhlf = Nmax/2 + 1
					upstream, dnstream = list(), list()
					for inline in inlines:
					
						# export header..
						if k == 0:
							print >>f_outfile, "dataset\t" + header.strip() + "\tsignal.ref\tsignal.nrm\tsignal.rng"
							print >>s_outfile, "dataset\t" + header.strip() + "\tsignal.ref\tsignal.nrm\tsignal.rng\torder\tratio"
							k += 1
						
						# export signal...
						if inline.strip() != "":
							initems = inline.strip().split("\t")
							nrmSignal = float(initems[headerDict["signal.avg"]])/refSignal
							rngSignal = float(float(initems[headerDict["signal.avg"]])-minSignal)/(refSignal-minSignal)
							output = datasetID + "\t" + inline.strip() + "\t" + str(refSignal) + "\t" + str(nrmSignal) + "\t" + str(rngSignal)
							print >>f_outfile, output
							
							if N < Nhlf - option.spacer:
								upstream.append(rngSignal)
							elif N >= Nhlf + option.spacer:
								dnstream.append(rngSignal)
							
							ratiosDict[datasetID] = numpy.log2(float(numpy.mean(upstream))/numpy.mean(dnstream))
							outputDict[datasetID][N] = output
							N += 1
		
		print "Exporting signal per dataset..."
		order = 1
		datasetIDs = general.valuesort(ratiosDict)
		for datasetID in datasetIDs:
			for N in sorted(outputDict[datasetID].keys()):
				print >>s_outfile, outputDict[datasetID][N] + "\t" + str(order) + "." + str(N).zfill(5) + "\t" + str(ratiosDict[datasetID])
			order += 1
		print
				
		# close output file:
		f_outfile.close()
		s_outfile.close()
	
	
	"""
	# extract signal densities from input matrix:
	if option.mode == "orders":
		
		# load correlation and quantile functions:
		from scipy.stats.stats import pearsonr
		from quantile import Quantile
		
		# define output path:
		profilepath = profilepath + option.name + "/"
		viewerpath = profilepath + "viewer/"
		orderspath = profilepath + "orders/"
		general.pathGenerator(orderspath)
		
		print
		print "Loading signal data per experiment..."
		regionsReps = 0
		profileDict, averageDict = dict(), dict()
		indata = open(viewerpath + option.infile)
		inline = indata.readline()
		while inline:
			initems = inline.strip().split("\t")
			region = initems.pop(0)
			initems = map(float, initems)
			if region in profileDict:
				regionsReps += 1
				#print region
				#pdb.set_trace()
				
			#if sum(initems) > 10
			profileDict[region] = initems
			averageDict[region] = sum(initems)
			inline = indata.readline()
		
			# Note: This could be done in a chunky mode:
			#flanks = len(initems) - 1
			#midValue, upValues, dnValues = initems[flanks], initems[:flanks], initems[flanks+1:]
			#upChunks, dnChunks = list(), list()
			#for chunk in general.chunks(upValues, 10):
			#	upChunks.append(numpy.mean(chunk))
			#for chunk in general.chunks(dnValues, 10):
			#	dnChunks.append(numpy.mean(chunk))
			#profileDict[region] = upChunks + [midValue] + dnChunks
			
		print "Loaded region profiles:", len(profileDict)
		
		matrix = dict()
		print "Mean signal sum", numpy.mean(averageDict.values())
		print "10% Quantile:", Quantile(averageDict.values(), 0.10)
		print "25% Quantile:", Quantile(averageDict.values(), 0.25)
		print "50% Quantile:", Quantile(averageDict.values(), 0.50)
		print "75% Quantile:", Quantile(averageDict.values(), 0.75)
		threshold = Quantile(averageDict.values(), 0.50)
		
		profileTargets = list()
		profileRegions = sorted(profileDict.keys())
		for profileRegion in profileRegions:
			if averageDict[profileRegion] > threshold:
				profileTargets.append(profileRegion)
		print "Targeted regions:", len(profileTargets)
		
		print "Calculating distances..."
		for i in profileRegions:
			if averageDict[i] > threshold:
				if not i in matrix:
					matrix[i] = dict()
				for j in profileRegions:
					if averageDict[j] > threshold:
						correlation, corPvalue = pearsonr(profileDict[i], profileDict[j])
						matrix[i][j] = correlation
						#print i, j, correlation
		print
	"""
	
	# examine changes (delta) in peak density as a function of distance:
	if option.mode == "deltas":
	
		# define output path:
		profilepath = profilepath + option.name + "/"
		matrixpath = profilepath + "matrix/"
		vectorpath = profilepath + "signal/"
		inputspath = profilepath + "inputs/"
		deltaspath = profilepath + "deltas/"
		general.pathGenerator(matrixpath)
		general.pathGenerator(vectorpath)
		general.pathGenerator(inputspath)
		general.pathGenerator(deltaspath)
		
		# load input matrix output file:
		print
		print "Loading signal data per experiment..."
		matrixfile = matrixpath + "mapprofile_matrix_" + option.name +  "_" + option.peaks + "_matrix.txt"
		matrixDict = general.build2(matrixfile, i="dataset", j="position", x=option.metric, mode="matrix")
		for datasetID in matrixDict:
			for position in matrixDict[datasetID]:
				matrixDict[datasetID][position] = float(matrixDict[datasetID][position])
		
		# parse factor information per dataset:
		contextDict = dict()
		for datasetID in matrixDict:
			parts = datasetID.split(".")
			factor = ".".join(parts[:len(parts)-1])
			context = parts[len(parts)-1]
			if not factor in contextDict:
				contextDict[factor] = dict()
			contextDict[factor][context] = datasetID
		
		# examine changes per factor:
		comboDict, factorDict = dict(), dict()
		print "Parsing signal data per factor..."
		for factor in contextDict:
			if len(contextDict[factor]) > 1:
				comboDict[factor] = dict()
				factorDict[factor] = dict()
				contexts = contextDict[factor].keys()
				combiCount, deltaCount = 0, 0
				for combination in itertools.combinations(contexts, 2):
					contextX, contextY = sorted(list(combination))
					datasetX, datasetY = factor + "." + contextX, factor + "." + contextY
					positionsX, positionsY = general.valuesort(matrixDict[datasetX]), general.valuesort(matrixDict[datasetY])
					positionsX.reverse()
					positionsY.reverse()
					maximalX, maximalY = positionsX[0], positionsY[0]
					if (int(maximalX) < 0 and int(maximalY) < 0) or (int(maximalX) >= 0 and int(maximalY) >= 0):
						deltaStatus = 0
					else:
						deltaStatus = 1
						deltaCount += 1
					combiCount += 1
					if not contextX in comboDict[factor]:
						comboDict[factor][contextX] = dict()
					comboDict[factor][contextX][contextY] = [int(maximalX), int(maximalY), matrixDict[datasetX][maximalX], matrixDict[datasetY][maximalY], deltaStatus]
				factorDict[factor] = [combiCount, deltaCount, float(deltaCount)/combiCount]
		
		# export data:
		f_output = open(deltaspath + "mapprofile_delta_" + option.name +  "_" + option.peaks + "_factor.txt", "w")
		c_output = open(deltaspath + "mapprofile_delta_" + option.name +  "_" + option.peaks + "_combos.txt", "w")
		print >>f_output, "\t".join(["factor", "context.count", "combo.count", "delta.count", "delta.fraction"])
		print >>c_output, "\t".join(["factor", "i", "j", "i.pos", "j.pos", "i.max", "j.max", "status"])
		print "Exporting signal changes per factor..."
		for factor in comboDict:
			print >>f_output, "\t".join(map(str, [factor, len(contextDict[factor].keys())] + factorDict[factor]))
			for contextX in sorted(comboDict[factor].keys()):
				for contextY in sorted(comboDict[factor][contextX].keys()):
					print >>c_output, "\t".join(map(str, [factor, contextX, contextY] + comboDict[factor][contextX][contextY]))
					
		# close output files:
		f_output.close()
		c_output.close()
					
	
	# examine relative (ratios) signal density up- and down-stream of centers (TSSs):
	if option.mode == "ratios":
	
		# define output path:
		profilepath = profilepath + option.name + "/"
		matrixpath = profilepath + "matrix/"
		vectorpath = profilepath + "signal/"
		inputspath = profilepath + "inputs/"
		deltaspath = profilepath + "deltas/"
		ratiospath = profilepath + "ratios/"
		general.pathGenerator(matrixpath)
		general.pathGenerator(vectorpath)
		general.pathGenerator(inputspath)
		general.pathGenerator(deltaspath)
		general.pathGenerator(ratiospath)
		
		# load input matrix output file:
		print
		print "Loading signal data per experiment..."
		matrixfile = matrixpath + "mapprofile_matrix_" + option.name +  "_" + option.peaks + "_matrix.txt"
		matrixDict = general.build2(matrixfile, i="dataset", j="position", x=option.metric, mode="matrix")
		for datasetID in matrixDict:
			for position in matrixDict[datasetID]:
				matrixDict[datasetID][position] = float(matrixDict[datasetID][position])
		
		# generate output file:
		#ratiosfile = ratiospath + "mapprofile_ratios_" + option.name +  "_" + option.peaks + "_ratios.txt"
		#f_output = open(ratiosfile, "w")
		#print >>f_output, "\t".join(["dataset"])
		
		# parse factor information per dataset:
		contextDict = dict()
		for datasetID in matrixDict:
			parts = datasetID.split(".")
			factor = ".".join(parts[:len(parts)-1])
			context = parts[len(parts)-1]
			if not factor in contextDict:
				contextDict[factor] = dict()
			contextDict[factor][context] = datasetID
		
		# define output file:
		f_output = open(ratiospath + "mapprofile_ratios_" + option.name + "_" + option.peaks + "_ratios.txt", "w")
		print >>f_output, "\t".join(["dataset", "factor", "context", "upstream", "dnstream", "ratio"])
		
		print "Calculating upstream/downstream signal ratios..."
		upmost, dnmost = 0, 0
		ratiosDict, outputDict, factorDict = dict(), dict(), dict()
		for datasetID in matrixDict:
			parts = datasetID.split(".")
			factor = ".".join(parts[:len(parts)-1])
			context = parts[len(parts)-1]
			upstream, dnstream = list(), list()
			for position in matrixDict[datasetID]:
				if int(position) <= -option.spacer:
					upstream.append(matrixDict[datasetID][position])
				elif int(position) >= option.spacer:
					dnstream.append(matrixDict[datasetID][position])
			ratiosDict[datasetID] = numpy.log2(numpy.mean(upstream)/numpy.mean(dnstream))
			outputDict[datasetID] = [datasetID, factor, context, numpy.mean(upstream), numpy.mean(dnstream), ratiosDict[datasetID]]
			if ratiosDict[datasetID] > 0:
				upmost += 1
			elif ratiosDict[datasetID] < 0:
				dnmost += 1
		
		# export ratios in order:
		for datasetID in general.valuesort(ratiosDict):
			print >>f_output, "\t".join(map(str, outputDict[datasetID]))
		
		# get standard deviation of ratios:
		if option.threshold == "OFF" or option.threshold == "SD":
			ratiosThreshold = numpy.std(ratiosDict.values())
		else:
			ratiosThreshold = float(option.threshold)
		
		print "Examining factors with changes in signal ratios..."
		stable, change = 0, 0
		for factor in contextDict:
			if len(contextDict[factor]) > 1:
				upbias, dnbias = False, False
				for context in contextDict[factor]:
					datasetID = contextDict[factor][context]
					if ratiosDict[datasetID] > ratiosThreshold:
						upbias = True
					elif ratiosDict[datasetID] < -ratiosThreshold:
						dnbias = True
				if upbias and dnbias:
					change += 1
				else:
					stable += 1
		
		print "Factors with stable preferences:", stable, "(" + str(round(100*float(stable)/(stable+change), 2)) + "%)"
		print "Factors that change preferences:", change, "(" + str(round(100*float(change)/(stable+change), 2)) + "%)"
		print
		
		# close output file:
		f_output.close()
			
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())


#samtools view -H mapprofile_depths_OP102_HAM-1_L1_yale_stn_1_DNA_a_100716_ROCKFORD_FC624GH_3_GTAT_q30.sort.bam | grep @SQ | cut -f2,3 | sed 's/SN://g' | sed 's/LN://g' > genome.txt
#ce_wormbased_GEN_gx
#ce_windowing_500_gx