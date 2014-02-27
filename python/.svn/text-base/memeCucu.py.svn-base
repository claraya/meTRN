#!/usr/bin/env python
# derive motifs from transcription factor binding data

import sys
import time
import optparse
import general
import numpy
import fasta
import metrn
import modencode
import bed
import os
import copy
import pdb
import re
import network

from Bio import Motif

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())

""" define classes and functions of internal use """


""" define a function to build a fasta file from a bed file... """
def fastaGenerator(inputfile, coordfile, fastafile, genomefile, window=0, top="OFF"):

	# load peak dictionary:
	peak_dict, signal_dict = dict(), dict()
	inlines = open(inputfile).read().split("\n")
	for inline in inlines:
		if not inline == "":
			chrm, start, end, peak, score, strand, signal, pvalue, qvalue, point = inline.strip().split("\t")
			middle = int(start) + (int(end)-int(start))/2
			point = int(start) + int(point)
			if window != 0:
				start, end = point-window, point+window
			peak_dict[peak], signal_dict[peak] = [chrm, start, end, peak, score, strand, signal], float(signal)
	
	# generate target peaks file:
	peaks = general.valuesort(signal_dict)
	peaks.reverse()
	c_output = open(coordfile, "w")
	if top != "OFF":
		peaks = peaks[:int(top)]
	for peak in peaks:
		print >>c_output, "\t".join(map(str, peak_dict[peak]))
	c_output.close()
	
	# create fasta file:
	command = "fastaFromBed -fi " + genomefile + " -bed " + coordfile + " -fo " + fastafile + " -name"
	os.system(command)
	

""" define a function to build a BioPython Motif object from input motif lines... """
def motifGenerator(motiflines, format="ACGT"):
	motif = list()
	for motifline in motiflines:
		index = 0
		position = dict()
		for score in motifline.strip().replace("  "," ").replace("  "," ").split(" "):
			if score != " " and score != "":
				position[format[index]] = float(score)
				index += 1
		motif.append(position)
	return motif
	

""" define a function to export a TRANSFAC motif file from input motif lines... """
def transfacGenerator(motiflines, name, species="", format="ACGT", zmax=2, start=0, mode="standard"):
	transfaclines = list()
	transfaclines.append("ID " + name)
	transfaclines.append("BF " + species)
	transfaclines.append("P0\t" + "\t".join(format))
	position = 1
	consensus = ""
	
	# filter motif lines to grab just score lines...
	scorelines, record = list(), False
	if mode == "standard":
		for motifline in motiflines:
			if record:
				scorelines.append(motifline)
			if "letter-probability matrix" in motifline:
				record = True
	else:
		scorelines = motiflines
	
	# record scores into TRANSFAC format:
	for scoreline in scorelines:
		if general.clean(scoreline) != list() and not "URL" in scoreline:
			index = 0
			scores = list()
			conbase = "N"
			for score in scoreline.strip().replace("  "," ").replace("  "," ").split(" ")[start:]:
				if score != " " and score != "":
					scores.append(float(score))
					if float(score) > 0.5:
						conbase = format[index]
					index += 1
			transfaclines.append(str(position).zfill(zmax) + "\t" + "\t".join(map(str, scores)) + "\t" + conbase)
			consensus += conbase
			position += 1
	transfaclines.append("//")
	return transfaclines


def motifScanner(scanDict, orthologDict, cutoff=0.25, size=3, exclusions="OFF"):

	scannedFactors, matchedFactors, successFactors = 0, 0, 0
	successDict, matchDict = dict(), dict()
	for factor in scanDict:
		if factor.lower() in orthologDict:
			processFactor = True
			if exclusions != "OFF":
				for exclusion in exclusions.split(","):
					if exclusion in factor or exclusion.lower() in factor.lower():
						processFactor = False
			
			if processFactor:
				matchScores = list()
				for ortholog in orthologDict[factor.lower()]["hs"]:
					for motif in scanDict[factor]:
						if ortholog.lower()[:size] == motif.lower()[:size]:
							
							processMotif = True
							if exclusions != "OFF":
								for exclusion in exclusions.split(","):
									if exclusion in motif or exclusion.lower() in motif.lower():
										processMotif = False
							
							if processMotif:
								
								if not factor in matchDict:
									matchDict[factor] = dict()
								if not ortholog in matchDict[factor]:
									matchDict[factor][ortholog] = dict()
								matchDict[factor][ortholog][motif] = float(scanDict[factor][motif])
								matchScores.append(float(scanDict[factor][motif]))
								
								if float(scanDict[factor][motif]) > cutoff:
									if not factor in successDict:
										successDict[factor] = dict()
									if not ortholog in successDict[factor]:
										successDict[factor][ortholog] = dict()
									successDict[factor][ortholog][motif] = float(scanDict[factor][motif])
				
				# tally results:
				if factor in matchDict:
					matchedFactors += 1
					if factor in successDict:
						successFactors += 1
				scannedFactors += 1
	
	# return results:
	return successDict, matchDict, scannedFactors, matchedFactors, successFactors



def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "Path from script to files")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "Analysis mode: Determines operations to execute.")
	parser.add_option("--infile", action = "store", type = "string", dest = "infile", help = "Input file for analysis...", default="OFF")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "Peaks to be used for analysis...")
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "Motif database name...")
	parser.add_option("--max", action = "store", type = "string", dest = "max", help = "Maximum number of peaks to consider", default="OFF")
	parser.add_option("--window", action = "store", type = "string", dest = "window", help = "Window surrounding peak within which the search will be performed", default="OFF")
	parser.add_option("--include", action = "store", type = "string", dest = "include", help = "Targets/regions to include!", default="OFF")
	parser.add_option("--exclude", action = "store", type = "string", dest = "exclude", help = "Targets/regions to exclude!", default="OFF")
	parser.add_option("--repeatMask", action = "store", type = "string", dest = "repeatMask", help = "Should repeat-masked genome be used?", default="OFF")
	parser.add_option("--nuclear", action = "store", type = "string", dest = "nuclear", help = "Peaks are only nuclear?", default="ON")
	parser.add_option("--target", action = "store", type = "string", dest = "target", help = "Define analysis targets", default="OFF")
	parser.add_option("--parameters", action = "store", type = "string", dest = "parameters", help = "MEME parameters", default="OFF")
	parser.add_option("--background", action = "store", type = "string", dest = "background", help = "Base background frequency", default="OFF")
	parser.add_option("--sequence", action = "store", type = "string", dest = "sequence", help = "Sequence to search for...", default="OFF")
	parser.add_option("--threads", action = "store", type = "string", dest = "threads", help = "Multiprocessing threads", default="1")
	parser.add_option("--chunks", action = "store", type = "int", dest = "chunks", help = "", default=100)
	parser.add_option("--qsub", action = "store", type = "string", dest = "qsub", help = "qsub configuration header", default="OFF")
	parser.add_option("--server", action = "store", type = "string", dest = "server", help = "Are we on the server?", default="OFF")
	parser.add_option("--job", action = "store", type = "string", dest = "job", help = "Job name for cluster", default="OFF")
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
	if not option.mode in ["database", "measures"]:
		peakspath = peakspath + option.peaks + "/"
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
		
	# define organisms:
	organismTags = ["hs","mm","ce","dm"]
	
	# specify genome size file:
	if option.nuclear == "ON":
		chromosomes = metrn.chromosomes[organismTag]["nuclear"]
		genome_file = option.path + "/input/" + metrn.reference[organismTag]["genome"]
		genome_size_file = option.path + "/input/" + metrn.reference[organismTag]["nuclear_sizes"]
		genome_size_dict = general.build_config(genome_size_file, mode="single", separator="\t", spaceReplace=True)
		chromosome_path = option.path + "/fasta/" + metrn.reference[organismTag]["chromosome_path"]
		
	else:
		chromosomes = metrn.chromosomes[organismTag]["complete"]
		genome_file = option.path + "/input/" + metrn.reference[organismTag]["genome"]
		genome_size_file = option.path + "/input/" + metrn.reference[organismTag]["complete_sizes"]
		genome_size_dict = general.build_config(genome_size_file, mode="single", separator="\t", spaceReplace=True)
		chromosome_path = option.path + "/fasta/" + metrn.reference[organismTag]["chromosome_path"]
	
	# define window and max peaks flags:
	if not option.mode in ["database", "measures"]:
		windowFlag = "w" + option.window
		maxpeaksFlag = "m" + option.max
	
	# prepare output folders:
	if not option.mode in ["database", "measures"]:
		coordpath = memepath + option.peaks + "/" + windowFlag + "/" + maxpeaksFlag + "/coord/"
		fastapath = memepath + option.peaks + "/" + windowFlag + "/" + maxpeaksFlag + "/fasta/"
		outputpath = memepath + option.peaks + "/" + windowFlag + "/" + maxpeaksFlag + "/output/"
		general.pathGenerator(coordpath)
		general.pathGenerator(fastapath)
		general.pathGenerator(outputpath)
	
	# prepare database-dependent folders:
	if option.mode in ["graphing","pairwise","scanning","ortholog"]:
		graphingpath = memepath + option.peaks + "/" + windowFlag + "/" + maxpeaksFlag + "/graphing/" + option.name + "/"
		pairwisepath = memepath + option.peaks + "/" + windowFlag + "/" + maxpeaksFlag + "/pairwise/" + option.name + "/"
		scanningpath = memepath + option.peaks + "/" + windowFlag + "/" + maxpeaksFlag + "/scanning/" + option.name + "/"
		orthologpath = memepath + option.peaks + "/" + windowFlag + "/" + maxpeaksFlag + "/ortholog/" + option.name + "/"
		general.pathGenerator(graphingpath)
		general.pathGenerator(pairwisepath)
		general.pathGenerator(scanningpath)
		general.pathGenerator(orthologpath)
	
	# frequency measures mode:
	if option.mode == "measures":
		
		# define inputs and outputs:
		genomefile = genome_file
		coordfile = annotationspath + option.infile
		tempsfile = annotationspath + "memecucu_" + option.infile + "_sequence.coord"
		fastafile = annotationspath + "memecucu_" + option.infile + "_sequence.fasta"
		
		print
		if option.exclude != "OFF":
			print "Filtering:", option.exclude
			coordlines = open(coordfile).readlines()
			coordfile = str(tempsfile)
			c_output = open(tempsfile, "w")
			for coordline in coordlines:
				process = True
				coorditems = coordline.strip().split("\t")
				for exclusion in option.exclude.split(","):
					if exclusion in coorditems:
						process = False
						break
				if process:
					print >>c_output, "\t".join(coorditems)
			c_output.close()
			
		print "Building region sequences..."
		command = "fastaFromBed -fi " + genomefile + " -bed " + coordfile + " -fo " + fastafile + " -name"
		os.system(command)
		
		print "Building nucleotide frequencies..."
		frequencyDict = {"A":0, "C":0, "G":0, "T":0}
		fastaDict = fasta.buildfile(fastafile)
		for sequence in fastaDict:
			for nucleotide in sorted(frequencyDict.keys()):
				frequencyDict[nucleotide] += fastaDict[sequence].count(nucleotide)
		for nucleotide in sorted(frequencyDict.keys()):
			print nucleotide, float(frequencyDict[nucleotide])/sum(frequencyDict.values())
		print
		
		# remove temporary files:
		command = "rm -rf " + fastafile
		os.system(command)
		if option.exclude != "OFF":
			command = "rm -rf " + tempsfile
			os.system(command)
		
	
	# sequence generation mode:
	elif option.mode == "sequence":
		
		# specify repeat mask flag (not implemented):
		if option.repeatMask == "ON":
			rm_handle = "T"
		elif option.repeatMask == "OFF":
			rm_handle = "F"
		
		# define max peaks to evaluate and window around peak:
		max_peaks = option.max
		window = 0
		if option.max != "OFF":
			max_peaks = int(option.max)
		if option.window != "OFF":
			window = int(option.window)
		
		# load peak files:
		peakfiles = os.listdir(peakspath)
		
		# prepare a temporary background exclusion bed file:
		if option.exclude != "OFF":
			
			m_infile = memepath + option.exclude
			m_tmpfile = memepath + option.exclude + ".tmp"
			command = 'grep -v "feature" ' + m_infile + ' > ' + m_tmpfile
			os.system(command)
		
		# extract sequences from peak calls:
		master_dict, factors, k = dict(), list(), 0
		print
		print "Preparing target sequence files (FASTA):"
		for peakfile in peakfiles:
			if "_peaks.bed":
			
				k += 1
				print k, peakfile
				dataset = peakfile.replace("_peaks.bed","")
				organism, strain, factor, context, institute, method = metrn.labelComponents(dataset)
				inputfile = peakspath + peakfile
				coordfile = coordpath + peakfile
				fastafile = fastapath + "memecucu_w" + option.window + "_m" + option.max + "_rm" + rm_handle + "_" + dataset + ".fasta"
				
				# exclude genomic regions if necessary: 
				#if option.exclude != "OFF":
				#	
				#	p_infile = temppath
				#	p_outfile = macspath + "filtered/" + tempfile.replace(".bed","") + "_" + option.exclude.replace(".bed","") + "_filtered.bed"
				#	
				#	command = "intersectBed -v -a " + p_infile + " -b " + m_tmpfile + " > " + p_outfile
				#	os.system(command)
				#	
				#	tempfile = tempfile.replace(".bed","") + "_" + option.exclude.replace(".bed","") + "_filtered.bed"
				#	temppath = macspath + "filtered/" + tempfile
				
				# generate FASTA file:
				fastaGenerator(inputfile, coordfile, fastafile, genomefile=genome_file, window=window, top=option.max)
		print
	
	
	# motif discovery mode:
	elif option.mode == "discover":
		
		# launch motif analysis:
		print
		print "Executing MEME analysis:"
		fastafiles = os.listdir(fastapath)
		for fastafile in fastafiles:
			if "memecucu" in fastafile and ".fasta" in fastafile:
				command = "meme " + fastapath + fastafile + " -oc " + outputpath + fastafile.replace(".fasta","") + " " + option.parameters 
				print command
				os.system(command)
				print
	
	
	# directed search mode:
	elif option.mode == "fraction":
		
		# launch motif analysis:
		print
		print "Searching for motif matches:"
		searchDict = dict() 
		fastafiles = os.listdir(fastapath)
		for fastafile in fastafiles:
			if "memecucu" in fastafile and ".fasta" in fastafile:
				k, i = 0, 0
				sequenceDict = fasta.buildfile(fastapath + fastafile, options=["take.first"], append="")
				for feature in sequenceDict:
					if option.sequence in sequenceDict[feature]:
						i += 1
					k += 1
				searchDict[fastafile] = float(i)/k
		
		
		searchHits = general.valuesort(searchDict)
		searchHits.reverse()
		for searchHit in searchHits[:10]:
			print searchHit, round(searchDict[searchHit], 2)
		print

	
	# motif database mode:
	elif option.mode == "database":
	
		# output files (standard, TRANSFAC, MEME):
		b_infile = extraspath + option.background
		s_utfile = extraspath + option.infile.replace(".txt",".stn")
		t_utfile = extraspath + option.infile.replace(".txt",".dat")
		m_utfile = extraspath + option.infile.replace(".txt",".meme")
		s_output = open(s_utfile, "w")
		t_output = open(t_utfile, "w")
		
		# define header:
		if option.parameters != "OFF":
			header = True
			for line in option.parameters.split(";"):
				print >>s_output, line
		
		# export TRANSFAC headers:
		transfacLines = ["VV memeCucu:database transfac file", "XX"]
		print >>t_output, "\n".join(transfacLines)		
		
		print
		print "Loading input motifs..."
		k = 0
		motifs = open(extraspath + option.infile).read().split(">")
		for motif in motifs:
			lines = motif.replace("\r","\n").split("\n")
			lines = general.clean(lines)
			if len(lines) > 5:
				title = lines.pop(0)
				names = "MOTIF XXX"
				names = names.replace("XXX", title)
				infos = "letter-probability matrix: alength= ALENGTH w= WWW nsites= NSITES E= EVALUE "
				infos = infos.replace("ALENGTH", "4")
				infos = infos.replace("WWW", str(len(lines)))
				infos = infos.replace("NSITES", str(option.max))
				infos = infos.replace("EVALUE", "0")
				
				# export standard format file:
				print >>s_output, names
				print >>s_output, infos
				motiflines = list()
				for line in lines:
					motifline = " " + "  ".join(line.split(" ")[1:]) + " "
					motiflines.append(motifline)
					print >>s_output, motifline
				print >>s_output, ""
				
				# build TRANSFAC format:
				transfacLines = transfacGenerator(motiflines, name=title, mode="other")
				print >>t_output, "\n".join(transfacLines)
				k += 1
		
		print "Processed", k, "motifs."
		
		# close output file:
		s_output.close()
		t_output.close()
		
		# convert TRANSFAC file to MEME format:
		if option.background == "OFF":
			command = "transfac2meme -logodds" + t_utfile + " > " + m_utfile
		else:
			command = "transfac2meme -logodds -bg " + b_infile + " " + t_utfile + " > " + m_utfile
		os.system(command)
		print
	
	
	# motif graphing mode:
	elif option.mode == "graphing":
		
		print
		print "Graphing motif database..."
		titles = list()
		motifs = open(extraspath + option.infile).read().split("MOTIF")
		header = motifs.pop(0)
		for motif in motifs:
			inlines = motif.split("\n")
			motif = inlines.pop(0)
			title = motif.strip().split(" ")[0]
			if title in titles:
				print "Error:", title, "already processed!"
				pdb.set_trace()
			
			# select motifs of interest to search for matches:
			if option.parameters in motif or option.parameters == "OFF":
				titles.append(title)
				
				# build TRANSFAC format:
				transfacFile = graphingpath + title.replace("/","-") + ".dat"
				transfacLines = ["VV memeCucu:" + title, "XX"]
				transfacLines.extend(transfacGenerator(inlines, name=motif))
				t_output = open(transfacFile, "w")
				print >>t_output, "\n".join(transfacLines)
				t_output.close()
				
				# draw motif logo:
				graphingFile = graphingpath + title.replace("/","-") + ".eps"
				command = "weblogo -f TRANSFACFILE -D transfac -o GRAPHINGFILE --errorbars NO --size medium --color-scheme classic --units bits --sequence-type dna --composition none --aspect-ratio 4 --show-xaxis NO --show-yaxis NO --fineprint ''"
				os.system(command.replace("TRANSFACFILE", transfacFile).replace("GRAPHINGFILE",graphingFile))
				
				reversedFile = graphingpath + title.replace("/","-") + "-revComp.eps"
				command = "weblogo -f TRANSFACFILE -D transfac -o GRAPHINGFILE --errorbars NO --size medium --color-scheme classic --units bits --sequence-type dna --composition none --aspect-ratio 4 --show-xaxis NO --show-yaxis NO --fineprint '' --reverse --complement"
				os.system(command.replace("TRANSFACFILE", transfacFile).replace("GRAPHINGFILE",reversedFile))
				
				# remove TRANSFAC file:
				command = "rm -rf " + transfacFile
				os.system(command)
				
		print "Motifs graphed:", len(titles)
		print
		
	
	# motif comparison mode:
	elif option.mode == "pairwise":

		print
		print "Loading motif database..."
		titles = list()
		motifs = open(extraspath + option.infile).read().split("MOTIF")
		header = motifs.pop(0)
		for motif in motifs:
			inlines = motif.split("\n")
			motif = inlines.pop(0)
			title = motif.strip().split(" ")[0]
			if title in titles:
				print "Error:", title, "already processed!"
				pdb.set_trace()
			
			# select motifs of interest to search for matches:
			if option.parameters in motif or option.parameters == "OFF":
				
				# initialize output path for target motif:
				targetpath = pairwisepath + title + "/"
				general.pathGenerator(targetpath)
				
				# build motif (meme) file for target motif:
				t_utfile = targetpath + "motif.meme"
				t_output = open(t_utfile, "w")
				print >>t_output, header
				print >>t_output, "MOTIF " + title
				for inline in inlines:
					print >>t_output, inline.strip()
				t_output.close()
				
				# execute motif comparison:
				print "Processing:", title
				command = "tomtom -o OUTPUT TARGETMEME DATABASEMEME"
				command = command.replace("OUTPUT", targetpath + "tomtom")
				command = command.replace("TARGETMEME", t_utfile)
				command = command.replace("DATABASEMEME", extraspath + option.infile)
				os.system(command)
				print
		
		print
	
	
	# motif scanning mode:
	elif option.mode == "scanning":
	
		# specify repeat mask flag (not implemented):
		if option.repeatMask == "ON":
			rm_handle = "T"
		elif option.repeatMask == "OFF":
			rm_handle = "F"
		
		# define max peaks to evaluate and window around peak:
		max_peaks = option.max
		window = 0
		if option.max != "OFF":
			max_peaks = int(option.max)
		if option.window != "OFF":
			window = int(option.window)
		
		# define target flag:
		if option.target == "OFF":
			targetFlag = "factor"
		elif option.target == "ALL":
			targetFlag = "motifs"
		else:
			targetFlag = option.target
		
		# load peak files:
		peakfiles = os.listdir(peakspath)
		
		# initate scanning report file:
		f_utfile = scanningpath + "memecucu_" + option.mode + "_" + targetFlag + "_w" + option.window + "_m" + option.max + "_rm" + rm_handle + ".txt"
		s_utfile = scanningpath + "memecucu_" + option.mode + "_" + targetFlag + "_w" + option.window + "_m" + option.max + "_rm" + rm_handle + ".sum"
		f_output = open(f_utfile, "w")
		s_output = open(s_utfile, "w")
		print >>f_output, "\t".join(["factor","motif","fraction"])
		print >>s_output, "\t".join(["factor","motif","fraction"])
		
		print
		print "Loading motif database..."
		results = dict()
		titles = list()
		motifs = open(extraspath + option.infile).read().split("MOTIF")
		header = motifs.pop(0)
		for motif in motifs:
			inlines = motif.split("\n")
			motif = inlines.pop(0)
			title = motif.strip().split(" ")[0]
			if title in titles:
				print "Error:", title, "already processed!"
				pdb.set_trace()
			
			# select motifs of interest to search for matches:
			if (option.parameters in motif and option.target != "OFF") or (option.parameters in motif and option.target =="OFF" and "disc" in title) or (option.parameters == "OFF"):
				
				print "Processing:", title
				
				# initialize output path for target motif:
				targetpath = scanningpath + title + "/"
				general.pathGenerator(targetpath)
				
				# build TRANSFAC format:
				transfacFile = targetpath + "motif.dat"
				transfacLines = ["VV memeCucu:scanning transfac file", "XX"]
				transfacLines.extend(transfacGenerator(inlines, name=motif))
				t_output = open(transfacFile, "w")
				print >>t_output, "\n".join(transfacLines)
				t_output.close()
			
				# convert TRANSFAC file to MEME format:
				memeFile = targetpath + "motif.meme"
				if option.background == "OFF":
					command = "transfac2meme -logodds " + transfacFile + " > " + memeFile
				else:
					command = "transfac2meme -logodds -bg " + extraspath + option.background + " " + transfacFile + " > " + memeFile
				os.system(command)
				
				# find peak files for factor:
				matchHit, matchDict, matchFiles = "OFF", dict(), list()
				for peakfile in peakfiles:
					matchFlag = False
					organism, strain, factor, context, institute, method = metrn.labelComponents(peakfile)
					if option.target == "OFF" and factor.lower() in motif.replace(" ","_").split("_"):
						matchFlag = True
					elif option.target != "OFF" and factor in option.target.split(","):
						matchFlag = True
					elif option.target == "ALL":
						matchFlag = True
						
					if matchFlag:
						matchFiles.append(peakfile)
						matchHit = str(factor)
						if not matchHit in matchDict:
							matchDict[matchHit] = list()
						matchDict[matchHit].append(peakfile)
				
				# check that matches have been found:
				if matchHit == "OFF":
					print "No matches found:", motif
					print
					pdb.set_trace()
				
				# process matches:
				if matchDict != dict() and (matchHit.lower() in title or option.target != "OFF"):
					for matchHit in matchDict:
						if not matchHit in results:
							results[matchHit] = dict()
						fractions = list()
						for matchFile in matchDict[matchHit]:
							
							print "Scanning:", matchHit, "(" + matchFile + ")"
							
							# specify inputs and outputs:
							dataset = matchFile.replace("_peaks.bed", "")
							inputfile = peakspath + matchFile
							coordfile = coordpath + matchFile
							fastafile = fastapath + "memecucu_w" + option.window + "_m" + option.max + "_rm" + rm_handle + "_" + dataset + ".fasta"
							
							# generate FASTA file:
							fastaGenerator(inputfile, coordfile, fastafile, genomefile=genome_file, window=window, top=option.max)
							
							# determine number of input sequences:
							sequenceCount = len(open(coordfile).readlines())
							threshold = str(sequenceCount/10)
							
							# detect motif in sequences:
							if option.background == "OFF":
								command = "mast -ev " + threshold + " " + memeFile + " " + fastafile
							else:
								command = "mast -bfile " + extraspath + option.background + " -ev " + threshold + " " + memeFile + " " + fastafile
							os.system(command)
							
							# collect motif matches:
							mastlines = open("mast_out/mast.txt").read().split("SECTION I: HIGH-SCORING SEQUENCES")[1].split("SECTION II")[0].split("-------- ------")[1].split("**************")[0]
							mastlines = "".join(mastlines).strip().split("\n")
							matches, fraction = len(mastlines), round(float(len(mastlines))/sequenceCount, 2)
							fractions.append(fraction)
							print "Matches:", matches, "(" + str(100*fraction) + "%)"
							print
							
							# remove mast output:
							command = "rm -rf mast_out"
							os.system(command)
						
						# export summary file per motif:
						if fractions != list():
							results[matchHit][title] = max(fractions)
							print title, matchHit, str(100*max(fractions)) + "%"
							print >>f_output, "\t".join(map(str, [matchHit, title, max(fractions)]))
							print
		
		# export summary file per factor:
		for factor in sorted(results.keys()):
			motifs = general.valuesort(results[factor])
			motifs.reverse()
			print >>s_output, "\t".join(map(str, [factor, motifs[0], results[factor][motifs[0]]]))
		
		# close output files:
		f_output.close()
		s_output.close()
		print

	
	# motif orthology mode:
	elif option.mode == "ortholog":
	
		# specify repeat mask flag (not implemented):
		if option.repeatMask == "ON":
			rm_handle = "T"
		elif option.repeatMask == "OFF":
			rm_handle = "F"
		
		# define max peaks to evaluate and window around peak:
		max_peaks = option.max
		window = 0
		if option.max != "OFF":
			max_peaks = int(option.max)
		if option.window != "OFF":
			window = int(option.window)
		
		# define target flag:
		if option.target == "OFF":
			targetFlag = "factor"
		elif option.target == "ALL":
			targetFlag = "motifs"
		else:
			targetFlag = option.target
		
		print
		print "Loading motif scan results..."
		scanfile = scanningpath + "memecucu_scanning_" + targetFlag + "_w" + option.window + "_m" + option.max + "_rm" + rm_handle + ".txt"
		bestfile = scanningpath + "memecucu_scanning_" + targetFlag + "_w" + option.window + "_m" + option.max + "_rm" + rm_handle + ".sum"
		scanDict = general.build2(scanfile, i="factor", j="motif", x="fraction", mode="matrix")
		bestDict = general.build2(bestfile, i="factor", j="motif", x="fraction", mode="matrix")
		
		# define output files:
		f_output = open(orthologpath + "memecucu_ortholog_" + targetFlag + "_w" + option.window + "_m" + option.max + "_rm" + rm_handle + ".txt", "w")
		c_output = open(orthologpath + "memecucu_ortholog_" + targetFlag + "_w" + option.window + "_m" + option.max + "_rm" + rm_handle + ".cut", "w")		
		print >>f_output, "\t".join(["factor","ortholog","motif","fraction"])
		print >>c_output, "\t".join(["fraction.cutoff","fraction.factors"])
		
		print "Loading orthology tree..."
		orthologDict = metrn.orthologMapper(extraspath + option.infile, species=option.organism, targets=option.parameters.split(","))
		
		print "Scanning orthologous factors..."
		for cutoff in general.drange(0.2, 0.9, 0.05):
			successDict, matchDict, scannedFactors, matchedFactors, successFactors = motifScanner(scanDict, orthologDict, cutoff=cutoff, size=3, exclusions=option.exclude)
			print >>c_output, "\t".join(map(str, [cutoff, round(float(successFactors)/matchedFactors, 3)]))
		
		# find highest overlap ortholog and motif (per factor):
		rankDict, infoDict = dict(), dict()
		for factor in matchDict:
			highScore, highMotif, highOrtholog = 0, False, False
			for ortholog in matchDict[factor]:
				for motif in matchDict[factor][ortholog]:
					if matchDict[factor][ortholog][motif] > highScore:
						highScore = float(matchDict[factor][ortholog][motif])
						highOrtholog = str(ortholog)
						highMotif = str(motif)
			
			if highOrtholog:
				rankDict[factor] = highScore
				infoDict[factor] = [highOrtholog, highMotif, highScore]
		
		# export highest overlap ortholog and motif (per factor):
		rankFactors = general.valuesort(rankDict)
		rankFactors.reverse()
		for factor in rankFactors:
			print >>f_output, "\t".join(map(str, [factor] + infoDict[factor]))
		
		# close output files:
		f_output.close()
		c_output.close()
		print
		
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

#weblogo -f HMBOX1_DBD.dat -D transfac -o test.eps --errorbars NO --size medium --color-scheme classic