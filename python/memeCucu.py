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

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())

""" define classes and functions of internal use """

def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "Path from script to files")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "Peaks to be used for analysis")
	parser.add_option("--max", action = "store", type = "string", dest = "max", help = "Maximum number of peaks to consider")
	parser.add_option("--window", action = "store", type = "string", dest = "window", help = "Window surrounding peak within which the search will be performed")
	parser.add_option("--exclude", action = "store", type = "string", dest = "exclude", help = "Regions to remove!", default="OFF")
	parser.add_option("--repeatMask", action = "store", type = "string", dest = "repeatMask", help = "Should repeat-masked genome be used?", default="OFF")
	parser.add_option("--nuclear", action = "store", type = "string", dest = "nuclear", help = "Peaks are only nuclear?", default="ON")
	parser.add_option("--parameters", action = "store", type = "string", dest = "parameters", help = "MEME parameters")
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
	
	# load sequence dictionary:
	if option.repeatMask == "ON":
		# Not yet implemented!
		#genome_dict = fasta.buildfile(genome_fastafile, options=["take.first"], append="")
		rm_handle = "T"
	elif option.repeatMask == "OFF":
		#genome_dict = metrn.genomeBuilder(genome_file, chromosome_path, organism=organismTag, options=["take.first"], append="")
		rm_handle = "F"
	
	# set some of the basic parameters:
	max_peaks = int(option.max)
	window = int(option.window)
	
	# prepare output folders:
	coordpath = memepath + option.peaks + "/coord/"
	fastapath = memepath + option.peaks + "/fasta/"
	outputpath = memepath + option.peaks + "/output/"
	general.pathGenerator(coordpath)
	general.pathGenerator(fastapath)
	general.pathGenerator(outputpath)
		
	# get peak call files:
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
		if "_peaks.bed": #and ("HLH-1" in peakfile or "PHA-4" in peakfile):
		
			k += 1
			print k, peakfile
			dataset = peakfile.replace("_peaks.bed","")
			organism, strain, factor, context, institute, method = metrn.labelComponents(dataset)
			inputfile = peakspath + peakfile
			coordfile = coordpath + peakfile
			
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
			for peak in peaks[:int(option.max)]:
				print >>c_output, "\t".join(map(str, peak_dict[peak]))
			c_output.close()
			
			# create fasta file:
			fastafile = fastapath + "memecucu_w" + option.window + "_m" + option.max + "_rm" + rm_handle + "_" + dataset + ".fasta"
			command = "fastaFromBed -fi " + genome_file + " -bed " + coordfile + " -fo " + fastafile + " -name"
			os.system(command)
	
	# remove temporary background exclusion bed file:
	if option.exclude != "OFF":
		command = "rm -rf " + memepath + "*.tmp"
		os.system(command)			
				
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
	
	
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())


#python memeCucu.py --path ~/meTRN --organism hs --peaks hs_standards_reg_cx_hct --max 250 --window 100 --repeatMask OFF --parameters "-maxsize 10000000 -nmotifs 3 -minsites 50 -maxw 12"
#python memeCucu.py --path ~/meTRN --organism ce --peaks ce_selection_reg_cx_hct --max 250 --window 100 --repeatMask OFF --parameters "-maxsize 10000000 -nmotifs 3 -minsites 50 -maxw 12"
#python memeCucu.py --path ~/meTRN --organism dm --peaks dm_selection_reg_cx_hct --max 250 --window 100 --repeatMask OFF --parameters "-maxsize 10000000 -nmotifs 3 -minsites 50 -maxw 12"
