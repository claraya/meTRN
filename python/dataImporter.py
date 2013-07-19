#!/usr/bin/env python
# Import peak data into the meTRN structures!

import sys
import time
import optparse
import general
import numpy
import pickle
import pdb
import metrn
import modencode
import copy
import os

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())

""" define functions of internal use """

def uppify(indict):
	output = dict()
	for key in indict:
		newkey = key.upper()
		output[newkey] = indict[key]
	return output

def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "Path from script to files")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "Type of operations to perform")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "Basename for target peaks", default="OFF")
	parser.add_option("--source", action = "store", type = "string", dest = "source", help = "Source path or file", default="OFF")
	parser.add_option("--infile", action = "store", type = "string", dest = "infile", help = "Input file for analysis", default="OFF")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--rank", action = "store", type = "string", dest = "rank", help = "Name peaks by dataset and rank?", default="ON")
	parser.add_option("--label", action = "store", type = "string", dest = "label", help = "Type of labels to generate", default="factor.context")
	parser.add_option("--method", action = "store", type = "string", dest = "method", help = "Keep method in dataset labels?", default="ON")
	parser.add_option("--target", action = "store", type = "string", dest = "target", help = "Targets...", default="OFF")
	parser.add_option("--include", action = "store", type = "string", dest = "include", help = "Targets to include", default="OFF")
	parser.add_option("--exclude", action = "store", type = "string", dest = "exclude", help = "Targets to exclude", default="OFF")
	parser.add_option("--rename", action = "store", type = "string", dest = "rename", help = "Targets to rename. Comma-separated list of 'target:replacement' pairs to search and replace.", default="OFF")
	parser.add_option("--nonredundant", action = "store", type = "string", dest = "nonredundant", help = "Filter redundants?", default="OFF")
	parser.add_option("--parameters", action = "store", type = "string", dest = "parameters", help = "Variable parameters...", default="OFF")
	parser.add_option("--server", action = "store", type = "string", dest = "server", help = "Server variables...", default="OFF")
	parser.add_option("--fixed", action = "store", type = "string", dest = "fixed", help = "Should a fixed files be used?", default="OFF")
	parser.add_option("--cutChr", action = "store", type = "string", dest = "cutChr", help = "Should first 3 letters (chr) be removed?", default="OFF")
	parser.add_option("--filterChr", action = "store", type = "string", dest = "filterChr", help = "Remove peaks in these chromosomes", default="OFF")
	parser.add_option("--idrSource", action = "store", type = "string", dest = "idrSource", help = "Take peaks from idr/final path?", default="ON")
	parser.add_option("--reformat", action = "store", type = "string", dest = "reformat", help = "Should peaks be re-formatted?", default="OFF")
	(option, args) = parser.parse_args()
	
	# import paths:
	path_dict = modencode.configBuild(option.path + "/input/" + "configure_path.txt")
	
	# specify input and output paths:
	inpath = path_dict["input"]
	extraspath = path_dict["extras"]
	pythonpath = path_dict["python"]
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
	bindingpath = path_dict["binding"]
	networkpath = path_dict["network"]
	peakspath = path_dict["peaks"]
	gopath = path_dict["go"]
	hotpath = path_dict["hot"]
	cellspath = path_dict["cells"]
	neuronspath = path_dict["neurons"]
	
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
	elif option.organism == "mm" or option.organism == "m.musculus":
		organismTag = "mm"
	elif option.organism == "ce" or option.organism == "c.elegans":
		organismTag = "ce"
	elif option.organism == "dm" or option.organism == "d.melanogaster":
		organismTag = "dm"

	# import peaks mode:
	elif option.mode == "download.peaks":
		print
		
		# define destination:
		sourcepath = idrpath + "final/" + option.source + "/"
		general.pathGenerator(sourcepath)
		
		# launch command:
		print "Downloading peaks..."
		command = "scp " + option.server + ":" + metrn.sharedSCG + option.parameters + " " + sourcepath
		os.system(command)
		print
	
	# import peaks mode:
	elif option.mode == "import.peaks":
		print
		
		# define input path:
		if option.idrSource == "ON":
			sourcepath = idrpath + "final/" + option.source + "/"
		else:
			sourcepath = option.path + "/" + option.source + "/"
			sourcepath = sourcepath.replace("//", "/")
		
		# define output path:
		peakspath = idrpath + "peaks/" + option.peaks + "/"
		general.pathGenerator(peakspath)
		
		# gather peak files and transfer to the IDR peaks folder:
		print "Gathering peaks into peak folder..."
		infiles = os.listdir(sourcepath)
		for infile in infiles:
			
			# generate outfile name:
			outfile = copy.deepcopy(infile)
			
			# rename elements if necessary:
			if option.rename != "OFF":
				for scheme in option.rename.split(","):
					target, replace = scheme.split(":")
					outfile = outfile.replace(target, replace)
				
			dataset = outfile.replace("_peaks.bed", "")
			organism, strain, factor, context, institute, method = metrn.labelComponents(dataset)
			print "\t".join([organism, strain, factor, context, institute, method])
			#print infile
			
			# rank-name peaks?
			if option.rank == "ON":
				
				i = 1
				f_output = open(peakspath + outfile,"w")
				inlines = open(sourcepath + infile).readlines()
				for inline in inlines:
					initems = inline.strip().split("\t")
					
					if option.fixed == "ON" and len(initems) == 3:
						chrm, start, stop = inline.strip().split("\t")
						peak, score, strand, signal, pvalue, qvalue, point = "Peak_" + str(i), "0", ".", "0", "0", "0", "0"
					elif option.reformat == "ON" and len(initems) == 5:
						chrm, start, stop, name, score = inline.strip().split("\t")
						peak, score, strand, signal, pvalue, qvalue, point = "Peak_" + str(i), score, ".", "0", "0", "0", "0"
					else:
						chrm, start, stop, peak, score, strand, signal, pvalue, qvalue, point = inline.strip().split("\t")
						
					if option.method == "ON":
						peak = dataset + ":P" + str(i)
					else:
						peak = "_".join([organism, strain, factor, context, institute]) + ":P" + str(i)
					
					if option.cutChr == "ON":
						chrm = chrm.lstrip("chr")
					
					if option.filterChr == "OFF" or not chrm in option.filterChr.split(","):
						print >>f_output, "\t".join([chrm, start, stop, peak, score, strand, signal, pvalue, qvalue, point])
						i += 1
						
				f_output.close()
			
			# simply copy file...
			else:
				print "Transferring:", infile
				command = "cp " + sourcepath + infile + " " + peakspath + outfile
				os.system(command)
			
		print			
		
	
	# select peaks mode:
	if option.mode == "select.peaks":
		print
		
		# define input path and load input file names:
		sourcepath = idrpath + "peaks/" + option.source + "/"
		infiles = os.listdir(sourcepath)
		
		# define output path:
		peakspath = idrpath + "peaks/" + option.peaks + "/"
		general.pathGenerator(peakspath)
		
		# select dataset files that match the desired criteria:
		i, j = 0, 0
		print "Selecting compliant datasets..."
		selections = list()
		for infile in infiles:
			inclusions, exclusions = list(), list()
			if option.include != "OFF":
				for target in option.include.split(","):
					if target in infile:
						inclusions.append(target)
				#for target in option.include.split(";"):
				#	if target in infile:
				#		inclusions.append(target)
			if option.exclude != "OFF":
				for exclusion in option.exclude.split(","):
					if exclusion in infile:
						exclusions.append(exclusion)
			if inclusions != list() or len(inclusions) == len(option.include.split(";")) or option.include == "OFF":
				if exclusions == list() or option.exclude == "OFF":
					selections.append(infile)
					i += 1
			
		# remove redundant datasets if necessary:
		m, n = 0, 0
		if option.nonredundant == "ON":
			print "Filtering redundant datasets..."
			dataset_dict = dict()
			for infile in selections:
				label = metrn.labelGenerator(option.label, mode="label", dataset=infile)
				if not label in dataset_dict:
					dataset_dict[label] = dict()
				dataset_dict[label][infile] = general.countLines(sourcepath + infile)
			selections = list()
			for label in dataset_dict:
				m += 1
				sources = general.valuesort(dataset_dict[label])
				sources.reverse()
				if len(sources) > 1:
					print label, len(sources)
					n += 1
				selections.append(sources[0])
			print
		
		# transfer peak files that match the desired targets:
		print "Transferring selected datasets..."
		for infile in selections:
			command = "cp " + sourcepath + infile + " " + peakspath + infile
			os.system(command)
			j += 1
		print "Transferred selections:", j
		print "Redundancies fixed:", n
		print	
	
	
	# transfer peaks mode (from idr/peaks to data/peaks):
	if option.mode == "transfer.peaks":
		print
		
		# define input path:
		sourcepath = idrpath + "peaks/" + option.source + "/"
		
		# define output path:
		general.pathGenerator(peakspath)
		
		# gather peak files and transfer to the IDR peaks folder:
		print "Transferring peaks to data folder..."
		infiles = os.listdir(sourcepath)
		for infile in infiles:
			print "Transferring:", infile
			command = "cp " + sourcepath + infile + " " + peakspath + infile
			os.system(command)
		print		
						
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())
	
