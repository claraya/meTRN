#!/usr/bin/env python
# extract HOT-region analysis from R script outputs!

import sys
import time
import optparse
import general
import numpy
import scipy
import pickle
import pdb
import metrn
import modencode
import fasta
import random
import shutil
import os
from runner import *

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())

# define a function to generate chromosome coordinates:
def genomeCoords(chrms, chrm_size_dict):
	coord_dict = dict()
	coord_base = 1
	for chrm in chrms:
		coord_dict[chrm] = [coord_base, coord_base + chrm_size_dict[chrm]]
		coord_base += chrm_size_dict[chrm]
		coord_base += 1
	coord_min, coord_max = 1, coord_base
	return coord_dict, coord_min, coord_max, coord_dict

"""" define a function to generate region coordinates """
def regionCoords(chrms, region_size_dict, splitTag=":"):
	coord_dict, edge_dict = dict(), dict()
	coord_base = 1
	for chrm in chrms:
		if chrm in region_size_dict:
			if not chrm in coord_dict:
				coord_dict[chrm] = dict()
				edgeStart = coord_base
			for start in sorted(region_size_dict[chrm]):
				for end in sorted(region_size_dict[chrm][start]):
					coordStart, coordStop = coord_base, coord_base + region_size_dict[chrm][start][end]
					if not coord_base in coord_dict[chrm]:
						coord_dict[chrm][coordStart] = dict()
					if not coordStop in coord_dict[chrm][coordStart]:
						coord_dict[chrm][coordStart][coordStop] = sorted([start, end])
						coord_base += region_size_dict[chrm][start][end]
			coord_base += 1
			edgeStop = coord_base
			edge_dict[chrm] = [edgeStart, edgeStop]
	coord_min, coord_max = 1, coord_base
	return coord_dict, coord_min, coord_max, edge_dict

""" define a function to randomly select bases in the genome based on a coordinate dictionary """
def randomBase(coord_dict, coord_min, coord_max, edge_dict, mode="random.genomic"):
	coord = random.randrange(coord_min, coord_max)
	
	if mode == "random.genomic":
		for chrm in coord_dict:
			start, stop = coord_dict[chrm]
			if coord >= start and coord <= stop:
				return chrm, coord-start + 1
	
	elif mode == "random.regions":
		for chrm in coord_dict:
			chrmStart, chrmStop = edge_dict[chrm]
			if coord >= chrmStart and coord <= chrmStop:
				for start in sorted(coord_dict[chrm]):
					if coord >= start and coord <= max(coord_dict[chrm][start].keys()):
						for stop in sorted(coord_dict[chrm][start]):
							if coord <= stop:
								start, stop = coord_dict[chrm][start][stop]
								return chrm, coord-start + 1
	
	elif mode == "random.focused":
		for chrm in coord_dict:
			chrmStart, chrmStop = edge_dict[chrm]
			if coord >= chrmStart and coord <= chrmStop:
				start = random.choice(coord_dict[chrm].keys())
				stop = random.choice(coord_dict[chrm][start].keys())
				start, stop = coord_dict[chrm][start][stop]
				return chrm, coord-start + 1

""" define a function to... """
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

""" define a function to count peaks in a gffkde output file """
def densityPeakCounter(indata, mode="file", gffkde="OFF"):
	if mode == "file":
		processed, inlines = list(), open(indata).readlines()
	elif mode == "list":
		processed, inlines = list(), indata
	for inline in inlines:
		if gffkde == "ON":
			for peakData in inline.strip().split("\t")[7].rstrip(";").split(";")[1:]:
				processed.append(peakData.split(",")[0])
		else:
			for peakData in inline.strip().split("\t")[10].rstrip(";").split(";")[1:]:
				processed.append(peakData.split(",")[0])
	return len(list(set(processed)))

""" define a function to scan regions for peaks from the same dataset """
def gffkdeDuplicateScanner(indata, mode="file"):
	if mode == "file":
		processed, inlines = list(), open(indata).readlines()
	elif mode == "list":
		processed, inlines = list(), indata
	r, k = 0, 0
	for inline in inlines:
		datasets = list()
		for peakData in inline.strip().split("\t")[7].rstrip(";").split(";")[1:]:
			processed.append(peakData.split(",")[0])
			datasets.append(peakData.split(",")[0].split("_peaks_")[0])
		if (len(datasets) - len(set(datasets))) > 0:
			duplicates = list()
			for dataset in sorted(list(set(datasets))):
				if datasets.count(dataset) > 1:
					duplicates.append(dataset)
			print len(duplicates), len(datasets) - len(set(datasets)), ", ".join(duplicates)
			k += 1
		r += 1
	return len(list(set(processed))), r, k, round(float(k)/r, 2)
	
""" define a function to scan regions for peaks from the same dataset """
def gffkde2bed(indata, mode="file", outfile="", chrmID="chrm"):
	if mode == "file":
		processed, inlines = list(), open(indata).readlines()
	elif mode == "list":
		processed, inlines = list(), indata
	if outfile:
		f_output = open(outfile, "w")
	else:
		f_output = open(indata + ".bed", "w")
	k = 1
	for inline in inlines:
		#IV	InferredHS300bw	HOTspot	59650	59950	1.00030548632118	.	.	1,1;OP64_HLH-1_EM_yale_stn_peaks_Rank_1068,59800,1;
		chrm, method, feature, start, stop, score, strand, period, info = inline.strip().split("\t")
		feature = feature + "." + str(k)
		if chrmID == "feature":
			chrm, start, stop = feature, str(1), str(int(stop)-int(start))
		rounded = str(round(numpy.ceil(float(score))))
		print >>f_output, "\t".join([chrm, start, stop, feature, rounded, "+", score, info])
		k += 1
	f_output.close()

""" define a function to scan regions for peaks from the same dataset """
def gffkde2sizes(indata, mode="file", outfile=""):
	if mode == "file":
		processed, inlines = list(), open(indata).readlines()
	elif mode == "list":
		processed, inlines = list(), indata
	if outfile:
		f_output = open(outfile, "w")
	else:
		f_output = open(indata + "_sizes.txt", "w")
	for inline in inlines:
		#IV	InferredHS300bw	HOTspot	59650	59950	1.00030548632118	.	.	1,1;OP64_HLH-1_EM_yale_stn_peaks_Rank_1068,59800,1;
		chrm, method, feature, start, stop, score, strand, period, info = inline.strip().split("\t")
		feature = feature + "." + str(k)
		chrm, size = feature, stop
		print >>f_output, "\t".join([chrm, size])
		k += 1
	f_output.close()

""" define a function that converts a local user path to SCG3 or GS server paths """
def serverPath(inpath, server="ON"):
	if server == "ON":
		return inpath.replace("/Users/claraya/", "/srv/gs1/projects/snyder/claraya/")
	elif server == "GS":
		return inpath.replace("/Users/claraya/", "/net/fields/vol1/home/araya/")
	
def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "Path from script to files")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "Type of operations to be performed: scan or filter")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "Basename for target peaks", default="OFF")
	parser.add_option("--infile", action = "store", type = "string", dest = "infile", help = "Input HOT-regions peak/density file")
	parser.add_option("--target", action = "store", type = "string", dest = "target", help = "How should the detail splitting be performed?", default="peaks")
	parser.add_option("--contribution", action = "store", type = "float", dest = "contribution", help = "Minimum contribution cutoff", default=0)
	parser.add_option("--significance", action = "store", type = "string", dest = "significance", help = "Significance density cutoff (mode:scan) or file of HOT regions to exclude (mode:filter)!", default="OFF")
	parser.add_option("--cutoff", action = "store", type = "string", dest = "cutoff", help = "Maximum density cutoff (mode:scan) or file of HOT regions to exclude (mode:filter)!", default="OFF")
	parser.add_option("--limits", action = "store", type = "string", dest = "limits", help = "Additional text cutoffs for detail filtering!", default="OFF")
	parser.add_option("--metric", action = "store", type = "string", dest = "metric", help = "Metric: occupancy, density, or complexity", default="occupancy")
	parser.add_option("--minOccupancy", action = "store", type = "string", dest = "minOccupancy", help = "Minimum occupancy required for filtering (mode:scan)", default=2)
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "Output file name-tag.", default="OFF")
	parser.add_option("--overlap", action = "store", type = "string", dest = "overlap", help = "Overlap comparison names", default="OFF")
	parser.add_option("--source", action = "store", type = "string", dest = "source", help = "Path of peaks to be filtered or other files...")
	parser.add_option("--contexts", action = "store", type = "string", dest = "contexts", help = "What contexts of development to track in 'temperature' mode.", default="total.extended")
	parser.add_option("--regions", action = "store", type = "string", dest = "regions", help = "Regions for Monte-Carlo simulations", default="OFF")
	parser.add_option("--start", action = "store", type = "int", dest = "start", help = "Start simulation index for Monte-Carlo simulations", default=1)
	parser.add_option("--stop", action = "store", type = "int", dest = "stop", help = "End simulation index for Monte-Carlo simulations", default=1000)
	parser.add_option("--total", action = "store", type = "int", dest = "total", help = "Total simulations (indexes) for Monte-Carlo simulations", default=1000)
	parser.add_option("--shuffle", action = "store", type = "string", dest = "shuffle", help = "Should region selection be size-based (OFF; slow) or shuffle (ON; faster)", default="ON")
	parser.add_option("--scramble", action = "store", type = "string", dest = "scramble", help = "Scramble regions in which to simulate peaks? Note that this can create overlapping regions so its best not to activate.", default="ON")
	parser.add_option("--adjust", action = "store", type = "string", dest = "adjust", help = "Increase binding regions for simulations by X bases?", default="OFF")
	parser.add_option("--overwrite", action = "store", type = "string", dest = "overwrite", help = "Overwrite stuff?", default="OFF")
	parser.add_option("--round", action = "store", type = "string", dest = "round", help = "Decimal numbers to which the 'density' metrics should be rounded.", default="2")
	parser.add_option("--nuclear", action = "store", type = "string", dest = "nuclear", help = "Peaks are only nuclear?", default="ON")
	parser.add_option("--gffkde", action = "store", type = "string", dest = "gffkde", help = "Should we perform GFFKDE-analysis?", default="OFF")
	parser.add_option("--bw", action = "store", type = "string", dest = "bw", help = "GFFKDE-analysis bandwidth", default="300")
	parser.add_option("--cs", action = "store", type = "string", dest = "cs", help = "GFFKDE-analysis cutoff-score", default="0.1")
	parser.add_option("--cp", action = "store", type = "string", dest = "cp", help = "GFFKDE-analysis cutoff-peak", default="0.00001")
	parser.add_option("--pl", action = "store", type = "string", dest = "pl", help = "GFFKDE-analysis peak local optima", default="30")
	parser.add_option("--parameters", action = "store", type = "string", dest = "parameters", help = "Variable parameters...", default="")
	parser.add_option("--threads", action = "store", type = "int", dest = "threads", help = "Parallel processing threads", default=1)
	parser.add_option("--chunks", action = "store", type = "int", dest = "chunks", help = "", default=100)
	parser.add_option("--module", action = "store", type = "string", dest = "module", help = "", default="md1")
	parser.add_option("--qsub", action = "store", type = "string", dest = "qsub", help = "Qsub configuration header", default="OFF")
	parser.add_option("--server", action = "store", type = "string", dest = "server", help = "Are we on the server?", default="OFF")
	parser.add_option("--job", action = "store", type = "string", dest = "job", help = "Job name for cluster", default="OFF")
	parser.add_option("--copy", action = "store", type = "string", dest = "copy", help = "Copy simulated peaks to analysis folder?", default="OFF")
	parser.add_option("--tag", action = "store", type = "string", dest = "tag", help = "Add tag to TFBS?", default="")
	(option, args) = parser.parse_args()
	
	# import paths:
	if option.server == "OFF":
		path_dict = modencode.configBuild(option.path + "/input/" + "configure_path.txt")
	elif option.server == "ON":
		path_dict = modencode.configBuild(option.path + "/input/" + "configure_server.txt")
	elif option.server == "GS":
		path_dict = modencode.configBuild(option.path + "/input/" + "configure_nexus.txt")
	
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
	
	# specify genome size file:
	if option.nuclear == "ON":
		chromosomes = metrn.chromosomes[organismTag]["nuclear"]
		genome_size_file = option.path + "/input/" + metrn.reference[organismTag]["nuclear_sizes"]
		genome_size_dict = general.build_config(genome_size_file, mode="single", separator="\t", spaceReplace=True)
	else:
		chromosomes = metrn.chromosomes[organismTag]["complete"]
		genome_size_file = option.path + "/input/" + metrn.reference[organismTag]["complete_sizes"]
		genome_size_dict = general.build_config(genome_size_file, mode="single", separator="\t", spaceReplace=True)
		
	# predefine and prebuild input/output paths:
	densitypath = hotpath + "density/"
	analysispath = hotpath + "analysis/"
	comparepath = hotpath + "compare/"
	overlappath = hotpath + "overlap/"
	simulationpath = hotpath + "simulation/"
	regionspath = hotpath + "regions/"
	temperaturepath = hotpath + "temperature/"
	general.pathGenerator(densitypath)
	general.pathGenerator(analysispath)
	general.pathGenerator(comparepath)
	general.pathGenerator(overlappath)
	general.pathGenerator(simulationpath)
	general.pathGenerator(regionspath)
	general.pathGenerator(temperaturepath)
	
	# generate gffkde handle:
	gffkde_handle = "_bw" + option.bw.replace("0.","") + "_cs" + option.cs.replace("0.","") + "_cp" + option.cp.replace("0.","") + "_pl" + option.pl.replace("0.","")
	
	# check that the index range is coherent:
	if option.stop > option.total:
		print
		print "Error: Range exceeded! Stop index is larger than total."
		print
		return
	
	# master mode:
	if "master:" in option.mode:
	
		# capture master mode:
		master, mode = option.mode.split(":")
		
		# determine whether peaks should cluster on a predefined number of regions:
		if option.shuffle == "OFF" and option.regions == "OFF":
			regions_flag = "random.genomic"
			regions = False
		elif option.shuffle == "OFF" and option.regions == "ON":
			regions_flag = "random.regions"
			regions = True
		elif option.shuffle == "ON" and option.regions == "OFF":
			regions_flag = "shuffle.genomic"
			regions = True
		elif option.shuffle == "ON" and option.regions != "OFF":
			regions_flag = "shuffle.regions"
			regions = True
		
		# prepare for qsub:
		bash_path = str(option.path + "/data/hot/simulation/runs/").replace("//","/")
		bash_base = "_".join([mode, regions_flag, option.peaks, option.name]) + "-M"
		qsub_base = "_".join([mode, regions_flag, option.peaks, option.name])
		general.pathGenerator(bash_path)
		if option.qsub != "OFF":
			qsub_header = open(qsubpath + option.qsub).read()
			qsub = True
		else:
			qsub_header = ""
			qsub = False
		if option.job == "QSUB":
			qsub_header = qsub_header.replace("qsubRunner", "qsub-" + qsub_base)
		elif option.job != "OFF":
			qsub_header = qsub_header.replace("qsubRunner", "qsub-" + option.job)
			bash_base = option.job + "-M"
			
		# update server path:
		if option.qsub != "OFF":
			option.path = serverPath(option.path, server=option.server)
			if option.server != "OFF":
				bash_mode = ""
			
		# prepare slave modules:
		m, modules, commands, sequences, chunks, start, complete = 1, list(), list(), list(), option.chunks, option.start, False
		for index in range(option.start, option.stop+1):
			run = "mc" + general.indexTag(index, option.total)
			
			# montecarlo peak simulation mode:
			if mode == "montecarlo":
			
				command = "python <<CODEPATH>>mapHOT.py --path <<PATH>> --organism <<ORGANISM>> --mode <<MODE>> --peaks <<PEAKS>> --start <<START>> --stop <<STOP>> --total <<TOTAL>> --contexts <<CONTEXTS>> --regions <<REGIONS>> --shuffle <<SHUFFLE>> --metric <<METRIC>> --scramble <<SCRAMBLE>> --adjust <<ADJUST>> --overwrite <<OVERWRITE>> --round <<ROUND>> --name <<NAME>> --qsub <<QSUB>> --server <<SERVER>> --copy <<COPY>> --module <<MODULE>>"
				command = command.replace("<<CODEPATH>>", option.path + "/python/")
				command = command.replace("<<PATH>>", option.path)
				command = command.replace("<<ORGANISM>>", option.organism)
				command = command.replace("<<MODE>>", mode)
				command = command.replace("<<PEAKS>>", option.peaks)
				command = command.replace("<<START>>", str(start))
				command = command.replace("<<STOP>>", str(index))
				command = command.replace("<<TOTAL>>", str(option.total))
				command = command.replace("<<CONTEXTS>>", option.contexts)
				command = command.replace("<<REGIONS>>", option.regions)
				command = command.replace("<<SHUFFLE>>", option.shuffle)
				command = command.replace("<<METRIC>>", option.metric)
				command = command.replace("<<SCRAMBLE>>", option.scramble)
				command = command.replace("<<ADJUST>>", option.adjust)
				command = command.replace("<<OVERWRITE>>", option.overwrite)
				command = command.replace("<<ROUND>>", option.round)
				command = command.replace("<<NAME>>", option.name)
				command = command.replace("<<QSUB>>", option.qsub)
				command = command.replace("<<SERVER>>", option.server)
				command = command.replace("<<COPY>>", option.copy)
				command = command.replace("<<MODULE>>", "md" + str(m))
				
			# simulation analysis mode:
			if mode == "simulation":
			
				command = "python <<CODEPATH>>mapHOT.py --path <<PATH>> --organism <<ORGANISM>> --mode <<MODE>> --peaks <<PEAKS>> --start <<START>> --stop <<STOP>> --total <<TOTAL>> --contexts <<CONTEXTS>> --regions <<REGIONS>> --shuffle <<SHUFFLE>> --metric <<METRIC>> --scramble <<SCRAMBLE>> --adjust <<ADJUST>> --overwrite <<OVERWRITE>> --round <<ROUND>> --name <<NAME>> --qsub <<QSUB>> --server <<SERVER>> --copy <<COPY>> --module <<MODULE>>"
				command = command.replace("<<CODEPATH>>", option.path + "/python/")
				command = command.replace("<<PATH>>", option.path)
				command = command.replace("<<ORGANISM>>", option.organism)
				command = command.replace("<<MODE>>", mode)
				command = command.replace("<<PEAKS>>", option.peaks)
				command = command.replace("<<START>>", str(start))
				command = command.replace("<<STOP>>", str(index))
				command = command.replace("<<TOTAL>>", str(option.total))
				command = command.replace("<<CONTEXTS>>", option.contexts)
				command = command.replace("<<REGIONS>>", option.regions)
				command = command.replace("<<SHUFFLE>>", option.shuffle)
				command = command.replace("<<METRIC>>", option.metric)
				command = command.replace("<<SCRAMBLE>>", option.scramble)
				command = command.replace("<<ADJUST>>", option.adjust)
				command = command.replace("<<OVERWRITE>>", option.overwrite)
				command = command.replace("<<ROUND>>", option.round)
				command = command.replace("<<NAME>>", option.name)
				command = command.replace("<<QSUB>>", option.qsub)
				command = command.replace("<<SERVER>>", option.server)
				command = command.replace("<<COPY>>", option.copy)
				command = command.replace("<<MODULE>>", "md" + str(m))
					
			# simulation collecting mode:
			if mode == "collecting":
			
				command = "python <<CODEPATH>>mapHOT.py --path <<PATH>> --organism <<ORGANISM>> --mode <<MODE>> --peaks <<PEAKS>> --start <<START>> --stop <<STOP>> --total <<TOTAL>> --contexts <<CONTEXTS>> --regions <<REGIONS>> --shuffle <<SHUFFLE>> --metric <<METRIC>> --scramble <<SCRAMBLE>> --adjust <<ADJUST>> --overwrite <<OVERWRITE>> --round <<ROUND>> --name <<NAME>> --qsub <<QSUB>> --server <<SERVER>> --copy <<COPY>> --module <<MODULE>>"
				command = command.replace("<<CODEPATH>>", option.path + "/python/")
				command = command.replace("<<PATH>>", option.path)
				command = command.replace("<<ORGANISM>>", option.organism)
				command = command.replace("<<MODE>>", mode)
				command = command.replace("<<PEAKS>>", option.peaks)
				command = command.replace("<<START>>", str(option.start))
				command = command.replace("<<STOP>>", str(option.stop))
				command = command.replace("<<TOTAL>>", str(option.total))
				command = command.replace("<<CONTEXTS>>", option.contexts)
				command = command.replace("<<REGIONS>>", option.regions)
				command = command.replace("<<SHUFFLE>>", option.shuffle)
				command = command.replace("<<METRIC>>", option.metric)
				command = command.replace("<<SCRAMBLE>>", option.scramble)
				command = command.replace("<<ADJUST>>", option.adjust)
				command = command.replace("<<OVERWRITE>>", option.overwrite)
				command = command.replace("<<ROUND>>", option.round)
				command = command.replace("<<NAME>>", option.name)
				command = command.replace("<<QSUB>>", option.qsub)
				command = command.replace("<<SERVER>>", option.server)
				command = command.replace("<<COPY>>", option.copy)
				command = command.replace("<<MODULE>>", "md" + str(m))
				break
					
			# is it time to export a chunk?
			if index-start+1 == chunks:
				
				# update start, modules, commands, and module count (m):
				start = index + 1
				commands.append(command)
				modules.append(commands)
				commands = list()
				complete = True
				m += 1
			
			# store whether the most recent index/command has been stored:
			else:
				complete = False
		
		# update if there are additional commands:
		if not complete:
			commands.append(command)
			modules.append(commands)
			m += 1
		
		# launch commands:
		print
		print "Launching comparisons:", len(modules)
		#for module in modules:
		#	for command in module:
		#		print command
		runCommands(modules, threads=option.threads, mode="module.run", run_mode="verbose", run_path=bash_path, run_base=bash_base, record=True, qsub_header=qsub_header, qsub=qsub)
		print "Comparisons performed:", len(modules)
		print
	
	# perform Monte-Carlo simulation of peaks:
	elif option.mode == "montecarlo":
	
		# prepare chromosome sizes:
		chrm_size_dict, chrm_base_dict = dict(), dict()
		for chrm in chromosomes:
			chrm_size_dict[chrm] = genome_size_dict[chrm]
			chrm_base_dict[chrm] = [1, genome_size_dict[chrm]]
			
		# get context code and target contexts:
		contextCode, contextTargets = metrn.contexts[organismTag]["list"][option.contexts]
		
		# determine whether peaks should cluster on a predefined number of regions:
		if option.shuffle == "OFF" and option.regions == "OFF":
			regions_flag = "random.genomic"
			regions = False
		elif option.shuffle == "OFF" and option.regions == "ON":
			regions_flag = "random.regions"
			regions = True
		elif option.shuffle == "ON" and option.regions == "OFF":
			regions_flag = "shuffle.genomic"
			regions = True
		elif option.shuffle == "ON" and option.regions != "OFF":
			regions_flag = "shuffle.regions"
			regions = True
		
		# shuffle-based methods:
		if option.shuffle == "ON":
			
			# define simulation setup path:
			carbonpath = hotpath + "simulation/montecarlo/" + option.peaks + "/" + regions_flag + "/"
			general.pathGenerator(carbonpath)
			
			# prebuild overlap regions if necessary:
			print
			print "Preparing configuration files..."
			if regions_flag == "shuffle.regions":
					
				# define density file (gffkde), regions file (bed), annexed file (bed), overlap file (bed):
				kdensityfile = densitypath + "gffkde2_" + option.regions + gffkde_handle + ".Peaks.Extended"
				regionalfile = carbonpath + "maphot_montecarlo_master_regional." + option.module + ".bed"
				completefile = carbonpath + "maphot_montecarlo_master_complete." + option.module + ".bed"
				collapsefile = carbonpath + "maphot_montecarlo_master_collapse." + option.module + ".bed"
				compiledfile = carbonpath + "maphot_montecarlo_master_compiled." + option.module + ".bed"
				
				#gffkde2bed(kdensityfile, mode="file", outfile=regionalfile)
				metrn.completeBed(peakspath + option.peaks + "/", outfile=completefile)
				
				# make overlap regions (collapse peaks):
				command = "mergeBed -i " + completefile + " -nms > " + collapsefile
				os.system(command)
				
				# generate a density-style analysis file:
				metrn.densityBed(collapsefile, compiledfile)
			
			# simulate peaks per factor:
			print "Launching simulations:", time.asctime(time.localtime())
			for i in range(option.start, option.stop+1):
				run = "mc" + general.indexTag(i, option.total)
				print
				print "\tSimulation:", run
				
				if option.overwrite == "ON" or not run in os.listdir(carbonpath):
				
					# define simulation paths:
					randompath = carbonpath + run + "/random/"
					configpath = carbonpath + run + "/config/"
					resultpath = carbonpath + run + "/result/"
					
					# clear prexisting simulations:
					if run in os.listdir(carbonpath):
						command = "rm -rf " + randompath
						os.system(command)
					
					# generate output paths:
					general.pathGenerator(randompath)
					general.pathGenerator(configpath)
					general.pathGenerator(resultpath)
					c = 0
					
					# generate exclusion regions if necessary:
					if regions_flag == "shuffle.regions":
						
						# define shuffle file (bed), exclude file (bed):
						shufflefile = configpath + "maphot_montecarlo_" + run + "_shuffle.bed"
						regionsfile = configpath + "maphot_montecarlo_" + run + "_regions.bed"
						excludefile = configpath + "maphot_montecarlo_" + run + "_exclude.bed"
						
						# make shuffled regions:
						if option.scramble == "ON":
							command = "shuffleBed -i " + collapsefile + " -g " + genome_size_file + " -chrom > " + shufflefile
							os.system(command)
						else:
							shufflefile = collapsefile
							
						# make adjusted regions:
						if option.adjust != "OFF":
							command = "slopBed -i " + shufflefile + " -g " + genome_size_file + " -b " + option.adjust + " > " + regionsfile
							os.system(command)
						else:
							regionsfile = shufflefile
						
						# make excluded regions:
						command = "complementBed -i " + regionsfile + " -g " + genome_size_file + " > " + excludefile
						os.system(command)
					
					# simulate the peaks:
					for peak_file in os.listdir(peakspath + option.peaks):
						organism, strain, factor, context, institute, method = metrn.labelComponents(peak_file)
						if option.contexts in ["all","any","XX","total.extended","total.condense"] or context in contextTargets and ".bed" in peak_file:
							c += 1
							dataset = peak_file.split("_peaks.")[0]
							strain, factor, context, institute, method = peak_file.split("_")[:5]
							inpeakfile = peakspath + option.peaks + "/" + peak_file
							randomfile = randompath + peak_file.replace("_peaks.bed", "_random.bed")
							print "\t", run, ":", peak_file
						
							# genomic shuffling:
							if regions_flag == "shuffle.genomic":
								command = "shuffleBed -i " + inpeakfile + " -g " + genome_size_file + " -chrom > " + randomfile
								os.system(command)	
								
							# regions shuffling:
							elif regions_flag == "shuffle.regions":
								command = "shuffleBed -i " + inpeakfile + " -g " + genome_size_file + " -excl " + excludefile + " > " + randomfile
								os.system(command)
				
					# define run complete and collapsed files:
					completerun = resultpath + "maphot_montecarlo_" + run + "_complete.bed"
					collapserun = resultpath + "maphot_montecarlo_" + run + "_collapse.bed"
					
					# make run complete and collapsed files:
					metrn.completeBed(randompath, outfile=completerun)
					command = "mergeBed -i " + completerun + " -nms > " + collapserun
					os.system(command)
		
		# non-shuffle methods:
		elif option.shuffle == "OFF":

			# prepare peak size and peak score dictionaries:
			maxSize = 0
			peak_size_dict, peak_score_dict = dict(), dict()
			for peak_file in os.listdir(peakspath + option.peaks):
				strain, factor, context, institute, method = peak_file.split("_")[:5]
				if option.contexts in ["all","any","XX"] or context in contextTargets and ".bed" in peak_file:
					inlines = open(peakspath + option.peaks + "/" + peak_file).readlines()
					for inline in inlines:
						chrm, start, end, feature, score = inline.strip().split("\t")[:5]
						if not peak_file in peak_size_dict:
							peak_size_dict[peak_file] = list()
							peak_score_dict[peak_file] = list()
						peak_size_dict[peak_file].append(int(end)-int(start))
						peak_score_dict[peak_file].append(int(float(score)))
						if int(end)-int(start) > maxSize:
							maxSize = int(end)-int(start)

			# determine whether peak should cluster on a predefined number of regions:
			if regions:
				
				# match density file:
				process = False
				for densityfile in os.listdir(densitypath):
					if "gffkde2_" + option.regions + "_bw" in densityfile and ".Peaks.Detail" in densityfile:
						process = True
						print "Matched density file:", densityfile
						break
				
				print "Preparing region dictionaries..."
				if process:
					region_size_dict, region_base_dict = dict(), dict()
					inlines = open(densitypath + densityfile).readlines()
					for inline in inlines:
						if inline[0] != "#":
							chrm, parameters, hotspot, start, end = inline.strip().split("\t")[:5]
							chrm, start, end = chrm.replace("MTDNA","MtDNA"), int(start), int(end)
							start, end = sorted([start, end])
							if not chrm in region_size_dict:
								region_size_dict[chrm] = dict()
								region_base_dict[chrm] = dict()
							if not start in region_size_dict[chrm]:
								region_size_dict[chrm][start] = dict()
								region_base_dict[chrm][start] = dict()
							region_size_dict[chrm][start][end] = abs(end-start)
							region_base_dict[chrm][start][end] = [start, end]
	 			else:
 					raise "Error: Density file not found!"
 				
			# generate coordinate system:
			print "Generating coordinate dictionaries..."
			if option.regions == "OFF":
				coord_dict, coord_min, coord_max, edge_dict = genomeCoords(chromosomes, chrm_size_dict)
				base_dict = chrm_base_dict
			else:
				coord_dict, coord_min, coord_max, edge_dict = regionCoords(chromosomes, region_size_dict)
				base_dict = region_base_dict
			
			# simulate peaks per factor:
			print "Starting simulations:", time.asctime(time.localtime())
			for i in range(option.start, option.stop+1):
				run = "mc" + general.indexTag(i, option.total)
				print "Simulation index:", run
				
				# define simulation paths:
				randompath = hotpath + "simulation/montecarlo/" + option.peaks + "/" + regions_flag + "/" + run + "/random/"
				configpath = hotpath + "simulation/montecarlo/" + option.peaks + "/" + regions_flag + "/" + run + "/config/"
				general.pathGenerator(randompath)
				general.pathGenerator(configpath)
				c = 0
				for peak_file in peak_size_dict:
					c += 1
					dataset = peak_file.split("_peaks.")[0]
					strain, factor, context, institute, method = peak_file.split("_")[:5]
					f_output = open(montecarlopath + peak_file, "w")
					k = 1
					#print "Processing:", peak_file, len(peak_file), "peaks"
					for size in peak_size_dict[peak_file]:
						#print "Size:", size
						chrm, base = randomBase(coord_dict, coord_min, coord_max, edge_dict, mode=regions_flag)
						#print "Random:", chrm, base
						while (size > maxSize and chrm == "MtDNA") or (base + size > chrm_size_dict[chrm]):
							chrm, base = randomBase(coord_dict, coord_min, coord_max, edge_dict, mode=regions_flag)
							#print "Repeat:", chrm, base
						score = peak_score_dict[peak_file][k-1]
						code = dataset + ".Random_peak_" + str(k)
						peak = "Random_peak_" + str(k)
						output = [chrm, base, base + size, code, score, "+", strain, factor, context, institute, method, peak]
						print >>f_output, "\t".join(map(str, output))
						k += 1
						#pdb.set_trace()
					f_output.close()
		print "Simulations complete!"
		print
		
		
	# perform simulation sampling peaks from pre-computed random peaks (monte-carlo):
	elif option.mode == "simulation":
	
		# prepare chromosome sizes:
		chrm_size_dict, chrm_base_dict = dict(), dict()
		for chrm in chromosomes:
			chrm_size_dict[chrm] = genome_size_dict[chrm]
			chrm_base_dict[chrm] = [1, genome_size_dict[chrm]]
			
		# get context code and target contexts:
		contextCode, contextTargets = metrn.contexts[organismTag]["list"][option.contexts]
		
		# determine whether peaks should cluster on a predefined number of regions:
		if option.shuffle == "OFF" and option.regions == "OFF":
			regions_flag = "random.genomic"
			regions = False
		elif option.shuffle == "OFF" and option.regions == "ON":
			regions_flag = "random.regions"
			regions = True
		elif option.shuffle == "ON" and option.regions == "OFF":
			regions_flag = "shuffle.genomic"
			regions = True
		elif option.shuffle == "ON" and option.regions != "OFF":
			regions_flag = "shuffle.regions"
			regions = True

		# setup montecarlo simulation (source) path:
		carbonpath = hotpath + "simulation/montecarlo/" + option.peaks + "/" + regions_flag + "/"
		
		# define configuration files:
		completefile = "maphot_montecarlo_master_complete.bed"
		collapsefile = "maphot_montecarlo_master_collapse.bed"
		compiledfile = "maphot_montecarlo_master_compiled.bed"
		
		# condense (remove replicate modules):
		if not completefile in os.listdir(carbonpath):
			command = "cp " + carbonpath + completefile.replace(".bed", ".md1.bed") + " " + carbonpath + completefile
			os.system(command)
		if not collapsefile in os.listdir(carbonpath):
			command = "cp " + carbonpath + collapsefile.replace(".bed", ".md1.bed") + " " + carbonpath + collapsefile
			os.system(command)
		if not compiledfile in os.listdir(carbonpath):
			command = "cp " + carbonpath + compiledfile.replace(".bed", ".md1.bed") + " " + carbonpath + compiledfile
			os.system(command)
		command = "rm -rf " + carbonpath + "maphot_montecarlo_master_complete.md*bed"
		os.system(command)
		command = "rm -rf " + carbonpath + "maphot_montecarlo_master_collapse.md*bed"
		os.system(command)
		command = "rm -rf " + carbonpath + "maphot_montecarlo_master_compiled.md*bed"
		os.system(command)
		
		# transfer simulated peaks and execute density analysis:
		print
		for i in range(option.start, option.stop+1):
			run = "mc" + general.indexTag(i, option.total)
			print "\tSimulation:", run
			
			# define simulation paths:
			randompath = carbonpath + run + "/random/"
			configpath = carbonpath + run + "/config/"
			resultpath = carbonpath + run + "/result/"
			#general.pathGenerator(randompath)
			#general.pathGenerator(configpath)
			#general.pathGenerator(resultpath)
			
			# define analysis paths:
			subsetpath = hotpath + "simulation/analysis/" + option.peaks + "/" + regions_flag + "/" + option.name + "/" + run + "/subset/"
			gffkdepath = hotpath + "simulation/analysis/" + option.peaks + "/" + regions_flag + "/" + option.name + "/" + run + "/gffkde/"
			mergedpath = hotpath + "simulation/analysis/" + option.peaks + "/" + regions_flag + "/" + option.name + "/" + run + "/merged/"
			general.pathGenerator(subsetpath)
			general.pathGenerator(gffkdepath)
			general.pathGenerator(mergedpath)
			
			# define complete peaks file and merged (collapse) peaks file:
			completefile = mergedpath + "maphot_montecarlo_" + run + "_complete.bed"
			collapsefile = mergedpath + "maphot_montecarlo_" + run + "_collapse.bed"
			compiledfile = mergedpath + "maphot_montecarlo_" + run + "_compiled.bed"
			
			c = 0
			# scan simulation peaks:
			for randomfile in os.listdir(randompath):
				strain, factor, context, institute, method = randomfile.split("_")[:5]
				dataset = "_".join([strain, factor, context, institute, method])
				if option.contexts in ["all","any","XX","total.extended","total.condense"] or context in contextTargets and ".bed" in randomfile:
					c += 1
					
					# transfer simulated peaks:
					#command = "cp " + randompath + randomfile + " " + subsetpath + randomfile
					#os.system(command)
			
			""" This section added to regenerate compiledfile... """
			
			# make complete simulated peaks file (complete peaks):
			#metrn.completeBed(subsetpath, outfile=completefile)
				
			# make overlap regions (collapse peaks):
			#command = "mergeBed -i " + completefile + " -nms > " + collapsefile
			#os.system(command)
				
			# generate a density-style analysis file:
			metrn.densityBed(collapsefile, compiledfile)
			
			"""
			# launch kernel-density analysis (if desired):
			if option.gffkde == "ON":
				command = "perl <<SCRIPTSPATH>>GFFKDE2.pl <<RANDOMPATH>> <<BW>> <<CS>> <<CP>> <<PL>> <<GFFKDEPATH>>gffkde2_<<NAME>><<GFFKDEHANDLE>>"
				command = command.replace("<<SCRIPTSPATH>>", scriptspath)
				command = command.replace("<<RANDOMPATH>>", randompath)
				command = command.replace("<<GFFKDEPATH>>", gffkdepath)
				command = command.replace("<<GFFKDEHANDLE>>", gffkde_handle)
				command = command.replace("<<BW>>", option.bw)
				command = command.replace("<<CS>>", option.cs)
				command = command.replace("<<CP>>", option.cp)
				command = command.replace("<<PL>>", option.pl)
				command = command.replace("<<NAME>>", option.name)
				print command
				os.system(command)
				print
			
				# remove large density file:
				command = "rm -rf <<GFFKDEPATH>>gffkde2_<<NAME>><<GFFKDEHANDLE>>.Density"
				command = command.replace("<<GFFKDEPATH>>", gffkdepath)
				command = command.replace("<<GFFKDEHANDLE>>", gffkde_handle)
				command = command.replace("<<NAME>>", option.name)
				os.system(command)
			
				# remove unwanted region file:
				command = "rm -rf <<GFFKDEPATH>>gffkde2_<<NAME>><<GFFKDEHANDLE>>.Peaks"
				command = command.replace("<<GFFKDEPATH>>", gffkdepath)
				command = command.replace("<<GFFKDEHANDLE>>", gffkde_handle)
				command = command.replace("<<NAME>>", option.name)
				os.system(command)
			
				# remove unwanted detail file:
				command = "rm -rf <<GFFKDEPATH>>gffkde2_<<NAME>><<GFFKDEHANDLE>>.Peaks.Detail"
				command = command.replace("<<GFFKDEPATH>>", gffkdepath)
				command = command.replace("<<GFFKDEHANDLE>>", gffkde_handle)
				command = command.replace("<<NAME>>", option.name)
				os.system(command)
			
			# copy (keep) or remove the used peaks?
			if option.copy == "OFF":
				command = "rm -rf " + subsetpath + "*.bed"
				os.system(command)
			"""
	
	# collect results from simulated binding regions with random peaks (monte-carlo, simulation):
	elif option.mode == "collecting":
		
		# set analysis metric and flag:
		if option.metric == "occupancy":
			metric_flag = "_occupan"
		elif option.metric == "density":
			metric_flag = "_density"
		elif option.metric == "complexity":
			metric_flag = "_complex"
	
		# determine whether peaks should cluster on a predefined number of regions:
		if option.shuffle == "OFF" and option.regions == "OFF":
			regions_flag = "random.genomic"
			regions = False
		elif option.shuffle == "OFF" and option.regions == "ON":
			regions_flag = "random.regions"
			regions = True
		elif option.shuffle == "ON" and option.regions == "OFF":
			regions_flag = "shuffle.genomic"
			regions = True
		elif option.shuffle == "ON" and option.regions != "OFF":
			regions_flag = "shuffle.regions"
			regions = True
		
		# setup montecarlo simulation (source) path:
		carbonpath = hotpath + "simulation/montecarlo/" + option.peaks + "/" + regions_flag + "/"
		
		# get context code and target contexts:
		contextCode, contextTargets = metrn.contexts[organismTag]["list"][option.contexts]
		
		# transfer analyzed mergedBed compiled (density) file to density path:
		compiledfile = carbonpath + "maphot_montecarlo_master_compiled.bed"
		densityfile = "mergeBed_" + option.peaks + "_compiled.bed"
		if not densityfile in os.listdir(densitypath):
			command = "cp " + compiledfile + " " + densitypath + densityfile
			os.system(command)
		
		# define analysis path:
		reportpath = hotpath + "simulation/analysis/" + option.peaks + "/" + regions_flag + "/" + option.name + "/"
			
		# prepare region reports:
		f_output = open(reportpath + "maphot_simulation_report" + metric_flag + "_region.txt", "w")
		
		# gather density frequencies from simulations:
		print
		value_dict, simulation_dict, regions_dict, dataset_dict, factor_dict, context_dict = dict(), dict(), dict(), dict(), dict(), dict()
		for i in range(option.start, option.stop+1):
			run = "mc" + general.indexTag(i, option.total)
			#print "\tSimulation:", run
			
			sys.stdout.write("\rCollecting densities from simulations: {0}".format(run))
			sys.stdout.flush()
			
			# define simulation paths:
			randompath = carbonpath + run + "/random/"
			configpath = carbonpath + run + "/config/"
			
			# define analysis paths:
			subsetpath = hotpath + "simulation/analysis/" + option.peaks + "/" + regions_flag + "/" + option.name + "/" + run + "/subset/"
			gffkdepath = hotpath + "simulation/analysis/" + option.peaks + "/" + regions_flag + "/" + option.name + "/" + run + "/gffkde/"
			mergedpath = hotpath + "simulation/analysis/" + option.peaks + "/" + regions_flag + "/" + option.name + "/" + run + "/merged/"
			general.pathGenerator(reportpath)
			
			# define complete peaks file, merged (collapse) peaks file, and GFFKDE density file:
			kdensityfile = gffkdepath + "gffkde2_<<NAME>><<GFFKDEHANDLE>>.Peaks.Extended".replace("<<NAME>>", option.name).replace("<<GFFKDEHANDLE>>", gffkde_handle)
			completefile = mergedpath + "maphot_montecarlo_" + run + "_complete.bed"
			collapsefile = mergedpath + "maphot_montecarlo_" + run + "_collapse.bed"
			compiledfile = mergedpath + "maphot_montecarlo_" + run + "_compiled.bed"
			
			# load kernel-density data (if desired):
			if option.gffkde == "ON":
			
				# load densities for individual simulation:
				inlines = open(kdensityfile).readlines()
				for inline in inlines:
					chrm, method, feature, start, stop, score, strand, period, info = inline.strip().split("\t")
					
					# set analysis metric and flag:
					if option.metric == "occupancy":
						value = int(round(float(score)))
					elif option.metric == "density":
						size = int(stop) - int(start) + 1
						value = 1000*float(int(round(float(score))))/size
						if option.round != "OFF":
							value = round(value, int(option.round))
					elif option.metric == "complexity":
						raise "Error: metric not implemented!"
	
					# store simulation value:
					if not run in simulation_dict:
						simulation_dict[run] = dict()
					if not value in simulation_dict[run]:
						simulation_dict[run][value] = 0
					simulation_dict[run][value] += 1
			
			# load mergedBed density data:
			else:
			
				# load densities for individual simulation:
				inlines = open(compiledfile).readlines()
				index = 1
				for inline in inlines:
					chrm, start, stop, feature, occupancy, strand, density, datasetCount, factorCount, contextCount, info = inline.strip().split("\t")
					
					# set analysis metric and flag:
					if option.metric == "occupancy":
						value = int(occupancy)
					elif option.metric == "density":
						value = float(density)
						if option.round != "OFF":
							value = round(value, int(option.round))
					elif option.metric == "complexity":
						value = int(datasetCount)
					
					# store simulation value:
					if not run in simulation_dict:
						simulation_dict[run] = dict()
					if not value in simulation_dict[run]:
						simulation_dict[run][value] = 0
					simulation_dict[run][value] += 1
					
					# export simulated region data:
					size = str(int(stop) - int(start) + 1)
					complexity = datasetCount
					regularity = str(1000*float(datasetCount)/int(size))
					print >>f_output, "\t".join([run + "." + str(index), occupancy, density, complexity, regularity, size])
					index += 1
					
			# load densities for simulation summary:
			for value in simulation_dict[run]:
				if not value in value_dict:
					value_dict[value] = list()
				value_dict[value].append(simulation_dict[run][value])
		
			# load region count for individual simulation summary:
			regions_dict[run] = len(inlines)
		
		# determine number of regions and simulations analyzed:
		regions = sum(regions_dict.values())
		simulations = len(simulation_dict)
		
		# prepare density reports:
		g_output = open(reportpath + "maphot_simulation_report" + metric_flag + "_global.txt", "w")
		s_output = open(reportpath + "maphot_simulation_report" + metric_flag + "_single.txt", "w")
		v_output = open(reportpath + "maphot_simulation_report" + metric_flag + "_values.txt", "w")
		
		# print headers:
		print >>g_output, "\t".join(["value", "regions.sum", "regions.mean", "regions.median", "regions.std", "regions.total", "simulations"])
		print >>s_output, "\t".join(["simulation", "value", "regions.sum", "regions.mean", "regions.median", "regions.std", "regions.total", "simulations"])
		print >>v_output, "\t".join(["simulation", "value"])
		
		# export global simulation density report:
		print "\nExporting global simulation data..."
		for value in sorted(value_dict.keys()):
			values = value_dict[value]
			print >>g_output, "\t".join(map(str, [value, sum(values), numpy.mean(values), int(numpy.median(values)), numpy.std(values), regions, simulations]))
		
		# export single simulation density report:
		print "Exporting individual simulation data..."
		for simulation in sorted(simulation_dict.keys()):
			values = list()
			for value in sorted(simulation_dict[simulation].keys()):
				observations = [value]*simulation_dict[simulation][value]
				values += observations
				for observation in observations:
					print >>v_output, "\t".join([simulation, str(observation)])
			for value in sorted(simulation_dict[simulation].keys()):
				print >>s_output, "\t".join(map(str, [simulation, value, simulation_dict[simulation][value], numpy.mean(values), int(numpy.median(values)), numpy.std(values), sum(simulation_dict[simulation].values()), simulations]))	
		
		# close output files:
		f_output.close()
		g_output.close()
		s_output.close()
		v_output.close()
		print
		
	
	# download results of simulated binding regions from cluster:
	elif option.mode == "downloader":
	
		# determine whether peaks should cluster on a predefined number of regions:
		if option.shuffle == "OFF" and option.regions == "OFF":
			regions_flag = "random.genomic"
			regions = False
		elif option.shuffle == "OFF" and option.regions == "ON":
			regions_flag = "random.regions"
			regions = True
		elif option.shuffle == "ON" and option.regions == "OFF":
			regions_flag = "shuffle.genomic"
			regions = True
		elif option.shuffle == "ON" and option.regions != "OFF":
			regions_flag = "shuffle.regions"
			regions = True
		
		# define server path dictionary:
		server_dict = modencode.configBuild(option.path + "/input/" + "configure_server.txt")
		
		# define simulation (source) paths, remote (serverpath) and local (reportpath):
		serverpath = server_dict["hot"] + "simulation/analysis/" + option.peaks + "/" + regions_flag + "/" + option.name + "/"
		reportpath = hotpath + "simulation/analysis/" + option.peaks + "/" + regions_flag + "/" + option.name + "/"
		general.pathGenerator(reportpath)
		
		# recover server login base:
		command = "rsync -avz " + option.parameters + ":" + serverpath + "maphot_* " + reportpath
		os.system(command)

		# recover density (source) paths, remote (serverpath) and local (reportpath):
		serverpath = server_dict["hot"] + "simulation/montecarlo/" + option.peaks + "/" + regions_flag + "/"
		reportpath = hotpath + "simulation/montecarlo/" + option.peaks + "/" + regions_flag + "/"
		general.pathGenerator(reportpath)
		
		# recover server login base:
		command = "rsync -avz " + option.parameters + ":" + serverpath + "maphot_* " + reportpath
		os.system(command)


	# perform scan to pull-out HOT regions:
	elif option.mode == "scan":
	
		# determine input density file:
		if option.gffkde == "ON":
			densityfile = "gffkde2_" + option.peaks + gffkde_handle + ".Peaks.Extended"
		else:
			densityfile = "mergeBed_" + option.peaks + "_compiled.bed"
		
		# set analysis metric and flag:
		if option.metric == "occupancy":
			metric_flag = "_occupan"
		elif option.metric == "density":
			metric_flag = "_density"
		elif option.metric == "complexity":
			metric_flag = "_complex"
		elif option.metric == "combined":
			metric_flag = "_combine"
			
		# set contribution behavior:
		if option.contribution == 0:
			contribution_handle = "_con" + "XX"
		else:
			contribution_handle = "_con" + "%.0e" % (float(option.contribution))
		
		# set cutoff behavior:
		if option.cutoff == "OFF" or option.cutoff == "X":
			cutoff_status = False
			cutoff_handle = "_cut" + "XX"
		elif option.significance != "OFF":
			cutoff_status = True
			cutoff_handle = "_sig" + option.significance
		else:
			cutoff_status = True
			cutoff_handle = "_cut" + str(option.cutoff).zfill(3)
		
		# set configuration output handles:
		if option.name == "OFF":
			configuration_name = "_simple"
			configuration_handle = "_simple" + cutoff_handle + contribution_handle
		else:
			configuration_name = "_" + option.name
			configuration_handle = "_" + option.name + cutoff_handle + contribution_handle
		
		# generate output files:
		firstPart, lastPart = densityfile.split(option.peaks)
		p_outfile = "maphot_" + option.target + configuration_name + "_" + firstPart + option.peaks + configuration_handle + lastPart.lower().replace(".bed","").replace(".Peaks.Extended","_compiled").replace("_compiled", metric_flag).replace(".","_") + "_pass.bed"
		f_outfile = "maphot_" + option.target + configuration_name + "_" + firstPart + option.peaks + configuration_handle + lastPart.lower().replace(".bed","").replace(".Peaks.Extended","_compiled").replace("_compiled", metric_flag).replace(".","_") + "_fail.bed"
		
		# open input and output files:
		p_output = open(analysispath + p_outfile, "w")
		f_output = open(analysispath + f_outfile, "w")
		#print >>p_output, "\t".join(["chrm", "start", "end", "feature", "score", "strand", "density", "details"])
		#print >>f_output, "\t".join(["chrm", "start", "end", "feature", "score", "strand", "density", "details"])
		
		# scan HOT-region analysis output:
		fi, pi = 1, 1
		inlines = open(densitypath + densityfile).readlines()
		for inline in inlines:
			
			if option.gffkde == "ON":
				chrm, runMethod, className, start, end, density, unknown, unspecified, details = inline.strip().split("\t")
			else:
				chrm, start, end, feature, occupancy, strand, density, datasetCount, factorCount, contextCount, details = inline.strip().split("\t")
		
			# parse details:
			identifiers, datasets, factors, contexts = metrn.parseDetails(details, target=option.target, gffkde=option.gffkde, contributionCutoff=option.contribution)
			
			# scan details:
			catch = False
			if option.limits != "OFF":
				for limit in option.limits.split(","):
					text, limit = limit.split(":")
					if details.count(text) > int(limit):
						catch = True
			
			# count total hits (density) and unique hits (complexity):
			size = int(end) - int(start) + 1
			occupancy = len(identifiers)
			density = 1000*float(occupancy)/size
			datasetCounts = len(set(datasets))
			factorCounts = len(set(factors))
			contextCounts = len(set(contexts))
			
			# set analysis metric and flag:
			if option.metric == "occupancy":
				value = occupancy
				if value > float(option.cutoff):
					catch = True
			elif option.metric == "density":
				value = density
				if value > float(option.cutoff) and occupancy >= option.minOccupancy:
					catch = True
			elif option.metric == "complexity":
				value = complexity
				if value > float(option.cutoff) and occupancy >= option.minOccupancy:
					catch = True
			elif option.metric == "combined":
				maxOccupancy, maxDensity = map(float, option.cutoff.split(","))
				if occupancy > maxOccupancy and density > maxDensity:
					catch = True
		
			# apply cutoffs:
			if catch:
				print >>f_output, "\t".join([chrm, start, end, "HOT." + str(fi), str(occupancy), "+", str(density), str(datasetCounts), str(factorCounts), str(contextCounts), details])
				fi += 1
			else:
				print >>p_output, "\t".join([chrm, start, end, "TFBS." + option.tag + str(pi), str(occupancy), "+", str(density), str(datasetCounts), str(factorCounts), str(contextCounts), details])
				pi += 1
		
		# close output file:
		f_output.close()
		
		print
		print "Regions analyzed:", pi + fi
		print "Regions passed:", pi, "(" + str(round(100*float(pi)/(pi + fi), 2)) + "%)"
		print "Regions failed:", fi, "(" + str(round(100*float(fi)/(pi + fi), 2)) + "%)"
		print

	
	# perform comparative analysis with other region files (such as the Yip 2012 HOT regions):
	if option.mode == "compare":
	
		# define query and target files:
		queryfile = analysispath + option.infile
		targetfile = option.source
		
		# define output files:
		abOutfile = comparepath + "maphot_compare_" + option.name + "_AB_collapse.bed"
		baOutfile = comparepath + "maphot_compare_" + option.name + "_BA_collapse.bed"
		
		print
		print "Query:", option.infile
		print "Target:", option.source.split("/")[len(option.source.split("/"))-1]
		
		# determine overlaps:
		command = "intersectBed -a " + queryfile + " -b " + targetfile + " -u >" + abOutfile
		os.system(command)

		command = "intersectBed -a " + targetfile + " -b " + queryfile + " -u >" + baOutfile
		os.system(command)
		
		# measure comparison:
		qyCount = general.countLines(queryfile)
		tgCount = general.countLines(targetfile)
		abCount = general.countLines(abOutfile)
		baCount = general.countLines(baOutfile)
		
		print
		print "Queries matched:", abCount, "(" + str(round(100*float(abCount)/(qyCount), 2)) + "%)"
		print "Targets matched:", baCount, "(" + str(round(100*float(baCount)/(tgCount), 2)) + "%)"
		print
		
	# perform overlap analysis amongst HOT region files:
	if option.mode == "regions":
		
		if option.source == "analysis":
		
			# define output files:
			hotfile = regionspath + "maphot_" + option.regions + "_hot.bed"
			rgbfile = regionspath + "maphot_" + option.regions + "_rgb.bed"
		
			print
			print "Finding HOT-region and RGB-region files..."
			for infile in sorted(os.listdir(analysispath)):
				if option.name in infile and "fail.bed" in infile:
					print "HOT regions:", infile
					command = "cp " + analysispath + infile + " " + hotfile
					os.system(command)
				if option.name in infile and "pass.bed" in infile:
					print "RGB regions:", infile
					command = "cp " + analysispath + infile + " " + rgbfile
					os.system(command)
			print "Copying files to region folders..."
			print
		
		if option.source == "overlap":
		
			# define output files:
			anyfile = regionspath + "maphot_" + option.regions + "_any.bed"
			allfile = regionspath + "maphot_" + option.regions + "_all.bed"
			
			print
			print "Finding overlap HOT-region files..."
			for infile in sorted(os.listdir(overlappath)):
				if option.name in infile and "_region_global.bed" in infile:
					print "Global HOT regions (context-dependent):", infile
					command = "cp " + overlappath + infile + " " + anyfile
					os.system(command)
				if option.name in infile and "_region_shared.bed" in infile:
					print "Shared HOT regions (ubiquitously HOT):", infile
					command = "cp " + overlappath + infile + " " + allfile
					os.system(command)
			print "Copying files to region folders..."
			print
						
	# perform overlap analysis amongst HOT region files:
	if option.mode == "overlap":
		
		# define output files:
		#j_outfile = "maphot_" + option.mode + "_" + option.overlap + "_" + option.name + "_report_chords.txt"
		f_outfile = "maphot_" + option.mode + "_" + option.overlap + "_" + option.name + "_report_matrix.txt"
		s_outfile = "maphot_" + option.mode + "_" + option.overlap + "_" + option.name + "_region_shared.bed"
		b_outfile = "maphot_" + option.mode + "_" + option.overlap + "_" + option.name + "_region_subset.bed"
		m_outfile = "maphot_" + option.mode + "_" + option.overlap + "_" + option.name + "_merged_subset.bed"
		c_outfile = "maphot_" + option.mode + "_" + option.overlap + "_" + option.name + "_region_global.bed"
		x_outfile = "maphot_" + option.mode + "_" + option.overlap + "_" + option.name + "_merged_global.bed"
		u_outfile = "maphot_" + option.mode + "_" + option.overlap + "_" + option.name + "_region_unique.bed"
		
		# prepare output file:
		f_output = open(overlappath + f_outfile, "w")
		print >>f_output, "\t".join(["region.a", "region.b", "count.a", "count.b", "overlap", "overlap.max", "overlap.sum"])
		
		# prepare key:target dictionary:
		key_dict = dict()
		for target in option.target.split(","):
			key, target = target.split(":")
			key_dict[key] = target
		
		print
		print "Print associating keys to HOT-region files..."
		hotfile_dict = dict()
		for hotfile in sorted(os.listdir(analysispath)):
			if option.name in hotfile and "fail.bed" in hotfile:
				for key in sorted(key_dict.keys()):
					if key_dict[key] in hotfile:
						hotfile_dict[key] = hotfile
						print key, ":", hotfile
		
		print
		print "Preparing HOT-overlap matrix..."
		matrix_dict, chords_dict = dict(), dict()
		for keyA in sorted(hotfile_dict.keys()):
			targetA = key_dict[keyA]
			if not keyA in matrix_dict:
				matrix_dict[keyA] = dict()
			for keyB in sorted(hotfile_dict.keys()):
				targetB = key_dict[keyB]
				
				command = 'grep -v "feature" ' + analysispath + hotfile_dict[keyA] + ' > ' + overlappath + "maphot_targetA.tmp"
				os.system(command)
		
				command = 'grep -v "feature" ' + analysispath + hotfile_dict[keyB] + ' > ' + overlappath + "maphot_targetB.tmp"
				os.system(command)
		
				command = "intersectBed -a " + overlappath + "maphot_targetA.tmp" + " -b " + overlappath + "maphot_targetB.tmp" + " -u > " + overlappath + "maphot_intersect.tmp"
				os.system(command)
				
				regions_a = general.countLines(overlappath + "maphot_targetA.tmp")
				regions_b = general.countLines(overlappath + "maphot_targetB.tmp")
				regions_o = general.countLines(overlappath + "maphot_intersect.tmp")
				
				overlap = regions_o
				union = regions_o + regions_a - regions_o + regions_b - regions_o
				if union == 0:
					overlap_max =0
					overlap_sum =0
				else:
					overlap_max = float(overlap)/max(regions_a, regions_b)
					overlap_sum = float(overlap)/union
				
				output = [keyA, keyB, regions_a, regions_b, overlap, overlap_max, overlap_sum]
				print >>f_output, "\t".join(map(str, output))
				
				# store overlap for chords-graphing:
				if not keyA in chords_dict:
					chords_dict[keyA] = dict()
				chords_dict[keyA][keyB] = regions_o
		
		# find HOT regions that are shared amongst all datasets:
		print "Finding shared HOT regions..."
		print
		processed = list()
		for key in sorted(hotfile_dict.keys()):
			
			command = 'grep -v "feature" ' + analysispath + hotfile_dict[key] + " > " + overlappath + "maphot_source.tmp"
			os.system(command)
		
			# create or update the shared HOT regions:
			if processed == list():
				command = "cp " + overlappath + "maphot_source.tmp" + " " + overlappath +  "maphot_shared.tmp"
				os.system(command)
			else:
				command = "intersectBed -a " + overlappath + s_outfile + " -b " + overlappath + "maphot_source.tmp" + " -u > " + overlappath + "maphot_shared.tmp"
				os.system(command)
			
			command = "cp " + overlappath + "maphot_shared.tmp" + " " + overlappath +  s_outfile
			os.system(command)
			
			# index/re-name HOT regions with the respective key:
			indexColumns(overlappath + "maphot_source.tmp", column=4, base="HOT-" + key + ".", header=False, extension="")
			
			print "Peaks in", key,"HOT regions:", densityPeakCounter(overlappath + "maphot_source.tmp")
			
			# aggregate context HOT regions:
			if processed == list():
				command = "cp " + overlappath + "maphot_source.tmp" + " " + overlappath + "maphot_aggregate.tmp"
				os.system(command)
			else:
				command = "cp " + overlappath + "maphot_aggregate.tmp" + " " + overlappath + "maphot_current.tmp"
				os.system(command)
				
				command = "cat " + overlappath + "maphot_current.tmp" + " " + overlappath + "maphot_source.tmp" + " > " + overlappath + "maphot_aggregate.tmp"
				os.system(command)
			
			# update processed keys:
			processed.append(key)
		
		# export shared HOT regions (also re-name HOT regions):
		indexColumns(overlappath + "maphot_shared.tmp", column=4, base="HOT.", header=False, extension="")
		command = "cp " + overlappath + "maphot_shared.tmp" + " " + overlappath + s_outfile
		os.system(command)
		
		# generate and export chords .json file:
		#sharedRegions = general.countLines(overlappath + "maphot_shared.tmp")
		#chords_dict["uHOT"] = dict()
		#chords_dict["uHOT"]["uHOT"] = 0
		#for keyX in chords_dict:
		#	chords_dict[keyX]["uHOT"] = sharedRegions
		#	chords_dict["uHOT"][keyX] = sharedRegions
		
		# export json file:
		#chords_list = list()
		#for keyA in chords_dict:
		#	keyA_list = list()
		#	for keyB in chords_dict:
		#		keyA_list.append(chords_dict[keyA][keyB])
		#	chords_list.append(keyA_list)
		#import simplejson as json
		#j_output = open(overlappath + j_outfile, "w")
		#json.dump(chords_list, j_output)
		#j_output.close()
		
		print
		print "Peaks in shared HOT regions:", densityPeakCounter(overlappath + s_outfile)
		
		# collapse context-wide subset HOT regions (not-shared):
		command = "mergeBed -i " + overlappath + "maphot_aggregate.tmp" + " -nms > " + overlappath + x_outfile
		os.system(command)
		
		# export aggregate context HOT regions:
		#indexColumns(overlappath + "maphot_aggregate.tmp", column=4, base="HOT.", header=False, extension="")
		command = "cp " + overlappath + "maphot_aggregate.tmp" + " " + overlappath + c_outfile
		os.system(command)
		
		
		print "Peaks across all HOT regions:", densityPeakCounter(overlappath + c_outfile)
		
		
		# find HOT regions that are only in a subset of the datasets (not-shared):
		print "Finding non-shared (subset) HOT regions..."
		processed = list()
		for key in sorted(hotfile_dict.keys()):
			
			# create dataset non-shared (subset) HOT regions:
			command = 'grep -v "feature" ' + analysispath + hotfile_dict[key] + " > " + overlappath + "maphot_source.tmp"
			os.system(command)
		
			command = "intersectBed -a " + overlappath + "maphot_source.tmp" + " -b " + overlappath + "maphot_shared.tmp" + " -v > " + overlappath + "maphot_subset.tmp"
			os.system(command)
			
			# export context-specific subset HOT regions (not-shared):
			k_outfile = b_outfile.replace("_cx", "_" + key.lower())
			command = "cp " + overlappath + "maphot_subset.tmp" + " " + overlappath + k_outfile
			os.system(command)
			
			# index/re-name HOT regions with the respective key:
			indexColumns(overlappath + "maphot_subset.tmp", column=4, base="HOT-" + key + ".", header=False, extension="")
			
			# aggregate subset HOT regions (for context-wide analysis):
			if processed == list():
				command = "cp " + overlappath + "maphot_subset.tmp" + " " + overlappath + "maphot_aggregate.tmp"
				os.system(command)
			else:
				command = "cp " + overlappath + "maphot_aggregate.tmp" + " " + overlappath + "maphot_current.tmp"
				os.system(command)
				
				command = "cat " + overlappath + "maphot_current.tmp" + " " + overlappath + "maphot_subset.tmp" + " > " + overlappath + "maphot_aggregate.tmp"
				os.system(command)
			
			# update processed keys:
			processed.append(key)
		
		# collapse context-wide subset HOT regions (not-shared):
		command = "mergeBed -i " + overlappath + "maphot_aggregate.tmp" + " -nms > " + overlappath + m_outfile
		os.system(command)
		
		# export context-wide subset HOT regions (not-shared); employ re-naming:
		#indexColumns(overlappath + "maphot_aggregate.tmp", column=4, base="HOT.", header=False, extension="")
		command = "cp " + overlappath + "maphot_aggregate.tmp" + " " + overlappath + b_outfile
		os.system(command)
		
		print "Peaks in subset HOT regions:", densityPeakCounter(overlappath + b_outfile)
		
		# find HOT regions that are unique to individual datasets:
		print "Finding dataset-specific (unique) HOT regions..."
		processed = list()
		for keyA in sorted(hotfile_dict.keys()):
			
			# prepare unique HOT regions for dataset A:
			command = 'grep -v "feature" ' + analysispath + hotfile_dict[keyA] + " > " + overlappath + "maphot_unique.tmp"
			os.system(command)
			
			# filter unique HOT regions for dataset A:
			for keyB in sorted(hotfile_dict.keys()):
				if keyA != keyB:
				
					# prepare HOT regions for dataset B:
					command = 'grep -v "feature" ' + analysispath + hotfile_dict[keyB] + " > " + overlappath + "maphot_subset.tmp"
					os.system(command)
					
					command = "cp " + overlappath + "maphot_unique.tmp" + " " + overlappath + "maphot_current.tmp"
					os.system(command)
					
					command = "intersectBed -a " + overlappath + "maphot_current.tmp" + " -b " + overlappath + "maphot_subset.tmp" + " -v > " + overlappath + "maphot_unique.tmp"
					os.system(command)
			
			# export unique HOT regions for dataset A; employ re-naming:
			k_outfile = u_outfile.replace("_cx", "_" + keyA.lower())
			indexColumns(overlappath + "maphot_unique.tmp", column=4, base="HOT.", header=False, extension="")
			command = "cp " + overlappath + "maphot_unique.tmp" + " " + overlappath + k_outfile
			os.system(command)
			
			# index/re-name HOT regions with the respective key:
			indexColumns(overlappath + "maphot_unique.tmp", column=4, base="HOT-" + keyA + ".", header=False, extension="")
			
			# aggregate unique HOT regions:
			if processed == list():
				command = "cp " + overlappath + "maphot_unique.tmp" + " " + overlappath + "maphot_aggregate.tmp"
				os.system(command)
			else:
				command = "cp " + overlappath + "maphot_aggregate.tmp" + " " + overlappath + "maphot_current.tmp"
				os.system(command)
				
				command = "cat " + overlappath + "maphot_current.tmp" + " " + overlappath + "maphot_unique.tmp" + " > " + overlappath + "maphot_aggregate.tmp"
				os.system(command)
			
			processed.append(key)
		
		# export aggregate unique HOT regions:
		#indexColumns(overlappath + "maphot_aggregate.tmp", column=4, base="HOT.", header=False, extension="")
		command = "cp " + overlappath + "maphot_aggregate.tmp" + " " + overlappath + u_outfile
		os.system(command)
		
		print "Peaks in unique HOT regions:", densityPeakCounter(overlappath + u_outfile)
		print	
		
		# remove temporary files:
		command = "rm -rf " + overlappath + "*tmp"
		os.system(command)
		

	# filter or select input HOT regions that overlap with a given set of HOT regions:
	if "filter" in option.mode:

		# define filter/select operation flag
		if option.mode == "filter:remove":
			mode_flag = option.mode.replace(":","_")
			operation_flag = " -v"
		elif option.mode == "filter:select":
			mode_flag = option.mode.replace(":","_")
			operation_flag = " -u"
		elif option.mode == "filter:invert":
			mode_flag = option.mode.replace(":","_")
			operation_flag = " -u"
		
		# filter peak files...
		if option.peaks != "OFF":
		
			# load peak files:
			option.source = str(option.source + "/").replace("//","/")
			infiles = os.listdir(option.source)
			
			# update peaks path:
			peakspath = peakspath + option.peaks + "/"
			general.pathGenerator(peakspath)
			
			# filter peaks:
			for infile in infiles:
				
				if option.mode == "filter:remove":
					command = "intersectBed -a " + option.source + infile + " -b " + regionspath + option.infile + operation_flag + " > " + peakspath + infile
					os.system(command)
			
				elif option.mode == "filter:select":
					command = "intersectBed -a " + option.source + infile + " -b " + regionspath + option.infile + operation_flag + " > " + peakspath + infile
					os.system(command)
			
		
		# filter HOT region files...
		if option.peaks == "OFF":
		
			# update source path:
			sourcepath = hotpath + option.source + "/"
			
			# designate output file:
			f_outfile = "maphot_" + mode_flag + "_" + option.name + "_" + option.infile.replace("maphot_","")
			
			# filter out HOT regions that overlap :
			print
			print "Filtering HOT regions that overlap..."
			command = 'grep -v "feature" ' + analysispath + option.infile + " > " + analysispath + "maphot_infile.tmp"
			os.system(command)
			
			command = 'grep -v "feature" ' + sourcepath + option.target + " > " + analysispath + "maphot_target.tmp"
			os.system(command)
			
			if option.mode == "filter:remove":
				command = "intersectBed -a " + analysispath + "maphot_infile.tmp" + " -b " + analysispath + "maphot_target.tmp" + operation_flag + " > " + analysispath + "maphot_filtered.tmp"
				os.system(command)
			
			elif option.mode == "filter:select":
				command = "intersectBed -a " + analysispath + "maphot_infile.tmp" + " -b " + analysispath + "maphot_target.tmp" + operation_flag + " > " + analysispath + "maphot_filtered.tmp"
				os.system(command)
			
			elif option.mode == "filter:invert":
				command = "intersectBed -a " + analysispath + "maphot_target.tmp" + " -b " + analysispath + "maphot_infile.tmp" + operation_flag + " > " + analysispath + "maphot_filtered.tmp"
				os.system(command)
			
			# export filtered HOT regions:
			print "Exporting filtered HOT regions..."
			indexColumns(analysispath + "maphot_filtered.tmp", column=4, base="HOT.", header=False, extension="")
			
			# remove temporary files:
			command = "rm -rf " + analysispath + "*tmp"
			os.system(command)
			print
			
	
	# measure temperature of regions across different contexts:
	if option.mode == "temperature":
		
		# determine input density file:
		if option.gffkde == "ON":
			densityfile = "gffkde2_" + option.peaks + gffkde_handle + ".Peaks.Extended"
		else:
			densityfile = "mergeBed_" + option.peaks + "_compiled.bed"
		g_infile = densitypath + densityfile
		
		# set analysis metric and flag:
		if option.metric == "occupancy":
			metric_flag = "_occupan"
		elif option.metric == "density":
			metric_flag = "_density"
		elif option.metric == "complexity":
			metric_flag = "_complex"
		elif option.metric == "combined":
			metric_flag = "_combine"
			
		# update source path:
		sourcepath = hotpath + option.source + "/"
		
		# load target contexts:
		codeContexts, targetContexts = metrn.options_dict["contexts.extended"][option.contexts]
		
		# designate and open output file:
		f_outfile = "maphot_" + option.mode + "_" + option.peaks + "_" + codeContexts + "_" + option.infile.replace("maphot_","").replace(".bed","") + "_counts"
		n_outfile = "maphot_" + option.mode + "_" + option.peaks + "_" + codeContexts + "_" + option.infile.replace("maphot_","").replace(".bed","") + "_normal"
		a_outfile = "maphot_" + option.mode + "_" + option.peaks + "_" + codeContexts + "_" + option.infile.replace("maphot_","").replace(".bed","") + "_average"
		f_output = open(temperaturepath + f_outfile, "w")
		n_output = open(temperaturepath + n_outfile, "w")
		a_output = open(temperaturepath + a_outfile, "w")
		
		# collapse contexts and export header:
		header = ["region"]
		collapsedContexts = list()
		for context in targetContexts:
			if context in ["EE", "EM", "LE"]:
				if not "EX" in collapsedContexts:
					collapsedContexts.append("EX")
			else:
				collapsedContexts.append(context)
		header += collapsedContexts 
		print >>f_output, "\t".join(header)
		print >>n_output, "\t".join(header)
		print >>a_output, "\t".join(["context", "mean", "st.dev", "percentile.75", "percentile.25"])
		
		print
		print "Targeting contexts (total):", ",".join(targetContexts)
		print "Targeting contexts (collapsed):", ",".join(collapsedContexts)
		
		# generate a copy of the target HOT regions:
		t_infile = temperaturepath + "maphot_targets.tmp"
		command = "cp " + g_infile + " " + t_infile
		os.system(command)
		
		# scan HOT-regions!
		print
		print "Processing HOT-region overlaps..."
		hot_dict, k = dict(), 1
		for inline in open(sourcepath + option.infile).readlines():
			
			# load HOT-region information (individual region):
			chrm, start, end, region, occupancy, strand, density, complexity, factorCount, contextCount, infos = inline.strip().split("\t")
			k += 1
			
			# make current region file:
			r_infile = temperaturepath + "maphot_region.tmp"
			r_outfile = open(r_infile, "w")
			print >>r_outfile, inline.strip()
			r_outfile.close()
			
			# determine overlap between query (genomic regions) and targeted HOT-region (i.e. extract overlap from background regions):
			#command = "cat " + g_infile  + "| awk 'BEGIN{OFS=\"\\t\"; FS=\"\\t\"} {print $1, $4, $5, $3, $6, \"+\", $9}' > " + temperaturepath + "maphot_tinfile.tmp" 
			#os.system(command)
			
			#command = "cat " + r_infile  + "| awk 'BEGIN{OFS=\"\\t\"; FS=\"\\t\"} {print $1, $2, $3, $4, $5, $6, $7}' > " + temperaturepath + "maphot_rinfile.tmp" 
			#os.system(command)
			
			#t_infile = temperaturepath + "maphot_tinfile.tmp"
			#r_infile = temperaturepath + "maphot_rinfile.tmp"
			
			command = "intersectBed -a " + t_infile + " -b " + r_infile + " -u > " + temperaturepath + "maphot_regions.tmp"
			os.system(command)
			
			# scan density in this region across contexts and load into dictionary:
			hot_dict[region] = dict()
			inlines = open(temperaturepath + "maphot_regions.tmp").readlines()
			for inline in inlines:
				items = inline.strip().split("\t")
				infos = items[10].split(";")
				infos.pop(0)
				for info in infos:
					if info != "":
						dataset, peak = info.split(":")
						organism, strain, factor, context, institute = dataset.split("_")
						for targetContext in targetContexts:
							if not targetContext in hot_dict[region]:
								hot_dict[region][targetContext] = 0
							if context == targetContext:
								hot_dict[region][targetContext] +=1
		
		
		# export counts and normalized data:
		#print
		print "Scoring HOT-region overlaps..."
		average_dict = dict()
		for inline in open(sourcepath + option.infile).readlines():
			
			# load HOT-region information (individual region):
			chrm, start, end, region, occupancy, strand, density, complexity, factorCount, contextCount, infos = inline.strip().split("\t")
			
			# prepare matrix entry and store into average dictionary:
			counts = [region]
			for context in collapsedContexts:
				if not context in average_dict:
					average_dict[context] = list()
				hits = 0
				if context == "EX":
					hitContexts = ["EE", "EM", "LE"]
				else:
					hitContexts = [context]
				for hitContext in hitContexts:
					if hitContext in hot_dict[region]:
						hits += hot_dict[region][hitContext]
				average_dict[context].append(hits)
				counts.append(hits)
			
			# normalize the counts:
			normal = [region]
			for count in counts[1:]:
				if counts[1] != 0:
					normal.append(float(count)/counts[1])
				else:
					normal.append(0)
			
			# export normalized and unnormalized counts:
			print >>f_output, "\t".join(map(str, counts))
			print >>n_output, "\t".join(map(str, normal))
		
		# export averages :
		print "Exporting average overlaps..."
		for context in collapsedContexts:
			print >>a_output, "\t".join([context, str(numpy.mean(average_dict[context])), str(numpy.std(average_dict[context])), str(scipy.percentile(average_dict[context], 75)), str(scipy.percentile(average_dict[context], 25))])
		print
		
		# close outputs:
		f_output.close()
		n_output.close()
		a_output.close()
		
		# remove temporary files:
		command = "rm -rf " + temperaturepath + "*tmp"
		os.system(command)
				
					
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

#montecarlo
#simulation
#collecting
#downloader