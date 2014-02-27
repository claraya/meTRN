#!/usr/bin/env python
# execute Co-association analysis with IntervalStats!

import sys
import time
import optparse
import general
import numpy
import pickle
import pdb
import metrn
import modencode
import multiprocessing
import random
import os

from runner import *

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())


""" define a function that converts a local user path to SCG3 or GS server paths """
def serverPath(inpath, server="ON"):
	if server == "ON":
		return inpath.replace("/Users/claraya/", "/srv/gs1/projects/snyder/claraya/")
	elif server == "GS":
		return inpath.replace("/Users/claraya/", "/net/fields/vol1/home/araya/")

def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "Path from script to files")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "Operations to perform")
	parser.add_option("--peakflag", action = "store", type = "string", dest = "peakflag", help = "Peak type: raw, xot, hot", default="raw")
	parser.add_option("--technique", action = "store", type = "string", dest = "technique", help = "Matrix type: binary or signal or other", default="binary")
	parser.add_option("--iterations", action = "store", type = "int", dest = "iterations", help = "Iterations to perform", default=100)
	parser.add_option("--input", action = "store", type = "string", dest = "input", help = "Inputs and names for multiple sample processing (mode: sample)", default="OFF")
	parser.add_option("--threads", action = "store", type = "string", dest = "threads", help = "Multiprocessing threads", default="1")
	parser.add_option("--chunks", action = "store", type = "int", dest = "chunks", help = "In how many chunks should these be processed?", default=1)
	parser.add_option("--qsub", action = "store", type = "string", dest = "qsub", help = "Qsub configuration header", default="OFF")
	parser.add_option("--server", action = "store", type = "string", dest = "server", help = "Are we on the server?", default="OFF")
	parser.add_option("--qsubreplace", action = "store", type = "string", dest = "qsubreplace", help = "Target replacement for execution code", default="OFF")
	parser.add_option("--job", action = "store", type = "string", dest = "job", help = "Job name for cluster", default="OFF")
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
	bindingpath = path_dict["binding"]
	coassociationspath = path_dict["coassociations"]
	neuronspath = path_dict["neurons"]
	cellspath = path_dict["cells"]
	
	# standardize paths for analysis:
	alignerpath = bwapath
	indexpath = alignerpath + "index/"
	alignmentpath = alignerpath + "alignment/"
	qcfilterpath = alignerpath + "qcfilter/"
	qcmergepath = alignerpath + "qcmerge/"
	
	# launch SOM iterations:
	if option.mode == "launch":
		
		# prepare for qsub:
		bash_path = str(neuronspath + "runs/").replace("//","/")
		bash_base = "_".join([option.peakflag, option.technique]) + "-M"
		qsub_base = "_".join([option.peakflag, option.technique])
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
		
		# generate seeds (random):
		population = range(1, 5000)
		seeds = random.sample(population, option.iterations)
			
		# prepare R commands:
		modules, commands, sequence = list(), list(), list()
		for seed in seeds:
			
			# generate R command:
			#Rscript mapSOM_codebase-matrix.r ~/meTRN/data/neurons/ raw binary 272
			command = "Rscript <<RSCRIPT>> <<NEURONSPATH>> <<PEAKFLAG>> <<TECHNIQUE>> <<SEED>>"
			command = command.replace("<<RSCRIPT>>", scriptspath + "mapSOM_codebase-matrix.r")
			command = command.replace("<<NEURONSPATH>>", neuronspath)
			command = command.replace("<<PEAKFLAG>>", option.peakflag)
			command = command.replace("<<TECHNIQUE>>", option.technique)
			command = command.replace("<<SEED>>", str(seed))
			commands.append(command)
			sequence.append(command)
			
			if len(commands) == option.chunks:
				modules.append(commands)
				commands = list()
		
		if commands != list():
			modules.append(commands)
		
		# launch R commands:
		print
		print "Launching comparisons:", len(modules)
		k = 0
		for module in modules:
			for command in module:
				k += 1
				#print command
				#print
		runCommands(modules, threads=option.threads, mode="module.run", run_mode="verbose", run_path=bash_path, run_base=bash_base, record=True, qsub_header=qsub_header, qsub=qsub, qsub_replace=option.qsubreplace)
		print "Comparisons performed:", len(modules), "(" + str(k) + " commands)"
		print		

	
	# recover SOM iterations:
	elif option.mode == "recover":
	
		print
		print "Recovering quantization errors from seeds..."
		matrix = dict()
		for peakset in os.listdir(neuronspath):
			if not ".py" in peakset and not peakset == "runs":
				for technique in os.listdir(neuronspath + peakset):
					for analysis in os.listdir(neuronspath + peakset + "/" + technique + "/reports/"):
						for seed in os.listdir(neuronspath + peakset + "/" + technique + "/reports/" + analysis):
							error = open(neuronspath + peakset + "/" + technique + "/reports/" + analysis + "/" + seed).read().strip()
							if not analysis in matrix:
								matrix[analysis] = dict()
							seed = seed.replace("seed.", "")
							matrix[analysis][seed] = float(error)
		
		"Exporting recommended seeds for each analysis..."
		for analysis in sorted(matrix.keys()):
			seeds = general.valuesort(matrix[analysis])
			print analysis, ":", seeds[0], "(" + str(matrix[analysis][seeds[0]]) + ")"
		print
	
	
	# examine SOM results:
	elif option.mode == "examine":
	
		print
		print "Examining SOM neuron contents..."
		matrix = dict()
		for peakset in os.listdir(neuronspath):
			if not ".py" in peakset and not peakset == "runs":
				for technique in os.listdir(neuronspath + peakset):
					for analysis in os.listdir(neuronspath + peakset + "/" + technique + "/results/"):
						if "retired" != analysis:
							counter, total = 0, 0
							for neuron in os.listdir(neuronspath + peakset + "/" + technique + "/results/" + analysis + "/regions/raw/"):
								regions = general.countLines(neuronspath + peakset + "/" + technique + "/results/" + analysis + "/regions/raw/" + neuron)
								total += 1
								if regions > 0:
									counter += 1
							print analysis, counter, "(" + str(total) + ")"
		print
	
	
	# launch SOM emulations for subset-modes:
	if option.mode in ["sample", "filter", "remove"]:
		
		# prepare for qsub:
		bash_path = str(neuronspath + "runs/").replace("//","/")
		bash_base = "_".join([option.peakflag, option.technique]) + "-M"
		qsub_base = "_".join([option.peakflag, option.technique])
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
		
		# generate seeds (random):
		population = range(1, 5000)
		seeds = random.sample(population, option.iterations)
		
		# generate index tag:
		targets = list()
		for index in range(1, int(option.iterations) + 1):
			indexTag = general.indexTag(index, int(option.iterations), minimum=3)
			targets.append(indexTag)
			
		# prepare R commands:
		modules, commands, sequence = list(), list(), list()
		for specification in option.input.split(","):
			for target in targets:
				
				# load specifications:
				inPut, name = specification.split(":")
			
				# define script path:
				if option.mode == "sample":
					script = "mapSOM_codebase-sample.r"
				elif option.mode == "filter":
					script = "mapSOM_codebase-filter.r"
				elif option.mode == "remove":
					script = "mapSOM_codebase-remove.r"
				
				# re-define name:
				name = name.replace("XXX", target)
				
				# generate R command:
				#Rscript mapSOM_codebase-runs.r ~/meTRN/data/neurons/ raw binary 272
				command = "Rscript <<RSCRIPT>> <<NEURONSPATH>> <<PEAKFLAG>> <<TECHNIQUE>> <<SEED>> <<TARGET>> <<INPUT>> <<NAME>>"
				command = command.replace("<<RSCRIPT>>", scriptspath + script)
				command = command.replace("<<NEURONSPATH>>", neuronspath)
				command = command.replace("<<PEAKFLAG>>", option.peakflag)
				command = command.replace("<<TECHNIQUE>>", option.technique)
				command = command.replace("<<SEED>>", str(100))
				command = command.replace("<<TARGET>>", target)
				command = command.replace("<<INPUT>>", inPut)
				command = command.replace("<<NAME>>", name)
				commands.append(command)
				sequence.append(command)
				
				if len(commands) == option.chunks:
					modules.append(commands)
					commands = list()
		
		if commands != list():
			modules.append(commands)
		
		# launch R commands:
		print
		print "Launching comparisons:", len(modules)
		k = 0
		for module in modules:
			for command in module:
				k += 1
				#print command
				#print
		runCommands(modules, threads=option.threads, mode="module.run", run_mode="verbose", run_path=bash_path, run_base=bash_base, record=True, qsub_header=qsub_header, qsub=qsub, qsub_replace=option.qsubreplace)
		print "Comparisons performed:", len(modules), "(" + str(k) + " commands)"
		print
		
	
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

