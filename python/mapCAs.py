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
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "path from script to files")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "launch")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "basename for target peaks", default="OFF")
	parser.add_option("--source", action = "store", type = "string", dest = "source", help = "Path of peaks to be filtered or other files...")
	parser.add_option("--domain", action = "store", type = "string", dest = "domain", help = "domain regions BED file", default="OFF")
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "name for domain", default="OFF")
	parser.add_option("--label", action = "store", type = "string", dest = "label", help = "Type of labels to generate", default="factor.context")
	parser.add_option("--cells", action = "store", type = "string", dest = "cells", help = "Use cellular-resolution peaks?", default="OFF")
	parser.add_option("--cutoff", action = "store", type = "float", dest = "cutoff", help = "p-value cutoff for fraction calculation", default=0.05)
	parser.add_option("--threads", action = "store", type = "string", dest = "threads", help = "multiprocessing threads", default="1")
	parser.add_option("--chunks", action = "store", type = "int", dest = "chunks", help = "", default=100)
	parser.add_option("--include", action = "store", type = "string", dest = "include", help = "Targets to include", default="OFF")
	parser.add_option("--exclude", action = "store", type = "string", dest = "exclude", help = "Targets to exclude", default="OFF")
	parser.add_option("--qsub", action = "store", type = "string", dest = "qsub", help = "qsub configuration header", default="OFF")
	parser.add_option("--server", action = "store", type = "string", dest = "server", help = "are we on the server?", default="OFF")
	parser.add_option("--job", action = "store", type = "string", dest = "job", help = "job name for cluster", default="OFF")
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
	orthologspath = path_dict["orthologs"]
	coassociationspath = path_dict["coassociations"]
	neuronspath = path_dict["neurons"]
	cellspath = path_dict["cells"]
	
	# standardize paths for analysis:
	alignerpath = bwapath
	indexpath = alignerpath + "index/"
	alignmentpath = alignerpath + "alignment/"
	qcfilterpath = alignerpath + "qcfilter/"
	qcmergepath = alignerpath + "qcmerge/"
	
	# define P-value cutoff handle:
	pvaluecutoff_handle = "%.0e" % (float(option.cutoff))
	
	# use cellular resolution peaks?
	if option.cells == "ON":
		peakspath = cellspath + "peaks/"
		annotationspath = cellspath + "annotations/"
	
	# launch IntervalStats:
	if option.mode == "launch":
		
		# prepare for qsub:
		bash_path = str(coassociationspath + "runs/").replace("//","/")
		bash_base = "_".join([option.mode, option.peaks, option.name]) + "-M"
		qsub_base = "_".join([option.mode, option.peaks, option.name])
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
				
		# create output paths:
		resultspath = coassociationspath + option.peaks + "/" + option.name + "/results"
		summarypath = coassociationspath + option.peaks + "/" + option.name + "/summary"
		general.pathGenerator(resultspath)
		general.pathGenerator(summarypath)
			
		# modify domain if necessary (for cellular-resolution analyses):
		domainpath = path_dict[option.source]
		if option.cells == "ON":
			peakspath = peakspath.replace("data/peaks/", "data/cells/peaks/")	
			domainpath = domainpath.replace("data/annotations/", "data/cells/annotations/")	
	
		# load peak files:
		peak_files = os.listdir(peakspath + option.peaks)
		
		# prepare IntervalStats commands:
		modules, commands, sequences = list(), list(), list()
		for PK1 in peak_files:
			for PK2 in peak_files:
				if PK1 != PK2:
				
					# check inclusion flags:
					process = False
					if option.include == "OFF":
						process = True
					else:
						inclusionFlag1, inclusionFlag2 = False, False
						for inclusion in option.include.split(","):
							if inclusion in PK1:
								inclusionFlag1 = True
							if inclusion in PK2:
								inclusionFlag2 = True
						if inclusionFlag1 and inclusionFlag2:
							process = True
					
					if process:
					
						# define output files:
						outfile = "mapcas_" + option.name + "_" + PK1.replace(".bed","") + "_vs_" + PK2
						outpath = coassociationspath + option.peaks + "/" + option.name + "/results/"
					
						# generate IntervalStats command:
						command = "IntervalStats -q <<FILE1>> -r <<FILE2>> -d <<DOMAINS>> -o <<OUTPUT>>"
						command = command.replace("<<FILE1>>", peakspath + option.peaks + "/" + PK1)
						command = command.replace("<<FILE2>>", peakspath + option.peaks + "/" + PK2)
						command = command.replace("<<DOMAINS>>", domainpath + option.domain)
						command = command.replace("<<OUTPUT>>", outpath + outfile)
						commands.append(command)
						sequences.append(command)
					
						if len(commands) == option.chunks:
							modules.append(commands)
							commands = list()
		
		if commands != list():
			modules.append(commands)
		
		# launch IntervalStats commands:
		print
		print "Launching comparisons:", len(modules)
		k = 0
		for module in modules:
			for command in module:
				k += 1
				#print command
				#print
		runCommands(modules, threads=option.threads, mode="module.run", run_mode="verbose", run_path=bash_path, run_base=bash_base, record=True, qsub_header=qsub_header, qsub=qsub)
		print "Comparisons performed:", len(modules), "(" + str(k) + " commands)"
		print

	
	# crunch IntervalStats results:
	if option.mode == "crunch":
		
		# define input and output paths:
		resultspath = coassociationspath + option.peaks + "/" + option.name + "/results/"
		summarypath = coassociationspath + option.peaks + "/" + option.name + "/summary/"
		passingpath = coassociationspath + option.peaks + "/" + option.name + "/passing/p" + pvaluecutoff_handle + "/"
		general.pathGenerator(passingpath)
		
		# load result files:
		result_files = os.listdir(resultspath)
		
		# count executed tests:
		test_count = len(result_files)
		
		# count executed tests:
		comp_count = int(numpy.sqrt(test_count))
		
		# gather fractions and passing lines from results files:
		matrix, datasets, r = dict(), list(), 1
		print
		for result_file in result_files:
			sys.stdout.write("\rImporting co-association results: {0}".format(r))
			sys.stdout.flush()
			r += 1
			basename = result_file.replace("mapcas_" + option.name + "_", "").replace("_peaks", "").replace(".bed", "")
			i, j = basename.split("_vs_")
			datasets.append(i)
			datasets.append(j)
			inlines = open(resultspath + result_file).readlines()
			line_count = len(inlines)
			if not i in matrix:
				matrix[i] = dict()
			p, c, x, pvalues = 0, 0, 0, list()
			if line_count > 0:
				for inline in inlines:
					pvalue = float(inline.strip().split("\t")[6])
					pvalues.append(pvalue)
					if pvalue < option.cutoff:
						p += 1
					if pvalue*comp_count < option.cutoff:
						c += 1
					if pvalue*test_count < option.cutoff:
						x += 1
				if len(inlines) != 0:
					fraction_passing = float(p)/len(inlines)
					fraction_correct = float(c)/len(inlines)
					fraction_extreme = float(x)/len(inlines)
				else:
					fraction_passing = 0
					fraction_correct = 0
					fraction_extreme = 0
				avgPvalue = numpy.mean(pvalues)
				medPvalue = numpy.median(pvalues)
			else:
				fraction_passing = 0
				fraction_correct = 0
				fraction_extreme = 0
				avgPvalue = 1
				medPvalue = 1
			matrix[i][j] = [fraction_passing, fraction_correct, fraction_extreme, avgPvalue, medPvalue]
			
			#mapcas_orthologous_TSSs_ce_YL445_EFL-1_YA_yale_stn_peaks_vs_ce_OP90_NHR-6_L4_yale_stn_peaks.bed
			#V:5582699:5585109	V:5592093:5592995	2411	6984	220375	1748292	0.126052
		
		# export fractions into matrix format:
		print "\nExporting co-association summary:", len(result_files)
		f_output = open(summarypath + "mapcas_crunch_" + option.name + "_p" + pvaluecutoff_handle + "_matrix.txt", "w")
		print >>f_output, "\t".join(["i","j","fraction.passing","fraction.corrected","fraction.extreme","pvalue.mean","pvalue.median"])
		datasets = list(set(datasets))
		for i in datasets:
			for j in datasets:
				if i == j:
					print >>f_output, "\t".join(map(str, [i, i] + [1, 1, 1, 0, 0]))
				elif i in matrix and j in matrix[i]:
					print >>f_output, "\t".join(map(str, [i, j] + matrix[i][j]))
				else:
					print >>f_output, "\t".join(map(str, [i, j] + [0, 0, 0, 1, 1]))
		f_output.close()
		print
	
	
	# report IntervalStats results:
	if option.mode == "report":
		
		# define input and output paths:
		resultspath = coassociationspath + option.peaks + "/" + option.name + "/results/"
		summarypath = coassociationspath + option.peaks + "/" + option.name + "/summary/"
		passingpath = coassociationspath + option.peaks + "/" + option.name + "/passing/p" + pvaluecutoff_handle + "/"
	
		# define output file:
		f_output = open(summarypath + "mapcas_report_" + option.name + "_p" + pvaluecutoff_handle + "_matrix.txt", "w")
		
		# load header dictionary:
		HD = general.build_header_dict(summarypath + "mapcas_crunch_" + option.name + "_p" + pvaluecutoff_handle + "_matrix.txt")
		
		# generate matrix for mirrored calculations:
		matrix_passing, matrix_correct, matrix_extreme = dict(), dict(), dict()
		inlines = open(summarypath + "mapcas_crunch_" + option.name + "_p" + pvaluecutoff_handle + "_matrix.txt").readlines()
		inlines.pop(0)
		for inline in inlines:
			i, j, fraction_passing, fraction_correct, fraction_extreme, pvalue_mean, pvalue_median = inline.strip().split("\t")
			if not i in matrix_passing:
				matrix_passing[i] = dict()
				matrix_correct[i] = dict()
				matrix_extreme[i] = dict()
			matrix_passing[i][j] = float(fraction_passing)
			matrix_correct[i][j] = float(fraction_correct)
			matrix_extreme[i][j] = float(fraction_extreme)
		
		# parse dataset summary file:
		HD = general.build_header_dict(summarypath + "mapcas_crunch_" + option.name + "_p" + pvaluecutoff_handle + "_matrix.txt")
		print >>f_output, "\t".join(["i","j","i.factor","j.factor","i.context","j.context","fraction.passing","fraction.corrected","fraction.extreme","pvalue.mean","pvalue.median", "mirror.passing", "mirror.corrected", "mirror.extreme", "maximo.passing", "maximo.corrected", "maximo.extreme"])
		for inline in inlines:
			i, j, fraction_passing, fraction_correct, fraction_extreme, pvalue_mean, pvalue_median = inline.strip().split("\t") 
			mirror_passing = numpy.mean([matrix_passing[i][j], matrix_passing[j][i]])
			mirror_correct = numpy.mean([matrix_correct[i][j], matrix_correct[j][i]])
			mirror_extreme = numpy.mean([matrix_extreme[i][j], matrix_extreme[j][i]])
			maximo_passing = max([matrix_passing[i][j], matrix_passing[j][i]])
			maximo_correct = max([matrix_correct[i][j], matrix_correct[j][i]])
			maximo_extreme = max([matrix_extreme[i][j], matrix_extreme[j][i]])
			iorganism, istrain, ifactor, icontext, iinstitute, imethod = metrn.labelComponents(i)
			jorganism, jstrain, jfactor, jcontext, jinstitute, jmethod = metrn.labelComponents(j)
			i = metrn.labelGenerator(option.label, mode="label", dataset=i).replace("-S3", "S3").replace("-hESC", "hesc")
			j = metrn.labelGenerator(option.label, mode="label", dataset=j).replace("-S3", "S3").replace("-hESC", "hesc")
			#i = ifactor + "." + icontext
			#j = jfactor + "." + jcontext
			print >>f_output, "\t".join(map(str, [i, j, ifactor, jfactor, icontext, jcontext, fraction_passing, fraction_correct, fraction_extreme, pvalue_mean, pvalue_median, mirror_passing, mirror_correct, mirror_extreme, maximo_passing, maximo_correct, maximo_extreme]))
		f_output.close()
		
	
	# classify paralogous and non-paralogous results:
	if option.mode == "family":
		
		# define input and output paths:
		resultspath = coassociationspath + option.peaks + "/" + option.name + "/results/"
		summarypath = coassociationspath + option.peaks + "/" + option.name + "/summary/"
		passingpath = coassociationspath + option.peaks + "/" + option.name + "/passing/p" + pvaluecutoff_handle + "/"
	
		# determine species:
		organismTag = option.peaks[:2]
		
		# load paralogous factors:
		paralogMatrix = general.build2(orthologspath + "paralogs/mapdata_paralogs_" + organismTag + "_matrix.txt", i="i", j="j", mode="matrix", counter=True)
		
		# define output file:
		f_output = open(summarypath + "mapcas_family_" + option.name + "_p" + pvaluecutoff_handle + "_matrix.txt", "w")
		
		# load header dictionary:
		HD = general.build_header_dict(summarypath + "mapcas_report_" + option.name + "_p" + pvaluecutoff_handle + "_matrix.txt")
		
		# generate matrix for mirrored calculations:
		matrix_passing, matrix_correct, matrix_extreme = dict(), dict(), dict()
		inlines = open(summarypath + "mapcas_report_" + option.name + "_p" + pvaluecutoff_handle + "_matrix.txt").readlines()
		header = inlines.pop(0)
		print >>f_output, header.strip() + "\t" + "paralog"
		for inline in inlines:
			initems = inline.strip().split("\t")
			i, j = initems[HD["i.factor"]], initems[HD["j.factor"]]
			if i in paralogMatrix and j in paralogMatrix[i]:
				paralogFlag = "1"
			else:
				paralogFlag = "0"
			print >>f_output, inline.strip() + "\t" + paralogFlag
		f_output.close()
		
		
	
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

#python mapCAs.py --path ~/ceTRN --mode launch --peaks optimal_standard_totals_sx_filter --domain in2shape_wbTrans_adjust_up2000_dn200.nh.bed --name wbTrans_adjust_up2000_dn200 --threads 3

#python mapCAs.py --path ~/ceTRN --mode crunch --peaks optimal_standard_totals_sx_filter --domain in2shape_wbTrans_adjust_up2000_dn200.nh.bed --name wbTrans_adjust_up2000_dn200 --threads 5 --cutoff 0.05
#python mapCAs.py --path ~/ceTRN --mode report --peaks optimal_standard_totals_sx_filter --domain in2shape_wbTrans_adjust_up2000_dn200.nh.bed --name wbTrans_adjust_up2000_dn200 --threads 5 --cutoff 0.05
#python mapCAs.py --path ~/ceTRN --mode crunch --peaks optimal_standard_totals_sx_filter --domain in2shape_wbTrans_adjust_up2000_dn200.nh.bed --name wbTrans_adjust_up2000_dn200 --threads 5 --cutoff 0.01
#python mapCAs.py --path ~/ceTRN --mode report --peaks optimal_standard_totals_sx_filter --domain in2shape_wbTrans_adjust_up2000_dn200.nh.bed --name wbTrans_adjust_up2000_dn200 --threads 5 --cutoff 0.01

#python mapCAs.py --path ~/ceTRN --mode launch --peaks optimal_standard_totals_sx_limits --domain in2shape_wbTrans_adjust_up2000_dn200.nh.bed --name wbTrans_adjust_up2000_dn200 --threads 5
#python mapCAs.py --path ~/ceTRN --mode crunch --peaks optimal_standard_totals_sx_limits --domain in2shape_wbTrans_adjust_up2000_dn200.nh.bed --name wbTrans_adjust_up2000_dn200 --threads 5
#python mapCAs.py --path ~/ceTRN --mode report --peaks optimal_standard_totals_sx_limits --domain in2shape_wbTrans_adjust_up2000_dn200.nh.bed --name wbTrans_adjust_up2000_dn200 --threads 5