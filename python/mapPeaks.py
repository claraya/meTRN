#!/usr/bin/env python
# generate peak set complete files, binding region files, and report files!

import sys
import time
import optparse
import general
import numpy
import pickle
import pdb
import metrn
import modencode
import os

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())


""" define a function to format a peak BED file as a genome file """
def regionFormat(peaksfile, bedfile, outfile, header="OFF"):
	peakDict = general.build2(peaksfile, id_column="feature", value_columns=["chrm","start","end","strand"], header_dict=metrn.bedHeader, header=False)
	f_output = open(outfile, "w")
	indata = open(bedfile)
	inline = indata.readline()
	if header == "ON":
		inline = indata.readline()
	while inline:
		chrm, start, stop, feature, occupancy, strand, density, datasetCounts, factorCounts, contextCounts, details,  = inline.strip().split("\t")
		for peak in details.split(";"):
			if peak in peakDict:
				output = [feature]
				output.append(peakDict[peak]["start"])
				output.append(peakDict[peak]["end"])
				output.append(peak)
				output.append(occupancy)
				output.append(peakDict[peak]["strand"])
				output.append(density)
				output.append(datasetCounts)
				output.append(factorCounts)
				output.append(contextCounts)
				output.append(peakDict[peak]["chrm"])
				print >>f_output, "\t".join(output)
		inline = indata.readline()
	f_output.close()
	

""" define a function to format a region BED file as a genome file """
def genomeFormat(bedfile, outfile, header="OFF"):
	f_output = open(outfile, "w")
	indata = open(bedfile)
	inline = indata.readline()
	if header == "ON":
		inline = indata.readline()
	while inline:
		chrm, start, stop, feature, score, strand = inline.strip().split("\t")[:6]
		print >>f_output, "\t".join([feature, str(max([int(start),int(stop)]))])
		inline = indata.readline()
	f_output.close()


def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "path from script to files")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "peaks to be used for analysis")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "operations to be performed")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--library", action = "store", type = "string", dest = "library", help = "library configuration final report card")
	parser.add_option("--rmdup", action = "store", type = "string", dest = "rmdup", help = "Remove duplicates library", default="OFF")
	parser.add_option("--overwrite", action = "store", type = "string", dest = "overwrite", help = "Overwrite files?", default="OFF")
	parser.add_option("--rename", action = "store", type = "string", dest = "rename", help = "Targets to rename. Comma-separated list of 'target:replacement' pairs to search and replace.", default="OFF")
	parser.add_option("--contexts", action = "store", type = "string", dest = "contexts", help = "Contexts to compare", default="OFF")
	parser.add_option("--species", action = "store", type = "string", dest = "species", help = "Species to compare", default="OFF")
	parser.add_option("--orthology", action = "store", type = "string", dest = "orthology", help = "Use 'direct', 'family' (Yong's), or 'group' (Pouya's) orthologs?", default="direct")
	parser.add_option("--nametag", action = "store", type = "string", dest = "nametag", help = "Orthology nametag: nametagHsCe", default="ortho")
	parser.add_option("--commonNames", action = "store", type = "string", dest = "commonNames", help = "Grab common names file?", default="ON")
	parser.add_option("--familyFiles", action = "store", type = "string", dest = "familyFiles", help = "Grab cleaned files?", default="formatted")
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
	orthologspath = path_dict["orthologs"]
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
	
	# define organisms:
	organismTags = ["hs","mm","ce","dm"]
	
	# build mode:
	if option.mode == "build":
	
		# scan peaks and signal tracks:
		print
		peakSets = os.listdir(peakspath)
		mx = 1
		for peakSet in sorted(peakSets):
			if not "mappeaks" in peakSet and not ".DS" in peakSet:
				
				# set output headers:
				#peakHeaders = ["chrm","start","end","feature","score","strand","signal","pvalue","qvalue","point","size","dataset","organism","strain","factor","context","institute","method"]
				dataHeaders = ["dataset","peaks","avg.size", "organism", "strain", "factor", "context", "institute", "method"]
				
				# specify output files:
				completefile = peakspath + "mappeaks_" + peakSet + "_complete.bed"
				collapsefile = peakspath + "mappeaks_" + peakSet + "_collapse.bed"
				compiledfile = peakspath + "mappeaks_" + peakSet + "_compiled.bed"
				reportfile = peakspath + "mappeaks_" + peakSet + "_report.txt"
				
				# determine if the file has already been processed:
				process = False
				if option.overwrite == "ON" or not "mappeaks_" + peakSet + "_compiled.bed" in os.listdir(peakspath):
					process = True
				
				# overwrite or process missing files:
				if process:
				
					print "Processing:", peakSet
					
					# define output files:
					c_output = open(completefile, "w")
					r_output = open(reportfile, "w")
					#print >>m_output, "\t".join(peakHeaders)
					print >>r_output, "\t".join(dataHeaders)
					
					# export and prebuild peaks dictionary:
					print "Generating complete peaks file..."
					peakfiles = os.listdir(peakspath + peakSet)
					for peakfile in peakfiles:
						if not ".DS" in peakfile:
							
							dataset = peakfile.replace("_peaks.bed","")
							organism, strain, factor, context, institute, method = metrn.labelComponents(dataset)
								
							peaks, sizes = 0, list()
							for inline in open(peakspath + peakSet + "/" + peakfile).readlines():
								chrm, start, end, feature, score, strand, signal, pvalue, qvalue, point = inline.strip().split("\t")[:10]
								start, end, score = int(start), int(end), float(score)
								size = end-start
								strand = "+"
								output = [chrm, start, end, feature, str(score), strand, signal, pvalue, qvalue, point, size, dataset, organism, strain, factor, context, institute, method]
								print >>c_output, "\t".join(map(str, output))
								sizes.append(size)
								peaks += 1
							
							output = [dataset, peaks, numpy.mean(sizes), organism, strain, factor, context, institute, method]
							print >>r_output, "\t".join(map(str, output))
					
					# close output file:
					c_output.close()
					r_output.close()
					
					# generate collapsed file:
					print "Generating collapsed peaks file..."
					command = "mergeBed -i " + completefile + " -nms > " + collapsefile
					os.system(command)
					
					# generate a density-style analysis file:
					print "Annotating collapsed peaks file..."
					metrn.densityBed(collapsefile, compiledfile)
					
					# delete collapse file (since it is redundant):
					command = "rm -rf " + collapsefile
					os.system(command)
					
					print
	
	
	# binding core mode:
	if option.mode == "core":
	
		# process peak sets:
		print
		print "Generating core-regions..."
		#for peakSet in ["dm_selection_reg_ee_raw"]:
		#for peakSet in ["hs_selection_reg_h1_raw"]:
		for peakSet in os.listdir(peakspath):
			
			# determine organism:
			organismTag = peakSet[:2]
			
			if organismTag in metrn.reference:
			
				# specify input files:
				completefile = peakspath + "mappeaks_" + peakSet + "_complete.bed"
				collapsefile = peakspath + "mappeaks_" + peakSet + "_collapse.bed"
				compiledfile = peakspath + "mappeaks_" + peakSet + "_compiled.bed"
				
				# specify operation files:
				genomicsinput = inpath + metrn.reference[organismTag]["complete_sizes"]
				completeinput = peakspath + "mappeaks_" + peakSet + "_complete.tmp.bed"
				compiledinput = peakspath + "mappeaks_" + peakSet + "_compiled.tmp.bed"
				#completegraph = peakspath + "mappeaks_" + peakSet + "_complete.bdg"
				#compiledgraph = peakspath + "mappeaks_" + peakSet + "_compiled.bdg"
				completegraph = peakspath + "mappeaks_" + peakSet + "_bedgraph.bed"
					
				# for processing with feature information retained:
				#regionFormat(completefile, compiledfile, completeinput)
				#genomeFormat(compiledfile, compiledinput)
				
				# should we process this peak set?
				if not "mappeaks_" + peakSet + "_bedgraph.bed" in os.listdir(peakspath) or option.overwrite == "ON":
				
					print "Processing:", peakSet
				
					# for processing with chromosome information only (standard):
					#command = "sort -k 1,1 -k2,2 -n " + completefile + " > " + completeinput
					command = "sortBed -i " + completefile + " > " + completeinput
					os.system(command)
					
					# generate bedGraph from peak set:
					command = "genomeCoverageBed -bg -i " + completeinput + " -g " + genomicsinput + " > " + completegraph
					os.system(command)
					
					# remove temporary inputs:
					command = "rm -rf " + completeinput
					os.system(command)
					
					command = "rm -rf " + compiledinput
					os.system(command)
		
		print
	
	
	# orthologs mode:
	if option.mode == "orthologs":
	
		# identify target specie and comparison species:
		targetTag = option.peaks.split("_")[0]
		speciesTags = option.species.split(",")
		
		# generate output peaks name:
		orthologTag = metrn.orthologLabel(targetTag, speciesTags)
		orthologPeaks = option.peaks.replace(option.rename, option.nametag + orthologTag)
		general.pathGenerator(peakspath + orthologPeaks)
		
		# target specie orthologs:
		if option.orthology == "direct":
			orthologs = metrn.orthologFinder(targetTag, speciesTags, path=orthologspath + "orthologs/", orthology=option.orthology, commonNames=option.commonNames)
		elif option.orthology == "family":
			orthologs = metrn.orthologFinder(targetTag, speciesTags, path=orthologspath + "families/", orthology=option.orthology, familyFiles=option.familyFiles)
		elif option.orthology == "groups":
			orthologs = metrn.orthologFinder(targetTag, speciesTags, path=orthologspath + "groups/", orthology=option.orthology, familyFiles=option.familyFiles)
		
		print
		print "Orthologs", "(" + str(len(orthologs)) + "):"
		print "\n".join(orthologs)
		print
		
		# gather peak files:
		k = 0
		peakfiles = os.listdir(peakspath + option.peaks)
		#print orthologs
		#print len(peakfiles)
		print "Generating peak set:", orthologPeaks
		for ortholog in orthologs:
			for peakfile in peakfiles:
				dataset = peakfile.replace("_peaks.bed","")
				organism, strain, factor, context, institute, method = metrn.labelComponents(dataset)
				if ortholog.upper() == factor.upper():
					command = "cp " + peakspath + option.peaks + "/" + peakfile + " " + peakspath + orthologPeaks + "/" + peakfile
					os.system(command)
					k += 1
		print "Orthologs found:", len(orthologs)
		print "Datasets found:", k
		print
	
	# context comparison mode:
	if option.mode == "comparison":
	
		# identify target context comparison:
		targetContexts = option.contexts.split(",")
		
		# generate output peaks name:
		print
		print "Generating comparison dictionary..."
		comparisonDict, alternateDict = dict(), dict()
		for context in targetContexts:
			opposed = option.contexts.replace(context, "").replace(",", "")
			contextPeaks = option.peaks.replace("_xx_", "_" + context + "_")
			opposedPeaks = option.peaks.replace("_xx_", "_" + opposed + "_")
			comparisonPeaks = option.peaks.replace(option.rename, option.nametag).replace("_xx_", "_" + context + "_")
			comparisonDict[context] = comparisonPeaks
			general.pathGenerator(peakspath + comparisonPeaks)
			comparisonDict[context] = comparisonPeaks, contextPeaks, opposedPeaks
			alternateDict[context] = opposed
			#print context, comparisonPeaks, contextPeaks, opposedPeaks, opposed
			
		# find comparisons:
		print "Transferring matched peak files..."
		for context in targetContexts:
			comparisonPeaks, contextPeaks, opposedPeaks = comparisonDict[context]
			contextReport = general.build2(peakspath + "mappeaks_" + contextPeaks + "_report.txt", id_column="factor", header="auto", skip=True, verbose=False)
			opposedReport = general.build2(peakspath + "mappeaks_" + opposedPeaks + "_report.txt", id_column="factor", header="auto", skip=True, verbose=False)
			overlapReport = sorted(list(set(contextReport.keys()).intersection(set(opposedReport.keys()))))
			general.pathGenerator(peakspath + comparisonPeaks)
			for factor in overlapReport:
				command = "cp " + peakspath + contextPeaks + "/" + contextReport[factor]["dataset"] + "_peaks.bed" + " " + peakspath + comparisonPeaks + "/"
				os.system(command)
		print "Matches (factor overlap):", len(overlapReport)
		print
	
						
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())
