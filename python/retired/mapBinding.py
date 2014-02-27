#!/usr/bin/env python
# map transcription factor binding (peak calls) into unique binding regions

import sys
import time
import optparse
import general
import numpy
import metrn
import modencode
import network
import bed
import os
import copy
import pdb
import re
import pickle

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())


""" define classes and functions of internal use """
def setorder(inlist):
	outlist = list()
	for item in inlist:
		if not item in outlist:
			outlist.append(item)
	return outlist


""" define a function to count peaks in a gffkde output file """
def gffkdePeakCounter(indata, mode="file"):
	if mode == "file":
		processed, inlines = list(), open(indata).readlines()
	elif mode == "list":
		processed, inlines = list(), indata
	for inline in inlines:
		for peakData in inline.strip().split("\t")[7].rstrip(";").split(";")[1:]:
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


""" define a function to annotate peaks based on the density files """
def annotatePeaks(infile, peak_dict, maphot=True, export_file=False):
	
	# open output merged file:
	if export_file:
		m_output = open(export_file, "wb")
		print >>m_output, "\t".join(["chrm", "start", "end", "feature", "score", "strand", "dataset.count", "factor.count", "stage.count", "signal.avg", "signal.med", "signal.sum", "signal.min", "signal.max", "collapsed.strains", "collapsed.factors", "collapsed.stages", "collapsed.institutes", "collapsed.methods", "collapsed.info"])
	
	# gather peaks from merged file:
	bind_dict, region_dict, index, peaks = dict(), dict(), 1, list()
	inlines = open(infile).readlines()
	if maphot:
		inlines.pop(0)
	for inline in inlines:
		if maphot:
			#V	InferredHS300bw	HOTspot	58437	58743	1.00084493975139	.	.	1,1;OP327_F23F12.9_EM_yale_stn_sample_optimal_narrowPeak_MACS_peak_9720,58593,0.999950001249979;
			chrm, start, end, feature, hits, strand, density, collapsed_details = inline.strip("\n").split("\t")
			signals, contributions, info = list(), list(), list()
			datasets, strains, factors, stages, institutes, methods = list(), list(), list(), list(), list(), list()
			collapsed_details = collapsed_details.split(";")
			collapsed_details.pop(0)
			for collapsed_detail in collapsed_details:
				if not collapsed_detail == "":
					collapsed_peak, midpoint, contribution = collapsed_detail.split(",")
					dataset, peak = collapsed_peak.split("_peaks_")
					strain, factor, stage, institute, method = dataset.split("_")[:5]
					pchrm, pstart, pend, pscore, pstrand, psignal, ppvalue, pqvalue, ppoint = peak_dict[strain][factor][stage][institute][method][peak]
					signals.append(float(psignal))
					contributions.append(str(contribution))
					datasets.append(dataset)
					strains.append(strain)
					factors.append(factor)
					stages.append(stage)
					institutes.append(institute)
					methods.append(method)
					info.append(",".join(map(str, [dataset,peak,pstart,pend,psignal,contribution])))
					if not strain in bind_dict:
						bind_dict[strain] = dict()
					if not factor in bind_dict[strain]:
						bind_dict[strain][factor] = dict()
					if not stage in bind_dict[strain][factor]:
						bind_dict[strain][factor][stage] = dict()
					if not institute in bind_dict[strain][factor][stage]:
						bind_dict[strain][factor][stage][institute] = dict()
					if not method in bind_dict[strain][factor][stage][institute]:
						bind_dict[strain][factor][stage][institute][method] = dict()
					bind_dict[strain][factor][stage][institute][method][peak] = [feature, chrm, start, end, contribution, hits, density, pstrand, psignal, ppvalue, pqvalue, ppoint]
			
			# store tfbs (region) data:
			number = len(signals)
			values = [numpy.mean(signals), numpy.median(signals), sum(signals), min(signals), max(signals)]
			datasets = setorder(datasets)
			strains = setorder(strains)
			factors = setorder(factors)
			stages = setorder(stages)
			institutes = setorder(institutes)
			methods = setorder(methods)
			
			dataset_count = len(datasets)
			factor_count = len(factors)
			stage_count = len(stages)
			
			collapsed_strains = ";".join(strains)
			collapsed_factors = ";".join(factors)
			collapsed_stages = ";".join(stages)
			collapsed_institutes = ";".join(institutes)
			collapsed_methods = ";".join(methods)
			collapsed_info = ";".join(info)
				
			#hits = int(hits)
			#if hits != number or number != dataset_count:
			#	print feature, hits, float(density), number, dataset_count, factor_count, stage_count
			#	pdb.set_trace()
			
			region_dict[index] = map(str, [chrm, start, end, feature, number, strand, dataset_count, factor_count, stage_count] + values + [collapsed_strains, collapsed_factors, collapsed_stages, collapsed_institutes, collapsed_methods, collapsed_info])
			index += 1
	
	# return the peak (bind) and region (tfbs) dictionaries, and the peak count!
	return bind_dict, region_dict, gffkdePeakCounter(inlines, mode="list")
	
	
def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "Path from script to files")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "Operations to be executed: build, map")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "Basename for target peaks", default="OFF")
	parser.add_option("--source", action = "store", type = "string", dest = "source", help = "Which peaks should be used as source for 'build' and 'expand' modes?", default=False)
	parser.add_option("--infile", action = "store", type = "string", dest = "infile", help = "input BED file of feature coordinates", default=False)
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "Output name for mapping, input name for map", default="OFF")
	parser.add_option("--header", action = "store", type = "string", dest = "header", help = "Does the annotation file have a header?", default="ON")
	parser.add_option("--target", action = "store", type = "string", dest = "target", help = "Feature column:type to target", default="feature")
	parser.add_option("--fraction", action = "store", type = "string", dest = "fraction", help = "Fractional overlap required", default="0.1")
	parser.add_option("--queries", action = "store", type = "string", dest = "queries", help = "Which types of features should be registered", default="standard")
	parser.add_option("--ids", action = "store", type = "string", dest = "ids", help = "How should feature IDs be keyed in?", default="standard")
	parser.add_option("--others", action = "store", type = "string", dest = "others", help = "Should 'other' category be created?", default="OFF")
	parser.add_option("--elsewhere", action = "store", type = "string", dest = "elsewhere", help = "Should features that do not overlap categories be counted?", default="OFF")
	parser.add_option("--start", action = "store", type = "int", dest = "start", help = "Start index for range", default=1)
	parser.add_option("--stop", action = "store", type = "int", dest = "stop", help = "End index for range", default=40)
	parser.add_option("--headerDict", action = "store", type = "string", dest = "headerDict", help = "Header dictionary...", default="bed")
	parser.add_option("--label", action = "store", type = "string", dest = "label", help = "Dataset labeling mode", default="factor.context")
	parser.add_option("--rename", action = "store", type = "string", dest = "rename", help = "Targets to rename. Comma-separated list of 'target:replacement' pairs to search and replace.", default="OFF")
	parser.add_option("--reference", action = "store", type = "string", dest = "reference", help = "Query to which ratios should be scaled...", default="OFF")
	parser.add_option("--order", action = "store", type = "string", dest = "order", help = "How should the inputs be ordered?", default="OFF")
	parser.add_option("--prioritize", action = "store", type = "string", dest = "prioritize", help = "Should we prioritize counts? That is, assign overlaps preferentially to classes.", default="OFF")
	parser.add_option("--complexity", action = "store", type = "string", dest = "complexity", help = "Max complexity (factor.count) allowed", default=False)
	parser.add_option("--method", action = "store", type = "string", dest = "method", help = "Filter for a specific method?", default=False)
	parser.add_option("--exclude", action = "store", type = "string", dest = "exclude", help = "Comma-separated list of factors to exclude", default="")
	parser.add_option("--window", action = "store", type = "string", dest = "window", help = "Window surrounding feature for promoter and downstream analysis", default="0")
	parser.add_option("--max", action = "store", type = "string", dest = "max", help = "Maximum number of peaks allowed per experiment", default=False)
	parser.add_option("--filter", action = "store", type = "string", dest = "filter", help = "Filter peaks based on density or overlap?", default="OFF")
	parser.add_option("--quality", action = "store", type = "string", dest = "quality", help = "Filter peaks based on quality?", default="OFF")
	parser.add_option("--nuclear", action = "store", type = "string", dest = "nuclear", help = "Peaks are only nuclear?", default="ON")
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
	neuronspath = path_dict["neurons"]
	cellspath = path_dict["cells"]
	
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
		
	# load gene ID dictionaries:
	#id2name_dict, name2id_dict = modencode.idBuild(inpath + metrn.reference[organismTag]["gene_ids"], "Sequence Name (Gene)", "Gene Public Name", mode="label", header=True, upper=True)
	
	# set some of the basic parameters:
	exclude = option.exclude.split(",")
	window = int(option.window)
	if option.max:
		max_peaks = int(option.max)
	else:
		max_peaks = "X"
	
	# define output name:
	outkey = "mapbinding_" + option.peaks
	
	# define relevant paths:
	bindingpath = bindingpath + option.peaks + "/" + option.name + "/"
	mappingpath = bindingpath + "mapping/"
	datasetpath = bindingpath + "dataset/"
	overlappath = bindingpath + "overlap/"
	picklespath = bindingpath + "pickles/"
	reportspath = bindingpath + "reports/"
	summarypath = bindingpath + "summary/"
	general.pathGenerator(mappingpath)
	general.pathGenerator(datasetpath)
	general.pathGenerator(overlappath)
	general.pathGenerator(picklespath)
	general.pathGenerator(reportspath)
	general.pathGenerator(summarypath)
	
	# Note: These are the retired queries:
	# queries = ['CDS', 'protein_coding_primary_transcript', 'ncRNA_primary_transcript', 'snoRNA', 'Pseudogene', 'tRNA', 'miRNA_primary_transcript', 'snRNA', 'rRNA_primary_transcript', 'intron', 'exon', 'TSS', 'SL1_acceptor_site', 'SL2_acceptor_site', 'transcription_end_site', 'polyA_site', 'polyA_signal_sequence', 'five_prime_UTR', 'three_prime_UTR', 'snlRNA', 'DNAse_I_hypersensitivity']
	
	# map overlap between peaks and genomic features:
	if option.mode == "map:overlap":
		
		# define input and output files:
		c_infile = peakspath + "mappeaks_" + option.peaks + "_complete.bed"
		c_outfile = overlappath + outkey + "_complete_overlap_"
		c_tmpfile = overlappath + outkey + "_complete.tmp"
		
		m_infile = peakspath + "mappeaks_" + option.peaks + "_compiled.bed"
		m_outfile = overlappath + outkey + "_compiled_overlap_"
		m_tmpfile = overlappath + outkey + "_compiled.tmp"
		
		i_infile = annotationspath + option.infile
		i_tmpfile = overlappath + option.infile + ".tmp"
		i_colfile = overlappath + option.infile + "_header.tmp"
		
		f_outfile = overlappath + outkey + "_compiled"
		s_outfile = summarypath + outkey + "_compiled"
		
		# define feature key setup:
		idComplex = ["chrm", "start", "end", "feature"]
		
		# define header presence
		if option.header == "ON":
			headerFlag = True
		else:
			headerFlag = False
		
		# define annotation headers:
		if option.headerDict == "auto":
			annotationHeader = general.build_header_dict(i_infile)
		elif option.headerDict == "bed":
			annotationHeader = metrn.bedHeader
		else:
			annotationHeader = dict()
			for entry in option.headerDict.split(","):
				key, value = entry.split(":")
				annotationHeader[key] = int(value)
			
		# load annotations:
		print
		print "Loading input annotations..."
		annotationDict = general.build2(i_infile, id_complex=idComplex, header_dict=annotationHeader, header=headerFlag, id_include=True, separator=":")
		
		#x = annotationDict.keys()[0]
		#y = annotationDict[x]
		#overlapBed = annotationDict
		#query = "7_Egn4"
		#queryColumn = "feature"
		#queryBed = bed.valueFilter(overlapBed, filterDict={queryColumn:query}, modeDict={queryColumn:"match"}, structure="tree")
		#print queryBed.keys()
		#pdb.set_trace()
		
		# define standard targets/queries:
		if option.queries == "standard":
			queries = ['tss', 'five_prime_utr', 'three_prime_utr', 'exon', 'intron', 'protein_coding_primary_transcript']
		elif option.queries == "auto":
			queries = list()
			for feature in annotationDict:
				if not feature in option.exclude.split(","):
					queries.append(annotationDict[feature][option.target])
			queries = sorted(list(set(queries)))
		elif option.queries != "OFF":
			queries = list()
			for query in option.queries.split(","):
				if not query in option.exclude.split(","):
					queries.append(query) 
		
		print "Preparing temporary files..."
		command = 'cp ' + c_infile + ' ' + c_tmpfile
		os.system(command)
		
		command = 'cp ' + m_infile + ' ' + m_tmpfile
		os.system(command)
		
		if option.header == "ON":
			command = 'grep -v "feature" ' + i_infile + ' > ' + i_tmpfile
			os.system(command)
		
			command = "head -n 1 " + i_infile + ' > ' + i_colfile
			os.system(command)
		
		else:
			command = "cp " + i_infile + " " + i_tmpfile
			os.system(command)
			
			
		print "Intersecting annotation features with peaks..."
		command = "intersectBed -u -f 0.1 -a " + i_tmpfile + " -b " + c_tmpfile + " > " + c_outfile
		os.system(command)
		
		print "Intersecting annotation features with regions..."
		command = "intersectBed -u -f 0.1 -a " + i_tmpfile + " -b " + m_tmpfile + " > " + m_outfile
		os.system(command)
		
		print "Gathering compiled peak regions..."
		mergedBed = general.build2(m_infile, id_column="feature", header_dict=metrn.regionHeader, header=False, id_include=True)
		
		print
		overlapDict = dict()
		technique_dict = {
			"ox":["equal","exactly"],
			"og":["greater.equal","at least"],
			"ol":["lesser.equal","at most"],
			"ob":["greater","more than"]
			}
		for technique in ["ox","og","ol","ob"]:
			overlapDict[technique] = dict()
			for factorCount in range(option.start, option.stop + 1):
				modeCall, modeText = technique_dict[technique]
				
				fm_outfile = f_outfile.replace("_compiled", "_compiled_" + technique + str(factorCount).zfill(2)) + ".bed"
				om_outfile = fm_outfile.replace(".bed","") + "_overlap"
				
				filteredBed = bed.valueFilter(mergedBed, filterDict={"occupancy":factorCount}, modeDict={"occupancy":modeCall}, structure="tree")
				bed.export(filteredBed, fm_outfile, headerMode="explicit", headerDict=metrn.regionHeader, structure="tree")
				print "...with", modeText, factorCount, "factor(s) bound:", len(filteredBed)
				
				command = "intersectBed -u -a " + i_tmpfile + " -b " + fm_outfile + " > " + om_outfile
				os.system(command)
				
				queryCount = 0
				overlapDict[technique][factorCount] = dict()
				#overlapBed = general.build2(om_outfile, id_column=idColumn, header_dict=annotationHeader, header=True, id_complex=idComplex, id_include=True, id_index=True, separator=":")
				overlapBed = general.build2(om_outfile, id_complex=idComplex, header_dict=annotationHeader, header=headerFlag, id_include=True, separator=":")
				for query in queries:
					#queryBed = bed.valueFilter(overlapBed, filterDict={queryColumn:query.lower()}, modeDict={queryColumn:"match"}, structure="tree")
					queryBed = bed.valueFilter(overlapBed, filterDict={option.target:query}, modeDict={option.target:"match"}, structure="tree")
					overlapDict[technique][factorCount][query] = len(queryBed)
					queryCount += len(queryBed)
				if option.others == "ON":
					overlapDict[technique][factorCount]["others"] = len(overlapBed) - queryCount
			print
				
		print
		print "Removing temporary files..."
		command = "rm -rf " + bindingpath + "*.tmp"
		os.system(command)
		
		command = "rm -rf " + overlappath + "*.tmp"
		os.system(command)
		
		command = "rm -rf " + annotationspath + "*.tmp"
		os.system(command)
		print
		
		# update queries to include "others" if necessary:
		if option.others == "ON":
			queries.append("others")
		
		# export per factor-count feature-type summary:
		for technique in overlapDict:
			sv_sumfile = s_outfile.replace("_compiled", "_compiled_" + technique + "XX") + "_overlap_summary_values"
			sn_sumfile = s_outfile.replace("_compiled", "_compiled_" + technique + "XX") + "_overlap_summary_normal"
			sv_output = open(sv_sumfile, "w")
			sn_output = open(sn_sumfile, "w")
			print >>sv_output, "\t".join(queries)
			print >>sn_output, "\t".join(queries)
			for factorCount in sorted(overlapDict[technique].keys()):
				values, normalized = list(), list()
				for query in queries:
					values.append(overlapDict[technique][factorCount][query])
				for value in values:
					if sum(values) > 0:
						normalized.append(float(value)/sum(values))
					else:
						normalized.append(0)
				print >>sv_output, "\t".join(map(str, values))
				print >>sn_output, "\t".join(map(str, normalized))
			sv_output.close()
			sn_output.close()
		
	# map overlap between peaks/regions and genomic features:
	if option.mode == "map:dataset" or option.mode == "map:regions":
		
		# define source path:
		if option.mode == "map:dataset":
			sourcepath = peakspath + option.peaks + "/"
		elif option.mode == "map:regions":
			sourcepath = path_dict["binding"] + option.peaks + "/input/"
			
		# define input and output files:
		i_infile = annotationspath + option.infile
		i_tmpfile = datasetpath + option.infile + ".tmp"
		i_colfile = datasetpath + option.infile + "_header.tmp"
		
		f_outfile = datasetpath + outkey + "_compiled"
		s_outfile = summarypath + outkey + "_compiled"
		
		# define header presence
		if option.header == "ON":
			headerFlag = True
		else:
			headerFlag = False
		
		# define annotation headers:
		if option.headerDict == "auto":
			annotationHeader = general.build_header_dict(i_infile)
		elif option.headerDict == "bed":
			annotationHeader = metrn.bedHeader
		else:
			annotationHeader = dict()
			for entry in option.headerDict.split(","):
				key, value = entry.split(":")
				annotationHeader[key] = int(value)
			
		# define feature key setup:
		if option.ids == "standard":
			idComplex = ["chrm", "start", "end", "feature"]
			idOverlap = len(open(i_infile).readline().split("\t")) + 3
		elif option.ids == "feature":
			idComplex = ["feature"]
			idOverlap = len(open(i_infile).readline().split("\t")) + 3
		
		# load annotations:
		print
		print "Loading input annotations..."
		annotationDict = general.build2(i_infile, i=option.target, j=idComplex, x="", mode="matrix", header_dict=annotationHeader, header=headerFlag, separator=":", counter=True)
		
		#k = annotationDict.keys()[0]
		#print k
		#print annotationDict[k]
		#pdb.set_trace()
				
		# define standard targets/queries:
		if option.queries == "standard":
			queries = ['tss', 'five_prime_utr', 'three_prime_utr', 'exon', 'intron', 'protein_coding_primary_transcript']
		elif option.queries == "auto":
			queries = list()
			for query in sorted(annotationDict.keys()):
				if not query in option.exclude.split(","):
					queries.append(query)
			queries = sorted(list(set(queries)))
		elif option.queries != "OFF":
			queries = list()
			for query in option.queries.split(","):
				if not query in option.exclude.split(","):
					queries.append(query) 
		
		if option.header == "ON":
			command = 'grep -v "feature" ' + i_infile + ' > ' + i_tmpfile
			os.system(command)
		
			command = "head -n 1 " + i_infile + ' > ' + i_colfile
			os.system(command)
		
		else:
			command = "cp " + i_infile + " " + i_tmpfile
			os.system(command)
		
		print "Intersecting annotation features with peaks..."
		print
		overlapDict = dict()
		for dataset in os.listdir(sourcepath):
			
			# generate dataset label:
			if option.mode == "map:dataset":
				datasetID = metrn.labelGenerator(option.label, dataset=dataset, mode="label")
			elif option.mode == "map:regions":
				datasetID = dataset.replace(".bed", "")
			
			# rename elements if necessary:
			if option.rename != "OFF":
				for scheme in option.rename.split(","):
					target, replace = scheme.split(":")
					datasetID = datasetID.replace(target, replace)
				
			# store dataset in dictionary:
			overlapDict[datasetID] = dict()
				
			# define peak source and output
			p_infile = sourcepath + dataset
			p_outfile = datasetpath + dataset.replace(".bed", "_intersect.bed")
			
			print "Processing:", datasetID
			if option.fraction == "OFF":
				command = "intersectBed -wo -a " + i_tmpfile + " -b " + p_infile + " > " + p_outfile
				os.system(command)
			else:
				command = "intersectBed -wo -f " + str(option.fraction) + " -a " + i_tmpfile + " -b " + p_infile + " > " + p_outfile
				os.system(command)
			
			# gather annotation overlap peak regions
			overlapBed = general.build2(p_outfile, i=option.target, j=idOverlap, x="", mode="matrix", header_dict=annotationHeader, header=headerFlag, separator=":", counter=True)
			
			#k = overlapBed.keys()[0]
			#print k
			#print overlapBed[k]
			#print p_outfile
			#pdb.set_trace()
			
			# count everything ...
			if option.prioritize == "OFF":
				queryCount = 0
				for query in queries:
					if query in overlapBed:
						overlapDict[datasetID][query] = sum(overlapBed[query].values())
					else:
						overlapDict[datasetID][query] = 0
					queryCount += overlapDict[datasetID][query]
						
				if option.others == "ON":
					totalCount = 0
					for entry in overlapBed:
						totalCount += sum(overlapBed[entry].values())
					overlapDict[datasetID]["others"] = totalCount - queryCount
			
			# ...or preferentially assign counts to classes?
			else:
				prioritize = sorted(overlapBed.keys())
				invertedBed = general.dictinvert(overlapBed, mode="matrix")
				queryCount = 0
				for feature in invertedBed:
					for query in prioritize:
						if not query in overlapDict[datasetID]:
							overlapDict[datasetID][query] = 0
						if query in invertedBed[feature]:
							overlapDict[datasetID][query] += 1
							break
				for query in prioritize:
					queryCount += overlapDict[datasetID][query]
				if option.others == "ON":
					totalCount = len(invertedBed)
					overlapDict[datasetID]["others"] = totalCount - queryCount
			
			# is necessary to determine feature counts?
			if option.elsewhere == "ON":
				featureCount = general.countLines(p_infile)
				overlapDict[datasetID]["others"] = featureCount - queryCount
				
			#if True:
			#	featureCount = general.countLines(p_infile)
			#	invertedBed = general.dictinvert(overlapBed, mode="matrix")
			#	feature = invertedBed.keys()[0]
			#	print featureCount, len(invertedBed), queryCount
			#	print feature
			#	print invertedBed[feature]
			#	pdb.set_trace()
				
		# update queries to include "others" if necessary:
		if option.others == "ON":
			queries.append("others")
			
		# score datasets with highest promoter ratio:
		if option.reference != "OFF":
			referenceDict = dict()
			for datasetID in sorted(overlapDict.keys()):
				values = list()
				for query in queries:
					values.append(overlapDict[datasetID][query])
					if query == option.reference:
						reference = overlapDict[datasetID][query]
				if sum(values) != 0:
					referenceDict[datasetID] = float(reference)/sum(values)
				else:
					referenceDict[datasetID] = 0
			datasetIDs = general.valuesort(referenceDict)
			datasetIDs.reverse()
		elif option.order != "OFF":
			datasetIDs = option.order.split(",")
		else:
			datasetIDs = sorted(overlapDict.keys())
		
		# export per factor-count feature-type summary:
		print
		print "Exporting datasets as ranked by overlap with reference:"
		print
		sv_sumfile = s_outfile.replace("_compiled", "_compiled_dataset_summary_values")
		sn_sumfile = s_outfile.replace("_compiled", "_compiled_dataset_summary_normal")
		sv_output = open(sv_sumfile, "w")
		sn_output = open(sn_sumfile, "w")
		print >>sv_output, "\t".join(["dataset"] + queries)
		print >>sn_output, "\t".join(["dataset"] + queries)
		rank = 1
		for datasetID in datasetIDs:
			print "#" + str(rank), ":", datasetID
			values, normalized = list(), list()
			for query in queries:
				values.append(overlapDict[datasetID][query])
			for value in values:
				if sum(values) > 0:
					normalized.append(float(value)/sum(values))
				else:
					normalized.append(0)
			print >>sv_output, "\t".join(map(str, [datasetID] + values))
			print >>sn_output, "\t".join(map(str, [datasetID] + normalized))
			rank += 1
		sv_output.close()
		sn_output.close()
		
		print
		print "Removing temporary files..."
		command = "rm -rf " + bindingpath + "*.tmp"
		os.system(command)
		
		command = "rm -rf " + datasetpath + "*.tmp"
		os.system(command)
		
		command = "rm -rf " + annotationspath + "*.tmp"
		os.system(command)
		print
		
			
	
	# map peak calls to genomic features:
	elif option.mode == "map:annotation":
		
		# get feature coordinates and names:
		print
		print "loading annotations..."
		feature_dict = bed.build2(annotationspath + option.infile, mode="tree", chrmTree=True, featureMode="feature", header="auto", headerDict="auto", featureExtension=False, filtering=True, filterDict={"feature.type":"cds"}, modeDict={"feature.type":"match"}, testPrint=False)
		
		# define output:	
		f_output = open(mappingpath + outkey + "_w" + str(option.window) + "_m" + str(max_peaks) + "_" + option.infile.replace(".bed","").replace(".txt","") + "_" + option.target,"w")
		print >>f_output, "\t".join(["chrm", "start", "end", "feature", "score", "strand", "strain", "factor", "stage", "institute", "method", "peak", "up", "in", "dn", "tfbs.mark", "tfbs.a", "tfbs.b", "tfbs.c", "tfbs.d"])
		
		# load peak dictionary:
		print "loading peaks..."
		complete_data = open(picklespath + outkey + "_complete.pickle")
		complete_dict = pickle.load(complete_data)
		
		# scan peak call hits to genomic features:
		check = False
		print "mapping peaks to features..."
		print
		network_dict, p = dict(), 0
		for strain in complete_dict:
			for factor in complete_dict[strain]:
				for stage in complete_dict[strain][factor]:

					if not stage in network_dict:
						network_dict[stage] = dict()
					
					for institute in complete_dict[strain][factor][stage]:
						for method in complete_dict[strain][factor][stage][institute]:
							
							if not method in network_dict[stage]:
								network_dict[stage][method] = dict()
							
							p += 1
							print p, "_".join([strain,factor,stage,institute,method]), "(" + str(len(complete_dict[strain][factor][stage][institute][method])) + " peaks)"
							
							for peak in complete_dict[strain][factor][stage][institute][method]:
								feature, chrm, pstart, pend, contribution, hits, density, strand, signal, pvalue, qvalue, point = complete_dict[strain][factor][stage][institute][method][peak]
								
								#if check:
								#	print [chrm, pstart, pend, pscore, strand, signal, pvalue, qvalue, point]
								#	pdb.set_trace()
								
								if chrm in cetrn.chrm2code_dict:
									chrm = cetrn.chrm2code_dict[chrm]
								
								for feature in feature_dict[chrm]:
									fstart = int(feature_dict[chrm][feature]["start"])
									fend = int(feature_dict[chrm][feature]["end"])
									fstrand = feature_dict[chrm][feature]["strand"]
									
									if (pstart >= fstart - window and pstart <= fend + window) or (pend >= fstart - window and pend <= fend + window) or (pstart <= fstart - window and pend >= fend + window):
									
										upstream, inside, dnstream = "-","-","-"
										
										if fstrand == "+":
											upL = fstart - window
											upR = fstart
											dnL = fend
											dnR = fend + window
										else:
											upL = fend
											upR = fend + window
											dnL = fstart - window
											dnR = fstart
										
										inL, inR = sorted([fstart, fend])
										
										if (pstart >= upL and pstart <= upR) or (pend >= upL and pend <= upR) or (pstart <= upL and pend >= upR):
											upstream = "+"
										if (pstart >= inL and pstart <= inR) or (pend >= inL and pend <= inR) or (pstart <= inL and pend >= inR):
											inside = "+"
										if (pstart >= dnL and pstart <= dnR) or (pend >= dnL and pend <= dnR) or (pstart <= dnL and pend >= dnR):
											dnstream = "+"
										
										#if feature == "T07D1.2" and peak == "RepAll_peak_16299":
										#	check = True
									
										print >>f_output, "\t".join(map(str, [chrm, pstart, pend, feature, signal, "+", strain, factor, stage, institute, method, peak, upstream, inside, dnstream, strand, signal, pvalue, qvalue, point]))
						
		# close output file:
		f_output.close()
	
		
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

