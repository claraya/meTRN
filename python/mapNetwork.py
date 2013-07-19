#!/usr/bin/env python
# convert expression data into expressed and not-expressed genes

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
import fasta
import os

from runner import *

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())


""" define a function for mapping IDs """
def geneConverter(organismTag, inpath, i, j, upper=True, nameUpper=True, idUpper=True):
	geneFile = metrn.reference[organismTag]["gene_ids"]
	idLabels = metrn.reference[organismTag][i]
	nameLabels = metrn.reference[organismTag][j]
	id2name_dict, name2id_dict = modencode.idBuild(inpath + geneFile, idLabels, nameLabels, mode="label", header=True, nameUpper=nameUpper, idUpper=idUpper)
	return id2name_dict, name2id_dict


""" define a function that loads the orthologs between species from file """
def orthologBuilder(aspecies, bspecies, infile):

	speciesDict = {
		"celegans" : "ce",
		"dmel" : "dm",
		"human" : "hs",
		"mouse" : "mm"
		}
	
	removalDict = {
		"ce" : "Cel-",
		"dm" : "dmel_",
		"hs" : "human_",
		"mm" : "mouse_"
		}
	
	indata = open(infile)
	inline = indata.readline()
	iDict, jDict = dict(), dict()
	while inline:
		index, i, j, iGene, jGene, iValue, jValue = inline.strip().split("\t")
		if i in speciesDict and j in speciesDict:
			i, j = speciesDict[i], speciesDict[j]
			if i in [aspecies, bspecies] and j in [aspecies, bspecies]:
				if i == aspecies and j == bspecies:
					iGene, jGene = iGene.lstrip(removalDict[i]), jGene.lstrip(removalDict[j])
				elif j == aspecies and i == bspecies:
					jGene, iGene = iGene.lstrip(removalDict[i]), jGene.lstrip(removalDict[j])
				if not iGene in iDict:
					iDict[iGene] = list()
				if not jGene in jDict:
					jDict[jGene] = list()
				iDict[iGene].append(jGene)
				jDict[jGene].append(iGene)
		inline = indata.readline()
	return iDict, jDict


""" define a function that assigns features to others and returns a dictionary... """
def featureMapper(featurefile, regionfile, outputfile, fraction, target, headerDict="auto", header="OFF", overwrite="ON"):
	
	# get output part names:
	outputname = outputfile.split("/")[len(outputfile.split("/"))-1]
	outputpath = "/".join(outputfile.split("/")[:len(outputfile.split("/"))-1]) + "/"
	
	# is there a header on the region file?
	if overwrite == "ON" and header == "ON":
		
		headfile = regionfile.replace(".bed", ".head")
		tempfile = regionfile.replace(".bed", ".tmp")
		command = "cp " + regionfile + " " +  tempfile
		os.system(command)
		
		command = "head -n 1 " + regionfile + ' > ' + headfile
		os.system(command)
		
		command = 'grep -v "feature" ' + regionfile + ' > ' + tempfile
		os.system(command)
	
	# define header presence
	if header == "ON":
		headerFlag = True
	else:
		headerFlag = False
		
	# intersect BED files:
	if overwrite == "ON" or not outputname in os.listdir(outputpath):
		
		if fraction == "OFF":
			command = "intersectBed -wo -a " + regionfile + " -b " + featurefile + " > " + outputfile
			os.system(command)
		else:
			command = "intersectBed -wo -f " + str(fraction) + " -a " + regionfile + " -b " + featurefile + " > " + outputfile
			os.system(command)
	
	# define annotation headers:
	if headerDict == "auto":
		annotationHeader = general.build_header_dict(regionfile)
	elif headerDict == "bed":
		annotationHeader = metrn.bedHeader
	else:
		annotationHeader = dict()
		for entry in headerDict.split(","):
			key, value = entry.split(":")
			annotationHeader[key] = int(value)
			
	# define feature key setup:
	idOverlap = len(open(regionfile).readline().split("\t")) + 3
	
	# gather annotation overlap peak regions:
	overlapBed = general.build2(outputfile, i=idOverlap, j=target, x="", mode="matrix", header_dict=annotationHeader, header=headerFlag, separator=":", counter=True)
	return overlapBed


""" define a function that converts a local user path to SCG3 or GS server paths """
def serverPath(inpath, server="ON"):
	if server == "ON":
		return inpath.replace("/Users/claraya/", "/srv/gs1/projects/snyder/claraya/")
	elif server == "GS":
		return inpath.replace("/Users/claraya/", "/net/fields/vol1/home/araya/")

				
def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "path from script to files")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "Basename for target peaks", default="OFF")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "launch")
	parser.add_option("--infile", action = "store", type = "string", dest = "infile", help = "input expression file", default="OFF")
	parser.add_option("--source", action = "store", type = "string", dest = "source", help = "folder location of input files", default="OFF")
	parser.add_option("--target", action = "store", type = "string", dest = "target", help = "How to format dataset labels...", default="factor.context")
	parser.add_option("--mapping", action = "store", type = "string", dest = "mapping", help = "Dataset mapping guide...", default="OFF")
	parser.add_option("--exclude", action = "store", type = "string", dest = "exclude", help = "Targets to exclude", default="OFF")
	parser.add_option("--rename", action = "store", type = "string", dest = "rename", help = "Targets to rename. Comma-separated list of 'target:replacement' pairs to search and replace.", default="OFF")
	parser.add_option("--tip", action = "store", type = "string", dest = "tip", help = "Are these TIP prediction files?", default="OFF")
	parser.add_option("--fraction", action = "store", type = "string", dest = "fraction", help = "Fractional overlap required", default="0.1")
	parser.add_option("--nuclear", action = "store", type = "string", dest = "nuclear", help = "Peaks are only nuclear?", default="ON")
	parser.add_option("--species", action = "store", type = "string", dest = "species", help = "Species to be compared; comma-separated", default="OFF")
	parser.add_option("--genes", action = "store", type = "string", dest = "genes", help = "reference gene file", default="OFF")
	parser.add_option("--transcripts", action = "store", type = "string", dest = "transcripts", help = "reference transcript file")
	parser.add_option("--orthology", action = "store", type = "string", dest = "orthology", help = "Use 'direct' or 'family' orthologs?", default="direct")
	parser.add_option("--nametag", action = "store", type = "string", dest = "nametag", help = "Orthology nametag: nametagHsCe", default="ortho")
	parser.add_option("--commonNames", action = "store", type = "string", dest = "commonNames", help = "Grab common names file?", default="ON")
	parser.add_option("--familyFiles", action = "store", type = "string", dest = "familyFiles", help = "Grab cleaned files?", default="formatted")
	parser.add_option("--coord", action = "store", type = "string", dest = "coord", help = "reference coordinates to export: RNA or TSS", default="RNA")
	parser.add_option("--cutoff", action = "store", type = "float", dest = "cutoff", help = "expression cutoff", default=0.05)
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "name for domain", default="OFF")
	parser.add_option("--strip", action = "store", type = "string", dest = "strip", help = "Remove RNA transcript periods (last)?", default="OFF")
	parser.add_option("--a", action = "store", type = "string", dest = "a", help = "Input A network", default="OFF")
	parser.add_option("--b", action = "store", type = "string", dest = "b", help = "Input B network", default="OFF")
	parser.add_option("--overwrite", action = "store", type = "string", dest = "overwrite", help = "Overwrite stuff?", default="OFF")
	parser.add_option("--threads", action = "store", type = "string", dest = "threads", help = "multiprocessing threads", default="1")
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
	orthologspath = path_dict["orthologs"]
	peakspath = path_dict["peaks"]
	gopath = path_dict["go"]
	hotpath = path_dict["hot"]
	networkpath = path_dict["network"]
	qsubpath = path_dict["qsub"]
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
	
	# load gene ID dictionaries:
	id2name_dict, name2id_dict = modencode.idBuild(inpath + metrn.reference[organismTag]["gene_ids"], "Sequence Name (Gene)", "Gene Public Name", mode="label", header=True, idUpper=True, nameUpper=True)
	
	# define expression cutoff handle:
	expression = "%.0e" % (float(option.cutoff))
	
	# define flagging dicts:
	flagDict = dict()
	flagDict["nameUpper"] = { "ce":True, "dm":True, "hs":True }
	flagDict["idUpper"] = { "ce":False, "dm":True, "hs":True }
	
	# build network mode:
	if option.mode == "build":
		
		# determine path of input files:
		networkpath = networkpath + option.peaks + "/" + option.name + "/"
		targetspath = networkpath + "targets/"
		mappingpath = networkpath + "mapping/"
		general.pathGenerator(targetspath)
		general.pathGenerator(mappingpath)
	
		# update peaks path:
		peakspath = peakspath + option.peaks + "/"
		
		# process TIP predictions:
		if option.tip == "ON":
		
			# load mapping:
			if option.mapping != "OFF":
				mappingDict = general.build2(annotationspath + option.mapping, id_column="filename")
			
			# load binding targets:
			print 
			print "Loading target predictions..."
			mapDict = fasta.buildfile(path_dict[option.source] + option.infile)
			skipped = list()
			
			# prebuild mapping if necessary:
			if option.mapping != "OFF":
				
				peak2dataDict, data2peakDict = dict(), dict()
				for dataset in sorted(mapDict.keys()):
					dataset  = dataset.replace(".bam", "").replace("AlnRep0", "Aln").replace("AlnRep1", "Aln")
				
					# rename datasets if indicated:
					if option.rename != "OFF":
						for scheme in option.rename.split(","):
							target, replace = scheme.split(":")
							dataset = dataset.replace(target, replace)
					
					if not dataset in mappingDict:
						skipped.append(dataset)
					else:
						strain, factor, context, institute, method = mappingDict[dataset]["strain"], mappingDict[dataset]["factor"], mappingDict[dataset]["context"], mappingDict[dataset]["institute"], mappingDict[dataset]["method"]
						filelabel = organismTag + "_" + "_".join([strain, factor, context, institute])
						for peakfile in os.listdir(peakspath):
							if filelabel in peakfile:
								if not peakfile in peak2dataDict:
									peak2dataDict[peakfile] = dict()
								if not dataset in data2peakDict:
									data2peakDict[dataset] = dict()
								peak2dataDict[peakfile][dataset] = 1
								data2peakDict[dataset][peakfile] = 1
				
				# check 1-to-1 mapping between peak files and datasets...
				for dataset in data2peakDict:
					if len(data2peakDict[dataset]) > 1:
						raise "Error: Multiple peak-files matched!", dataset
		
		# assign targets by distance:
		else:
			
			# specify the target regions file:
			i_infile = path_dict[option.source] + option.infile
			
			# process peaks:
			mapDict = dict()
			skipped = list()
			for peakfile in os.listdir(peakspath):
				
				# generate dataset key:
				dataset = peakfile.replace("_peaks.bed", "")
				print "Processing:", dataset
				
				# define peak source and output
				p_infile = peakspath + peakfile
				p_outfile = mappingpath + peakfile.replace(".bed", ".tmp")
				
				# assign peaks to features:
				peakDict = featureMapper(p_infile, i_infile, p_outfile, fraction=option.fraction, target="feature", headerDict="bed", header="OFF", overwrite=option.overwrite)
				
				# map dataset to features:
				mapDict[dataset] = list()
				for peak in peakDict:
					for feature in peakDict[peak]:
						mapDict[dataset].append(feature)
				mapDict[dataset] = sorted(list(set(mapDict[dataset])))
			
			#key = mapDict.keys()[0]
			#print key 
			#print mapDict[key]
			#pdb.set_trace()
			
		# transfer binding targets:
		print "Transfering target predictions..."
		networkDict = dict()
		s, t = 0, 0
		for dataset in sorted(mapDict.keys()):
		
			# get targets:
			if option.tip == "ON":
				targets = general.clean(mapDict[dataset].split("\t"), "")
				dataset  = dataset.replace(".bam", "").replace("AlnRep0", "Aln").replace("AlnRep1", "Aln")
			else:
				targets = general.clean(mapDict[dataset], "")
			
			# rename datasets if indicated:
			if option.rename != "OFF":
				for scheme in option.rename.split(","):
					target, replace = scheme.split(":")
					dataset = dataset.replace(target, replace)
			
			# substitute signal filename for peaks:
			process = False
			if option.tip == "ON" and option.mapping == "OFF":
				dataset = organismTag + "_" + dataset
				peakfile = dataset.replace(".fc.signal", "_peaks.bed")
				dataset = dataset.replace(".fc.signal", "")
				process = True
				t += 1
			
			# generate peak file names from mapping:
			elif option.tip == "ON" and option.mapping != "OFF":
				if dataset in data2peakDict:
					peakfile = data2peakDict[dataset].keys()[0]
					process = True
					t += 1
					
			# store TIP target assignments:
			if option.tip == "ON" and peakfile in os.listdir(peakspath) and targets != list() and process:
				networkDict[dataset] = targets
				s += 1
			
			# store overlap target assignments:
			elif option.tip == "OFF":
				networkDict[dataset] = targets
				s += 1
	
		
		# define output files:
		n_output = open(targetspath + "mapnetwork_build_" + option.peaks + "_" + option.name + ".txt", "w")
		print >>n_output, "\t".join(["dataset", "target", "organism", "strain", "factor", "context", "institute", "method"])
		
		# export network targets:
		k = 0
		print "Exporting target predictions..."
		for dataset in sorted(networkDict.keys()):
			
			if option.mapping == "OFF":
				datasetID = metrn.labelGenerator(target=option.target, mode="label", dataset=dataset)
				organism, strain, factor, context, institute, method = metrn.labelComponents(dataset)
			else:
				organism, strain, factor, context, institute, method = organismTag, mappingDict[dataset]["strain"], mappingDict[dataset]["factor"], mappingDict[dataset]["context"], mappingDict[dataset]["institute"], mappingDict[dataset]["method"]
				datasetID = "_".join([organism, strain, factor, context, institute, method])
				datasetID = metrn.labelGenerator(target=option.target, mode="label", dataset=datasetID)
				
			for target in sorted(networkDict[dataset]):
				print >>n_output, "\t".join([datasetID, target, organism, strain, factor, context, institute, method])
				k += 1
		
		# close output files:
		n_output.close()
		
		print "Loaded regulator-target interactions:", k
		print "Input peak call files:", len(os.listdir(peakspath))
		print "Input signal sources:", t
		print "Stored signal sources:", s
		print "Skipped sources:", len(skipped)
		print
		
	# resolve network mode:
	elif option.mode == "resolve":
	
		# determine path of input files:
		networkpath = networkpath + option.peaks + "/" + option.name + "/"
		targetspath = networkpath + "targets/"
		mappingpath = networkpath + "mapping/"
		general.pathGenerator(targetspath)
		general.pathGenerator(mappingpath)
		
		# load gene ids...
		rna2gen_dict, gen2rna_dict = modencode.idBuild(inpath + metrn.reference[organismTag]["gene_ids"], metrn.reference[organismTag]["rna.link"], metrn.reference[organismTag]["gene.link"], mode="label", header=True, idUpper=flagDict["idUpper"], nameUpper=flagDict["nameUpper"])
		rna2pro_dict, pro2rna_dict = modencode.idBuild(inpath + metrn.reference[organismTag]["gene_ids"], metrn.reference[organismTag]["rna.link"], metrn.reference[organismTag]["protein.link"], mode="label", header=True, idUpper=flagDict["idUpper"], nameUpper=flagDict["nameUpper"])
		rna2sym_dict, sym2rna_dict = modencode.idBuild(inpath + metrn.reference[organismTag]["gene_ids"], metrn.reference[organismTag]["rna.link"], metrn.reference[organismTag]["symbol.link"], mode="label", header=True, idUpper=flagDict["idUpper"], nameUpper=flagDict["nameUpper"])
		
		# load and process network:
		networkfile = targetspath + "mapnetwork_build_" + option.peaks + "_" + option.name + ".txt"
		n_output = open(targetspath + "mapnetwork_resolve_" + option.peaks + "_" + option.name + ".txt", "w")
		print >>n_output, "\t".join(["dataset", "target", "organism", "strain", "factor", "context", "institute", "method", "symbol", "gene", "protein"])
		inlines = open(networkfile).readlines()
		inlines.pop(0)
		for inline in inlines:
			if inline.strip() != "":
				dataset, target = inline.strip().split("\t")[:2]
				if option.strip == "ON":
					length = len(target.split("."))
					query = ".".join(target.split(".")[:length-1])
				else:
					query = str(target)
				
				if query in rna2gen_dict and query in rna2pro_dict:
					gene = rna2gen_dict[query]
					protein = rna2pro_dict[query]
					symbol = rna2sym_dict[query]
					print >>n_output, inline.strip() + "\t" + symbol + "\t" + gene + "\t" + protein
		n_output.close()
		
	
	# cascade network analysis:
	elif option.mode == "cascade":
	
		# determine path of input files:
		networkpath = networkpath + option.peaks + "/" + option.name + "/"
		targetspath = networkpath + "targets/"
		mappingpath = networkpath + "mapping/"
		cascadepath = networkpath + "cascade/"
		general.pathGenerator(targetspath)
		general.pathGenerator(mappingpath)
		general.pathGenerator(cascadepath)
		
		print 
		print "Loading network data..."
		networkfile = targetspath + "mapnetwork_resolve_" + option.peaks + "_" + option.name + ".txt"
		networkDict = general.build2(networkfile, i="factor", j="symbol", x="context", mode="matrix", skip=True, verbose=False)
		networkHead = open(networkfile).readline()
		
		# generate line-matching dictionary:
		linesDict = dict()
		headerDict = general.build_header_dict(networkfile)
		inlines = open(networkfile).readlines()
		inlines.pop(0)
		for inline in inlines:
			initems = inline.strip().split("\t")
			factor, gene = initems[headerDict["factor"]], initems[headerDict["symbol"]]
			if not factor in linesDict:
				linesDict[factor] = dict()
			linesDict[factor][gene] = inline.strip()
			
		# initiate output files:
		k_output = open(cascadepath + "mapnetwork_cascade_" + option.peaks + "_" + option.name + "_all.txt", "w")
		m_output = open(cascadepath + "mapnetwork_cascade_" + option.peaks + "_" + option.name + "_com.txt", "w")
		z_output = open(cascadepath + "mapnetwork_cascade_" + option.peaks + "_" + option.name + "_reg.txt", "w")
		print >>k_output, networkHead.strip()
		print >>m_output, networkHead.strip()
		print >>z_output, networkHead.strip()
		
		print "Scanning overlaps..."
		k, m, z = 0, 0, 0
		cascadeDict, regularDict = dict(), dict()
		cellSets = os.listdir(cellspath + "cellset/" + option.mapping)
		for factor in networkDict:
			if factor in cellSets:
				factorCells = open(cellspath + "cellset/" + option.mapping + "/" + factor).read().split("\n")
				for gene in networkDict[factor]:
					if gene in cellSets:
						geneCells = open(cellspath + "cellset/" + option.mapping + "/" + gene).read().split("\n")
						overlapCells = set(factorCells).intersection(set(geneCells))
						
						if len(overlapCells) >= float(option.fraction)*len(geneCells) and len(factorCells) >= len(geneCells):
							#print factor, gene, len(networkDict[factor])
							if not factor in cascadeDict:
								cascadeDict[factor] = dict()
							cascadeDict[factor][gene] = networkDict[factor][gene]
							print >>m_output, linesDict[factor][gene]
							m += 1
							
							if gene in networkDict:
								#print factor, gene, len(networkDict[factor])
								if not factor in regularDict:
									regularDict[factor] = dict()
								regularDict[factor][gene] = networkDict[factor][gene]
								print >>z_output, linesDict[factor][gene]
								z += 1
						
						#print factor, gene, len(networkDict[factor])
						print >>k_output, linesDict[factor][gene]
						k += 1
		
		# close output files:
		k_output.close()
		m_output.close()
		z_output.close()
		
		print
		print "Possibles cascade interactions:", k 
		print "Supported cascade interactions:", m
		print "Regulator cascade interactions:", z
		print
		
	
	# find common regulatory interactions mode:
	elif option.mode == "commons":
	
		# Note: This method searches for each factor, say UNC-62, all of the name-related targets across contexts. 
		# As such, this method will find that UNC-62 in the L4 stage binds to VIT-1, VIT-3, VIT-4 and VIT-5.
	
		# determine path of input files:
		networkpath = networkpath + option.peaks + "/" + option.name + "/"
		targetspath = networkpath + "targets/"
		mappingpath = networkpath + "mapping/"
		summarypath = networkpath + "commons/summary/"
		datasetpath = networkpath + "commons/dataset/"
		general.pathGenerator(targetspath)
		general.pathGenerator(mappingpath)
		general.pathGenerator(summarypath)
		general.pathGenerator(datasetpath)
		
		# load input network:
		print
		print "Loading regulatory network..."
		networkfile = targetspath + "mapnetwork_resolve_" + option.peaks + "_" + option.name + ".txt"
		networkDict = general.build2(networkfile, i="dataset", j="symbol", x="context", mode="matrix", skip=True, verbose=False)
		
		# load protein interaction network:
		print "Loading protein interaction network..."
		if option.a != "OFF":
			i, j, x = option.a, option.b, option.b
			ppiDict = general.build2(extraspath + option.mapping, i=i, j=j, x=x, mode="matrix")
			
			# establish cutoff handle:
			cutoffHandle = "_cut" + str(int(option.cutoff))
			
			# prepare integrated network output:
			integratedfile = summarypath + "mapnetwork_commons_" + option.peaks + "_" + option.name + cutoffHandle + "_network.txt"
			i_output = open(integratedfile, "w")
			print >>i_output, "\t".join(["source","type","target","factor","stage"])
			
			# transfer interaction network:
			proteins, missing = list(), 0
			for source in ppiDict:
				for target in ppiDict[source]:
					proteins.append(source)
					proteins.append(target)
					if source in id2name_dict and target in id2name_dict:
						source, target = id2name_dict[source], id2name_dict[target]
						print >>i_output, "\t".join([source, "pp", target, source, "interaction"])
			
			# tally missing proteins:
			proteins = sorted(list(set(proteins)))
			for protein in proteins:
				if not protein in id2name_dict:
					missing += 1
			print "Protein missing:", round(100*float(missing)/len(proteins), 2), "%"
			print
		
		# find name-related targets network:
		print "Finding name-related targets..."
		commonsDict, sizeDict = dict(), dict()
		for dataset in networkDict:
			for target in networkDict[dataset]:
				basename = target.split("-")[0]
				stage = networkDict[dataset][target]
				factor = dataset.replace("." + stage, "")
				if not dataset in commonsDict:
					commonsDict[dataset] = dict()
					sizeDict[dataset] = dict()
				if not basename in commonsDict[dataset]:
					commonsDict[dataset][basename] = list()
				commonsDict[dataset][basename].append(target)
				sizeDict[dataset][basename] = len(commonsDict[dataset][basename])
		
		# transfer regulatory network, applying basename cutoffs...
		for dataset in sorted(commonsDict.keys()):
			for basename in commonsDict[dataset]:
				if sizeDict[dataset][basename] >= int(option.cutoff):
					for target in commonsDict[dataset][basename]:
						stage = networkDict[dataset][target]
						factor = dataset.replace("." + stage, "")
						print >>i_output, "\t".join([dataset, "pd", target, factor, stage])
		
		# close integrated network output:
		if option.a != "OFF":
			i_output.close()
		
		# find & export interesting candidates:
		summaryfile = summarypath + "mapnetwork_commons_" + option.peaks + "_" + option.name + cutoffHandle + "_summary.txt"
		f_output = open(summaryfile, "w")
		print >>f_output, "\t".join(["dataset", "target.class", "target.count", "target.ids"])
		print "Exporting summary and aggregates..."
		for dataset in sorted(commonsDict.keys()):
			
			#print "Processing:", dataset
			datasetfile = datasetpath + "mapnetwork_commons_" + option.peaks + "_" + option.name + cutoffHandle + "_" + dataset + ".txt"
			d_output = open(datasetfile, "w")
			print >>d_output, "\t".join(["source", "target"])
			for basename in sizeDict[dataset]:
				if sizeDict[dataset][basename] >= int(option.cutoff):
					#print dataset, sizeDict[dataset][basename], ",".join(commonsDict[dataset][basename])
					print >>f_output, "\t".join(map(str, [dataset, basename, sizeDict[dataset][basename], ",".join(commonsDict[dataset][basename])]))
					print >>d_output, "\t".join([dataset, basename])
					for target in commonsDict[dataset][basename]:
						print >>d_output, "\t".join([basename, target])
			d_output.close()
			#print
		
		# close output summary file:
		f_output.close()
		print
		
	
	# hybridize networks mode:
	elif option.mode == "hybrid":
	
		# determine path of input files:
		comparisonpath = networkpath + "/comparison/"
		general.pathGenerator(comparisonpath)
		
		# collect species names:
		aspecies, bspecies = option.species.split(",")
		
		# define input network files:
		ainfile = networkpath + option.a + "/targets/mapnetwork_resolve_" + option.a.replace("/", "_") + ".txt"
		binfile = networkpath + option.b + "/targets/mapnetwork_resolve_" + option.b.replace("/", "_") + ".txt"
		
		# identify target specie and comparison species:
		speciesTags = option.species.split(",")
		
		# define orthology path:
		if option.orthology == "direct":
			orthologypath = orthologspath + "orthologs/"
		elif option.orthology == "family":
			orthologypath = orthologspath + "families/"
		
		# generate ortholog tag name:
		orthologTag = metrn.orthologLabel(option.organism, speciesTags)
		
		# generate orthology dictionary:
		ortholog_dict = metrn.orthologBuilder(speciesTags, path=orthologypath, orthology=option.orthology, commonNames=option.commonNames, familyFiles=option.familyFiles, verbose="OFF")
		
		# target specie orthologs:
		if option.orthology == "direct":
			aOrthologs = metrn.orthologFinder(aspecies, speciesTags, path=orthologspath + "orthologs/", orthology=option.orthology, commonNames=option.commonNames)
			bOrthologs = metrn.orthologFinder(bspecies, speciesTags, path=orthologspath + "orthologs/", orthology=option.orthology, commonNames=option.commonNames)
			
		elif option.orthology == "family":
			aOrthologs = metrn.orthologFinder(aspecies, speciesTags, path=orthologspath + "families/", orthology=option.orthology, familyFiles=option.familyFiles)
			bOrthologs = metrn.orthologFinder(bspecies, speciesTags, path=orthologspath + "families/", orthology=option.orthology, familyFiles=option.familyFiles)
		
		print
		print "Orthologs for " + aspecies.upper() + ":", str(len(aOrthologs))
		print "Orthologs for " + bspecies.upper() + ":", str(len(bOrthologs))
		print
		
		# load ortholog and gene IDs:
		print "Loading ortholog identifiers..."
		aid2name_dict, aname2id_dict = geneConverter(aspecies, inpath=inpath, i="orth.link", j="symbol.link", nameUpper=flagDict["nameUpper"][aspecies], idUpper=flagDict["idUpper"][aspecies])
		bid2name_dict, bname2id_dict = geneConverter(bspecies, inpath=inpath, i="orth.link", j="symbol.link", nameUpper=flagDict["nameUpper"][bspecies], idUpper=flagDict["idUpper"][bspecies])
		
		#print btarget2id_dict.keys()[:5]
		#print bid2target_dict.keys()[:5]
		#pdb.set_trace()
		
		# load network files:
		print "Loading species networks..."
		aNetwork = general.build2(ainfile, i="dataset", j="symbol", x="factor", mode="matrix", skip=True, verbose=False)
		bNetwork = general.build2(binfile, i="dataset", j="symbol", x="factor", mode="matrix", skip=True, verbose=False)
		
		# load ortholog mappings:
		print "Mapping orthologs between species..."
		aMapping, bMapping = orthologBuilder(aspecies, bspecies, extraspath + option.infile)
		
		# scan network overlaps (make reverse networks):
		# note: these network have gene-name keys but the targets are in list format...
		print "Identifying network factors..."
		aRevwork = dict()
		for aDataset in aNetwork:
			for aTarget in aNetwork[aDataset]:
				aGene = aNetwork[aDataset][aTarget]
				if not aGene in aRevwork:
					aRevwork[aGene] = dict()
				if not aDataset in aRevwork[aGene]:
					aRevwork[aGene][aDataset] = list()
				aRevwork[aGene][aDataset].append(aTarget)
		
		bRevwork = dict()
		for bDataset in bNetwork:
			for bTarget in bNetwork[bDataset]:
				bGene = bNetwork[bDataset][bTarget]
				if not bGene in bRevwork:
					bRevwork[bGene] = dict()
				if not bDataset in bRevwork[bGene]:
					bRevwork[bGene][bDataset] = list()
				bRevwork[bGene][bDataset].append(bTarget)
		
		t1, l1, l2, l3, l4 = list(), list(), list(), list(), list()
		for aGene in aOrthologs:
			aFactor = aname2id_dict[aGene]
			if aGene in aRevwork:
				t1.append(aGene)
				
				# fix cases where the correct (mapping) ID isn't found:
				if aFactor in aMapping:
					l1.append(aGene)
				else:
					aMatches = list()
					for aQuery in aMapping:
						if aQuery in aid2name_dict and aGene == aid2name_dict[aQuery]:
							aMatches.append(aQuery)
					if aMatches != list():
						aFactor = aMatches[0]	
						if len(aMatches) > 1:
							print "Caution: More than matching ID found for", aGene, aFactor, aMatches
							pdb.set_trace()
						l1.append(aGene)
					else:
						l2.append(aGene)
				
				# capture dataset global targets:
				aGlobal = list()
				for aDataset in aRevwork[aGene]:
					aGlobal.extend(aRevwork[aGene][aDataset])
				aGlobal = sorted(list(set(aGlobal)))
				
				# check ortholog networks:
				for bGene in ortholog_dict[aspecies][aGene][bspecies]:
					if bGene in bRevwork:
						
						# capture ortholog global targets:
						bGlobal = list()
						for bDataset in bRevwork[bGene]:
							bGlobal.extend(bRevwork[bGene][bDataset])
						bGlobal = sorted(list(set(bGlobal)))
				
						# process dataset targets:
						for aDataset in aRevwork[aGene]:
							aTargets = aRevwork[aGene][aDataset]
							aMatches = list()
							aQueries = list()
							for aTarget in aTargets:
								if aTarget in aname2id_dict:
									aMatches.append(aname2id_dict[aTarget])
									if aname2id_dict[aTarget] in aMapping:
										aQueries.extend(aMapping[aname2id_dict[aTarget]])
							aQueries = sorted(list(set(aQueries)))
							
							# process ortholog targets:
							for bDataset in bRevwork[bGene]:
								bTargets = bRevwork[bGene][bDataset]
								bMatches = list()
								bQueries = list()
								for bTarget in bTargets:
									if bTarget in bname2id_dict:
										bMatches.append(bname2id_dict[bTarget])
										if bname2id_dict[bTarget] in bMapping:
											bQueries.extend(bMapping[bname2id_dict[bTarget]])
								bQueries = sorted(list(set(bQueries)))
								
								# determine overlaps:
								aOverlap = set(aQueries).intersection(set(bMatches))
								bOverlap = set(bQueries).intersection(set(aMatches))
								
								if len(aQueries) != 0:
									aQueryFraction = float(len(aOverlap))/len(aQueries)
									aUnionFraction = float(len(aOverlap))/len(set(aQueries).union(set(bMatches)))
								else:
									aQueryFraction = 0
									aUnionFraction = 0
								
								if len(bQueries) != 0:
									bQueryFraction = float(len(bOverlap))/len(bQueries)
									bUnionFraction = float(len(bOverlap))/len(set(bQueries).union(set(aMatches)))
								else:
									bQueryFraction = 0
									bUnionFraction = 0
								
								#print aGene, len(aTargets)
								#print bGene, len(bTargets)
								#print len(aQueries), len(bQueries)
								#print len(aOverlap), len(bOverlap)
								#pdb.set_trace()
								
								# export data:
								output = [aDataset, bDataset, aGene, bGene, len(aMatches), len(bMatches), len(aQueries), len(bQueries), len(aOverlap), round(aQueryFraction, 3), round(aUnionFraction, 3), len(bOverlap), round(bQueryFraction, 3), round(bUnionFraction, 3)]
								print "\t".join(map(str, output))
											
						#	print aGene, aDataset, aFactor, bFactor
						#	pdb.set_trace()
					
					
				"""
				for bFactor in aMapping[aFactor]:
					if bFactor in bid2name_dict:
						bGene = bid2name_dict[bFactor]
						if bGene in bOrthologs:
								
								# now we have two orthologs, in their networks:
								if bGene in bRevwork:
									
									#print aGene, aFactor, bGene, bFactor
									#pdb.set_trace()
									
									for aDataset in aRevwork[aGene]:
										aTargets = aRevwork[aGene][aDataset]
										aTargets = sorted(list(set(aTargets)))
										aMatchs, aMissed = list(), list()
										for aTarget in aTargets:
											if aTarget in atarget2id_dict:
												aMatchs.append(atarget2id_dict[aTarget])
											else:
												aMissed.append(aTarget)
										aMapped, aUnmaps = list(), list()
										for aMatch in aMatchs:
											if aMatch in aMapping:
												aMapped.extend(aMapping[aMatch])
											else:
												aUnmaps.append(aMatch)
										aMatchs = sorted(list(set(aMatchs)))
										aMissed = sorted(list(set(aMissed)))
										aMapped = sorted(list(set(aMapped)))
										aUnmaps = sorted(list(set(aUnmaps)))
										
										for bDataset in bRevwork[bGene]:
											bTargets = bRevwork[bGene][bDataset]
											bTargets = sorted(list(set(bTargets)))
											bMatchs, bMissed = list(), list()
											for bTarget in bTargets:
												if bTarget in btarget2id_dict:
													bMatchs.append(btarget2id_dict[bTarget])
												else:
													bMissed.append(bTarget)
											bMatchs = sorted(list(set(bMatchs)))
											bMissed = sorted(list(set(bMissed)))
											
											Overlap = sorted(list(set(aMapped).intersection(set(bMatchs))))
											Fraction = float(len(Overlap))/len(aMapped)
											Normaled = float(len(Overlap))/len(set(aMapped).union(set(bMatchs)))
											#print aDataset, len(aTargets), len(aMatchs), len(aMapped), len(aUnmaps)
											#print bDataset, len(bTargets), len(bMatchs), len(bMissed)
											#print Overlap
											#pdb.set_trace()
											
											print aDataset, bDataset, aGene, bGene, len(aMatchs), len(bMatchs), len(aMapped), len(Overlap), round(Fraction, 3), round(Normaled, 3)
				"""
					
		"""									
							else:
								l4.append(bGene)
						else:
							l3.append(bFactor)
				else:
					l2.append(aFactor)
			else:
				l1.append(aGene)
		"""
		
		print "Target orthologs:", t1
		print "Orthologs in network:", len(l1)
		print "Orthologs not mapped:", len(l2)
		print "...match name missing:", len(l3)
		print "...match not assayed:", len(l4)
		
		#for aDataset in aRevwork[aGene]:
			#for aTarget in aNetwork[aDataset]:
			#	aFactor = aNetwork[aDataset][aTarget]
			#	if aFactor in aname2id_dict:
			#		aGene = aname2id_dict[aFactor]
			#		aTarget = aid2target_dict[aGene]
			#		if aGene in aOrthologs:
			#			for bGene in aOrthologs[aGene]:
			#				print aDataset, aFactor, aGene, bGene
			#				pdb.set_trace()
		
		#a1 = aNetwork.keys()[0]
		#a2 = aNetwork[a1].keys()[0]
		#print a1, a2, aNetwork[a1][a2]
		#pdb.set_trace()
		
	
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

