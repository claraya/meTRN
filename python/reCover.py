#!/usr/bin/env python
# grab subfolder pictures!

import sys
import time
import optparse
import general
import os
import pdb
import re

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())

def grabFiletype(filename, cutoff=5):
	fileparts = filename.split(".")
	if len(fileparts) == 1:
		return "NA"
	else:
		filetype = fileparts[len(fileparts)-1]
		if (re.match("^[A-Za-z0-9_-]*$", filetype)) and len(filetype) <= cutoff:
			return "." + filetype
		else:
			return "NA"

def listFiles(dir, output=list()):
	basedir = dir
	#print "Files in ", os.path.abspath(dir), ": "
	subdirlist = []
	for item in os.listdir(dir):
		if os.path.isfile(os.path.join(basedir,item)):
			output.append(os.path.join(basedir,item))
		else:
			subdirlist.append(os.path.join(basedir, item))
	for subdir in subdirlist:
		listFiles(subdir, output)
	return sorted(list(set(output)))
        
def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "path from script to files")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "processing mode")
	parser.add_option("--output", action = "store", type = "string", dest = "output", help = "output path for files")
	parser.add_option("--include", action = "store", type = "string", dest = "include", help = "include files with these flag in name/path; comma-separated", default="OFF")
	parser.add_option("--exclude", action = "store", type = "string", dest = "exclude", help = "exclude files with these flag in name/path; comma-separated", default="OFF")
	parser.add_option("--filetype", action = "store", type = "string", dest = "filetype", help = "target filetypes; comma-separated", default="OFF")
	parser.add_option("--subfolders", action = "store", type = "string", dest = "subfolders", help = "search subfolders only?", default="OFF")
	parser.add_option("--place", action = "store", type = "string", dest = "place", help = "subfolder, subpath, or OFF", default="OFF")
	parser.add_option("--compare", action = "store", type = "string", dest = "compare", help = "comparison path for files")
	parser.add_option("--extendName", action = "store", type = "string", dest = "extendName", help = "encode path into filename?", default="OFF")
	parser.add_option("--slashReplace", action = "store", type = "string", dest = "slashReplace", help = "symbol to replace dashes in filename extension", default="-")
	parser.add_option("--spaceReplace", action = "store", type = "string", dest = "spaceReplace", help = "symbol to replace spaces in filename extension", default="_")
	(option, args) = parser.parse_args()
	
	# image search mode:
	if option.mode == "images":
		
		# load folders:
		contents, filtered, filetypes =  listFiles(option.compare, output=list()), list(), list()
		temppath = "/Users/claraya/Desktop/TMP/"
		general.pathGenerator(temppath)
		
		for infile in contents:
			if option.include in infile and "gallery" in infile.lower():
				command = "cp " + infile + " " + temppath
				os.system(command)
				print infile
	
	
	# difference scan mode:
	if option.mode == "scan":
	
		# load folders:
		contents, filtered, filetypes =  listFiles(option.path, output=list()), list(), list()
		backedup, biltered, biletypes =  listFiles(option.compare, output=list()), list(), list()
		
		for infile in contents:
			filtered.append(infile.replace(option.path,""))
		for infile in backedup:
			biltered.append(infile.replace(option.compare,""))
		
		bmissing = set(biltered).difference(set(filtered))
		print sorted(list(bmissing))
		print
		print "Files scanned:", len(contents)
		print "Files backedup:", len(backedup)
		print "Files in backup (only):", len(bmissing)
		print
		pdb.set_trace()
		
		
		"""
		# examine folder contents:
		index, collected = 1, list()
		for content in contents:
				#print index, ":", content
				index += 1
				
				filetype = grabFiletype(content)
					
				include_test = False
				if option.include == "OFF":
					include_test = True
				else:
					for include_flag in option.include.split(","):
						if include_flag in content:
							include_test = True
							
				exclude_test = True
				if option.exclude == "OFF":
					exclude_test = True
				else:
					for exclude_flag in option.exclude.split(","):
						if exclude_flag in content:
							exclude_test = False
					
				filetype_test = False
				if option.filetype == "OFF":
					filetype_test = True
				else:
					for filetype_flag in option.filetype.split(","):
						if filetype_flag == filetype.lower():
							filetype_test = True
						
				subfolder_test = False
				if option.subfolders == "OFF":
					subfolder_test = True
				else:
					folder = content.split("/")[0]
					if os.path.isdir(option.path + "/" + folder):
						subfolder_test = True
				
				if include_test and exclude_test and filetype_test and subfolder_test:
					filtered.append(content)
					collected.append(content)
					filetypes.append(filetype.lower())
		"""
	
			
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

#python reCover.py --path /Users/claraya/meTRN/python/ --mode scan --compare "/Volumes/Drobo/Projects/meTRN - v3/python/"
#python reCover.py --path /Users/claraya/meTRN/extras/ --mode scan --compare "/Volumes/Drobo/Projects/meTRN - v3/extras/"
#python reCover.py --path /Users/claraya/meTRN/fasta/ --mode scan --compare "/Volumes/Drobo/Projects/meTRN - v3/fasta/"
#python reCover.py --path /Users/claraya/meTRN/idr/ --mode scan --compare "/Volumes/Drobo/Projects/meTRN - v3/idr/"
#python reCover.py --path /Users/claraya/meTRN/input/ --mode scan --compare "/Volumes/Drobo/Projects/meTRN - v3/input/"
#python reCover.py --path /Users/claraya/meTRN/meme/ --mode scan --compare "/Volumes/Drobo/Projects/meTRN - v3/meme/"
#python reCover.py --path /Users/claraya/meTRN/qsub/ --mode scan --compare "/Volumes/Drobo/Projects/meTRN - v3/qsub/"
#python reCover.py --path /Users/claraya/meTRN/scripts/ --mode scan --compare "/Volumes/Drobo/Projects/meTRN - v3/scripts/"
#python reCover.py --path /Users/claraya/meTRN/setup/ --mode scan --compare "/Volumes/Drobo/Projects/meTRN - v3/setup/"
#python reCover.py --path /Users/claraya/meTRN/signal/ --mode scan --compare "/Volumes/Drobo/Projects/meTRN - v3/signal/"
#python reCover.py --path /Users/claraya/meTRN/data/ --mode scan --compare "/Volumes/Drobo/Projects/meTRN - v3/data/"

#python reCover.py --path /Users/claraya/meTRN/data/ --mode images --compare "/Users/claraya/meTRN/extras/modencode-images/"

