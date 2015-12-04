#!/usr/bin/python
import numpy as np
import itertools
import textwrap
import sys
import os
import re
import timeit
import glob
from collections import OrderedDict
from collections import defaultdict
import fileinput

start = timeit.default_timer()

Usage = """
Usage:  ./reconstruct_gbk_from_get_homologue_output.py OUTPUT.faa\n

The script MUST be run in the same directory as where your original .gbk files that were run through "get_homologues" are located.

The script takes one arguments:
i) the output from "compile_pangenome_output.py" which is a large fasta file containing unwieldly headings which represent condensed GBK information from RAST annotations

Note: spaces at the end of your file can cause this script to flip.

"""

if len(sys.argv)<2:
        print Usage

else:
        FAA = sys.argv[1]

##PHASE I: Extract all information from .faa output from get_homologues and parse into dictionary

#Open .faa file
InFile = open(FAA, 'r')

#Bring in the data and separate header information from sequence
Fasta_header = []
foo = []
AA = []

for Line in InFile:
        Fasta_header.append(Line.strip('\n\r'))
	foo = InFile.next()
	foo = foo.strip('\n\r')
	AA.append(foo)

InFile.close()

#Split and capture information of interested from header
#Store original Fasta_header and AA sequence in dictionary with key1= Genome_ID

Genome_ID = []

for Datum in Fasta_header:
	Diced = Datum.split("|")
	Genome_ID.append(Diced[1])

#Remove "[""]" brackets from Genome_ID
foo=[]

for Push in Genome_ID:
	foo1 = re.sub(r'\[|\]','', Push)
	foo.append(foo1)

Genome_ID = foo	

#Make dictionary with Genome_ID as key and Fasta_header and AA as values for each key
count = 0
BigDict = defaultdict(list)

for count in range(len(Genome_ID)):
	BigDict[Genome_ID[count]].append([Fasta_header[count], AA[count]])

##From HERE on out, I will iterate through each Genome_ID and go ALL the way to the written file
#For each Genome_ID in Dictionary unpack its values and create second dictionary with contig as key
r = re.compile(r'.*Contig ([\d]+)')

#Open file to write to named according to your input
FileName = re.sub('output_','',FAA)
FileName = re.sub('_list.txt','',FileName)
output = open(FileName + '.temp', 'w')

##For each genome create a dictionary for each contig
for key in BigDict.iterkeys():
	ContigBin = defaultdict(list)	#NOTE: CDS will be overwritten for each new Genome_ID
	Contig = []
	CDS = []
	
	#Page contains all information for a single genome
	for Page in BigDict[key]:
		CDS = [Page[0],Page[1]]
		#Create a dictionary with all contgs fround in a single genome
		Diced = []
		#Extract only first element from CDS (which contains the Fasta_header
		for Cut in Page:
			Diced = Cut.split("|")
			if re.match(r'>', Diced[0]):
				foo = Diced[5]
				match = r.match(foo)
				ContigBin[match.group(1)].append(CDS)

#print ContigBin ##ContigBin is a dictionary with values being a list of a list

	##HERE You have a library for each Genome_ID which consists of a key for each contig and a value for
	#all the genes found on that contig
	GBK_list = []

	for all_gbk_in_directory in glob.glob("*.gbk"):
		foo = re.sub('.gbk', '', all_gbk_in_directory)
		GBK_list.append(foo)

	
	#Pick out file name in which the current contig appears
	#This assumes that part of the file ID will correspond with the ID used in the .gbk file
	Current_Genome = []
	Contig_List = []

	for Contig in ContigBin.iterkeys():
		foo = []
		for Contig_CDS in ContigBin[Contig][0]:
			Diced = Contig_CDS.split("|")
			if re.match(r'>', Diced[0]):
				Current_Genome = Diced[1]
				Current_Genome = re.sub(r'\[|\]','',Current_Genome)
				
				foo1 = Diced[5]
				match = r.match(foo1)
				foo = match.group(1)
		Contig_List.append(foo)

	#Pick out the exact location within the file
	Contig_Header_dict = {}

	for Current_Contig in Contig_List:
		Contig_Header_start = []
		Contig_Header_end = []
		Contig_Origin_start = []
		Contig_Origin_end = []

	        for Partial_Genome_ID in Current_Genome.split(" "):
        	        if Partial_Genome_ID in GBK_list:

				##Locate the start and end of the header and origin portion of contig
				lookup_start = "Contig " + Current_Contig
				lookup_end = "CDS"
				lookup_origin_start = "ORIGIN"
				lookup_origin_end = "//"

				#Find lines which match key pieces of information
				with open(Partial_Genome_ID + ".gbk") as myFile:
				    for num, line in enumerate(myFile, 1):
			        	if lookup_start in line:
				            Contig_Header_start = num-1

				with open(Partial_Genome_ID + ".gbk") as myFile:
				    for num, line in enumerate(myFile, 1):
				        if lookup_end in line:
					    Contig_Header_end.append(num)

				with open(Partial_Genome_ID + ".gbk") as myFile:
				    for num, line in enumerate(myFile, 1):
				        if lookup_origin_start in line:
				            Contig_Origin_start.append(num)

				with open(Partial_Genome_ID + ".gbk") as myFile:
				    for num, line in enumerate(myFile, 1):
				        if lookup_origin_end in line:
					    Contig_Origin_end.append(num)
			
				#Using the fact that "Contig + Contig Number" is unique to each Contig entry
				#I am appending the line number that matches this into each large list containing the 
				#multiple entries for "CDS", "ORIGIN", and "//"
				#Then, I sort the list and pick out the element of the list beside the known Contig_Header_start

				Contig_Header_end.append(Contig_Header_start)		
				Contig_Header_end.sort()
			
				Contig_Origin_start.append(Contig_Header_start)		
				Contig_Origin_start.sort()

				Contig_Origin_end.append(Contig_Header_start)		
				Contig_Origin_end.sort()
			
				Target = [i for i,x in enumerate(Contig_Header_end) if x == Contig_Header_start]
				Contig_Header_end = Contig_Header_end[Target[0]+1]-1
			
				Target = [i for i,x in enumerate(Contig_Origin_start) if x == Contig_Header_start]
				Contig_Origin_start = Contig_Origin_start[Target[0]+1]-2

				Target = [i for i,x in enumerate(Contig_Origin_end) if x == Contig_Header_start]
				Contig_Origin_end = Contig_Origin_end[Target[0]+1]

				Contig_Header_start = Contig_Header_start - 1
				Contig_Header_dict[Current_Contig] = [Partial_Genome_ID, Contig_Header_start, Contig_Header_end, Contig_Origin_start, Contig_Origin_end]

#Contig_Header_dict now contains all the start and stop information for each contig present

##Now begins the process of actually writing the information to file
#Writer Header to File for the Current Contig

	for ContigPrint in ContigBin:

		#Extract pertinent CDS info from Contig dictionary
		Current_contig_Gene_name_and_EC = []
		Current_contig_Strand = []
		Current_contig_FAA = []

		for ContigData in ContigBin[ContigPrint]:
			for ContigCDS in ContigData:
				Diced = ContigCDS.split("|")
				if re.match(r'>', Diced[0]):
					#Capture All Relevant Information
					Current_contig_Gene_name_and_EC.append(Diced[3])
					Current_contig_Strand.append(Diced[5])
				if re.match(r'^\w', Diced[0]):
					Current_contig_FAA.append(Diced)
		

		Strand = []
		locus = []
		for Find in Current_contig_Strand:
		        Diced = Find.split(':')
		        match = re.match(r'(^[\S]+) ', Diced[2])
			Strand.append(str.strip(match.group()))
			locus.append(re.sub(r'-','..',Diced[1]))

		GO = []
		for Find in Current_contig_Strand:
		        Diced = Find.split('^')
			if re.search("GO:[\S]+,", Diced[1]):
				Double_diced = Diced[1].split(',')
				GO.append(Double_diced)

			else:
				GO.append('')

		#Split Gene_name and EC
		EC_number = []
		Gene_name = []
		multiple_EC = []
		multiple_names = []

##Deal with the fact many records are missing EC numbers and then that many records have multiple EC numbers
#Note: missing EC_numbers are given place holders "", which are removed in the final step
		for Ream in Current_contig_Gene_name_and_EC:
			if re.search("EC \d", Ream):
				iter = re.finditer(r"\(EC ", Ream)
				indices = [m.start(0) for m in iter]

				if len(indices) > 1:
					cnt = 0
					for ICI in indices:
						foo1 = re.sub(r'/',"",Ream[cnt:(ICI-1)]) 
						foo1 = str.strip(foo1)
						multiple_EC.append(Ream[ICI:(ICI+12)])
						multiple_names.append(foo1)
						cnt = cnt + (ICI+13)
					
					EC_number.append(multiple_EC)
					Gene_name.append(multiple_names)
							
					multiple_EC = []
					multiple_names = []

				else:
					strt = re.search("\(EC ", Ream).start()
					Gene_name.append(str.strip(Ream[0:(strt-1)]))
					EC_number.append(Ream[strt:(strt+12)])
	
			else:
				EC_number.append('')
				Gene_name.append(Ream)

		preformed_EC = []
		preformed_name = []
		preformed_product = []
		preformed_translation = []
		preformed_GO = []

	#Gonna have to deal with multiple EC_numbers
		EC_holder = []
		count = 0

		#Clean EC numbers to appropriate format (subtract brackets and EC)
		for count in range(len(EC_number)):
			if np.size(EC_number[count]) > 1:
				foo = []
				for EC_clean in EC_number[count]:
					foo1 = re.sub(r'\)|\(|EC ','', EC_clean)
					foo.append(foo1)
				EC_holder.append(foo)

			else:
				if re.search("EC \d", EC_number[count]):
					EC_holder.append(re.sub(r'\)|\(|EC ','', EC_number[count]))

				else:
					EC_holder.append('')

		#Print EC data to preformed state
		for count in range(len(EC_holder)):
			if (np.size(EC_holder[count]) > 1):  
	
				foo = []
				for EC in EC_holder[count]:
					foo1 = "/EC_number=\""+EC+"\""+"\n"
					foo.append(foo1)
				preformed_EC.append(foo)
		
			else:
				preformed_EC.append("/EC_number=\""+EC_holder[count]+"\""+"\n")

		#Print GO data to preformed state
		for count in range(len(GO)):
			if (np.size(GO[count]) > 1):  
				foo = []
				for GO_get in GO[count]:
					foo1 = "/db_xref=\""+GO_get+"\""+"\n"
					foo.append(foo1)
				preformed_GO.append(foo)
		
			else:
				preformed_GO.append("/db_xref=\""+GO[count]+"\""+"\n")

		#Print Gene name preformed state
		for count in range(len(Gene_name)):
			if (np.size(Gene_name[count]) > 1):  
	
				foo = []
				foof = []
				for Gene in Gene_name[count]:
					foo1 = "/gene=\""+Gene+"\""+"\n"
					foo.append(foo1)
					foo2 = "/product=\""+Gene+"\""+"\n"
					foof.append(foo2)
				preformed_name.append(foo)
				preformed_product.append(foof)
		
			else:
				preformed_name.append("/gene=\""+Gene_name[count]+"\""+"\n")
				preformed_product.append("/product=\""+Gene_name[count]+"\""+"\n")


	#Make translated sequence one block	
		foo = []
		for Translated_AA in Current_contig_FAA:
			preformed_translation.append("/translation=\""+Translated_AA[0]+"\""+"\n")

#Begin writing the meat of the .gbk file starting with a fork for CDS which are complementary (-1) or not
#This is really ugly b/c I've included the correct number of whitespaces to have everything aligned
#The gist is that I fork for writing different CDS (complementary or not), then I account for having multiple
#gene names associated with each entry as well as EC numbers. Then I write the Gene_comment (which is hte 
#genome from which the gene originates), followed by the tranlsated protein sequence and then finally a
#repeate of the gene name, called the product (I had to include "gene" because pathway-tools requires it

#		print Current_Contig		
#		print Strand
#		print locus
#		print preformed_product
#		print preformed_name
#		print preformed_EC
#		print preformed_GO
#		print preformed_translation

		for line in open(Contig_Header_dict[ContigPrint][0] + ".gbk", "r").readlines()[Contig_Header_dict[ContigPrint][1]:Contig_Header_dict[ContigPrint][2]]:
			output.write(line.rstrip('\n')+"\n")
		
		for count in range(len(Strand)):

			if (Strand[count] == "-1"):
				output.write(
				"     CDS             complement("+locus[count]+")\n")
				if np.size(preformed_name[count]) > 1:
					for Gene in preformed_name[count]:
						output.write("                     "+Gene)
				else:
					output.write("                     "+preformed_name[count])

				if np.size(preformed_EC[count]) > 1:
					for EC in preformed_EC[count]: 
						output.write("                     "+EC)
				else:
			 		output.write("                     "+preformed_EC[count])

				if np.size(preformed_GO[count]) > 1:
					for GO_get in preformed_GO[count]: 
						output.write("                     "+GO_get)
				else:
			 		output.write("                     "+preformed_GO[count])

				output.write(textwrap.fill(preformed_translation[count],width=65,initial_indent='\t\t     ',subsequent_indent='\t\t     ')+"\n")
				if np.size(preformed_product[count]) > 1:
					for Product in preformed_product[count]:
						output.write("                     "+Product)
				else:
					output.write("                     "+preformed_product[count])

#Repeat above steps for (+1) CDS
			else:
				output.write(
				"     CDS             "+locus[count]+"\n")
				if np.size(preformed_name[count]) > 1:
					for Gene in preformed_name[count]:
						output.write("                     "+Gene)
				else:
					output.write("                     "+preformed_name[count])

				if np.size(preformed_EC[count]) > 1:
					for EC in preformed_EC[count]: 
						output.write("                     "+EC)
				else:
			 		output.write("                     "+preformed_EC[count])

				if np.size(preformed_GO[count]) > 1:
					for GO_get in preformed_GO[count]: 
						output.write("                     "+GO_get)
				else:
			 		output.write("                     "+preformed_GO[count])

				output.write(textwrap.fill(preformed_translation[count],width=65,initial_indent='\t\t     ',subsequent_indent='\t\t     ')+"\n")
				if np.size(preformed_product[count]) > 1:
					for Product in preformed_product[count]:
						output.write("                     "+Product)
				else:
					output.write("                     "+preformed_product[count])

	#Add Origin

		for line in open(Contig_Header_dict[ContigPrint][0] + ".gbk", "r").readlines()[Contig_Header_dict[ContigPrint][3]:Contig_Header_dict[ContigPrint][4]]:
			output.write(line.rstrip('\n')+"\n")

output.close()

#Use sed to replace all lines with empty values
#I went this route b/c it was easiest to keep "" as place holders until the very end, rather than write conditional statements to selectively print real content

os.system(' '.join([
	"sed -i '/\/EC_number=\"\"/d'",
	FileName+'.temp'
	]))

os.system(' '.join([
	"sed -i '/\/product=\"\"/d'",
	FileName+'.temp'
	]))

os.system(' '.join([
	"sed -i '/\/gene=\"\"/d'",
	FileName+'.temp'
	]))

os.system(' '.join([
	"sed -i '/\/db_xref=\"\"/d'",
	FileName+'.temp'
	]))

###Addendum to script after realizing that contig names *obviously* need to be unique
#Determine the occurrence of all contig names
r = re.compile(r'.*Contig ([\d]+)')

Contig_names = {}
seen = {}

ContigFind = "Contig "

o = open(FileName+".gbk", "w")

#Create dictionary of all occurrences of contig name
for line in open(FileName+".temp"):
        if ContigFind in line:
                match = r.match(line)

                if not Contig_names.has_key(match.group(1)):
                        Contig_names[match.group(1)] = 0

                Contig_names[match.group(1)] = Contig_names[match.group(1)]+1

#Find contig names with >1 occurrence and overwrite line with "contigname.x"
for line in open(FileName+".temp"):
        match1 = re.match(r'LOCUS[\s]+([\d]+)', line)   #LOCUS comes first in the ordering
        match2 = re.match(r'DEFINITION  Contig[\s]([\d]+)', line)

        if re.search(r'LOCUS[\s]+([\d]+)', line):
                if Contig_names[match1.group(1)] > 1:

                        #Make dictionary to count the number of occurrences of each contig name to create appropriate ".x"
                        if not match1.group(1) in seen:
                                seen[match1.group(1)] = 0

                        seen[match1.group(1)] = seen[match1.group(1)] + 1

                        if seen[match1.group(1)] > 1:
                                o.write(re.sub(match1.group(0)+" ",match1.group(0)+"."+str(seen[match1.group(1)])+" ", line))

                        else:
                                o.write(line)
                else:
                        o.write(line)
	elif re.search(r'DEFINITION  Contig[\s]([\d]+)', line):
        	if Contig_names[match2.group(1)] > 1:

                	if seen[match2.group(1)] > 1:
                        	o.write(re.sub(match2.group(0)+" ",match2.group(0)+"."+str(seen[match2.group(1)])+" ", line))

			else:
                        	o.write(line)
		else:
			o.write(line)
	else:
		o.write(line)

o.close()

os.system(' '.join([
	"rm",
	FileName+'.temp'
	]))


stop = timeit.default_timer()

print "This operation took " + str(stop - start) + " seconds."

print "Your input has been parsed and a new pseudo-genbank file has been created entitled: "+FileName+".gbk"



