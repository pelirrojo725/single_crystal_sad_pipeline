#!/usr/bin/env python

# Automated single crystal SAD/SSAD pipeline that uses
# fast_dp, fast_ep, cad, dm, arpwarp


# Revisions

# v2.4 Adding other atoms besides 'S' and 'Se' in the S-SAD determination
# step. This needs to be IMPORVED!

# v2.3 No need to manually edit arp.sh

# v2.2 Can determine number of molecules in asu from S-SAD & sequence
# taken from seq; no need to manually enter in arp.sh file. 

import sys
import os
import time
import copy
import exceptions 
import traceback
import subprocess
import os.path
import xml.etree.ElementTree as ET
import xmltodict
import math
import logging

from xml.dom import minidom
from fast_ep_log_reader import fast_ep_log_hand, fast_ep_log_nsites
from sequence_reader import seq_first_residue
from determine_asu import determine_asu
from determine_solvent_content import determine_solvent_content

# Search for new data. Will leave this for now. Need to talk to John & Matt

# Run fast_dp; enter fast_dp like you would in terminal
data_input = raw_input("Enter fast_dp options where/data/is/img: ")
print ("Your data is in %s" % data_input)
os.system("fast_dp %s" % data_input)

# Run fast_ep. Input file is fast_dp.mtz		
if(os.path.isfile("fast_dp.mtz")):
	os.system("fast_ep sad=fast_dp.mtz")
	
# Run dm. Input file is sad.mtz. Look at dm receipes (google dm ccp4)
# to see how improve this further. Need more if statements depending
# on inputs (e.g. does protein have NCS, then add AVER to MODE).
if(os.path.isfile('sad.mtz')):
	
	# Parsing fast_dp.xml
	xmldoc_fast_dp = ET.parse("fast_dp.xml")
	root = xmldoc_fast_dp.getroot()
	spacegroup = root[1].find("spaceGroup").text
	SG = spacegroup.replace(" ", "") # removes the white spaces
	a_cell = root[1].find("refinedCell_a").text
	b_cell = root[1].find("refinedCell_b").text
	c_cell = root[1].find("refinedCell_c").text
	alpha = root[1].find("refinedCell_alpha").text
	beta = root[1].find("refinedCell_beta").text
	gamma = root[1].find("refinedCell_gamma").text
	
	# Determine hand
	fast_ep_log_hand('fast_ep.log')
	
	# Parse nsites from fast_ep.log
	fast_ep_log_reader('fast_ep.log')
	
	# Parsing seq file, determine first residue
	seq_first_residue('m1sea.pir') # Needs to come from database
	
	# Determining if S-SAD
	atom = ['S', 'Se', 'Zn']
	atom = atom[2]
	if atom[0]:
		# Count nsites from seq
		if seq_first_char == ['M']: 
			nsitesSeq = seq.count('M') - 1
			print nsitesSeq
		else:
			nsitesSeq = seq.count('M')
			print nsitesSeq
	elif atom[1]:
		# Count nsites from seq
		if seq_first_char == ['M']:
			nsitesSeq = seq.count('M') -1 + seq.count('C')
		else:
			nsitesSeq = seq.count('M') + seq.count('C')
	else:
		# If atom is something else nsitesSeq=nsitesFastep
		nsitesSeq = nsitesFastep
			
	# Compare seq sites and fastep sites
	def asu(nsitesSeq, nsitesFastep):
		return round(nsitesSeq / nsitesFastep)
		
	if atom[0]:
		asu_fast_ep = asu(nsitesSeq, nsitesFastep) # Calling the function
	elif atom[1]:
		asu_fast_ep = asu(nsitesSeq, nsitesFastep)
	else:
		asu_fast_ep = asu(nsitesSeq, nsitesFastep)
	#asu_fast_ep = 1
	#logging.basicConfig(filename = 'asu.log', level = logging.DEBUG)
	#logging.debug('%s %s %s' % (nsitesFastep, nsitesSeq, asuFastep))	
	#determine_asu(1)
		
	# Parsing Matts coef.
	#xmldoc_MATT = ET.parse("MATTHEWS_COEF.xml")
	#MattRoot = xmldoc_MATT.getroot()
	#for elem in MattRoot.findall("result"):
		#if asuFastep == float(elem.get("nmol_in_asu")):
			#print asuFastep
			#print float(elem.get("nmol_in_asu"))
			#solvent = float(elem.get("percent_solvent"))/100
			#print solvent
			#break
	#matthews_coef_xml_file = 'MATTHEWS_COEF_m1sea.xml'
	determine_solvent_content('MATTHEWS_COEF_m1sea.xml')	

	# Solvent and hand check
	if solvent < 0.275:
		shelxe = open(os.path.join("0.25", "shelxe.sh"), "w")
		if "original" in fast_ep_string:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		else:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF\n"]
			shelxe.writelines(shelxe_sh)
		shelxe.close()
		os.system("sh 0.25/shelxe.sh")
	elif solvent >= 0.275 and solvent < 0.325:
		shelxe = open(os.path.join("0.30", "shelxe.sh"), "w")
		if "original" in fast_ep_string:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		else:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		shelxe.close()
		os.system("sh 0.30/shelxe.sh")
	elif solvent >= 0.325 and solvent < 0.375:
		shelxe = open(os.path.join("0.35", "shelxe.sh"), "w")
		if "original" in fast_ep_string:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		else:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		shelxe.close()
		os.system("sh 0.35/shelxe.sh")
	elif solvent >= 0.375 and solvent < 0.425:
		shelxe = open(os.path.join("0.40", "shelxe.sh"), "w")
		if "original" in fast_ep_string:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		else:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		shelxe.close()
		os.system("sh 0.40/shelxe.sh")
	elif solvent >= 0.425 and solvent < 0.475:
		shelxe = open(os.path.join("0.45", "shelxe.sh"), "w")
		if "original" in fast_ep_string:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		else:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		shelxe.close()
		os.system("sh 0.45/shelxe.sh")
	elif solvent >= 0.475 and solvent < 0.525:
		shelxe = open(os.path.join("0.50", "shelxe.sh"), "w")
		if "original" in fast_ep_string:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		else:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		shelxe.close()
		os.system("sh 0.50/shelxe.sh")
	elif solvent >= 0.525 and solvent < 0.575:
		shelxe = open(os.path.join("0.55", "shelxe.sh"), "w")
		if "original" in fast_ep_string:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		else:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		shelxe.close()
		os.system("sh 0.55/shelxe.sh")
	elif solvent >= 0.575 and solvent < 0.625:
		shelxe = open(os.path.join("0.60", "shelxe.sh"), "w")
		if "original" in fast_ep_string:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		else:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
			os.system("sh 0.60shelxe.sh")
		shelxe.close()
		os.system("sh 0.60/shelxe.sh")
	elif solvent >= 0.625 and solvent < 0.675:
		shelxe = open(os.path.join("0.65", "shelxe.sh"), "w")
		if "original" in fast_ep_string:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		else:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		shelxe.close()
		os.system("sh 0.65/shelxe.sh")
	elif solvent >= 0.675 and solvent < 0.725:
		shelxe = open(os.path.join("0.70", "shelxe.sh"), "w")
		if "original" in fast_ep_string:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		else:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		shelxe.close()
		os.system("sh 0.70/shelxe.sh")
	elif solvent >= 0.725:
		shelxe = open(os.path.join("0.75", "shelxe.sh"), "w")
		if "original" in fast_ep_string:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		else:
			shelxe_sh = ["shelxe sad sad_fa -s%.2f -a2 -m50 -b -i\n" % solvent,
						 "f2mtz HKLIN sad.phs HKLOUT phs.mtz <<EOF\n",
						 "CELL %s %s %s %s %s %s\n" % (a_cell, b_cell, 
						 c_cell, alpha, beta, gamma),
						 "SYMM %s\n" % SG,
						 "LABOUT H K L F FOM PHI SIGF\n",
						 "CTYP H H H F W P Q\n",
						 "END\n" "EOF"]
			shelxe.writelines(shelxe_sh)
		shelxe.close()
		os.system("sh 0.75/shelxe.sh")
#else:
	#print "Matthew's coefficent check failed"
 
# Adding PHI and FOM from phs.mtz to sad.mtz
if (os.path.isfile("phs.mtz")):
	cad = open("cad.sh", "w")
	cad_sh = ["cad HKLIN1 sad.mtz HKLIN2 phs.mtz HKLOUT cad.mtz <<EOF\n",
			  "LABIN  FILE 1 E1=F   E2=SIGF E3=FreeF_flag\n",
			  "LABOUT FILE 1 E1=F   E2=SIGF E3=FreeF_flag\n",
			  "CTYP   FILE 1 E1=F   E2=Q    E3=I\n",
			  "LABIN  FILE 2 E1=FOM E2=PHI\n",
			  "LABOUT FILE 2 E1=FOM E2=PHI\n",
			  "CTYP   FILE 2 E1=W   E2=P\n" "END\n" "EOF"]
	cad.writelines(cad_sh)
	cad.close()
	os.system("sh cad.sh")
#else:
#	print "CAD failed, couldn't add PHI and FOM to sad.mtz"

# Writes the sh file to run dm
if (os.path.isfile("cad.mtz")):
	dm = open("dm.sh", "w")
	dm_sh = ["dm HKLIN cad.mtz HKLOUT cad_dm.mtz << endd\n",
				"SOLC %s\n" % solvent,
				"MODE SOLV HIST MULT\n",
				"NCYCLE 40\n",
				"SCHEME ALL\n",
				"COMBINE PERT FREE 2\n",
				"LABIN FP=F SIGFP=SIGF PHIO=PHI FOMO=FOM FREE=FreeF_flag\n",
				"LABOUT PHIDM=PHIDM FOMDM=FOMDM\n",
				"END\n"
				"endd"]
	dm.writelines(dm_sh)	
	dm.close()
	os.system("sh dm.sh")
#else:
#	print "Density modification failed"
	
# Run arpWarp
if (os.path.isfile("cad_dm.mtz")):
	arp = open("arp.sh", "w")
	#Matt_out = open("MATTHEWS_COEF.xml", "r")
	#Matt_r = Matt_out.read()
	#MattDic = xmltodict.parse(Matt_r)
	nmol_in_asu = asu(nsitesSeq, nsitesFastep)
	#asu = nmol_in_asu.replace(" ", "")
	arp_sh = ["auto_tracing.sh datafile $PWD/cad_dm.mtz ",
				"residues %s cgr %s seqin $PWD/%s " % (NumResInSeq, nmol_in_asu, seqfile),
				"modelin $PWD/sad.pdb workdir $PWD fp FP ",
				"sigfp SIGFP phibest PHI fom FOM ",
				"freelabin FreeF_flag buildingcycles 2"]
	arp.writelines(arp_sh)
	arp.close()
	arp = open("arp.sh", "r")
	print arp.read()
	os.system("sh arp.sh")
#else:
#	print "Model buidling failed"
	









