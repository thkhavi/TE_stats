#!/usr/bin/python
# 06/01/2015
# Written by Sandra Truong and Ryan McCormick
#
# This python script takes in .gff3 file of Transposable Element (TE) Annotations 
#	constructed using the TEdenovo & TEannot pipelines (Flutre et al., 2011)
#	and gives summary statistics of TE coverage given the Wicker classification 
#	(Wicker et al., 2007).
#
# Example use:
# python classification_te_coverage.py ../super_1883.gff3

import sys
import HTSeq
import itertools
from matplotlib import pyplot

###
# Global variables
###
INDEX_COVERAGE = 0
INDEX_CHANGE = 1

""" Class I (retrotransposons)"""
TE_RXX = [int(0), False]
""" Order LTR """
TE_RLX = [int(0), False]
""" Superfamily Copia """
TE_RLC = [int(0), False]
""" Superfamily Gypsy """
TE_RLG = [int(0), False]
""" Superfamily Bel-Pao """
TE_RLB = [int(0), False]
""" Superfamily Retrovirus """
TE_RLR = [int(0), False]
""" Superfamily ERV """
TE_RLE = [int(0), False]
""" Order DIRS """
TE_RYX = [int(0), False]
""" Superfamily DIRS """
TE_RYD = [int(0), False]
""" Superfamily Ngaro """
TE_RYN = [int(0), False]
""" Superfamily VIPER """
TE_RYV = [int(0), False]
""" Order PLE """
TE_RPX = [int(0), False]
""" Superfamily Penelope """
TE_RPP = [int(0), False]
""" Order LINE """
TE_RIX = [int(0), False]
""" Superfamily R2 """
TE_RIR = [int(0), False]
""" Superfamily RTE """
TE_RIT = [int(0), False]
""" Superfamily Jockey """
TE_RIJ = [int(0), False]
""" Superfamily L1 """
TE_RIL = [int(0), False]
""" Superfamily I """
TE_RII = [int(0), False]
""" Order SINE """
TE_RSX = [int(0), False]
""" Superfamily tRNA """
TE_RST = [int(0), False]
""" Superfamily 7SL """
TE_RSL = [int(0), False]
""" Superfamily 5S """
TE_RSS = [int(0), False]

""" Class II (DNA transposons)"""
TE_DXX = [int(0), False]
""" Subclass I: Order TIR """
TE_DTX = [int(0), False]
""" Superfamily Tc1-Mariner """
TE_DTT = [int(0), False]
""" Superfamily hAT """
TE_DTA = [int(0), False]
""" Superfamily Mutator """
TE_DTM = [int(0), False]
""" Superfamily Merlin """
TE_DTE = [int(0), False]
""" Superfamily Transib """
TE_DTR = [int(0), False]
""" Superfamily P """
TE_DTP = [int(0), False]
""" Superfamily PiggyBac """
TE_DTB = [int(0), False]
""" Superfamily PIF-Harbinger """
TE_DTH = [int(0), False]
""" Superfamily CACTA """
TE_DTC = [int(0), False]
""" Subclass I: Order Crypton """
TE_DYX = [int(0), False]
""" Superfamily Crypton """
TE_DYC = [int(0), False]
""" Subclass II: Order Helitron """
TE_DHX = [int(0), False]
""" Superfamily Helitron """
TE_DHH = [int(0), False]
""" Subclass II: Order Maverick """
TE_DMX = [int(0), False]
""" Superfamily Maverick """
TE_DMM = [int(0), False]

""" SSRs"""
TE_SSR = [int(0), False]

###
# End global variables
###


###
# Utility functions
###
def usage():
     sys.stderr.write("\nclassification_te_coverage.py expects\
                       \n\t(1) path to GFF3 file for contig of interest\
                       \n\t(2) length of contig (bp)\
                       \nExample usage:\
                       \n\tpython classification_te_coverage.py ../super_1883.gff3 \n\n")
     sys.stderr.flush()
     sys.exit()

###
# End utility functions
###


###
# Begin main()
###
if len(sys.argv) <= 1:
     usage()
if sys.argv[1] == "--help" or sys.argv[1] == "-h":
     usage()

# Read input gff3 file.
try:
     GFF3_FILE = HTSeq.GFF_Reader(sys.argv[1], end_included=True)
     NUMBER_TE = len([line.strip() for line in open(sys.argv[1])]) - 2
     CONTIG_ID = ((sys.argv[1].split(".gff3"))[0]).split("/")[-1]
     LILI_TABLE = [line.strip() for line in open(sys.argv[1])]
     CONTIG_LENGTH = int(LILI_TABLE[1].split(" ")[-1]) 
except IOError:
     sys.stderr.write("\nCannot open target gff3 file. Please check your input:\n")
     usage()

# Populate genomic array of sets
print("Populating genomic array of sets.")
GAS = HTSeq.GenomicArrayOfSets( [CONTIG_ID], stranded=False )
for t_element in itertools.islice(GFF3_FILE,0,None):
   if t_element.source == "sbi1_REPET_TEs":
      GAS[t_element.iv] += "TE@" + (t_element.attr['Target'])[:3]
   elif t_element.source == "sbi1_REPET_tblastx" or t_element.source == "sbi1_REPET_blastx":
      GAS[t_element.iv] += "blast@" + t_element.attr['ID']
   elif t_element.source == "sbi1_REPET_SSRs":
      GAS[t_element.iv] += "SSR@" + t_element.attr['ID']
print("Finished populating genomic array of sets.")

###
# Begin Loop
###

for bp_position in range(1, CONTIG_LENGTH):
   for t_element in list(GAS[HTSeq.GenomicPosition(CONTIG_ID, bp_position)]):
      if t_element.split("@")[0] == "TE":
         wickers_class = t_element.split("@")[1]
         if list(wickers_class)[0] == "R":
            """ Class I (retrotransposons)"""
            if TE_RXX[INDEX_CHANGE] == False:
               TE_RXX[INDEX_COVERAGE] = int(1) + TE_RXX[INDEX_COVERAGE]
               TE_RXX[INDEX_CHANGE] = True
            if list(wickers_class)[1] == "L":
               """ Order LTR """
               if TE_RLX[INDEX_CHANGE] == False:
                  TE_RLX[INDEX_COVERAGE] = int(1) + TE_RLX[INDEX_COVERAGE]
                  TE_RLX[INDEX_CHANGE] = True
               if list(wickers_class)[2] == "C":
                  """ Superfamily Copia """
                  if TE_RLC[INDEX_CHANGE] == False:
                     TE_RLC[INDEX_COVERAGE] = int(1) + TE_RLC[INDEX_COVERAGE]
                     TE_RLC[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "G":
                  """ Superfamily Gypsy """
                  if TE_RLG[INDEX_CHANGE] == False:
                     TE_RLG[INDEX_COVERAGE] = int(1) + TE_RLG[INDEX_COVERAGE]
                     TE_RLG[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "G":
                  """ Superfamily Bel-Pao """
                  if TE_RLB[INDEX_CHANGE] == False:
                     TE_RLB[INDEX_COVERAGE] = int(1) + TE_RLB[INDEX_COVERAGE]
                     TE_RLB[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "R":
                  """ Superfamily Retrovirus """
                  if TE_RLR[INDEX_CHANGE] == False:
                     TE_RLR[INDEX_COVERAGE] = int(1) + TE_RLR[INDEX_COVERAGE]
                     TE_RLR[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "E":
                  """ Superfamily ERV """
                  if TE_RLE[INDEX_CHANGE] == False:
                     TE_RLE[INDEX_COVERAGE] = int(1) + TE_RLE[INDEX_COVERAGE]
                     TE_RLE[INDEX_CHANGE] = True
            elif list(wickers_class)[1] == "Y":
               """ Order DIRS """
               if TE_RYX[INDEX_CHANGE] == False:
                  TE_RYX[INDEX_COVERAGE] = int(1) + TE_RYX[INDEX_COVERAGE]
                  TE_RYX[INDEX_CHANGE] = True
               if list(wickers_class)[2] == "D":
                  """ Superfamily DIRS """
                  if TE_RYD[INDEX_CHANGE] == False:
                     TE_RYD[INDEX_COVERAGE] = int(1) + TE_RYD[INDEX_COVERAGE]
                     TE_RYD[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "N":
                  """ Superfamily Ngaro """
                  if TE_RYN[INDEX_CHANGE] == False:
                     TE_RYN[INDEX_COVERAGE] = int(1) + TE_RYN[INDEX_COVERAGE]
                     TE_RYN[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "E":
                  """ Superfamily VIPER """
                  if TE_RYV[INDEX_CHANGE] == False:
                     TE_RYV[INDEX_COVERAGE] = int(1) + TE_RYV[INDEX_COVERAGE]
                     TE_RYV[INDEX_CHANGE] = True
            elif list(wickers_class)[1] == "P":
               """ Order PLE """
               if TE_RPX[INDEX_CHANGE] == False:
                  TE_RPX[INDEX_COVERAGE] = int(1) + TE_RPX[INDEX_COVERAGE]
                  TE_RPX[INDEX_CHANGE] = True
               if list(wickers_class)[2] == "P":
                  """ Superfamily Penelope """
                  if TE_RPP[INDEX_CHANGE] == False:
                     TE_RPP[INDEX_COVERAGE] = int(1) + TE_RPP[INDEX_COVERAGE]
                     TE_RPP[INDEX_CHANGE] = True
            elif list(wickers_class)[1] == "I":
               """ Order LINE """
               if TE_RIX[INDEX_CHANGE] == False:
                  TE_RIX[INDEX_COVERAGE] = int(1) + TE_RIX[INDEX_COVERAGE]
                  TE_RIX[INDEX_CHANGE] = True
               if list(wickers_class)[2] == "R":
                  """ Superfamily R2 """
                  if TE_RIR[INDEX_CHANGE] == False:
                     TE_RIR[INDEX_COVERAGE] = int(1) + TE_RIR[INDEX_COVERAGE]
                     TE_RIR[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "T":
                  """ Superfamily RTE """
                  if TE_RIT[INDEX_CHANGE] == False:
                     TE_RIT[INDEX_COVERAGE] = int(1) + TE_RIT[INDEX_COVERAGE]
                     TE_RIT[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "J":
                  """ Superfamily Jockey """
                  if TE_RIJ[INDEX_CHANGE] == False:
                     TE_RIJ[INDEX_COVERAGE] = int(1) + TE_RIJ[INDEX_COVERAGE]
                     TE_RIJ[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "L":
                  """ Superfamily L1 """
                  if TE_RIL[INDEX_CHANGE] == False:
                     TE_RIL[INDEX_COVERAGE] = int(1) + TE_RIL[INDEX_COVERAGE]
                     TE_RIL[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "I":
                  """ Superfamily I """
                  if TE_RII[INDEX_CHANGE] == False:
                     TE_RII[INDEX_COVERAGE] = int(1) + TE_RII[INDEX_COVERAGE]
                     TE_RII[INDEX_CHANGE] = True
            elif list(wickers_class)[1] == "S":
               """ Order SINE """
               if TE_RSX[INDEX_CHANGE] == False:
                  TE_RSX[INDEX_COVERAGE] = int(1) + TE_RSX[INDEX_COVERAGE]
                  TE_RSX[INDEX_CHANGE] = True
               if list(wickers_class)[2] == "T":
                  """ Superfamily tRNA """
                  if TE_RST[INDEX_CHANGE] == False:
                     TE_RST[INDEX_COVERAGE] = int(1) + TE_RST[INDEX_COVERAGE]
                     TE_RST[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "L":
                  """ Superfamily 7SL """
                  if TE_RSL[INDEX_CHANGE] == False:
                     TE_RSL[INDEX_COVERAGE] = int(1) + TE_RSL[INDEX_COVERAGE]
                     TE_RSL[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "S":
                  """ Superfamily 5S """
                  if TE_RSS[INDEX_CHANGE] == False:
                     TE_RSS[INDEX_COVERAGE] = int(1) + TE_RSS[INDEX_COVERAGE]
                     TE_RSS[INDEX_CHANGE] = True
         elif list(wickers_class)[0] == "D":
            """ Class II (DNA transposons)"""
            if TE_DXX[INDEX_CHANGE] == False: 
               TE_DXX[INDEX_COVERAGE] = int(1) + TE_DXX[INDEX_COVERAGE]
               TE_DXX[INDEX_CHANGE] = True
               wickers_class = t_element.split("@")[1]
            if list(wickers_class)[1] == "T":
               """ Subclass I: Order TIR """
               if TE_DTX[INDEX_CHANGE] == False:
                  TE_DTX[INDEX_COVERAGE] = int(1) + TE_DTX[INDEX_COVERAGE]
                  TE_DTX[INDEX_CHANGE] = True
               if list(wickers_class)[2] == "T":
                  """ Superfamily Tc1-Mariner """
                  if TE_DTT[INDEX_CHANGE] == False:
                     TE_DTT[INDEX_COVERAGE] = int(1) + TE_DTT[INDEX_COVERAGE]
                     TE_DTT[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "A":
                  """ Superfamily hAT """
                  if TE_DTA[INDEX_CHANGE] == False:
                     TE_DTA[INDEX_COVERAGE] = int(1) + TE_DTA[INDEX_COVERAGE]
                     TE_DTA[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "M":
                  """ Superfamily Mutator """
                  if TE_DTM[INDEX_CHANGE] == False:
                     TE_DTM[INDEX_COVERAGE] = int(1) + TE_DTM[INDEX_COVERAGE]
                     TE_DTM[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "E":
                  """ Superfamily Merlin """
                  if TE_DTE[INDEX_CHANGE] == False:
                     TE_DTE[INDEX_COVERAGE] = int(1) + TE_DTE[INDEX_COVERAGE]
                     TE_DTE[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "R":
                  """ Superfamily Transib """
                  if TE_DTR[INDEX_CHANGE] == False:
                     TE_DTR[INDEX_COVERAGE] = int(1) + TE_DTR[INDEX_COVERAGE]
                     TE_DTR[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "P":
                  """ Superfamily P """
                  if TE_DTP[INDEX_CHANGE] == False:
                     TE_DTP[INDEX_COVERAGE] = int(1) + TE_DTP[INDEX_COVERAGE]
                     TE_DTP[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "B":
                  """ Superfamily PiggyBac """
                  if TE_DTB[INDEX_CHANGE] == False:
                     TE_DTB[INDEX_COVERAGE] = int(1) + TE_DTB[INDEX_COVERAGE]
                     TE_DTB[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "H":
                  """ Superfamily F-Harbinger """
                  if TE_DTH[INDEX_CHANGE] == False:
                     TE_DTH[INDEX_COVERAGE] = int(1) + TE_DTH[INDEX_COVERAGE]
                     TE_DTH[INDEX_CHANGE] = True
               elif list(wickers_class)[2] == "C":
                  """ Superfamily CACTA """
                  if TE_DTC[INDEX_CHANGE] == False:
                     TE_DTC[INDEX_COVERAGE] = int(1) + TE_DTC[INDEX_COVERAGE]
                     TE_DTC[INDEX_CHANGE] = True
            if list(wickers_class)[1] == "T":
               """ Subclass I: Order Crypton """
               if TE_DYX[INDEX_CHANGE] == False:
                  TE_DYX[INDEX_COVERAGE] = int(1) + TE_DYX[INDEX_COVERAGE]
                  TE_DYX[INDEX_CHANGE] = True
               if list(wickers_class)[2] == "T":
                  """ Superfamily Crypton """
                  if TE_DYC[INDEX_CHANGE] == False:
                     TE_DYC[INDEX_COVERAGE] = int(1) + TE_DYC[INDEX_COVERAGE]
                     TE_DYC[INDEX_CHANGE] = True
            elif list(wickers_class)[1] == "H":
               """ Subclass I: Order Helitron """
               if TE_DHX[INDEX_CHANGE] == False:
                  TE_DHX[INDEX_COVERAGE] = int(1) + TE_DHX[INDEX_COVERAGE]
                  TE_DHX[INDEX_CHANGE] = True
               if list(wickers_class)[2] == "H":
                  """ Superfamily Helitron """
                  if TE_DHH[INDEX_CHANGE] == False:
                     TE_DHH[INDEX_COVERAGE] = int(1) + TE_DHH[INDEX_COVERAGE]
                     TE_DHH[INDEX_CHANGE] = True
            elif list(wickers_class)[1] == "H":
               """ Subclass II: Order Maverick """
               if TE_DMX[INDEX_CHANGE] == False:
                  TE_DMX[INDEX_COVERAGE] = int(1) + TE_DMX[INDEX_COVERAGE]
                  TE_DMX[INDEX_CHANGE] = True
               if list(wickers_class)[2] == "H":
                  """ Superfamily Maverick """
                  if TE_DMM[INDEX_CHANGE] == False:
                     TE_DMM[INDEX_COVERAGE] = int(1) + TE_DMM[INDEX_COVERAGE]
                     TE_DMM[INDEX_CHANGE] = True
      elif t_element.split("@")[0] == "blast":
         wickers_class = t_element.split("@")[1]
         wickers_class = wickers_class.split(":")
         if (wickers_class)[1] == "ClassI":
            """ Class I (retrotransposons)"""
            if TE_RXX[INDEX_CHANGE] == False: 
               TE_RXX[INDEX_COVERAGE] = int(1) + TE_RXX[INDEX_COVERAGE]
               TE_RXX[INDEX_CHANGE] = True       
            if (wickers_class)[2] == "LTR":
               """ Order LTR """
               if TE_RLX[INDEX_CHANGE] == False:
                  TE_RLX[INDEX_COVERAGE] = int(1) + TE_RLX[INDEX_COVERAGE]
                  TE_RLX[INDEX_CHANGE] = True
               if (wickers_class)[3] == "Copia":
                  """ Superfamily Copia """
                  if TE_RLC[INDEX_CHANGE] == False:
                     TE_RLC[INDEX_COVERAGE] = int(1) + TE_RLC[INDEX_COVERAGE]
                     TE_RLC[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "Gypsy":
                  """ Superfamily Gypsy """
                  if TE_RLG[INDEX_CHANGE] == False:
                     TE_RLG[INDEX_COVERAGE] = int(1) + TE_RLG[INDEX_COVERAGE]
                     TE_RLG[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "Bel-Pao":
                  """ Superfamily Bel-Pao """
                  if TE_RLB[INDEX_CHANGE] == False:
                     TE_RLB[INDEX_COVERAGE] = int(1) + TE_RLB[INDEX_COVERAGE]
                     TE_RLB[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "Retrovirus":
                  """ Superfamily Retrovirus """
                  if TE_RLR[INDEX_CHANGE] == False:
                     TE_RLR[INDEX_COVERAGE] = int(1) + TE_RLR[INDEX_COVERAGE]
                     TE_RLR[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "ERV":
                  """ Superfamily ERV """
                  if TE_RLE[INDEX_CHANGE] == False:
                     TE_RLE[INDEX_COVERAGE] = int(1) + TE_RLE[INDEX_COVERAGE]
                     TE_RLE[INDEX_CHANGE] = True
            elif (wickers_class)[2] == "DIRS":
               """ Order DIRS """
               if TE_RYX[INDEX_CHANGE] == False:
                  TE_RYX[INDEX_COVERAGE] = int(1) + TE_RYX[INDEX_COVERAGE]
                  TE_RYX[INDEX_CHANGE] = True
               if (wickers_class)[3] == "DIRS":
                  """ Superfamily DIRS """
                  if TE_RYD[INDEX_CHANGE] == False:
                     TE_RYD[INDEX_COVERAGE] = int(1) + TE_RYD[INDEX_COVERAGE]
                     TE_RYD[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "Ngaro":
                  """ Superfamily Ngaro """
                  if TE_RYN[INDEX_CHANGE] == False:
                     TE_RYN[INDEX_COVERAGE] = int(1) + TE_RYN[INDEX_COVERAGE]
                     TE_RYN[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "VIPER":
                  """ Superfamily VIPER """
                  if TE_RYV[INDEX_CHANGE] == False:
                     TE_RYV[INDEX_COVERAGE] = int(1) + TE_RYV[INDEX_COVERAGE]
                     TE_RYV[INDEX_CHANGE] = True
            elif (wickers_class)[2] == "PLE":
               """ Order PLE """
               if TE_RPX[INDEX_CHANGE] == False:
                  TE_RPX[INDEX_COVERAGE] = int(1) + TE_RPX[INDEX_COVERAGE]
                  TE_RPX[INDEX_CHANGE] = True
               if (wickers_class)[3] == "Penelope":
                  """ Superfamily Penelope """
                  if TE_RPP[INDEX_CHANGE] == False:
                     TE_RPP[INDEX_COVERAGE] = int(1) + TE_RPP[INDEX_COVERAGE]
                     TE_RPP[INDEX_CHANGE] = True
            elif (wickers_class)[2] == "LINE":
               """ Order LINE """
               if TE_RIX[INDEX_CHANGE] == False:
                  TE_RIX[INDEX_COVERAGE] = int(1) + TE_RIX[INDEX_COVERAGE]
                  TE_RIX[INDEX_CHANGE] = True
               if (wickers_class)[3] == "R2":
                  """ Superfamily R2 """
                  if TE_RIR[INDEX_CHANGE] == False:
                     TE_RIR[INDEX_COVERAGE] = int(1) + TE_RIR[INDEX_COVERAGE]
                     TE_RIR[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "RTE":
                  """ Superfamily RTE """
                  if TE_RIT[INDEX_CHANGE] == False:
                     TE_RIT[INDEX_COVERAGE] = int(1) + TE_RIT[INDEX_COVERAGE]
                     TE_RIT[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "Jockey":
                  """ Superfamily Jockey """
                  if TE_RIJ[INDEX_CHANGE] == False:
                     TE_RIJ[INDEX_COVERAGE] = int(1) + TE_RIJ[INDEX_COVERAGE]
                     TE_RIJ[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "L1":
                  """ Superfamily L1 """
                  if TE_RIL[INDEX_CHANGE] == False:
                     TE_RIL[INDEX_COVERAGE] = int(1) + TE_RIL[INDEX_COVERAGE]
                     TE_RIL[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "I":
                  """ Superfamily I """
                  if TE_RII[INDEX_CHANGE] == False:
                     TE_RII[INDEX_COVERAGE] = int(1) + TE_RII[INDEX_COVERAGE]
                     TE_RII[INDEX_CHANGE] = True
            elif (wickers_class)[2] == "SINE":
               """ Order SINE """
               if TE_RSX[INDEX_CHANGE] == False:
                  TE_RSX[INDEX_COVERAGE] = int(1) + TE_RSX[INDEX_COVERAGE]
                  TE_RSX[INDEX_CHANGE] = True
               if (wickers_class)[3] == "tRNA":
                  """ Superfamily tRNA """
                  if TE_RST[INDEX_CHANGE] == False:
                     TE_RST[INDEX_COVERAGE] = int(1) + TE_RST[INDEX_COVERAGE]
                     TE_RST[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "7SL":
                  """ Superfamily 7SL """
                  if TE_RSL[INDEX_CHANGE] == False:
                     TE_RSL[INDEX_COVERAGE] = int(1) + TE_RSL[INDEX_COVERAGE]
                     TE_RSL[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "5S":
                  """ Superfamily 5S """
                  if TE_RSS[INDEX_CHANGE] == False:
                     TE_RSS[INDEX_COVERAGE] = int(1) + TE_RSS[INDEX_COVERAGE]
                     TE_RSS[INDEX_CHANGE] = True
         if (wickers_class)[1] == "ClassII":
            """ Class II (DNA transposons)"""
            if TE_DXX[INDEX_CHANGE] == False: 
               TE_DXX[INDEX_COVERAGE] = int(1) + TE_DXX[INDEX_COVERAGE]
               TE_DXX[INDEX_CHANGE] = True
            if (wickers_class)[2] == "TIR":
               """ Subclass I: Order TIR """
               if TE_DTX[INDEX_CHANGE] == False:
                  TE_DTX[INDEX_COVERAGE] = int(1) + TE_DTX[INDEX_COVERAGE]
                  TE_DTX[INDEX_CHANGE] = True
               if (wickers_class)[3] == "Tc1-Mariner":
                  """ Superfamily Tc1-Mariner """
                  if TE_DTT[INDEX_CHANGE] == False:
                     TE_DTT[INDEX_COVERAGE] = int(1) + TE_DTT[INDEX_COVERAGE]
                     TE_DTT[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "hAT":
                  """ Superfamily hAT """
                  if TE_DTA[INDEX_CHANGE] == False:
                     TE_DTA[INDEX_COVERAGE] = int(1) + TE_DTA[INDEX_COVERAGE]
                     TE_DTA[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "Mutator":
                  """ Superfamily Mutator """
                  if TE_DTM[INDEX_CHANGE] == False:
                     TE_DTM[INDEX_COVERAGE] = int(1) + TE_DTM[INDEX_COVERAGE]
                     TE_DTM[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "Merlin":
                  """ Superfamily Merlin """
                  if TE_DTE[INDEX_CHANGE] == False:
                     TE_DTE[INDEX_COVERAGE] = int(1) + TE_DTE[INDEX_COVERAGE]
                     TE_DTE[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "Transib":
                  """ Superfamily Transib """
                  if TE_DTR[INDEX_CHANGE] == False:
                     TE_DTR[INDEX_COVERAGE] = int(1) + TE_DTR[INDEX_COVERAGE]
                     TE_DTR[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "P":
                  """ Superfamily P """
                  if TE_DTP[INDEX_CHANGE] == False:
                     TE_DTP[INDEX_COVERAGE] = int(1) + TE_DTP[INDEX_COVERAGE]
                     TE_DTP[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "PiggyBac":
                  """ Superfamily PiggyBac """
                  if TE_DTB[INDEX_CHANGE] == False:
                     TE_DTB[INDEX_COVERAGE] = int(1) + TE_DTB[INDEX_COVERAGE]
                     TE_DTB[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "F-Harbinger":
                  """ Superfamily F-Harbinger """
                  if TE_DTH[INDEX_CHANGE] == False:
                     TE_DTH[INDEX_COVERAGE] = int(1) + TE_DTH[INDEX_COVERAGE]
                     TE_DTH[INDEX_CHANGE] = True
               elif (wickers_class)[3] == "CACTA":
                  """ Superfamily CACTA """
                  if TE_DTC[INDEX_CHANGE] == False:
                     TE_DTC[INDEX_COVERAGE] = int(1) + TE_DTC[INDEX_COVERAGE]
                     TE_DTC[INDEX_CHANGE] = True
            if (wickers_class)[2] == "Crypton":
               """ Subclass I: Order Crypton """
               if TE_DYX[INDEX_CHANGE] == False:
                  TE_DYX[INDEX_COVERAGE] = int(1) + TE_DYX[INDEX_COVERAGE]
                  TE_DYX[INDEX_CHANGE] = True
               if (wickers_class)[3] == "Crypton":
                  """ Superfamily Crypton """
                  if TE_DYC[INDEX_CHANGE] == False:
                     TE_DYC[INDEX_COVERAGE] = int(1) + TE_DYC[INDEX_COVERAGE]
                     TE_DYC[INDEX_CHANGE] = True
            elif (wickers_class)[2] == "Helitron":
               """ Subclass I: Order Helitron """
               if TE_DHX[INDEX_CHANGE] == False:
                  TE_DHX[INDEX_COVERAGE] = int(1) + TE_DHX[INDEX_COVERAGE]
                  TE_DHX[INDEX_CHANGE] = True
               if (wickers_class)[3] == "Helitron":
                  """ Superfamily Helitron """
                  if TE_DHH[INDEX_CHANGE] == False:
                     TE_DHH[INDEX_COVERAGE] = int(1) + TE_DHH[INDEX_COVERAGE]
                     TE_DHH[INDEX_CHANGE] = True
            elif (wickers_class)[2] == "Maverick ":
               """ Subclass II: Order Maverick """
               if TE_DMX[INDEX_CHANGE] == False:
                  TE_DMX[INDEX_COVERAGE] = int(1) + TE_DMX[INDEX_COVERAGE]
                  TE_DMX[INDEX_CHANGE] = True
               if (wickers_class)[3] == "Maverick":
                  """ Superfamily Maverick """
                  if TE_DMM[INDEX_CHANGE] == False:
                     TE_DMM[INDEX_COVERAGE] = int(1) + TE_DMM[INDEX_COVERAGE]
                     TE_DMM[INDEX_CHANGE] = True
      elif t_element.split("@")[0] == "SSR":
         """ SSRs """
         if TE_SSR[INDEX_CHANGE] == False: 
            TE_SSR[INDEX_COVERAGE] = int(1) + TE_SSR[INDEX_COVERAGE]
            TE_SSR[INDEX_CHANGE] = True
   TE_RXX[INDEX_CHANGE] = False
   TE_RLX[INDEX_CHANGE] = False
   TE_RLC[INDEX_CHANGE] = False
   TE_RLG[INDEX_CHANGE] = False
   TE_RLB[INDEX_CHANGE] = False
   TE_RLR[INDEX_CHANGE] = False
   TE_RLE[INDEX_CHANGE] = False
   TE_RYX[INDEX_CHANGE] = False
   TE_RYD[INDEX_CHANGE] = False
   TE_RYN[INDEX_CHANGE] = False
   TE_RYV[INDEX_CHANGE] = False
   TE_RPX[INDEX_CHANGE] = False
   TE_RPP[INDEX_CHANGE] = False
   TE_RIX[INDEX_CHANGE] = False
   TE_RIR[INDEX_CHANGE] = False
   TE_RIT[INDEX_CHANGE] = False
   TE_RIJ[INDEX_CHANGE] = False
   TE_RIL[INDEX_CHANGE] = False
   TE_RII[INDEX_CHANGE] = False
   TE_RSX[INDEX_CHANGE] = False
   TE_RST[INDEX_CHANGE] = False
   TE_RSL[INDEX_CHANGE] = False
   TE_RSS[INDEX_CHANGE] = False
   TE_DXX[INDEX_CHANGE] = False
   TE_DTX[INDEX_CHANGE] = False
   TE_DTT[INDEX_CHANGE] = False
   TE_DTA[INDEX_CHANGE] = False
   TE_DTM[INDEX_CHANGE] = False
   TE_DTE[INDEX_CHANGE] = False
   TE_DTR[INDEX_CHANGE] = False
   TE_DTP[INDEX_CHANGE] = False
   TE_DTB[INDEX_CHANGE] = False
   TE_DTH[INDEX_CHANGE] = False
   TE_DTC[INDEX_CHANGE] = False
   TE_DYX[INDEX_CHANGE] = False
   TE_DYC[INDEX_CHANGE] = False
   TE_DHX[INDEX_CHANGE] = False
   TE_DHH[INDEX_CHANGE] = False
   TE_DMX[INDEX_CHANGE] = False
   TE_DMM[INDEX_CHANGE] = False
   TE_SSR[INDEX_CHANGE] = False

###
# End Contig Loop
###

# Write to comma separated value file
ALL_TE_BP_COVERAGE_DATA_FILE_NAME = CONTIG_ID + '_all_te_bp_coverage_data.txt'
with open(ALL_TE_BP_COVERAGE_DATA_FILE_NAME, 'w') as FILE:
   FILE.write("Class, Order, Superfamily, # of bp covered\n")
   FILE.write("Class I (retrotransposons),,, " + str(TE_RXX[INDEX_COVERAGE]) + "\n")
   FILE.write(",Order LTR,, " + str(TE_RLX[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily Copia, " + str(TE_RLC[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily Gypsy, " + str(TE_RLG[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily Bel-Pao, " + str(TE_RLB[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily Retrovirus, " + str(TE_RLR[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily ERV, " + str(TE_RLE[INDEX_COVERAGE]) + "\n")
   FILE.write(",Order DIRS,, " + str(TE_RYX[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily DIRS, " + str(TE_RYD[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily Ngaro, " + str(TE_RYN[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily VIPER, " + str(TE_RYV[INDEX_COVERAGE])  + "\n")
   FILE.write(",Order PLE,, " + str(TE_RPX[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily Penelope, " + str(TE_RPP[INDEX_COVERAGE])  + "\n")
   FILE.write(",Order LINE,, " + str(TE_RIX[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily R2, " + str(TE_RIR[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily RTE, " + str(TE_RIT[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily Jockey, " + str(TE_RIJ[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily L1, " + str(TE_RIL[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily I, " + str(TE_RII[INDEX_COVERAGE])  + "\n")
   FILE.write(",Order SINE,, " + str(TE_RSX[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily tRNA, " + str(TE_RST[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily 7SL, " + str(TE_RSL[INDEX_COVERAGE])  + "\n")
   FILE.write(",,Superfamily 5S, " + str(TE_RSS[INDEX_COVERAGE])  + "\n")
   FILE.write("Class II (DNA transposons),,, " + str(TE_DXX[INDEX_COVERAGE]) + "\n")
   FILE.write(",Subclass I: Order TIR,, " + str(TE_DTX[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily Tc1-Mariner, " + str(TE_DTT[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily hAT, " + str(TE_DTA[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily Mutator, " + str(TE_DTM[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily Merlin, " + str(TE_DTE[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily Transib, " + str(TE_DTR[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily P, " + str(TE_DTP[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily PiggyBac, " + str(TE_DTB[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily PIF-Harbinger, " + str(TE_DTH[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily CACTA, " + str(TE_DTC[INDEX_COVERAGE]) + "\n")
   FILE.write(",Subclass I: Order Crypton,, " + str(TE_DYX[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily Crypton, "  + str(TE_DYC[INDEX_COVERAGE]) + "\n")
   FILE.write(",Subclass I: Order Helitron,, " + str(TE_DHX[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily Helitron, " + str(TE_DHH[INDEX_COVERAGE]) + "\n")
   FILE.write(",Subclass I: Order Maverick,, " + str(TE_DMX[INDEX_COVERAGE]) + "\n")
   FILE.write(",,Superfamily Maverick, " + str(TE_DMM[INDEX_COVERAGE]) + "\n")
   FILE.write("SSRs,,, " + str(TE_SSR[INDEX_COVERAGE]) + "\n")
###
# End main()
###
