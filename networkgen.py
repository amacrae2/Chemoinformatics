'''
Created on Dec 6, 2013

@author: alecmacrae
'''
# project 4: chemoinformatics

# establishes the imports and arguments that should be used to run the script
import sys
import random
import time

#constants
N = 100
p = 0.5
pCutoff = 0.05

drugs_file = sys.argv[1]
targets_file = sys.argv[2]
protein_nodes_file = sys.argv[3]
drugs = open(drugs_file)
targets = open(targets_file)
protein_nodes = open(protein_nodes_file)
network_output_file = "network.sif"
name_output_file = "name.nodeAttr"
indication_output_file = "indication.nodeAttr"
networkOF = open(network_output_file, "w")
nameOF = open(name_output_file, "w")
indicationOF = open(indication_output_file, "w")

#--FUNCTIONS-------------------------------------

def createproteinToDrugsMap(targets):
    """creates a map mapping from protein -> drugs"""
    result = {}
    lineNumber = 1
    for line in targets:
        if lineNumber != 1:
            line = line.strip().split(",")
            if not result.has_key(line[1]):
                result[line[1]] = []
            result[line[1]].append(line[0])
        lineNumber+=1
    return result

def createMap(data):
    """creates a map mapping the first column of the input data to a list of the third-end pieces of data"""
    result = {}
    lineNumber = 1
    for line in data:
        if lineNumber != 1:
            line = line.strip().split(",")
            result[line[0]] = line[2].split()
        lineNumber+=1
    return result

def createProteinMaps(protein_nodes):
    """creates 2 maps mapping from protein -> name and protein -> indications"""
    proteinToNameMap = {}
    proteinToIndicationsMap = {}
    lineNumber = 1
    for line in protein_nodes:
        if lineNumber != 1:
            line = line.strip().split(",")
            proteinToNameMap[line[0]] = line[1]
            proteinToIndicationsMap[line[0]] = line[2]
        lineNumber+=1
    return (proteinToNameMap,proteinToIndicationsMap)

def calculateTanimotoValue(drugOne,drugTwo):
    """find the Tanimoto constant value by calculating it from the number of shared bits"""
    both = len(set(drugOne).intersection(drugTwo))
    one = len(drugOne)
    two = len(drugTwo)
    either = one + two - both
    Tc = round(float(both)/float(either),6)
    return Tc

def computeTsummary(drugMap,drugsA,drugsB):
    """calculates the Tsummary value based on the equation in the project description"""
    Tsumm = 0.0
    for drug1 in drugsA:
        for drug2 in drugsB:
            Tc = calculateTanimotoValue(drugMap[drug1],drugMap[drug2])
            if Tc > p:
                Tsumm += Tc
    return Tsumm

def findTsummary(proteinToDrugsMap,drugMap,proteinA,proteinB):
    """finds the Tsummary value after creating a list of drugs that bind to each of the designated proteins"""
    drugsA = proteinToDrugsMap[proteinA]
    drugsB = proteinToDrugsMap[proteinB]
    Tsumm = computeTsummary(drugMap,drugsA,drugsB)
    return Tsumm

def findPBootstrap(proteinToDrugsMap, drugMap, Tsumm, proteinA, proteinB):
    """finds the value for P-bootstrap based on the calculation in the project discription"""
    nA = len(proteinToDrugsMap[proteinA])
    nB = len(proteinToDrugsMap[proteinB])
    countRandTsummGreaterThanTsumm = 0
    for i in xrange(N):
        randDrugsA = random.sample(drugMap.keys(),nA)
        randDrugsB = random.sample(drugMap.keys(),nB)
        TsummTemp = computeTsummary(drugMap,randDrugsA,randDrugsB)
        if TsummTemp >= Tsumm:
            countRandTsummGreaterThanTsumm += 1
    pBootstrap = round(float(countRandTsummGreaterThanTsumm)/float(N),6)
    return pBootstrap

def writeNetworkOF(proteinToDrugsMap,drugMap,proteinToNameMap):
    """writes to the network.sif output file as described in the project description"""
    for proteinA in sorted(proteinToNameMap.keys()):
        for proteinB in sorted(proteinToNameMap.keys()):
            if proteinA != proteinB and proteinA == min(proteinA,proteinB):
                Tsumm = findTsummary(proteinToDrugsMap,drugMap,proteinA,proteinB)
                pBootstrap = findPBootstrap(proteinToDrugsMap,drugMap,Tsumm,proteinA,proteinB)
                if pBootstrap <= pCutoff:
                    networkOF.write("%s edge %s\n" % (proteinA,proteinB))

def writeNodeAttr(proteinMap, header, fileName):
    fileName.write(header + "\n")
    for protein in proteinMap:
        fileName.write("%s = %s\n" % (protein,proteinMap[protein]))

def writeFiles(proteinToDrugsMap,drugMap,proteinToNameMap,proteinToIndicationsMap):
    writeNetworkOF(proteinToDrugsMap,drugMap,proteinToNameMap)
    writeNodeAttr(proteinToNameMap,"name",nameOF)
    writeNodeAttr(proteinToIndicationsMap,"indication",indicationOF)

def run():
    """runs the program"""
    proteinToDrugsMap = createproteinToDrugsMap(targets)
    drugMap = createMap(drugs)
    proteinMaps = createProteinMaps(protein_nodes)
    writeFiles(proteinToDrugsMap,drugMap,proteinMaps[0],proteinMaps[1])
    
    
#--END-FUNCTIONS----------------------------------

run()
drugs.close()
targets.close()
protein_nodes.close()
networkOF.close()
nameOF.close()
indicationOF.close()