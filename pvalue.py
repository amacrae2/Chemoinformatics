'''
Created on Dec 5, 2013

@author: alecmacrae
'''
# project 4: chemoinformatics

# establishes the imports and arguments that should be used to run the script
import sys
import random
import time

#constants
p = 0.5

if sys.argv[1] == "-n":
    N = int(sys.argv[2])
    drugs_file = sys.argv[3]
    targets_file = sys.argv[4]
    proteinA = sys.argv[5]
    proteinB = sys.argv[6]
else:
    N = 100
    drugs_file = sys.argv[1]
    targets_file = sys.argv[2]
    proteinA = sys.argv[3]
    proteinB = sys.argv[4]
        
drugs = open(drugs_file)
targets = open(targets_file)

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

def calculateTanimotoValue(drugOne,drugTwo):
    """find the Tanimoto constant value by calculating it from the number of shared bits"""
    both = len(set(drugOne).intersection(drugTwo))
    one = len(drugOne)
    two = len(drugTwo)
    either = one + two - both
    Tc = round(float(both)/float(either),6)
    return Tc

def computeTsummary(proteinToDrugsMap,drugMap,drugsA,drugsB):
    """calculates the Tsummary value based on the equation in the project description"""
    Tsumm = 0
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
    Tsumm = computeTsummary(proteinToDrugsMap,drugMap,drugsA,drugsB)
    return Tsumm

def findPBootstrap(proteinToDrugsMap, drugMap, Tsumm, proteinA, proteinB):
    """finds the value for P-bootstrap based on the calculation in the project description"""
    nA = len(proteinToDrugsMap[proteinA])
    nB = len(proteinToDrugsMap[proteinB])
    countRandTsummGreaterThanTsumm = 0
    for i in xrange(N):
        randDrugsA = random.sample(drugMap.keys(),nA)
        randDrugsB = random.sample(drugMap.keys(),nB)
        TsummTemp = computeTsummary(proteinToDrugsMap,drugMap,randDrugsA,randDrugsB)
        if TsummTemp >= Tsumm:
            countRandTsummGreaterThanTsumm += 1
    pBootstrap = round(float(countRandTsummGreaterThanTsumm)/float(N),6)
    return pBootstrap

def run():
    """runs the program"""
    proteinToDrugsMap = createproteinToDrugsMap(targets)
    drugMap = createMap(drugs)
    Tsumm = findTsummary(proteinToDrugsMap,drugMap,proteinA,proteinB)
    pBootstrap = findPBootstrap(proteinToDrugsMap,drugMap,Tsumm,proteinA,proteinB)
    print pBootstrap
    
#--END-FUNCTIONS----------------------------------

run()
drugs.close()
targets.close()