'''
Created on Dec 4, 2013

@author: alecmacrae
'''
# project 4: chemoinformatics

# establishes the imports and arguments that should be used to run the script
import sys
import time

drugs_file = sys.argv[1]
targets_file = sys.argv[2]
output_file = sys.argv[3]
drugs = open(drugs_file)
targets = open(targets_file)
output = open(output_file, "w")

#--FUNCTIONS-------------------------------------

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

def findTanimotoValues(drugMap):
    """finds the tanimoto value for each set of pairs of drugs and creates a mapping of it: drug1 -> drug2 -> Tc"""
    Tmap = {}
    for key1 in sorted(drugMap):
        for key2 in sorted(drugMap):
            if int(key1.strip("DB")) < int(key2.strip("DB")):  # removing characters can also be done using filter(type(x).isdigit, x) where x is the key in this case
                if not Tmap.has_key(key1):
                    Tmap[key1] = {}
                Tc = calculateTanimotoValue(drugMap[key1],drugMap[key2])
                Tmap[key1][key2] = Tc
    return Tmap

def determineIfMatch(targetMap,key1,key2):
    """based on 2 different keys find out if two drugs share the same target or not"""
    if targetMap.has_key(key1) and targetMap.has_key(key2):
        if targetMap[key1] == targetMap[key2]:
            return 1
        else:
            return 0
    else:
        return 0

def writeOutput(tanimotoValues,targetMap):
    """writes the desired output to a file line by line"""
    for key1 in sorted(tanimotoValues):
        for key2 in sorted(tanimotoValues[key1]):
            match = determineIfMatch(targetMap,key1,key2)
            output.write("%s,%s,%.6f,%i\n" % (key1,key2,tanimotoValues[key1][key2],match))

def run():
    """runs the program"""
    drugMap = createMap(drugs)
    targetMap = createMap(targets)
    tanimotoValues = findTanimotoValues(drugMap)
    writeOutput(tanimotoValues,targetMap)
    
#--END-FUNCTIONS----------------------------------

run()
drugs.close()
targets.close()
output.close()