import os
import sys
import math
import random
from sklearn import svm
import csv
import numpy as np
from TEMPy.MapParser import MapParser
from scipy import ndimage

aAList = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU',
          'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

def learn(workingDir):
    '''
    Learn after preparations
    :param workingDir: working directory, where protein_pdbs directory is found, e.g. D:/university/biology
    :return: nothing
    '''

    generate_csv_file(workingDir)

    #read the CSV file and generate the training set

    test_ARG(workingDir)
    test_all_categories(workingDir)


def test_ARG(workingDir):
    csvFilePath = workingDir+"/simulated/training_set.csv"
    featureList = []
    labelList = []
    testingSet = []

    totalArgTraining = 0.0
    totalNonArgTraining = 0.0

    print 'Generate training and testing sets...'
    with open(csvFilePath, "r") as csvFile:
        csvFileReader = csv.reader(csvFile, delimiter=',', quotechar='"')
        for emDataRow in csvFileReader:
            if (emDataRow[0] == "ARG" and random.random() > 0.4) or (random.random() > 0.96):
               featureList.append(emDataRow[1:])
               #labelList.append(emDataRow[0])

               if emDataRow[0] == "ARG":
                   labelList.append(1)
                   totalArgTraining += 1
               else:
                   labelList.append(0)
                   totalNonArgTraining += 1
            else:
                if emDataRow[0] == "ARG":
                    testingSet.append(emDataRow)
                else:
                    if random.random() > 0.97:
                        testingSet.append(emDataRow)

    X = np.array(featureList)
    y = labelList

    print 'Learning...'
    # kernel must be one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed' or a callable
    clf = svm.SVC(kernel='linear', C = 1.0)
    clf.fit(X,y)

    falsePositive = 0.0
    rightPositive = 0.0
    totalArgTesting = 0.0
    totalNonArgTesting = 0.0

    print 'Testing...'

    for testingRow in testingSet:
        prediction = clf.predict([testingRow[1:]])[0]
        #print "AA is {0}, predicted {1}".format(testingRow[0],prediction)
        if testingRow[0] == "ARG":
            totalArgTesting += 1
        else:
            totalNonArgTesting += 1
        if testingRow[0] == "ARG" and prediction==1:
            rightPositive += 1
        if testingRow[0] != "ARG" and prediction==1:
            falsePositive += 1

    falsePositivePercent = falsePositive/totalNonArgTesting
    rightPositivePercent = rightPositive/totalArgTesting

    print "ARG in training set: {0}".format(totalArgTraining)
    print "Non ARG in training set: {0}".format(totalNonArgTraining)

    print "ARG in testing set: {0}".format(totalArgTesting)
    print "Non ARG in testing set: {0}".format(totalNonArgTesting)

    print "False Positive Percent: {0:8.8f}".format(falsePositivePercent)
    print "Right Positive Percent: {0:8.8f}".format(rightPositivePercent)


def test_all_categories(workingDir):
    csvFilePath = workingDir+"/simulated/training_set.csv"
    featureList = []
    labelList = []
    testingSet = []

    print 'Generate training and testing sets...'
    with open(csvFilePath, "r") as csvFile:
        csvFileReader = csv.reader(csvFile, delimiter=',', quotechar='"')
        for emDataRow in csvFileReader:
            if random.random() > 0.7:
                featureList.append(emDataRow[1:])
                labelList.append(emDataRow[0])
            else:
                testingSet.append(emDataRow)
    X = np.array(featureList)
    y = labelList

    print 'Learning...'
    # kernel must be one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed' or a callable
    clf = svm.SVC(kernel='linear', C = 1.0)
    clf.fit(X,y)

    num_success = 0

    print 'Testing...'

    for testingRow in testingSet:
        prediction = clf.predict([testingRow[1:]])[0]
        if prediction == testingRow[0]:
            num_success += 1
    print 'Success in {} of {} ({}%)'.format(num_success, len(testingSet), num_success*100./len(testingSet) )


def features(em):
    em[ em < 1 ] = 0 # threshold. do not allow negative/small values
    zoomed_em = ndimage.interpolation.zoom(em, 0.2)
    zoomed_em[zoomed_em < 0] = 0 # threshold again (this is necessary, otherwise negatives sneak in)
    com = np.array(ndimage.measurements.center_of_mass( zoomed_em ))
    return [avg_dist_power(zoomed_em, com, exp) for exp in xrange(1,8)]


def avg_dist_power(em, ref_point, exp):
    ret = 0
    for x in xrange(len(em)):
        for y in xrange(len(em[x])):
            for z in xrange(len(em[x][y])):
                ret += np.linalg.norm( np.array([x,y,z]) - ref_point, exp) * em[x][y][z]
    ret /= len(em)*len(em[0])*len(em[0][0])
    return pow(ret , 1./exp)


def findMaxEMMapDimensions(workingDir):
    '''
    Find the maximal EM map dimensions
    :param workingDir:
    :return:
    '''
    print 'Finding the maximal EM map dimensions...'
    emDirectory = workingDir+"/simulated/EM"
    maxEmMapSize = [0,0,0]
    for aaDirName in os.listdir(emDirectory):
        if aaDirName in aAList:
            emAaDir = "{0}/{1}".format(emDirectory, aaDirName)
            for emfileName in os.listdir(emAaDir):
                emFilePath = "{0}/{1}".format(emAaDir,emfileName)
                emMap=MapParser.readMRC(emFilePath)
                emMapSize = emMap.box_size()
                if emMapSize[0]>maxEmMapSize[0]:
                    maxEmMapSize[0] = emMapSize[0]
                if emMapSize[1]>maxEmMapSize[1]:
                    maxEmMapSize[1] = emMapSize[1]
                if emMapSize[2]>maxEmMapSize[2]:
                    maxEmMapSize[2] = emMapSize[2]

    return maxEmMapSize


def numOfFilesSubdir(path):
    return sum([len(files) for r, d, files in os.walk(path)])


def generate_csv_file(workingDir):
    emDirectory = workingDir+"/simulated/EM"
    numFiles = numOfFilesSubdir(emDirectory)
    currFileNum = 0
    #generate CSV file of training set
    print 'Generate CSV file of training set...'
    csv_file_path = workingDir+"/simulated/training_set.csv"
    with open(csv_file_path, "w") as csvFile:
        csvFileWriter = csv.writer(csvFile, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for aaDirName in os.listdir(emDirectory):
            if aaDirName in aAList:
                emAaDir = "{0}/{1}".format(emDirectory, aaDirName)
                for emfileName in os.listdir(emAaDir):
                    emFilePath = "{0}/{1}".format(emAaDir,emfileName)
                    emMap = MapParser.readMRC(emFilePath)
                    emMapList = emMap.getMap()

                    csvFileWriter.writerow([aaDirName] + features(emMapList))

                    currFileNum += 1
                    sys.stdout.write('\r{:4}/{:4} ({:5.4}%), current file: {}'.format(currFileNum, numFiles, currFileNum*100./numFiles, emFilePath))# comma to suppress the newline
                    sys.stdout.flush()
    print '\n',