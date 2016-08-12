import os
import sys
import random
from sklearn import svm
import csv
import numpy as np
from TEMPy.MapParser import MapParser
from scipy import ndimage

aAList = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU',
          'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

def learn(workingDir, generateCSVFile = True):
    """
    Learn after preparations
    :param workingDir: working directory, where protein_pdbs directory is found, e.g. D:/university/biology
    :return: nothing
    """
    if generateCSVFile:
        generate_csv_file(workingDir)
    test_all_categories(workingDir)
    test_all_AAs_separately(workingDir)


def test_AA(workingDir, AA_name):
    print 'Testing AA: {0}'.format(AA_name)
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
            if (emDataRow[0] == AA_name and random.random() < 0.6) or random.random() < 0.05:
               featureList.append(emDataRow[1:])

               if emDataRow[0] == AA_name:
                   labelList.append(1)
                   totalArgTraining += 1
               else:
                   labelList.append(0)
                   totalNonArgTraining += 1
            else:
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
        if testingRow[0] == AA_name:
            totalArgTesting += 1
        else:
            totalNonArgTesting += 1
        if testingRow[0] == AA_name and prediction==1:
            rightPositive += 1
        if testingRow[0] != AA_name and prediction==1:
            falsePositive += 1

    falsePositivePercent = falsePositive/totalNonArgTesting
    rightPositivePercent = rightPositive/totalArgTesting

    print "{0} in training set: {1}".format(AA_name, int(totalArgTraining))
    print "Non {0} in training set: {1}".format(AA_name, int(totalNonArgTraining))

    print "{0} in testing set: {1}".format(AA_name, int(totalArgTesting))
    print "Non {0} in testing set: {1}".format(AA_name, int(totalNonArgTesting))

    print "False Positive Percent: {0:4.2f}%".format(100.*falsePositivePercent)
    print "Right Positive Percent: {0:4.2f}%".format(100.*rightPositivePercent)


def test_all_categories(workingDir):
    print '\nTesting all categories at once:'
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
    print 'Success in {} of {} ({:4.2f}%)'.format(num_success, len(testingSet), num_success*100./len(testingSet) )


def features(emMap):
    em = emMap.copy().normalise().getMap()
    em[ em < 1 ] = 0 # threshold. do not allow negative/small values
    em = ndimage.interpolation.zoom(em, 0.2)
    em[em < 0] = 0 # threshold again (this is necessary, otherwise negatives sneak in)
    return list(em.flatten())
    #com = np.array(ndimage.measurements.center_of_mass( em ))
    #return [avg_dist_power(em, com, exp) for exp in xrange(1,8)]# + list(zoomed_em.flatten())


def avg_dist_power(em, ref_point, exp):
    ret = 0
    for x in xrange(len(em)):
        for y in xrange(len(em[x])):
            for z in xrange(len(em[x][y])):
                ret += np.linalg.norm( np.array([x,y,z]) - ref_point, exp) * em[x][y][z]
    ret /= len(em)*len(em[0])*len(em[0][0])
    return pow(ret , 1./exp)


def findMaxEMMapDimensions(workingDir):
    """
    Find the maximal EM map dimensions
    :param workingDir:
    :return:
    """
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
    print 'Generate CSV file of features...'
    csv_file_path = workingDir+"/simulated/training_set.csv"
    with open(csv_file_path, "w") as csvFile:
        csvFileWriter = csv.writer(csvFile, delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for aaDirName in os.listdir(emDirectory):
            if aaDirName in aAList:
                emAaDir = "{0}/{1}".format(emDirectory, aaDirName)
                for emfileName in os.listdir(emAaDir):
                    emFilePath = "{0}/{1}".format(emAaDir,emfileName)
                    emMap = MapParser.readMRC(emFilePath)

                    csvFileWriter.writerow([aaDirName] + features(emMap))

                    currFileNum += 1
                    sys.stdout.write('\r{:4}/{:4} ({:5.4}%), current file: {}'.format(currFileNum, numFiles, currFileNum*100./numFiles, emFilePath))# comma to suppress the newline
                    sys.stdout.flush()
    print '\n',

def test_all_AAs_separately(workingDir):
    print '\nTesting all AAs separately...\n'
    for AA in aAList:
        test_AA(workingDir, AA)
        print '\n'