import re
import os
import sys
from TEMPy.StructureParser import PDBParser
from TEMPy.StructureBlurrer import StructureBlurrer

aAList = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU',
              'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

def prepare(workingDir):
    '''
    Prepare amino acids files
    :param workingDir: working directory, where protein_pdbs directory is found, e.g. D:/university/biology
    :return: nothing
    '''
    if workingDir[-1] == '/':
        workingDir = workingDir[:-1]
    generate_AA_PDBs_from_protein_PDBs(workingDir)
    print 'Generating normalized AA PDBs out of AA PDBs...'
    generate_normalized_AA_PDBs_from_AA_PDBs(workingDir)
    print 'generate EMs out of normalized PDBs'
    generate_EMs_from_normalized_AA_PDBs(workingDir)



def generate_AA_PDBs_from_protein_PDBs(workingDir):
    curAaSeq = 0
    curAaCode = 'UNKNOWN'
    numFiles = len(os.listdir(workingDir+'/protein_pdbs'))
    currFileNum = 0
    print 'Generating AA PDBs out of protein PDBS...'
    for proteinFileName in os.listdir(workingDir+'/protein_pdbs'):
        proteinFullFileName = "{0}/protein_pdbs/{1}".format(workingDir, proteinFileName)
        if os.path.isfile(proteinFullFileName):
            #print (proteinFullFileName)
            with open(proteinFullFileName) as proteinFile:
                proteinName = 'unknown'
                content = proteinFile.readlines()
                for line in content:
                    if line[:5]=='DBREF':
                        matchColumns = re.match( r'(\S+)(\s+)(\S+)', line, re.I)
                        proteinName = matchColumns.group(3)
                    if line[:4]=='ATOM':
                        matchColumns = re.match( r'(\S+)(\s+)(\S+)(\s+)(\S+)(\s+)(\S+)(\s+)(\S+)(\s+)(\S+)', line, re.I)
                        prevAaSeq = curAaSeq
                        prevAaCode = curAaCode
                        curAaSeq = matchColumns.group(11)
                        curAaCode = matchColumns.group(7)

                        #print "\rAA is %s" % matchColumns.group(7),
                        fileMode = 'a'
                        if prevAaSeq != curAaSeq:
                            fileMode = 'w'
                        if curAaCode in aAList:
                            pdbDdirectory = "%s/simulated/PDB/%s" % (workingDir, curAaCode)
                            if not os.path.exists(pdbDdirectory):
                                os.makedirs(pdbDdirectory)
                            aAPdbfileName = "{3}/simulated/PDB/{0}/{1}-{2}.pdb".format(curAaCode, proteinName, curAaSeq, workingDir)
                            with open(aAPdbfileName, fileMode) as aAfile:
                                aAfile.write(line)
        currFileNum += 1
        sys.stdout.write('\r{:4}/{:4} ({:5.4}%), current protein file: {}'.format(currFileNum, numFiles, currFileNum*100./numFiles, proteinFullFileName))
        sys.stdout.flush()
    print '\n',

def generate_normalized_AA_PDBs_from_AA_PDBs(workingDir):
    # generate normalized AA PDBs out of AA PDBs
    # (currently only normalized by moving c-alpha to 0,0,0)
    numFiles = numOfFilesSubdir(workingDir+'/simulated/PDB')
    currFileNum = 0
    pdbNormalizedDir = workingDir+"/simulated/PDB_normalized"
    for aaDirName in os.listdir(workingDir+'/simulated/PDB'):
        if aaDirName in aAList:
            for pdbFileName in os.listdir('{1}/simulated/PDB/{0}'.format(aaDirName, workingDir)):
                pdbFilePath = '{2}/simulated/PDB/{0}/{1}'.format(aaDirName,pdbFileName,workingDir)

                with open(pdbFilePath) as pdbFile:
                    content = pdbFile.readlines()
                    cAlphaPosition = False
                    for line in content:
                        if line[:4]=='ATOM':
                            atomName = line[12:15].strip()
                            if atomName=='CA':
                                cAlphaPosition = [float(line[30:37].strip()),float(line[38:45].strip()),float(line[46:53].strip())]

                    if cAlphaPosition:
                        aaNormalizedDir = "{0}/{1}".format(pdbNormalizedDir, aaDirName)
                        if not os.path.exists(aaNormalizedDir):
                            os.makedirs(aaNormalizedDir)
                        pdbNormalizedFilePath = "{0}/{1}".format(aaNormalizedDir,pdbFileName)

                        with open(pdbNormalizedFilePath, "w") as pdbNormalizedFile:
                            for line in content:
                                if line[:4]=='ATOM':
                                    posX = float(line[30:37].strip())
                                    posY = float(line[38:45].strip())
                                    posZ = float(line[46:53].strip())
                                    normalizedLine = line[:29]+"{0:8.3f}".format(posX-cAlphaPosition[0])+"{0:8.3f}".format(posY-cAlphaPosition[1])+"{0:8.3f}".format(posZ-cAlphaPosition[2])+line[54:]
                                    pdbNormalizedFile.write(normalizedLine)
                currFileNum += 1
                sys.stdout.write('\r{:4}/{:4} ({:5.4}%), current pdb file: {}'.format(currFileNum, numFiles, currFileNum*100./numFiles, pdbFilePath))
                sys.stdout.flush()
    print '\n',


def generate_EMs_from_normalized_AA_PDBs(workingDir):
    structureBlurrer = StructureBlurrer()
    emDirectory = workingDir+"/simulated/EM"
    pdbNormalizedDir = workingDir+ "/simulated/PDB_normalized"
    numFiles = numOfFilesSubdir(emDirectory)
    currFileNum = 0

    for aaDirName in os.listdir(pdbNormalizedDir):
        if aaDirName in aAList:
            emAaDir = "{0}/{1}".format(emDirectory, aaDirName)
            if not os.path.exists(emAaDir):
                os.makedirs(emAaDir)
            pdbPath = "{0}/{1}".format(pdbNormalizedDir,aaDirName)
            for pdbFileName in os.listdir(pdbPath):
                pdbFileNameMatch = re.match( r'(\S+)\.pdb', pdbFileName, re.I)
                pdbFileNameWoExtension = pdbFileNameMatch.group(1)
                pdbFilePath = "{0}/{1}".format(pdbPath,pdbFileName)
                aaStruture=PDBParser.read_PDB_file(pdbFileNameWoExtension,pdbFilePath)
                #aaSimMap = structureBlurrer.gaussian_blur(aaStruture, 2)
                aaSimMap = structureBlurrer.gaussian_blur_box(aaStruture, 2, 50, 50, 50)

                aAEmFileName = "{0}/{1}.map".format(emAaDir,pdbFileNameWoExtension)
                aaSimMap.write_to_MRC_file(aAEmFileName)
                currFileNum += 1
                sys.stdout.write('\r{:4}/{:4} ({:5.4}%), current file: {}'.format(currFileNum, numFiles, currFileNum*100./numFiles, aAEmFileName))
                sys.stdout.flush()
    print '\n',

def numOfFilesSubdir(path):
    return sum([len(files) for r, d, files in os.walk(path)])