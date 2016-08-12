

from time import time
start = time()

print 'Importing...',
import prepare
import learn
print 'done'

aAList = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU',
          'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
working_directory = 'C:/university/biology'
prepare.prepare(working_directory)
learn.learn(working_directory, True)
print 'Total time taken: {} seconds'.format(time()-start)

