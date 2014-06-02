import threading, Queue, random
import os, fnmatch
from subprocess import call

N_THREADS = 10 # define size of thread pool
GGVERSION = 'gg_97_otus_4feb2011'
GREENGENES = '/home/pah11816/Yatsunenko/Greengenes/gg_otus_4feb2011'
REF_SEQS = GREENGENES + '/rep_set/' + GGVERSION + '.fasta'
REF_TREE = GREENGENES + '/trees/' + GGVERSION + '.tre'
REF_TAX  = GREENGENES + '/taxonomies/greengenes_tax.txt'
PARAMFILE = '/home/pah11816/utils/qiime_params.conf'

# search for fasta files in $PWD
values = Queue.Queue()
for root, dirs, files in os.walk('.'):
    for basename in files:
        if fnmatch.fnmatch(basename, '*.fna'):
            filename = os.path.join(root, basename)
            values.put( filename ) # store found file in thread-safe queue

# define tasks to be executed by threads
def processQueue(threadId):
    while not values.empty():
        try:
            taskfile = values.get_nowait()
            taskdir = taskfile.rsplit('/',2)[0] + '/'
            print ''.join(['picking ', taskfile])
            call(['pick_closed_reference_otus.py', '-i', taskfile, '-r', REF_SEQS, '-p', PARAMFILE, '-t', REF_TAX, '-o', taskdir + 'picking'])
        except: # queue empty - done
            pass
    print ''.join(["Queue empty, thread ", str(threadId), " exiting..."])

# create and launch thread pool
for i in range(N_THREADS):
    t = threading.Thread(target=processQueue, args=[i])
    t.start()
