#!/usr/bin/python  
'''
The window_size means the window that get the fasta, but for the insertion and deletion, expecially deletion, we make the window size 100bp lager than offered window size.
'''

from Bio import SeqIO
import numpy as np
import sys
import h5py
import math
from os.path import basename,dirname

writeFasta=True
writevcf=True

'''
#if vcf_original_allele_check is False, always use the ref and alt alleles user specified. 
#Otherwise ref allele must match the reference genome
'''
vcf_original_allele_check=False
#==================================Load seq===================================

inputwindow = int(sys.argv[2])

#else:
#    inputwindow = 300 
mutpos = inputwindow/2-1
offeredwindow = inputwindow - 100
fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
np.random.seed(1)

#=================================Read seq======================
seqs=[str(fasta.seq) for fasta in fasta_sequences]

oris=[]
muts=[]
chrs=[]
poss=[]
annos=[]
names=[]
fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
for fasta in fasta_sequences:
    anno = fasta.name.split('_')
    annos.append(fasta.name)
    oris.append(anno[0])
    muts.append(anno[1])
    chrs.append(anno[2])
    poss.append(int(anno[3])+mutpos)
    if len(anno)>5:
        names.append('_'.join(anno[5:]))

oris=np.asarray(oris)
muts=np.asarray(muts)
chrs=np.asarray(chrs)
seqs=np.asarray(seqs)
poss=np.asarray(poss)
annos=np.asarray(annos)
names=np.asarray(names)

seqsmut=[]
inds=[]
for i in range(len(seqs)):
    if vcf_original_allele_check:
        if type(oris[i]) is not np.string_ or seqs[i][mutpos:(mutpos+len(oris[i]))]!=oris[i]:
            continue
    else:
        if type(oris[i]) is not np.string_:
            continue
    inds.append(i)
#Get the offeredwindow length
    mutseq = seqs[i][:mutpos]+ muts[i]+seqs[i][(mutpos+len(oris[i])):] 
    seqsmut.append(mutseq[int(math.floor(((len(mutseq)-offeredwindow)/2.0))):int(math.floor(len(mutseq)-(len(mutseq)-offeredwindow)/2.0))])
    oriseq = seqs[i][:mutpos]+ oris[i]+seqs[i][(mutpos+len(oris[i])):]
    seqs[i] = oriseq[int(math.floor(((len(oriseq)-offeredwindow)/2.0))):int(math.floor(len(oriseq)-(len(oriseq)-offeredwindow)/2.0))] 

print "Number of input variants:"
print(len(inds))
print "Number of valid variants:"
print(len(seqs))
inds = np.asarray(inds)
seqsmut = np.asarray(seqsmut)
seqs = seqs[inds]
oris = oris[inds]
muts = muts[inds]
chrs = chrs[inds]
poss = poss[inds]
annos = annos[inds]
if len(names)>0:
    names = names[inds]

#=================================Write out file================================================

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

if writeFasta:
    allrecwt=[]
    allrecmut=[]
    fo = open(sys.argv[1]+'.ref.fasta','a')
    for i in range(len(seqsmut)):
        fo.write(seqs[i]+'\n')
    fo.close()
    fo = open(sys.argv[1].replace('wt','mut')+'.ref.fasta','a')
    for i in range(len(seqsmut)):
        fo.write(seqsmut[i]+'\n')
    fo.close()

if writevcf:
    myfile=open(sys.argv[1]+'.ref.vcf','w')
    if len(names)>0:
        for i in range(len(seqs)):
            myfile.write(chrs[i]+'\t'+str(poss[i])+'\t'+names[i]+'\t'+oris[i]+'\t'+muts[i]+'\n')
    else:
        for i in range(len(seqs)):
            myfile.write(chrs[i]+'\t'+str(poss[i])+'\t1\t'+oris[i]+'\t'+muts[i]+'\n')
    myfile.close()

