import os, subprocess, re
from Bio import SeqIO
from Protein import Protein



#The purpose of this file is to write/manage the output stream to output files for the proteins of interest
#This includes managing the movement of data into tmp files that will be piped into executible programs. Ie. spidey.linux


def writeOutput(path, gene, CDS):

    cdsFilePath = os.path.join(path, 'CDS')
    genomeFilePath = os.path.join(path, 'Genome')

    P = Protein()

    with open(genomeFilePath, 'a') as genomeFileWrite:
        genomeFileWrite.write(gene.format('fasta'))

    with open(cdsFilePath, 'a') as cdsFileWrite:
        cdsFileWrite.write(CDS.format('fasta'))

def clearPreviousOutput(path):

    cdsFilePath = os.path.join(path, 'CDS')
    genomeFilePath = os.path.join(path, 'Genome')

    with open(cdsFilePath, 'w') as cdsFileClear:
        cdsFileClear.write('')
        print 'CDS output file has been cleared, ready to write new records...'

    with open(genomeFilePath, 'w') as genomeFileClear:
        genomeFileClear.write('')
        print 'Genome output file has been cleared, ready to write new records...'


