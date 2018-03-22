import os, re
from Bio import Entrez, SeqIO

from Protein import Protein
from FileManager import writeOutput, clearPreviousOutput


Entrez.email = 'bdighera@csu.fullerton.edu'
Entrez.tool = 'Ortholog Finder'

def Main():

    #Initializes the protein class which contains all main functions for manipulation and storage of protein data
    P = Protein()

    #Set input/output paths
    input_file_path = os.path.join('Input','ProteinInput')
    output_file_path = os.path.join('Output')

    #Clears previous CDS and genomic output from last run
    clearPreviousOutput(output_file_path)

    #Functions from the Protein class that use the input protein accession numbers to determine the corresponding protein ID, CDS, and genomic data
    P.Entrez_Protein_ID_Fetch(input_file_path)
    P.Entrez_Genome_Fetch()
    P.Entrez_CDS_Fetch()

    compiled_data = []

    #builds internal datastructure for each individual record containing all of each records corresponding data
    for i in range(len(P.input_protein_sequence)):
        compiled_data.append(['Protein Accession Number: ' + P.input_protein_accession_number[i],
                              'Protein Sequence: ' + P.input_protein_sequence[i],
                              'Protein ID Number: ' + P.retrieved_protein_ids[i],
                              'CDS: ' + P.retrieved_full_cds[i],
                              'Genome: ' + P.gene[i]])

        writeOutput(output_file_path, P.gene[i], P.retrieved_full_cds[i])

    #Function of the protein class that determines the intron phase and location of intron/exon boundry
    P.intronCalculator()

    print P.intron_phase
    print P.exon_lengths

    P.multiple_sequencing_alignment()
    P.renderingTreeImage()


if __name__ == '__main__':

    Main()


