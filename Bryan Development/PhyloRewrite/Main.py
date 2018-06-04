import os, sys
from Bio import Entrez

from Protein import Protein
from FileManager import writeOutput, clearPreviousOutput, repair_script, nameOfRun, timer, duplicate_management, welcomeStatement, file_merger

Entrez.email = 'bdighera@csu.fullerton.edu'
Entrez.tool = 'PhyloEpsilon: Protein Ortholog Finder'

def Main():

    #Initializes the protein class which contains all main functions for manipulation and storage of protein data
    P = Protein()
    # Gets user input for the name of the desired output file. Output graphic displays in the subfolder of completed trees
    name_of_run = nameOfRun()
    P.run_name = name_of_run
    # timer that determines how much pause is placed between calls to NCBI servers, reccomend 0.5 however 0 works when running NOT during peak hours
    time = timer()
    P.timer = time
    

    #Set input/output paths
    input_file_path = os.path.join('Input','ProteinInput')
    output_file_path = os.path.join('Output')

    #Clears previous CDS and genomic output from last run
    clearPreviousOutput(output_file_path)

    #Functions from the Protein class that use the input protein accession numbers to determine the corresponding protein ID, CDS, and genomic data
    P.Entrez_Protein_ID_Fetch(input_file_path)

    P.Entrez_Genome_Fetch()
    # Writes the fasta seq for the genome corresponding the protein. Establishing a parallel list. Ouput file path: PhyloRewrite/Output/Genome
    for i in range(len(P.gene)):
        data_type = 'gene'
        writeOutput(data_type, P.gene[i])
    print 'Genomic Sequences written to Output File... '

    P.Entrez_CDS_Fetch()
    #Writes the fasta seq for the CDS corresponding the protein. Establishing a parallel list. Ouput file path: PhyloRewrite/Output/CDS
    for i in range(len(P.retrieved_full_cds)):
        data_type = 'CDS'
        writeOutput(data_type, P.retrieved_full_cds[i])
    print 'CDS Sequences written to Output File... '

    #Checks as to whether the lists are in parallel, if they are not error will rise
    print '\n Checking the output for errors... \n'

    if len(P.gene) != len(P.retrieved_full_cds):
        repair_script()
        sys.exit()

    print 'No errors found, continuing... \n'


    #Function of the protein class that determines the intron phase and location of intron/exon boundry
    P.intronCalculator()

    #Contained within the lists is each
    print 'Intron Phases: ' + str(P.intron_phase)
    print 'Length of the Exons: ' + str(P.exon_lengths)

    #Uses clustal X linux executible to format a multiple sequencing alignment for the input protein sequences. Files found in execs/tmp
    P.multiple_sequencing_alignment()
    #Takes the multiple sequencing alignment output from clustal X and inputs into fasttree to generate an unrooted tree, then piped into ete2 to root. Files found in execs/tmp
    P.rootedTreeConstruction()
    #Builds the tree graphic
    P.renderingTreeImage()


if __name__ == '__main__':

    response = welcomeStatement()


    if response == '1':

        #Function that runs the entire program
        Main()

    elif response == '2':

        #Use following code if you have the Protein, CDS, and genomic sequences
        #Code produces the Phylogenetic tree and intron mapping
        P = Protein()
        P.multiple_sequencing_alignment()
        P.intronCalculator()
        P.rootedTreeConstruction()
        P.renderingTreeImage()

    elif response == '3':

        repair_script()

    elif response == '4':

        duplicate_management()

    elif response == '5':
        #Enter experimental code here:

        P = Protein()
        id_list = P.Entrez_Protein_ID_Fetch(os.path.join('Input', 'ProteinInput'))
        P.Genomic_Context(id_list)

    elif response == '6':
        file_merger()

    elif response != '1' or '2' or '3' or '4' or '5' or '6':

        print 'Sorry your response was not recognized! '
        sys.exit()








