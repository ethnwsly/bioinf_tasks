import os, subprocess, re, sys
from Bio import SeqIO
from time import sleep
from itertools import groupby


#The purpose of this file is to write/manage the output stream to output files for the proteins of interest

def writeOutput(path, input_sequence):

    if path == 'gene':
        genomeFilePath = os.path.join('Output', 'Genome')
        with open(genomeFilePath, 'a') as genomeFileWrite:
            for i in range(len(input_sequence)):
                genomeFileWrite.write(input_sequence[i].format('fasta'))

    if path == 'CDS':
        cdsFilePath = os.path.join('Output', 'CDS')
        with open(cdsFilePath, 'a') as cdsFileWrite:
            for i in range(len(input_sequence)):
                cdsFileWrite.write(input_sequence[i].format('fasta'))

def clearPreviousOutput(path):

    cdsFilePath = os.path.join(path, 'CDS')
    genomeFilePath = os.path.join(path, 'Genome')

    with open(cdsFilePath, 'w') as cdsFileClear:
        cdsFileClear.write('')
        print 'CDS output file has been cleared, ready to write new records...'

    with open(genomeFilePath, 'w') as genomeFileClear:
        genomeFileClear.write('')
        print 'Genome output file has been cleared, ready to write new records...'

def repair_script():
    protein_fetch = SeqIO.parse(os.path.join('Input', 'ProteinInput'), 'fasta')
    protein_list = []

    genome_fetch = SeqIO.parse(os.path.join('Output', 'Genome'), 'fasta')
    genome_list = []

    cds_fetch = SeqIO.parse(os.path.join('Output', 'CDS'), 'fasta')
    cds_list = []

    for record in protein_fetch:
        genus_name = str(record.description.split('[')[1].split(']')[0].split(' ')[0])
        species_name = str(record.description.split('[')[1].split(']')[0].split(' ')[1])

        protein_list.append(genus_name + ' ' + species_name)

    for record in cds_fetch:

        if record.description.split(' ')[1] == 'PREDICTED:':
            predicted_split = record.description.split(' ')
            predicted_genus_species_list = record.description.split(' ')[2:4]
            cds_list.append(str(predicted_genus_species_list[0] + ' ' + predicted_genus_species_list[1]))

        else:
            genus_species_list = record.description.split(' ')[1:3]
            cds_list.append(str(genus_species_list[0] + ' ' + genus_species_list[1]))

    for record in genome_fetch:
        genus_species_list = record.description.split(' ')[1:3]

        genome_list.append(str(genus_species_list[0] + ' ' + genus_species_list[1]))

    # print 'Protein List: ' + '%s' % len(protein_list)
    # print 'Genome List: ' + '%s' % len(genome_list)
    # print 'CDS List: ' + '%s' %  len(cds_list)
    # print '\n'


    # Checks as to whether the record in the protein list is in the genome list. Prints if not.
    for record in protein_list:
        if record not in genome_list:
            print record + ' is not located in the CDS output file. Please remove from input protein script and rerun script'
        if record not in cds_list:
            print record + ' CDS record is not found in the genome file. Please remove from input protein script and rerun script'

    for i in range(len(protein_list)):
        try:
            print 'List Number: ' '%s' % i
            print 'Protein List: ' '%s' % protein_list[i]
            print 'CDS List: ' '%s' % cds_list[i]
            print 'Genome List: ' '%s' % genome_list[i]
            print '\n'

        except IndexError:
            pass

def nameOfRun():

    file_output_name = str(raw_input('Enter the desired name of the run: '))

    return file_output_name

def timer():

    run_timer = float(raw_input('Enter timer value (in seconds): '))

    return run_timer

def msa_FileCorrection():

    #the purpose of this function is to convert the MSA protein file into one that can include the description


    path = os.path.join('execs', 'tmp', 'msaprotein_aligned')
    msa_input_file = SeqIO.parse(path, 'fasta')

    genus_species_pattern = re.compile('\[\w+ \w+\]')
    accession_pattern = re.compile('[XNP_]+\d*.\d')


    sequence_list = []

    for item in msa_input_file:

        accession_number = item.description.replace(' ', '_').replace(':', '_').replace('[', '_').replace(']', '_').replace('LOW_QUALITY_PROTEIN', '_')

        # accession_search = re.search(accession_pattern, str(item))
        # accession_number = accession_search.group()
        #
        # genus_species_search = re.search(genus_species_pattern, str(item))
        # genus_species = genus_species_search.group()
        #
        #new_msa_input = str('>' + accession_number + '_' + genus_species + '\n' + item.seq + '\n')
        new_msa_input = str('>' + accession_number + '\n' + item.seq + '\n')
        sequence_list.append(new_msa_input)
        # print new_msa_input
        # sequence_list.append(new_msa_input)

    with open(path, 'w') as clear_file:
            clear_file.write('')
            clear_file.close()
    for item in sequence_list:
        print item
        with open(path, 'a') as msa_file_rewrite:
            msa_file_rewrite.write(item)
            msa_file_rewrite.close()

def dendroPy_File_Correction():

    #THIS MODULE HAS BEEN DEPRECATED - CODE IS NO LONGER IN USE FOR MOST RECENT DEPLOYMENT

    path = os.path.join('execs', 'tmp', 'msaprotein_aligned')
    msa_input_file = SeqIO.parse(path, 'fasta')

    genus_species_pattern = re.compile('\[\w+ \w+\]')
    accession_pattern = re.compile('[XNP_]+\d*.\d')

    genus_species_list = []
    accession_list = []

    for item in msa_input_file:

        accession_search = re.search(accession_pattern, str(item))
        accession_number = accession_search.group()
        accession_list.append(accession_number)

        genus_species_search = re.search(genus_species_pattern, str(item))
        genus_species = genus_species_search.group()
        genus_species_list.append(genus_species)

    with open(os.path.join('execs', 'tmp', 'unrooted_tree.nwk'), 'r') as nwk_in_file:

        read_nwk = nwk_in_file.read()
        nwk_in_file.close()

    with open(os.path.join('execs', 'tmp', 'unrooted_tree.nwk'), 'w') as nwk_out_file:
        nwk_out_file.write('')

        for i in range(len(accession_list)):
            updated_nwk = read_nwk.replace(accession_list[i], str(accession_list[i] + '_' + genus_species_list[i]))
            nwk_out_file.write(updated_nwk) #

def duplicate_management():

    #All the paths used to read in and write the duplicate stripped files to
    protein_path = 'Input/ProteinInput'
    filtered_protein_path = 'Input/Filtered_Protein'
    genome_path = 'Output/Genome'
    filtered_genome_path = 'Output/Filtered_Genome'
    CDS_path = 'Output/CDS'
    filtered_CDS_path = 'Output/Filtered_CDS'

    # protein duplication stripping and writing to updated file
    # protein_accession = [sequence.name for sequence in SeqIO.parse(protein_path, 'fasta')]
    # filtered_protein_accession = list(set(protein_accession))
    #
    # filtered_protein_list = []
    # for i in range(len(filtered_protein_accession)):
    #     for record in SeqIO.parse(protein_path, 'fasta'):
    #         if str(filtered_protein_accession[i]) == record.name:
    #             filtered_protein_list.append(record)
    #
    #
    # SeqIO.write(filtered_protein_list, filtered_protein_path, 'fasta')


    ishead = lambda x: x.startswith('>')
    all_seqs = set()
    with open(protein_path) as handle:
        with open(filtered_protein_path, 'w') as outhandle:
            head = None
            for h, lines in groupby(handle, ishead):
                if h:
                    head = lines.next()
                else:
                    seq = ''.join(lines)
                    if seq not in all_seqs:
                        all_seqs.add(seq)
                        outhandle.write('%s%s\n' % (head, seq))

    # ishead = lambda x: x.startswith('>')
    # all_seqs = set()
    # with open(genome_path) as handle:
    #     with open(filtered_genome_path, 'w') as outhandle:
    #         head = None
    #         for h, lines in groupby(handle, ishead):
    #             if h:
    #                 head = lines.next()
    #             else:
    #                 eq = ''.join(lines)
    #                 if seq not in all_seqs:
    #                     all_seqs.add(seq)
    #                     outhandle.write('%s%s\n' % (head, seq))
    #
    # ishead = lambda x: x.startswith('>')
    # all_seqs = set()
    # with open(CDS_path) as handle:
    #     with open(filtered_CDS_path, 'w') as outhandle:
    #         head = None
    #         for h, lines in groupby(handle, ishead):
    #             if h:
    #                 head = lines.next()
    #             else:
    #                 eq = ''.join(lines)
    #                 if seq not in all_seqs:
    #                     all_seqs.add(seq)
    #                     outhandle.write('%s%s\n' % (head, seq))

    # genome_accession = [sequence.name for sequence in SeqIO.parse(genome_path, 'fasta')]
    # filtered_genome_accession = list(set(genome_accession))
    # print len(filtered_genome_accession)
    # for i in range(len(filtered_genome_accession)):
    #     for record in SeqIO.parse(genome_path, 'fasta'):
    #         if str(filtered_genome_accession[i]) == record.name:
    #             pass

def welcomeStatement():

    print 'Welcome to PhyloEpsilon! A Protein Ortholog Finding Tool.'
    print 'Created and maintained by Bryan Dighera with acknowledgements to the California State University, Fullerton Nikolaidis Lab!'
    print 'Questions or bugs can be directed to the developer at: bdighera@csu.fullerton.edu\n\n'


    print 'Please Choose Your Option: '
    print 'Full Run(1): Collects Genome and CDS, makes phylogenetic tree with intron and domain figure'
    print 'Half Run(2): Makes phylogenetic tree with intron and domain figure (Protein, CDS, and Genomic sequences must already collected and added to proper folders'
    print 'Repair Script(3): Displays the order of the list'
    print 'Duplicate Script(4): Attempts to Remove Duplicates (BETA)'
    print 'Experimental Code(5): [Currently Genomic Context Script]'
    print 'File Merger(6): Merges all of the ortholog files together to make one giant tree consisting of all of them\n\n'

    response = raw_input('Please Enter a Number Listed Above: ')

    return response

def file_merger():

    path_extension = raw_input('What is the name of the protein you wish to merge?: ')

    birds_path_extension = 'HSP%s_Birds' % (path_extension)
    mammals_path_extension = 'HSP%s_Mammals' % (path_extension)
    reptiles_path_extension = 'HSP%s_Reptiles' % (path_extension)
    amphibians_path_extension = 'HSP%s_Amphibians' % (path_extension)

    Birds_CDS = [item for item in SeqIO.parse(os.path.join('CompletedTrees', 'Orthologs', path_extension, birds_path_extension, 'CDS'), 'fasta')]
    Birds_Genome = [item for item in SeqIO.parse(os.path.join('CompletedTrees', 'Orthologs', path_extension, birds_path_extension, 'Genome'), 'fasta')]
    Birds_Protein = [item for item in SeqIO.parse(os.path.join('CompletedTrees', 'Orthologs', path_extension, birds_path_extension,'ProteinInput'), 'fasta')]

    Mammals_CDS = [item for item in SeqIO.parse(os.path.join('CompletedTrees', 'Orthologs', path_extension, mammals_path_extension,  'CDS'),'fasta')]
    Mammals_Genome = [item for item in SeqIO.parse(os.path.join('CompletedTrees', 'Orthologs', path_extension, mammals_path_extension, 'Genome'), 'fasta')]
    Mammals_Protein = [item for item in SeqIO.parse(os.path.join('CompletedTrees', 'Orthologs', path_extension, mammals_path_extension, 'ProteinInput'), 'fasta')]

    Reptiles_CDS = [item for item in SeqIO.parse(os.path.join('CompletedTrees', 'Orthologs', path_extension, reptiles_path_extension, 'CDS'), 'fasta')]
    Reptiles_Genome = [item for item in SeqIO.parse(os.path.join('CompletedTrees', 'Orthologs', path_extension, reptiles_path_extension, 'Genome'), 'fasta')]
    Reptiles_Protein = [item for item in SeqIO.parse(os.path.join('CompletedTrees', 'Orthologs', path_extension, reptiles_path_extension, 'ProteinInput'), 'fasta')]

    Amphibians_CDS  = [item for item in SeqIO.parse(os.path.join('CompletedTrees', 'Orthologs', path_extension, amphibians_path_extension, 'CDS'), 'fasta')]
    Amphibians_Genome  = [item for item in SeqIO.parse(os.path.join('CompletedTrees', 'Orthologs', path_extension, amphibians_path_extension, 'Genome'), 'fasta')]
    Amphibians_Protein = [item for item in SeqIO.parse(os.path.join('CompletedTrees', 'Orthologs', path_extension, amphibians_path_extension, 'ProteinInput'), 'fasta')]

    Complete_CDS = Birds_CDS + Mammals_CDS + Reptiles_CDS + Amphibians_CDS
    Complete_Genome = Birds_Genome + Mammals_Genome + Reptiles_Genome + Amphibians_Genome
    Complete_Protein = Birds_Protein + Mammals_Protein + Reptiles_Protein + Amphibians_Protein

    SeqIO.write(Complete_CDS, os.path.join('CompletedTrees', 'Orthologs', 'Complete%s_CDS' % (path_extension)), 'fasta')
    SeqIO.write(Complete_Genome, os.path.join('CompletedTrees', 'Orthologs', 'Complete%s_Genome' % (path_extension)), 'fasta')
    SeqIO.write(Complete_Protein, os.path.join('CompletedTrees', 'Orthologs', 'Complete%s_Protein' % (path_extension)), 'fasta')

    print 'Files have been successfully merged!'


    protein_output_file = 'cleaned_protein'
    protein_holder = []
    with open(os.path.join('CompletedTrees', 'Orthologs', 'Complete%s_Protein' % (path_extension)), 'r') as file:
        rec = file.read().split('>')[1:]
        rec = ['>' + i.strip() + '\n' for i in rec]
        protein_holder.extend(rec)
    protein_total = '\n'.join(list(set(protein_holder)))
    with open(protein_output_file, 'w') as out:
        out.write(protein_total)

    genome_output_file = 'cleaned_genome'
    genome_holder = []
    with open(os.path.join('CompletedTrees', 'Orthologs', 'Complete%s_Genome' % (path_extension)), 'r') as file:
        rec = file.read().split('>')[1:]
        rec = ['>' + i.strip() + '\n' for i in rec]
        genome_holder.extend(rec)
    genome_total = '\n'.join(list(set(genome_holder)))
    with open(genome_output_file, 'w') as out:
        out.write(genome_total)

    CDS_output_file = 'cleaned_CDS'
    CDS_holder = []
    with open(os.path.join('CompletedTrees', 'Orthologs', 'Complete%s_CDS' % (path_extension)), 'r') as file:
        rec = file.read().split('>')[1:]
        rec = ['>' + i.strip() + '\n' for i in rec]
        CDS_holder.extend(rec)
    CDS_total = '\n'.join(list(set(CDS_holder)))
    with open(CDS_output_file, 'w') as out:
        out.write(CDS_total)
