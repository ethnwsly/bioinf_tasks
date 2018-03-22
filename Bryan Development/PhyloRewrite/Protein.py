from Bio import Entrez, SeqIO
import os, subprocess, random
import re
from time import sleep
from Bio.Align.Applications import ClustalOmegaCommandline
import dendropy
from ete2 import Tree, TreeStyle

#The purpose of this file is to build a class for all of the protein of interests data
#This data includes holding its corresponding Protein,CDS,and Genomic data in addition to domains and introns

class Protein(object):


    def __init__(self):

        self.input_protein_accession_number = []
        self.input_protein_sequence = []
        self.retrieved_protein_ids = []
        self.retrieved_full_cds = []
        self.parsed_cds_acession = []
        self.gene_accession_list = []
        self.gene_start_sequence_list = []
        self.gene_end_sequence_list = []
        self.gene = []
        self.intron_phase = []
        self.exon_lengths = []

    def __repr__(self):
        return self.input_protein_accession_number, self.input_protein_sequence, self.retrieved_protein_ids

    #Takes the accession number from the input file and feeds it into Entrez biopython function
    #ID is then parsed from the Entrez output and the ID is piped into following CDS and Genome functions
    #ID data is saved to the Protein object
    def Entrez_Protein_ID_Fetch(self,path):

        print 'Acquiring Protein IDs...'

        #This module will take input of the path and take the protein accession numbers from the input file and convert them to their corresponding ID number
        #The protein ID number is necessary for retrieving its corresponding CDS and genomic data

        seq_records = SeqIO.parse(path, 'fasta')

        for record in seq_records:
            self.input_protein_accession_number.append(record.id)
            self.input_protein_sequence.append(record.seq)

        for item in self.input_protein_accession_number:
            e_search = Entrez.esearch(db='protein', term="%s" % item, retmax=1000)
            e_search_results = Entrez.read(e_search, validate=False)

            for id in e_search_results['IdList']:
                self.retrieved_protein_ids.append(id)

    #Takes the protein ID from Protein_ID_Fetch and determines the corresponding CDS
    #CDS is stripped from Entrez output using python regular expression
    #CDS data is saved to Protein object
    def Entrez_CDS_Fetch(self):

        print "Acquiring CDS corresponding to protein IDs..."

        crude_cds_list = []
        cds_accession_list = []
        cds_start_nucleotide_list = []
        cds_end_nucleotide_list = []
        cds_record_dict = {'cds_accession': cds_accession_list, 'cds_start_nucleotide': cds_start_nucleotide_list,
                           'cds_end_nucleotide': cds_end_nucleotide_list}

        for id in self.retrieved_protein_ids:

            cds_efetch_search = SeqIO.parse(Entrez.efetch(db='protein', id='%s' % (id), rettype='gb', retmode='text'),
                                            'genbank')
            for seq_record in cds_efetch_search:
                if seq_record.features:
                    for feature in seq_record.features:
                        if feature.type == "CDS":
                            # print feature.location
                            crude_cds_list.append(feature.qualifiers['coded_by'])


        cds_accession_pattern = re.compile('[NXM_]+\d*.\d|[BC]+\d*.\d')
        cds_start_nucleotide_pattern = re.compile('(?<=:)\d+|(?<=:)\<\d+')
        cds_end_nucleotide_pattern = re.compile('(?<=\..)\d+')

        for cds_accession_number in crude_cds_list:

            self.parsed_cds_acession.append(cds_accession_number)

            accession_expression = re.search(cds_accession_pattern, '%s' % (cds_accession_number))
            accession_expression_out = accession_expression.group()
            cds_accession_list.append(accession_expression_out)

            start_nucleotide_expression = re.search(cds_start_nucleotide_pattern, '%s' % (cds_accession_number))
            start_nucleotide_expression_out = start_nucleotide_expression.group()
            cds_start_nucleotide_list.append(start_nucleotide_expression_out)

            end_nucleotide_expression = re.search(cds_end_nucleotide_pattern, '%s' % (cds_accession_number))
            end_nucleotide_expression_out = end_nucleotide_expression.group()
            cds_end_nucleotide_list.append(end_nucleotide_expression_out)


        for i in range(0, len(cds_accession_list)):
            cds_hits = Entrez.efetch(db='nuccore', id=cds_record_dict['cds_accession'][i],
                                     seq_start=cds_record_dict['cds_start_nucleotide'][i],
                                     seq_stop=cds_record_dict['cds_end_nucleotide'][i], rettype='fasta', retmax=1000)

            for cds_record in SeqIO.parse(cds_hits, 'fasta'):
                self.retrieved_full_cds.append(cds_record.format('fasta'))

    #Takes the protein ID from Protein_ID_Fetch and determines corresponding Genomic data
    #Todo: Error prone, fail safes have not been properly tested
    def Entrez_Genome_Fetch(self):

        print 'Acquiring Genome corresponding to protein IDs...'

        i = 0

        for id in self.retrieved_protein_ids:

            try:
                genome_elink_search = Entrez.elink(db='gene', dbfrom='protein', id=id)
                genome_elink_results = Entrez.read(genome_elink_search)

                for gene_record in genome_elink_results:
                    for gene_link in gene_record['LinkSetDb']:
                        for identification in gene_link['Link']:
                            genome_id = identification['Id']
                            #print "Genome IDs: " + '%s' % genome_id
                            genome_fetch = Entrez.efetch(db='gene', id=genome_id, retmax=1000, rettype='fasta',
                                                         retmode='xml')
                            genome_fetch_parse = Entrez.parse(genome_fetch)

                            for gene_record in genome_fetch_parse:
                                # print gene_record.keys()

                                gene_accession = gene_record['Entrezgene_locus'][0]['Gene-commentary_accession'] + "." + \
                                                 gene_record['Entrezgene_locus'][0]['Gene-commentary_version']
                                gene_sequence_start_region = \
                                gene_record['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int'][
                                    'Seq-interval']['Seq-interval_from']
                                gene_sequence_end_reqion = \
                                gene_record['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int'][
                                    'Seq-interval']['Seq-interval_to']

                                #print "Sequence Numer: " + '%s' % i + '/' + '%s' % len(self.retrieved_protein_ids)
                                #print gene_record['Entrezgene_locus'][0]['Gene-commentary_accession'] + "." + \
                                      #gene_record['Entrezgene_locus'][0]['Gene-commentary_version']

                                self.gene_accession_list.append(gene_accession)
                                self.gene_start_sequence_list.append(gene_sequence_start_region)
                                self.gene_end_sequence_list.append(gene_sequence_end_reqion)

                                sleep(1)

                                if self.gene_accession_list[i] == '':
                                    print 'ERROR: Gene Accession is Empty!'
                                    break
                                if self.gene_start_sequence_list[i] == '':
                                    print "ERROR: Gene's Beginning Sequence is Empty!"
                                    break
                                if self.gene_end_sequence_list[i] == '':
                                    print "ERROR: Gene's End Sequence is Empty!"
                                    break

                                genome_handle = Entrez.efetch(db='nuccore', id=gene_accession,
                                                              seq_start=int(gene_sequence_start_region),
                                                              seq_stop=int(gene_sequence_end_reqion), rettype='fasta')

                                self.gene.append(genome_handle.read())

                                i += 1


            except KeyError:
                print "Error in protein ID: " '%s' % id
                print "Entire record being removed..."
                self.retrieved_protein_ids.remove(id)
                del self.input_protein_sequence[i]
                pass

    #Makes a formatted temporary CDS/Genomic file that is directly fed into spidey.linux.64
    #tmp file is cleared for each protein, output (in the form of intron phase and position) is funneled into Protein object
    #WARNING: if recieving error from spidey.linux.64 - ensure that you are running a 64 bit processor.You may need to update with 32 bit version
    #WARNING: if recieving permission denied error - make sure that properties allowing execution of the file are turned on in advanced settings
    def intronCalculator(self):

        def intronOutputParser(intronProcess):

            throwaway_lines = intronProcess[3:]
            number_of_introns = throwaway_lines[1].split()[3]
            introns = throwaway_lines[2:2 + int(number_of_introns)]

            if int(number_of_introns) > 1:

                intron_phase_list = []
                exon_lengths_list = []

                for exon_bit in introns[:-1]:
                    intron_location = exon_bit.split()[4].split()[0].split('-')
                    intron_start = int(intron_location[0])
                    intron_end = int(intron_location[1])
                    intron_phase = (intron_end - intron_start) % 3

                    intron_phase_list.append(intron_phase)
                    exon_lengths_list.append(exon_bit.split()[4].split())

                self.intron_phase.append(intron_phase_list)
                self.exon_lengths.append(exon_lengths_list)

            else:
                print 'No Introns! Moving on...'

        input_cds_path = os.path.join('Output', 'CDS')
        cds_seq = [sequence.seq for sequence in SeqIO.parse(input_cds_path, 'fasta')]

        input_genome_path = os.path.join('Output', 'Genome')
        genome_seq = [sequence for sequence in SeqIO.parse(input_genome_path, 'fasta')]

        for i in range(len(cds_seq)):
            spidey_cds_path = os.path.join('execs', 'tmp', 'spidey_cdsinput.txt')
            with open(spidey_cds_path, 'w') as cdsFileWrite:
                cdsFileWrite.write(str('>CDS' + '\n' + cds_seq[i]))
                cdsFileWrite.close()

            spidey_genome_path = os.path.join('execs', 'tmp', 'spidey_genomeinput.txt')
            with open(spidey_genome_path, 'w') as genomeFileWrite:
                genomeFileWrite.write(str('>' + genome_seq[i].description + '\n' + genome_seq[i].seq))
                genomeFileWrite.close()

            spidey_executible_path = os.path.join('execs', 'spidey.linux.64')
            print subprocess.list2cmdline([spidey_executible_path,
                                           "-i", spidey_genome_path,
                                           "-m", spidey_cds_path,
                                           "-p", "1"])

            intron_process = subprocess.Popen([spidey_executible_path,
                                               "-i", spidey_genome_path,
                                               "-m", spidey_cds_path,
                                               "-p", "1"], stdout=subprocess.PIPE)


            intronOutputParser(intron_process.stdout.readlines())

    #This function will place the input proteins into a tmp file that will be piped into clustal commandline tool
    #Output and input will be in the tmp directory - output will be a file containing aligned sequences
    def multiple_sequencing_alignment(self):

        in_file = os.path.join('execs','tmp','msaprotein_infile')
        out_file = os.path.join('execs','tmp', 'msaprotein_aligned')

        with open(in_file, 'w') as clearFile:
            clearFile.write('')
            print 'Multiple Sequencing Alignment File input file has been cleared. New run alignment can now be loaded...'

        protein_accession = [record.description for record in SeqIO.parse((os.path.join('Input', 'ProteinInput')),'fasta')]
        protein_sequence = [record.seq for record in SeqIO.parse((os.path.join('Input', 'ProteinInput')),'fasta')]

        with open(in_file, 'a') as msaprotein_writeFile:
            for i in range(len(protein_accession)):
                msaprotein_writeFile.write('\n'+ '>' + str(protein_accession[i]) + '\n' + str(protein_sequence[i]))


        print '\n' + 'CREATING MULTIPLE SEQUENCING ALIGNMENT...'

        clustalomega_cline = ClustalOmegaCommandline(cmd=os.path.join('execs', "clustalo-1.2.0"),
                                                     infile=in_file,
                                                     outfile=out_file, verbose=True, auto=True, force=True)

        clustalomega_cline()

        print '\n' + 'MULTIPLE SEQUENCING ALIGNMENT HAS BEEN CREATED...'

    #Takes the tmp multiple sequencing alignment file and pipes into Fast Tree standalone executable.
    #Output is a rooted tree that will be used by ete to form the visual phylogenetic tree
    def rootedTreeConstruction(self):

        in_file = os.path.join('execs','tmp', 'msaprotein_aligned')
        out_file = os.path.join('execs','tmp', 'unrooted_tree.nwk')

        subprocess.call(["./execs/FastTree", "-out", out_file, in_file])
        print('\n' + subprocess.list2cmdline(["./execs/FastTree", in_file, ">", out_file]))

        rooted_tree = dendropy.Tree.get_from_path(out_file, schema='newick')
        rooted_tree.reroot_at_midpoint()
        rooted_tree.write_to_path(os.path.join('execs','tmp', "rooted_tree.nwk"), schema='newick')

        print '\n' + 'ROOTED TREE HAS BEEN CONSTRUCTED...'

    #Function is full display functionality for the entire application
    #Todo: Basic tree implementation has been achieved, full graphics to come - intron phase, domain mapping
    def renderingTreeImage(self):

        with open(os.path.join('execs','tmp', "rooted_tree.nwk"), 'r') as nwk_tree_handle:
            nwk_tree = nwk_tree_handle.read()
            t = Tree(nwk_tree)
            print(t)

        ts = TreeStyle()
        ts.allow_face_overlap = True
        ts.show_leaf_name = True
        ts.show_branch_support = True
        domain_colors = []
        golden_ratio = 0.618033988749895
        h = random.random()
        h += golden_ratio
        h %= 1



        t.show(tree_style=ts)




















