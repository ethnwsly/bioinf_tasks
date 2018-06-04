from Bio import Entrez, SeqIO
import os, subprocess, sys, math, socket
import re
from time import sleep
from Bio.Align.Applications import ClustalOmegaCommandline
import dendropy
from randomcolor import RandomColor
from FileManager import nameOfRun, msa_FileCorrection, timer


from ete3 import Tree, SeqMotifFace, TreeStyle, TextFace, random_color
#The purpose of this file is to build a class for all of the protein of interests data
#This data includes holding its corresponding Protein,CDS,and Genomic data in addition to domains and introns

class Protein(object):



    def __init__(self):

        self.timer = 0
        self.run_name = ''
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
        self.accession_dict_with_introns = {}
        self.number_protein_IDs = 0
        self.number_gene_IDs = 0
        self.number_CDS_IDs = 0
        self.msa_aligned_protein = []

    def __repr__(self):
        return self.input_protein_accession_number, self.input_protein_sequence, self.retrieved_protein_ids

    #Takes the accession number from the input file and feeds it into Entrez biopython function
    #ID is then parsed from the Entrez output and the ID is piped into following CDS and Genome functions
    #ID data is saved to the Protein object
    def Entrez_Protein_ID_Fetch(self, path):

        print 'Acquiring Protein IDs...'

        #This module will take input of the path and take the protein accession numbers from the input file and convert them to their corresponding ID number
        #The protein ID number is necessary for retrieving its corresponding CDS and genomic data

        seq_records = SeqIO.parse(path, 'fasta')

        for record in seq_records:
            self.input_protein_accession_number.append(record.id)
            self.input_protein_sequence.append(record.seq)

        for item in self.input_protein_accession_number:
            sleep(self.timer)
            e_search = Entrez.esearch(db='protein', term="%s" % item, retmax=1000)
            e_search_results = Entrez.read(e_search, validate=False)

            for id in e_search_results['IdList']:
                print id
                self.retrieved_protein_ids.append(id)

        self.number_protein_IDs = len(self.retrieved_protein_ids)
        print 'Number of Protein IDs found: ' + str(self.number_protein_IDs)

        print 'Protein IDs: ' + str(self.retrieved_protein_ids)


        return self.retrieved_protein_ids

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

            sleep(self.timer)
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
            sleep(self.timer)
            cds_hits = Entrez.efetch(db='nuccore', id=cds_record_dict['cds_accession'][i],
                                     seq_start=cds_record_dict['cds_start_nucleotide'][i],
                                     seq_stop=cds_record_dict['cds_end_nucleotide'][i], rettype='fasta', retmax=1000)

            for cds_record in SeqIO.parse(cds_hits, 'fasta'):
                print str(cds_record.description) + '\n' + '%s' % cds_record.seq
                self.retrieved_full_cds.append(cds_record.format('fasta'))

        print 'Number of CDS records found: ' + str(len(self.retrieved_full_cds))

    #Takes the protein ID from Protein_ID_Fetch and determines corresponding Genomic data
    #Todo: Need to error handle if there is not a genome for a corresponding CDS
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
                            sleep(self.timer)
                            genome_fetch = Entrez.efetch(db='gene', id=genome_id, retmax=1000, rettype='fasta',
                                                         retmode='xml')
                            genome_fetch_parse = Entrez.parse(genome_fetch)

                            for gene_record in genome_fetch_parse:
                                #print gene_record.keys()

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

                                sleep(self.timer)

                                if self.gene_accession_list[i] == '':
                                    print 'ERROR: Gene Accession is Empty!'
                                    break
                                if self.gene_start_sequence_list[i] == '':
                                    print "ERROR: Gene's Beginning Sequence is Empty!"
                                    break
                                if self.gene_end_sequence_list[i] == '':
                                    print "ERROR: Gene's End Sequence is Empty!"
                                    break

                                sleep(self.timer)
                                genome_handle = Entrez.efetch(db='nuccore', id=gene_accession,
                                                              seq_start=int(gene_sequence_start_region),
                                                              seq_stop=int(gene_sequence_end_reqion), rettype='fasta')


                                self.gene.append(genome_handle.read())

                                print 'Genomic Sequence Run #: ' + '%s' % i
                                i += 1




            except KeyError:
                print "Error in protein ID: " '%s' % id
                print "Entire record being removed..."
                self.retrieved_protein_ids.remove(id)
                del self.input_protein_sequence[i]
                pass

        print 'Number of Gene Records found: ' + str(self.number_gene_IDs)
        self.number_gene_IDs = len(self.gene)

    #Makes a formatted temporary CDS/Genomic file that is directly fed into spidey.linux.64
    #tmp file is cleared for each protein, output (in the form of intron phase and position) is funneled into Protein object
    #WARNING: if recieving error from spidey.linux.64 - ensure that you are running a 64 bit processor.You may need to update with 32 bit version
    #WARNING: if recieving permission denied error - make sure that properties allowing execution of the file are turned on in advanced settings
    def intronCalculator(self):

        #Parses the output from spidey.linux to determine the intron phase and the location of the introns with respect to the exons
        def intronOutputParser(intronProcess):
            throwaway_lines = intronProcess[3:]
            number_of_introns = throwaway_lines[1].split()[3]
            introns = throwaway_lines[2:2 + int(number_of_introns)]
            intron_phase_list = []
            exon_lengths_list = []

            if int(number_of_introns) > 1:

                for exon_bit in introns[:-1]:
                    intron_location = exon_bit.split()[4].split()[0].split('-')
                    intron_start = int(intron_location[0])
                    intron_end = int(intron_location[1])
                    intron_phase = (intron_end - intron_start) % 3

                    intron_phase_list.append(intron_phase)
                    #print 'Intron Phase(s): ' + str(intron_phase)
                    exon_lengths_list.append(exon_bit.split()[4].split())
                    #print 'Exon Length(s): ' + str((exon_bit.split()[4].split()))


            else:
                print 'No Introns! Only One Exon! Moving on...'
                intron_phase_list.append('NONE')
                exon_lengths_list.append('NONE')

            self.intron_phase.append(intron_phase_list)
            self.exon_lengths.append(exon_lengths_list)

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
            subprocess.list2cmdline([spidey_executible_path,
                                           "-i", spidey_genome_path,
                                           "-m", spidey_cds_path,
                                           "-p", "1"])

            intron_process = subprocess.Popen([spidey_executible_path,
                                               "-i", spidey_genome_path,
                                               "-m", spidey_cds_path,
                                               "-p", "1"], stdout=subprocess.PIPE)

            #pipes the std output from the linux executible into the parser
            intronOutputParser(intron_process.stdout.readlines())

    #This function will place the input proteins into a tmp file that will be piped into clustal commandline tool
    #Output and input will be in the tmp directory - output will be a file containing aligned sequences
    def multiple_sequencing_alignment(self):

        in_file = os.path.join('execs','tmp','msaprotein_infile')
        out_file = os.path.join('execs','tmp', 'msaprotein_aligned')

        if self.run_name == '':
            nameOfRun()

        with open(in_file, 'w') as clearFile:
            clearFile.write('')
            print '\nMultiple Sequencing Alignment File input file has been cleared. New run alignment can now be loaded...'

        protein_accession = [record.description for record in SeqIO.parse((os.path.join('Input', 'ProteinInput')),'fasta')]
        protein_sequence = [record.seq for record in SeqIO.parse((os.path.join('Input', 'ProteinInput')),'fasta')]

        with open(in_file, 'a') as msaprotein_writeFile:
            for i in range(len(protein_accession)):
                protein = '\n'+ '>' + str(protein_accession[i]) + '\n' + str(protein_sequence[i])
                msaprotein_writeFile.write(protein)


        print '\n' + 'CREATING MULTIPLE SEQUENCING ALIGNMENT...'

        clustalomega_cline = ClustalOmegaCommandline(cmd=os.path.join('execs', "clustalo-1.2.0"),
                                                     infile=in_file,
                                                     outfile=out_file, verbose=True, auto=True, force=True)

        clustalomega_cline()
        msa_FileCorrection()

        self.msa_aligned_protein = [item.seq for item in SeqIO.parse(out_file, 'fasta')]

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

        newick_rooted_file = open(os.path.join('execs','tmp', "rooted_tree.nwk"),'r')
        read_edit_newick = newick_rooted_file.read()

        #The tree generated sometimes has the [&R] region - if not stripped it will throw error. Try except handles if the [&R] is not generated
        try:

            stripped_tree = read_edit_newick.strip('\[&R\] ')
            with open(os.path.join('execs', 'tmp', "rooted_tree.nwk"), 'w') as writeStrippedTree:
                writeStrippedTree.write('')
                writeStrippedTree.write(stripped_tree)

                with open(os.path.join('execs', 'tmp', "rooted_tree.nwk"), 'w') as writeStrippedTree:
                    writeStrippedTree.write('')
                    writeStrippedTree.write(stripped_tree)

        except AttributeError:
            pass


        print '\n' + 'ROOTED TREE HAS BEEN CONSTRUCTED...'

    #Function is full display functionality for the entire application
    #Todo: Basic tree implementation has been achieved with introns/intron phases, full graphics to come - domain mapping
    def renderingTreeImage(self):

        path = os.path.join('Input', 'ProteinInput')

        seq_records = SeqIO.parse(path, 'fasta')

        for record in seq_records:
            self.input_protein_accession_number.append(record.id)
            self.input_protein_sequence.append(record.seq)

        with open(os.path.join('execs','tmp', "rooted_tree.nwk")) as nwk_tree_handle:
            nwk_tree = nwk_tree_handle.read()
            t = Tree(nwk_tree)
            print(t)
            print '\n'

        ts = TreeStyle()
        ts.title.add_face(TextFace('PhyloEpsilon - Protein Ortholog Finding Tool by Bryan Dighera', fsize= 16,), column= 0)
        ts.allow_face_overlap = True
        ts.show_leaf_name = True
        ts.show_branch_support = True

        leaf_names = []
        for leaf in t.get_leaf_names():

            np_xp_pattern = re.compile('N[P]|X[P]')
            digits_pattern = re.compile('\d+.\d')

            np_xp_search_obj = re.search(np_xp_pattern, leaf)
            digits_search_obj = re.search(digits_pattern, leaf)

            np_xp_name = np_xp_search_obj.group()
            digits_name = digits_search_obj.group()
            final_accession = str(np_xp_name + '_' + digits_name)
            print final_accession
            leaf_names.append(final_accession)


        #print 'leaf names: ' + '%s' % leaf_names

        P = Protein()
        protein_domains, domain_colors, unrepeated_domains = P.Domains()
        print domain_colors

        #Creates a dictionary that corresponds the protein accession number to its corresponding introns
        for i in range(len(leaf_names)):
            self.accession_dict_with_introns[self.input_protein_accession_number[i]] = self.exon_lengths[i]

        i = 0

        print 'protein accession number: ' + '%s' % self.input_protein_accession_number
        print 'Accession dict: ' + '%s' % self.accession_dict_with_introns + '\n'

        #Iterates through the accession numbers that correspond the the order of the leaves of the phylogenetic tree to retrieve introns and build fig
        for accession_number in leaf_names:
            intron_motifs = [[0, 0, "[]", None, 12, "White", "White", None]]

            #Checks the accession number against the dictionary and retrieves the corresponding introns, if no introns then doesn't append any
            if accession_number in self.accession_dict_with_introns:
                print accession_number, self.accession_dict_with_introns[accession_number]
                exon_list = self.accession_dict_with_introns[accession_number]
                print exon_list

                for exon_length in exon_list:
                    if str(exon_length) != 'NONE':

                        for location in exon_length:
                            split_exon_location = str(location).split('-')
                            protein_seq_exon_location = int(math.floor(int(split_exon_location[1])/3))

                            #Calculates the intron phase and then checks the phase to append appropriate color indicating phase on diagram
                            intron_phase = (int(split_exon_location[1]) - int(split_exon_location[0])) % 3

                            if intron_phase == 0:
                                intron_motifs.append([protein_seq_exon_location - 2,
                                                      protein_seq_exon_location + 2,
                                                      "[]", None, 5, "Grey", "Grey", None])
                            elif intron_phase == 1:
                                intron_motifs.append([protein_seq_exon_location - 2,
                                                      protein_seq_exon_location + 2,
                                                      "[]", None, 5, "Black", "Black", None])

                            elif intron_phase == 2:
                                intron_motifs.append([protein_seq_exon_location - 2,
                                                      protein_seq_exon_location + 2,
                                                      "[]", None, 5, "Blue", "Blue", None])
                    else:
                        print 'NO INTRONS FOUND FOR RECORD'

                print str(intron_motifs) + '\n'
                msa_protein_seq = self.msa_aligned_protein[i].strip('-')

                #ete3 module that adds the introns(motifs) to the phylogenetic tree
                seqFace = SeqMotifFace(str(msa_protein_seq), gapcolor="black", seq_format= 'line', scale_factor=1, motifs=intron_motifs)
                (t & t.get_leaf_names()[i]).add_face(seqFace, 0, "aligned")

                i += 1

        n = 0


        # Iterates through the accession numbers that correspond to the order of the leaves of the phylogenetic tree and compare to domain dict values
        # TODO: Add the legend and possibly give a number to each of the domains so they can be easily identified in the legend
        for accession_number in leaf_names:

            domain_motifs = [[0, 0, "[]", None, 12, "White", "White", None]]

            for domain in protein_domains:

                if accession_number in domain:

                    print 'leaf accession #: ' + '%s' % accession_number
                    print 'domains accession: ' + '%s' % domain.keys()[0]
                    print domain.values()[0]

                    for each_domain in domain.values()[0]:

                        try:

                            domain_motif_color = domain_colors[each_domain[0]]
                            start_domain_loc = int(each_domain[1].split(':')[0])

                            end_domain_loc = int(each_domain[1].split(':')[1])
                            domain_name = str(each_domain[0])

                            domain_motifs.append([start_domain_loc,
                                                  end_domain_loc,
                                                  "<>", 20, 20, 'Black',
                                                  domain_motif_color,
                                                  'arial|8|black|'])
                        except ValueError:

                            domain_motif_color = domain_colors[each_domain[0]]

                            start_pattern = re.compile('(?<!=\W)\d+')
                            start_pattern_search = re.search(start_pattern, str(each_domain[1].split(':')[0]))
                            start_domain_loc = int(start_pattern_search.group())

                            end_pattern = re.compile('(?<!=\W)\d+')
                            end_pattern_search = re.search(end_pattern, str(each_domain[1].split(':')[1]))
                            end_domain_loc = int(end_pattern_search.group())

                            domain_motifs.append([start_domain_loc,
                                                  end_domain_loc,
                                                  "<>", 20, 20, 'Black',
                                                  domain_motif_color,
                                                  'arial|8|black|'])


            print domain_motifs

            msa_protein_seq = self.msa_aligned_protein[n].strip('-')
            print msa_protein_seq
            print len(msa_protein_seq)
            print '*' * 100

            domainFace = SeqMotifFace(str(msa_protein_seq), gapcolor="black", seq_format='line', scale_factor=1,
                                   motifs=domain_motifs)
            (t & t.get_leaf_names()[n]).add_face(domainFace, 0, "aligned")

            n += 1


        #Creating the legend

        print protein_domains
        for single_unrepeat, colors in domain_colors.iteritems():

            ts.legend.add_face(TextFace(single_unrepeat), column=0)
            ts.legend.add_face(SeqMotifFace("A" * 45, [[0, 80, "[]", None, 8, "Black", colors , None]]), column= 1)
            ts.legend_position = 1


        #name_of_run = nameOfRun()
        file_name = self.run_name
        t.show(tree_style=ts)
        t.render(os.path.join('CompletedTrees',  file_name + '.pdf'), tree_style=ts)

    #Function collects the domains of the input proteins
    #Function returns dictionary containing parsed genbank content from NCBI
    #Output is funneled into renderingTreeImage to build motifs for the domains
    def Domains(self):

        seq_records = SeqIO.parse(os.path.join('Input', 'ProteinInput'), 'fasta')

        input_protein_list = []
        for record in seq_records:
            input_protein_list.append(record.id)

        Complete_Domains = []
        domainNameList = []
        for item in input_protein_list:
            sleep(self.timer)
            e_fetch = Entrez.efetch(db='protein', id="%s" % item, retmax=1000, rettype='gb', retmode='fasta')

            for seq_record in SeqIO.parse(e_fetch, 'gb'):
                domain_list = []
                accession_number = seq_record.id

                for i in range(len(seq_record.features)):
                    if seq_record.features[i].type == 'Region':
                        domain_location = str(seq_record.features[i].location).split('[')[1].split(']')[0]
                        domain_name = str(seq_record.features[i].qualifiers['region_name'][0])
                        domainNameList.append(domain_name)

                        domain_list.append([domain_name, domain_location])

                Complete_Domains.append(dict([(accession_number, domain_list)]))


        rand_color = RandomColor()
        domains_dict_colors ={domain:rand_color.generate()[0] for domain in set(domainNameList)}
        Domains = [domain for domain in set(domainNameList)]

        return Complete_Domains, domains_dict_colors, Domains

    def Genomic_Context(self, id_list):

        gene_accession_list = []
        genome_id_list = []
        for record in id_list:
            genome_elink_search = Entrez.elink(db='gene', dbfrom='protein', id=record)
            genome_elink_results = Entrez.parse(genome_elink_search)

            for gene_record in genome_elink_results:
                for gene_link in gene_record['LinkSetDb']:
                    for identification in gene_link['Link']:
                        genome_id = identification['Id']
                        genome_id_list.append(genome_id)


        path = os.path.join('Output', 'Genome')
        genome_record_list = SeqIO.parse(path, 'fasta')

        coords_list = []
        for record in genome_record_list:

            genome_accession = record.description.split(':')[0]
            genome_cords = record.description.split(' ')[0].split(':')[1].split('-')
            genome_start_coord = int(genome_cords[0]) - 100000
            genome_end_coord = int(genome_cords[1]) - 100000
            complete_coords = (genome_start_coord,genome_end_coord)
            coords_list.append(complete_coords)
            gene_accession_list.append(genome_accession)

        print genome_id_list
        print coords_list

        for i in range(len(gene_accession_list)):
            genome_id = gene_accession_list[i]
            genome_start_coord = coords_list[i][1]
            genome_end_coord = coords_list[i][1]
            print '\n%s:%s-%s' % (genome_id, genome_start_coord, genome_end_coord)

            genomic_context_fetch = Entrez.parse(Entrez.efetch(db="nuccore",
                                                                       id= genome_id,
                                                                       seq_start=genome_start_coord,
                                                                       seq_stop=genome_end_coord,
                                                                       rettype="gb",
                                                                       retmode="xml"), validate=False)


            for gene_record in genomic_context_fetch:
                print gene_record['GBSeq_locus']

                # break
            # for gene_record in genomic_context_fetch:
            #     print gene_record.keys()
            #     print gene_record['Entrezgene_gene']
            #     print '**************' * 10
            #     print gene_record['Entrezgene_location']
            #     print '*******' * 10
            #     print gene_record['Entrezgene_locus']
            #
            #     gene_accession_1 = gene_record['Entrezgene_locus'][0]['Gene-commentary_accession'] + "." + \
            #                      gene_record['Entrezgene_locus'][0]['Gene-commentary_version']
            #     gene_accession_2 = gene_record['Entrezgene_locus'][1]['Gene-commentary_accession']
            #     gene_sequence_start_region = \
            #         gene_record['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int'][
            #             'Seq-interval']['Seq-interval_from']
            #     gene_sequence_end_reqion = \
            #         gene_record['Entrezgene_locus'][0]['Gene-commentary_seqs'][0]['Seq-loc_int'][
            #             'Seq-interval']['Seq-interval_to']
            #
            #     print '%s: %s' % (gene_accession_1, gene_accession_2)



        #
        #     for count, gb_fetch in enumerate(genomic_context_fetch):
        #
        #         genes = []
        #
        #         for gb_gene in gb_fetch['GBSeq_feature-table']:
        #             if gb_gene['GBFeature_key'] == 'gene':
        #                 for gb_info in gb_gene['GBFeature_quals']:
        #                     if gb_info['GBQualifier_name'] == 'gene':
        #                             name = gb_info['GBQualifier_value']
        #                             print name
        # #
                            #     try:
                            #         if 'GeneID:' in gb_info['GBQualifier_value']:
                            #             gene_id = int(gb_info['GBQualifier_value'].replace('GeneID:', ''))
                            #     except:
                            #         print("No GBQualifier_value")
                            #
                            #     start = int(gb_gene['GBFeature_intervals'][0]['GBInterval_from'])
                            # end = int(gb_gene['GBFeature_intervals'][0]['GBInterval_to'])
                            # if genome_accession == gene_id:
                            #     target_shift = end
                            # if start < end:
                            #     direction = "+"
                            # else:
                            #     direction = "-"
                            # genes.append(dict(input_accession = genome_accession, megene_name=name, gene_id=gene_id, gene_start=start, gene_end=end,
                            #                   gene_direction=direction))
                            # print genes
        #
        #
        #     except socket.error, se:
        #         'Error in Collecting GC'
        #         print se




















