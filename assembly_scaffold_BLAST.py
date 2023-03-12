import sys
import os
import subprocess
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
from random import shuffle


class FastAreader:
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
    print (head,seq)
    '''

    def __init__(self, fname):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            # header = ''
            # sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
                header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


class tree_maker:
    """

    """

    def __init__(self):

        # for assembly
        self.reverse_complement_dict = {}

        self.null_model_read_dict = {}

        # for baits
        self.baits_by_head_dict = {}

        self.seq_dict_by_head = {}
        self.head_dict_by_seq = {}

        self.duplicate_reads_by_file_set = set()
        self.all_reads_by_file_set = set()

        self.duplicate_sequence_dict = {}
        self.read_bait_correspondence_dict = {}

        self.unique_sequence_dict = {}
        self.duplicate_header_dict = {}

        # for fourmer frequency values of assembled paired end reads
        self.assembled_fourmer_freq_dict = {}

        # for file fourmer freq sums
        self.fourmer_frequency_dict = {}
        self.null_model_freq_dict = {}
        self.fourmer_frequency_by_file_dict = {}
        self.fourmer_sum_dict = {}

        self.forward_correspondence_dict = {}

        self.correspondence_node_dict = {}

        self.sink_dict = {}
        self.start_list = []
        self.end_list = []

        self.assembled_paired_end_read_dict = {}

        self.path_dict = {}
        self.direct_path_dict = {}

        self.assembly_seq_dict = {}

        self.Z_score_dict_by_read = {}

        self.agreement_seq_len_sum = 0
        self.total_num_agreements = 0

    def file_converter(self):
        """Designed to open a folder on desktop containing .fq.gz data from sample reads.

         These .fq.gz file names used to develop .fa files via bash subprocesses.

         Input: files in path
         Output: converted .fa files for use
        """

        path = "C://Users//Quin The Conquoror!//Desktop//Bobs Files fq.gz/Bobs files .fa"
        dir_list = os.listdir(path)

        for infile_name in dir_list:
            print(infile_name)
            outfile_name = infile_name.split('.')[0] + '.fa'
            print(outfile_name)

            # run the seqtk bash command which converts .fq.gz files to .fa for use
            # the .fa files are placed in a folder on the desktop for later use.
            subprocess_command = "seqtk seq -a " + infile_name + " > " + outfile_name
            print(subprocess_command)
            subprocess.run([subprocess_command], cwd=path, shell=True)

    def file_parse(self):
        """In parsing the fasta files, the bait sequences must first be assessed.

        Afterwards, the read files are unpacked sequentially.
        """

        path = "C://Users//Quin The Conquoror!//Desktop//Bobs Files fq.gz/Bobs files .fa"
        dir_list = os.listdir(path)

        for read_file_name in dir_list:
            read_file_data = FastAreader(read_file_name)
            print("file name: ", read_file_name)
            seq_number = 0
            self.fourmer_frequency_by_file_dict[read_file_name] = {}

            for head, seq in read_file_data.readFasta():
                # self.fourmer_frequency(read_file_name, head, seq, seq_number)
                self.seq_dict_by_head[head] = seq
                self.head_dict_by_seq[seq] = head
                seq_number += 1
                # if seq_number == 2500:
                #    break

    def reverse_complement(self, seq):
        """Develops the reverse complement of a parameter sequence."""

        return seq.lower().replace('a', 'T').replace('t', 'A').replace('g', 'C').replace('c', 'G').replace(' ', '')[::-1]

    def paired_end_assembly(self):
        """Assembles Illumina Hiseq paired end reads.

        3 cases:
        1) R1 ends in R2
        2) R1 begins in R2
        3) No overlap
        """

        for read_header in self.seq_dict_by_head.keys():
            if "2:N:0" in read_header:
                # provides a dictionary for the R2 reads
                self.reverse_complement_dict[read_header] = self.reverse_complement(self.seq_dict_by_head[read_header])

        # searching for overlap between paired sequences
        for read_header_r1, read_header_r2 in zip(self.seq_dict_by_head.keys(), self.reverse_complement_dict.keys()):
            if read_header_r1[:-15] == read_header_r2[:-15]:
                for i in range(2, len(self.seq_dict_by_head[read_header_r1])-2, 2):

                    if self.seq_dict_by_head[read_header_r1][:-i] == self.reverse_complement_dict[read_header_r2][i:]:
                        print("case 1 agreement: ", self.seq_dict_by_head[read_header_r1][:-i])
                        print("Read headers: ", read_header_r1, read_header_r2)
                        print("Original seqs: ", self.seq_dict_by_head[read_header_r1], self.reverse_complement_dict[read_header_r2])
                        print("Seq:", self.reverse_complement_dict[read_header_r2][:i] + self.seq_dict_by_head[read_header_r1])

                        case1_seq = self.reverse_complement_dict[read_header_r2][:i] + self.seq_dict_by_head[read_header_r1]

                        print("paired seq len: ", len(case1_seq))
                        self.assembled_paired_end_read_dict[read_header_r1] = case1_seq
                        print(" ")

                    if self.seq_dict_by_head[read_header_r1][i:] == self.reverse_complement_dict[read_header_r2][:-i]:
                        print("case 2 agreement:", self.seq_dict_by_head[read_header_r1][i:])
                        print("Read headers: ", read_header_r1, read_header_r2)
                        print("Original seqs: ", self.seq_dict_by_head[read_header_r1], self.reverse_complement_dict[read_header_r2])

                        case2_seq = self.seq_dict_by_head[read_header_r1][:i] + self.reverse_complement_dict[read_header_r2]

                        print("paired seq len: ", len(case2_seq))
                        print("Assembled seq: ", case2_seq)
                        self.assembled_paired_end_read_dict[read_header_r1] = case2_seq
                        print(" ")
                    """
                    if self.seq_dict_by_head[read_header_r1][i:] != self.reverse_complement_dict[read_header_r2][:-i] and\
                            self.seq_dict_by_head[read_header_r1][:-i] != self.reverse_complement_dict[read_header_r2][i:]:

                        print("case 3")
                        print("Read headers: ", read_header_r1, read_header_r2)
                        print("Original seqs: ", self.seq_dict_by_head[read_header_r1], self.reverse_complement_dict[read_header_r2])
                        print("Assembled seq: ", self.seq_dict_by_head[read_header_r1][:i] + self.reverse_complement_dict[read_header_r2])

                        case3_seq = self.seq_dict_by_head[read_header_r1] + "_" + self.reverse_complement_dict[read_header_r2]
                        print("Case 3 Seq:", case3_seq)
                        self.assembled_paired_end_read_dict[read_header_r1] = case3_seq
                        print(" ")
                    """

        self.assembled_fourmer_freq_dict = self.fourmer_frequency_of_paired_end_reads(self.assembled_paired_end_read_dict)

    def null_model(self):
        """The null model is provided as a randomness and error threshold.
        Each assembled paired end read will be scrambled and its fourmer frequencies analyzed for
        comparison to later data extracted from the exact read sequences.

        Input: paired end sequences
        Output: scrambled paired end sequences for fourmer analysis
        """

        for read_header in self.assembled_paired_end_read_dict.keys():
            exact_seq = self.assembled_paired_end_read_dict[read_header]

            char_list = list(exact_seq)
            shuffle(char_list)

            randomized_seq = ''.join(char_list)

            self.null_model_read_dict[read_header] = randomized_seq

        self.null_model_freq_dict = self.fourmer_frequency_of_paired_end_reads(self.null_model_read_dict)

    def fourmer_frequency_of_paired_end_reads(self, analysis_dict):
        """ This member function is designed """

        fourmer_list_list = []
        for read_header in analysis_dict.keys():
            seq = analysis_dict[read_header]

            assembly_fourmers = [seq[i:i + 4] for i in range(0, len(seq), 4)]
            fourmer_list_list.append(assembly_fourmers)

        # generate all possible fourmers as key to dictionary, save scores per fourmer by fourmer as key
        #print(fourmer_list_list)
        unique_fourmer_list = []
        all_fourmer_list = []
        for fourmer_list in fourmer_list_list:
            for fourmer in fourmer_list:
                all_fourmer_list.append(fourmer)
                if fourmer not in unique_fourmer_list and fourmer:
                    unique_fourmer_list.append(fourmer)

        fourmer_freq_dict = {}
        for fourmer in unique_fourmer_list:
            fourmer_freq_dict[fourmer] = all_fourmer_list.count(fourmer)/256

        return fourmer_freq_dict

    def fourmer_freq_fetcher(self, seq):
        """This function is designed to fragment the paired end read sequence agreement into fourmers and then access the fourmer frequencies from the fourmer frequency dictionary.

        Input: fourmer frequency dictionary

        Output: A list of fourmer frequencies as correspond to the fourmers of the agreement sequence.
        """
        read_agreement_fourmers = [seq[i:i + 4] for i in range(0, len(seq), 4)]

        read_agreement_fourmer_frequencies = [self.assembled_fourmer_freq_dict[fourmer] for fourmer in read_agreement_fourmers]

        print("seq: ", seq, "read agreement fourmer frequencies", read_agreement_fourmer_frequencies)
        if read_agreement_fourmer_frequencies:
            return read_agreement_fourmer_frequencies
        else:
            return seq


    def max_agreement_finder(self, agreement_seq, seq_1, seq_2):
        """ """

        self.agreement_seq_len_sum += len(agreement_seq)
        self.total_num_agreements += 1

        running_agreement_seq_avg = self.agreement_seq_len_sum / self.total_num_agreements

        print("running_agreement_seq_avg: ", running_agreement_seq_avg)

        agreement_position_1 = seq_1.find(agreement_seq)
        agreement_position_2 = seq_2.find(agreement_seq)

        print("agreement positions", agreement_position_1, agreement_position_2)

        if len(agreement_seq) > 4:
            longer_seq = agreement_seq
            print("longer_seq~~~~~~~~~~~", longer_seq)


    def null_model_comparison(self):
        """"""




    def paired_end_sequence_correspondence_graphing(self):
        """Develops a forwards correspondence graph used in acyclic graphing alignment per taxa."""

        for read_head in self.assembled_paired_end_read_dict.keys():
            forward_correspondence_list = []
            backwards_correspondence_list = []
            self.forward_correspondence_dict[read_head] = []
            self.correspondence_node_dict[read_head] = []
            seq_1 = self.assembled_paired_end_read_dict[read_head]

            for other_read_head in self.assembled_paired_end_read_dict.keys():
                if read_head != other_read_head:
                    seq_2 = self.assembled_paired_end_read_dict[other_read_head]

                    for j in range(4, len(self.seq_dict_by_head[read_head])-4, 2):

                        # end of sequence one is the start of sequence 2
                        if seq_1[j:] == seq_2[:len(seq_1[j:])]:
                            forward_correspondence_list.append(other_read_head)
                            agreement_seq = seq_1[j:]

                            print("end of sequence one is the start of sequence 2, agreement_seq", agreement_seq)

                            max_agreement_seq = self.max_agreement_finder(agreement_seq, seq_1, seq_2)

                            print("max_agreement_seq", max_agreement_seq)

                            # a node is developed,  holding a series of attributes
                            self.correspondence_node_dict[read_head].append((other_read_head, agreement_seq, seq_2, self.fourmer_freq_fetcher(agreement_seq)))

                        # start of sequence one is the end of sequence 2
                        if seq_1[:j] == seq_2[len(seq_1[j:]):]:
                            backwards_correspondence_list.append(other_read_head)
                            agreement_seq = seq_1[:j]


                            print("start of sequence one is the end of sequence 2, agreement_seq", agreement_seq)

                            self.correspondence_node_dict[read_head].append((other_read_head, agreement_seq, seq_2, self.fourmer_freq_fetcher(agreement_seq)))

                            max_agreement_seq = self.max_agreement_finder(agreement_seq, seq_1, seq_2)

                            print("max_agreement_seq", max_agreement_seq)

                        # start of sequence one is the end of reversed sequence 2
                        if seq_1[:j] == seq_2[len(seq_1[j:]):][::-1]:
                            backwards_correspondence_list.append(other_read_head)
                            agreement_seq = seq_1[:j]

                            print("start of sequence one is the end of reversed sequence 2, agreement_seq", agreement_seq)

                            self.correspondence_node_dict[read_head].append((other_read_head, agreement_seq, seq_2, self.fourmer_freq_fetcher(agreement_seq)))

                            max_agreement_seq = self.max_agreement_finder(agreement_seq, seq_1, seq_2)

                            print("max_agreement_seq", max_agreement_seq)

                        # end of sequence 1 is the beginning of reversed sequence 2
                        if seq_1[j:] == seq_2[:len(seq_1[j:])][::-1]:
                            forward_correspondence_list.append(other_read_head)
                            agreement_seq = seq_1[j:]

                            print("end of sequence 1 is the beginning of reversed sequence 2, agreement_seq", agreement_seq)

                            self.correspondence_node_dict[read_head].append((other_read_head, agreement_seq, seq_2, self.fourmer_freq_fetcher(agreement_seq)))

                            max_agreement_seq = self.max_agreement_finder(agreement_seq, seq_1, seq_2)

                            print("max_agreement_seq", max_agreement_seq)

        print(self.correspondence_node_dict)

    def driver(self):

        # self.file_converter()

        self.file_parse()
        self.paired_end_assembly()
        self.null_model()
        self.paired_end_sequence_correspondence_graphing()


def main():
    """Access class and call driver member function."""

    class_access = tree_maker()
    class_access.driver()


if __name__ == '__main__':
    main()