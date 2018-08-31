from Bio import SeqIO, SeqRecord, Alphabet
from subprocess import call
from collections import Counter
import os, sys, argparse, logging

#Argparser stuff
argParser = argparse.ArgumentParser(description='concatenateSeqs joins a collection of contigs together using a '
                                    'complete sequence as reference. Requires blast 2.5+ and the Biopython library.\n'
                                    'WARNING: The algorithm used for the assembly is not perfect - if the output looks '
                                                'wrong, it\'s probably the assembler\'s fault.')
argParser.add_argument('-i', '--input', nargs = 1, required = True, type=str,
                       help='fasta file containing the contigs to be merged')
argParser.add_argument('-r', '--reference', nargs= 1, required = True, type=str,
                       help='fasta file containing the reference genome')
argParser.add_argument('-o', '--output', nargs=1, required = True, type=str,
                       help='prefix for the output files')
argParser.add_argument('--join_threshold', nargs='?', default=5000, type=int,
                       help='maximum distance between contigs to be joined together (default: 5000). A negative value '
                            'will join all contigs in order (same result for all negative values)')
argParser.add_argument('--cut_threshold', nargs='?', default=0.99, type=float,
                       help= 'The script will attempt to cut contigs in order to align them to the reference. Value '
                             'given will mark the sequence percentage that a gap needs in order to be cut. If value is'
                             'zero program wont cut. Values between 0-1 (Default: 0.99)')
argParser.add_argument('--outfmt', nargs='?', default='fasta', type=str,
                       help='output format, either \'genbank\' or \'fasta\' (Default: fasta)')


#Initialize variables
#blastpath = '/opt/ncbi-blast-2.6.0+/bin/'
blastpath = '/home/rohit/alberto/ncbi-blast-2.7.1+/bin/'
ref_fasta = argParser.parse_args(sys.argv[1:]).reference[0]
contig_fasta = argParser.parse_args(sys.argv[1:]).input[0]
output_prefix = argParser.parse_args(sys.argv[1:]).output[0]
join_threshold = argParser.parse_args(sys.argv[1:]).join_threshold
cut_threshold = argParser.parse_args(sys.argv[1:]).cut_threshold
output_type = argParser.parse_args(sys.argv[1:]).outfmt

#Initalize log
def setup_logger(logger_name, log_file, level=logging.DEBUG):
    l = logging.getLogger(logger_name)
    formatter = logging.Formatter('%(message)s')
    file_handler = logging.FileHandler(log_file, mode='a')
    file_handler.setFormatter(formatter)
    l.setLevel(level)
    l.addHandler(file_handler)


setup_logger('concatenate', 'Concatenate_7_{}.log'.format(output_prefix))
logger = logging.getLogger('concatenate')

logger.info('Parameters passed to script: ref_fasta: {}, contig_fasta: {}, output_prefix: {}, join_threshold: {}, '
            'cut_threshold: {}, output_type: {}'.format(ref_fasta, contig_fasta,  output_prefix, join_threshold,
                                                        cut_threshold, output_type))

contig_list = {}
removed_contig_list = []
ref_list = []
similarity_threshold = None


def assemble_contigs(contig_list, reference, join_threshold):
    sol = []
    raw_list = []
    for contig in contig_list:
        contig.set_borders(reference)
        if contig.ref_start is None:
            sol.append(contig)
        else:
            raw_list.append(contig)
    contigs_by_start = sorted(raw_list, key=lambda contig: contig.ref_start)
    contigs_by_end = sorted(raw_list, key=lambda contig: contig.ref_end)
    if join_threshold < 0:
        logger.info('\tnegative join_threshold. Must join all the contigs')
    while len(contigs_by_end) > 0:
        actual_contig = contigs_by_end[0]
        name = actual_contig.id
        logger.info('\tLooking for contigs near {}'.format(name))
        sequence = actual_contig.seq
        start = actual_contig.ref_start
        end = actual_contig.ref_end
        dict = actual_contig.blasts_dict
        contigs_by_start.remove(actual_contig)
        used_contigs = []
        j=1
        for contig in contigs_by_start:
            if actual_contig.ref_seq == contig.ref_seq and (join_threshold < 0 or
                                                            end + join_threshold > contig.ref_start):
                name += contig.id
                sequence += contig.seq
                start = min(start, contig.ref_start)
                end = max(end, contig.ref_end)
                dict[reference] += contig.blasts_dict[reference]
                used_contigs.append(contig)
                if join_threshold >= 0:
                    logger.info('contig {} is near enough: end + join_threshold > contig.ref_start --> {} + {} = {} > '
                                '{}'.format(contig.id, end, join_threshold, end + join_threshold, contig.ref_start))
                j+=1
            else:
                break
        for elem in used_contigs:
            contigs_by_start.remove(elem)
            contigs_by_end.remove(elem)
        new_contig = Sequence(SeqRecord.SeqRecord(sequence, id=name), ref_seq=actual_contig.ref_seq)
        new_contig.ref_start = start
        new_contig.ref_end = end
        new_contig.blasts_dict = dict
        sol.append(new_contig)
        contigs_by_end.remove(actual_contig)
        if j == 1:
            logger.info('No contigs near enough to join')
        else:
            logger.info('Search finished. {} contigs have been joined\n'.format(j))
    logger.info('Finish assembling contigs with reference {}\n'.format(reference))
    return sol


class Sequence:
    def __init__(self, record, ref_seq = None):
        self.ref_seq = ref_seq
        self.id = record.id
        self.seq = record.seq
        self.seq.description = 'None'
        self.seq.alphabet = Alphabet.generic_dna
        self.blasts_dict = {}
        self.ref_start = None
        self.ref_end = None
        self.strand = None
        self.length = len(self.seq)
        self.removed_reason = None

    def add_blast_hit(self, blast_hit, ref_seq):
        if ref_seq.id not in self.blasts_dict.keys():
            self.blasts_dict[ref_seq.id] = [blast_hit]
        else:
            self.blasts_dict[ref_seq.id].append(blast_hit)

    def choose_ref_seq(self):
        logger.info('starting choose_ref_seq for contig {}'.format(self.id))
        #If there is only blast for one  ref sequence, then there is no need for processing
        if len(self.blasts_dict.keys()) < 1:
            logger.info('{} got no references, so it will be discarded'.format(self.id))
            self.removed_reason = 'No blast_hits assigned'
            self.ref_seq = 'None'
            removed_contig_list.append(self)
        elif len(self.blasts_dict.keys()) == 1:
            for key in self.blasts_dict.keys():
                logger.info('only one reference, assigning {} as ref_seq'.format(key))
                self.ref_seq = key

    def cut_sequence(self, ref_seq):
        logger.info('cutting sequence {}'.format(self.id))
        raw_blast_list = self.blasts_dict[ref_seq.id]
        blast_list = []
        for blast in raw_blast_list:
            if blast.similarity > similarity_threshold:
                blast_list.append(blast)
        blast_list.sort(key = lambda match:match.ref_start)
        max_diff = ref_seq.length * 0.01
        max_diff_cand = None
        for i in range(0, len(blast_list)-1):
            actual_blast = blast_list[i]
            next_blast = blast_list[i + 1]
            start_next_blast = next_blast.ref_start
            end_actual_blast = actual_blast.ref_end
            diff = abs(start_next_blast - end_actual_blast)
            if diff > max_diff + ref_seq.length * 0.01:
                max_diff = diff
                max_diff_cand = i
        if max_diff_cand is None:
            return None
        #A gap has been found
        new_dict = {}
        new_dict[ref_seq.id] = []
        cut_dtce = 0
        for j in range(0, max_diff_cand+1):
            if cut_dtce < blast_list[j].query_end < self.length:
                cut_dtce = blast_list[j].query_end
            new_dict[ref_seq.id].append(blast_list[j])
        if cut_dtce == 0:
            return [self, None]
        new_seq_1 = Sequence(SeqRecord.SeqRecord(self.seq[:cut_dtce], id=self.id + '_1'), ref_seq=self.ref_seq)
        new_seq_2 = Sequence(SeqRecord.SeqRecord(self.seq[cut_dtce:], id=self.id + '_2'), ref_seq=self.ref_seq)
        new_seq_1.ref_start = blast_list[0].ref_start
        new_seq_1.ref_end = blast_list[max_diff_cand].ref_end
        new_seq_1.blasts_dict = new_dict

        other_dict = {}
        other_dict[ref_seq.id] = []
        new_seq_2.ref_start = blast_list[max_diff_cand + 1].ref_start
        for j in range(max_diff_cand + 1, len(blast_list)):
            other_dict[ref_seq.id].append(blast_list[j])
        new_seq_2.ref_end = blast_list[-1].ref_end
        new_seq_2.blasts_dict = other_dict
        logger.info('2 new sequences: {} ({} - {}) w/ length {} and {} ({} - {}) w/ length {}'.format(
            new_seq_1.id, new_seq_1.ref_start, new_seq_1.ref_end, new_seq_1.length,
            new_seq_2.id, new_seq_2.ref_start, new_seq_2.ref_end, new_seq_2.length))
        return [new_seq_1, new_seq_2]

    def is_duplicated(self, final_dict):
        for key in final_dict:
            for contig in final_dict[key]:
                if self.seq in contig.seq:
                    logger.info('{} seq is included in {} seq, so will be discarded'.format(self.id, contig.id))
                    self.removed_reason = 'Duplicated contig'
                    removed_contig_list.append(self)
                    return True
        logger.info('Contig {} seems fine. No need to delete'.format(self.id))
        return False

    def place_contig_on_reference(self, ref_seq):
        min_start = self.ref_start
        max_end = self.ref_end
        if min_start is not None and max_end is not None:
            if max_end - min_start > ref_seq.length * cut_threshold:
                logger.warning('{} is considered for cutting: max_end - min_start > ref_seq.length * cut_threshold --> '
                               '{} - {} = {} > {} = {} * {}'.format(self.id, max_end, min_start, max_end - min_start,
                                                    ref_seq.length * cut_threshold, ref_seq.length, cut_threshold))
                return 'CUT'
            else:
                self.ref_start = min_start
                self.ref_end = max_end
                logger.info('{} has found a good match in {}: both max_end and min_start are not None but max_end - '
                            'min_start < ref_seq.length * cut_threshold --> {} - {} = {} < {} = {} * {}'.format(self.id,
                                                        ref_seq.id, max_end, min_start, max_end - min_start,
                                                        ref_seq.length * cut_threshold, ref_seq.length, cut_threshold))
                return 'OK'
        else:
            logger.warning('{} has no operating size, so it will be discarted (max_end and min_start are None)'.format(self.id))
            return 'NO'

    def reverse(self, key):
        blast_list = sorted(self.blasts_dict[key], key=lambda match: match.length, reverse=True)
        if blast_list[0].strand == 1:
            self.seq = self.seq.reverse_complement()
            logger.info('contig {} has been reversed'.format(self.id))

    def set_borders(self, key):
        #This function is an extension from get_operating_size
        min_start = None
        max_end = None
        for match in self.blasts_dict[key]:
            if len(self.blasts_dict[key]) == 1:
                min_start = match.ref_start
                max_end = match.ref_end
            elif match.similarity > similarity_threshold and match.length > 500:
                min_candidate = match.ref_start
                max_candidate = match.ref_end
                if min_start is None or min_start > min_candidate:
                    min_start = min_candidate
                if max_end is None or max_end < max_candidate:
                    max_end = max_candidate
            elif join_threshold < 0:
                min_candidate = match.ref_start
                max_candidate = match.ref_end
                if min_start is None or min_start > min_candidate:
                    min_start = min_candidate
                if max_end is None or max_end < max_candidate:
                    max_end = max_candidate
        self.ref_start = min_start
        self.ref_end = max_end

    def split_to_choose(self):
        logger.info('\tsplit_to_choose: now evaluating contig {}'.format(self.id))
        sorted_full_list = []
        for key in self.blasts_dict.keys():
            sorted_list = sorted(self.blasts_dict[key], key=lambda match: match.query_start)
            for match in sorted_list:
                if match.similarity > similarity_threshold and match.length > 500:
                    sorted_full_list.append(match)
        # Cases where there is 0 blasts need to go to choose_ref_seq() so them will be discarted
        if len(sorted_full_list) == 0:
            logger.info('{} got no useful blasts, so it will be discarded --> No blasts with similarity above {} and '
                        'length above 500 bp\n'.format(self.id, similarity_threshold))
            self.removed_reason = 'No blast_hits assigned'
            self.ref_seq = 'None'
            removed_contig_list.append(self)
            return []
        i = 0
        j = 0
        sol = []
        while i < len(sorted_full_list):
            actual_blast = sorted_full_list[i]
            key = actual_blast.parent_seq
            if i == 0:
                query_starting_pos = 0
            else:
                query_starting_pos = actual_blast.query_start
            if i == len(sorted_full_list)-1:
                query_ending_pos = self.length
            else:
                query_ending_pos = sorted_full_list[i+1].query_start
            reference_starting_pos = actual_blast.ref_start
            reference_ending_pos = actual_blast.ref_end
            dict = {}
            dict[key] = [actual_blast]
            new_seq = Sequence(SeqRecord.SeqRecord(self.seq[query_starting_pos:query_ending_pos],
                                                   id=self.id + '_{}'.format(j)), ref_seq=self.ref_seq)
            j += 1
            new_seq.ref_start = reference_starting_pos
            new_seq.ref_end = reference_ending_pos
            new_seq.blasts_dict = dict
            i += 1
            if new_seq.length == 0:
                new_seq.removed_reason = 'No real contig (sequence length = 0)'
                new_seq.ref_seq = 'None'
                removed_contig_list.append(new_seq)
                logger.info('contig {} removed because its sequence length was 0'.format(new_seq.id))
            else:
                new_seq.choose_ref_seq()
                sol.append(new_seq)
        logger.info('contig {} has been split into {} new contigs\n'.format(self.id, j))
        return sol

    def write_fasta(self):
        with open('ref_seq_temp.fasta', 'w') as fasta_file:
            fasta_file.write('>' + self.id + '\n')
            fasta_file.write(str(self.seq))


class BlastHit:
    def __init__(self, parent, length, similarity,  ref_start_pos, ref_end_pos, query_start_pos, query_end_pos,
                 bit_score):
        self.parent_seq = parent
        self.length = int(length)
        self.similarity = float(similarity)
        self.ref_start = int(ref_start_pos)
        self.ref_end = int(ref_end_pos)
        self.query_start = int(query_start_pos)
        self.query_end = int(query_end_pos)
        self.strand = 0
        self.bit_score = bit_score

        #check the strand of the fragment, fix the start/end positions if it is from the complementary strand
        if self.ref_start > self.ref_end:
            self.strand = 1
            self.ref_start = int(ref_end_pos)
            self.ref_end = int(ref_start_pos)

        #Fix bitscores in scientific notation
        if 'e' in self.bit_score:
            bit_score_split = self.bit_score.split('e')
            if bit_score_split[1][0] == '+':
                self.bit_score = float(bit_score_split[0]) * pow(10,  int(bit_score_split[-1]))

            elif bit_score_split[1][0] == '-':
                self.bit_score = float(bit_score_split[0]) * pow(10,  (int(bit_score_split[-1])* -1))

        else:
            self.bit_score = float(bit_score)

# Main body. Program starts here

#First step: Load up references and contig sequences
logger.info('loading ref file')
with open(ref_fasta, 'r') as reference_file:
    parsed_references = SeqIO.parse(reference_file, 'fasta')
    for record in parsed_references:
        ref_list.append(Sequence(record))
        logger.debug('\tref sequence {} loaded, with length {}'.format(record.id, len(record.seq)))

logger.info('loading fasta file')
with open(contig_fasta, 'r') as contig_file:
    parsed_contigs = SeqIO.parse(contig_file, 'fasta')
    for record in parsed_contigs:
        contig_list[record.id] = (Sequence(record))

#Perform blasts
logger.info('doing blasts')
good_contigs = {}
similarity_list = []
for ref_seq in ref_list:
    good_contigs[ref_seq.id] = []
    ref_seq.write_fasta()
    call([blastpath + 'makeblastdb', '-in', 'ref_seq_temp.fasta', '-out', 'dbTest', '-dbtype', 'nucl'])
    call([blastpath + 'blastn', '-query', contig_fasta, '-db', 'dbTest', '-out', 'blast_results.blastn', '-num_threads',
          '4', '-outfmt', '6'])
    with open('blast_results.blastn', 'r') as blast_file:
        for line in blast_file:
            split_line = line.split('\t')
            if split_line[0] in contig_list.keys() and split_line[1] == ref_seq.id:
                similarity_list.append(float(split_line[2]))
                contig_list[split_line[0]].add_blast_hit(BlastHit(split_line[1], split_line[3], split_line[2],
                                split_line[8], split_line[9], split_line[6], split_line[7], split_line[11]), ref_seq)
    os.remove('blast_results.blastn')
    os.remove('ref_seq_temp.fasta')
os.remove('dbTest.nsq')
os.remove('dbTest.nin')
os.remove('dbTest.nhr')
res = Counter(similarity_list).most_common(5)
#Look for the first four non-100.0 positions. Notice that only one of those can be 100.0 and it should be always at top
#five positions, but it wont be always the first one
first = res[1][0]
second = res[2][0]
third = res[3][0]
fourth = res[4][0]
if first == 100.0:
    first = res[0][0]
elif second == 100.0:
    second = res[0][0]
elif third == 100.0:
    third = res[0][0]
elif fourth == 100.0:
    fourth = res[0][0]
similarity_threshold = (first + second + third + fourth)/4
print('\nStep 1 done: References and contigs loaded. Blasts performed\n')

#Step 2: in order to get better matches split the contigs into smaller ones based on its blasts. Then assign each contig
#to a reference based on some values
logger.info('number of reference sequences: {}'.format(len(ref_list)))
logger.info('number of contig sequences: {}'.format(len(contig_list.keys())))
logger.info('choosing ref_seqs for each contig')
contig_split_dict = {}
for contig in contig_list:
    list = contig_list[contig].split_to_choose()
    for elem in list:
        #Each elem is a new contig
        contig_split_dict[elem.id] = elem
print('Step 2 done: Contigs split and reference chosen\n')

logger.info('\n\n')
logger.info('choosing ref_seq results:')
#Get some logger data
count_contig_dict = {}
ref_seq = None
for key in contig_split_dict:
    try:
        ref_seq = contig_split_dict[key].ref_seq
    except AttributeError:
        ref_seq = 'None'
    if ref_seq not in count_contig_dict.keys():
        count_contig_dict[ref_seq] = 1
    else:
        count_contig_dict[ref_seq] += 1
for key in count_contig_dict.keys():
    logger.info('{}: {}'.format(key, count_contig_dict[key]))

#Keep only non-removed contigs (obviously)
clean_contig_list = []
for key in contig_split_dict:
    if contig_split_dict[key] not in removed_contig_list:
        clean_contig_list.append(contig_split_dict[key])

logger.info('\n\n')
logger.info('starting place contig on reference')
#Step 3: evaluate each contig into its reference: check whose are good, bad and those who need to be cut
#Notice that contigs were cut becore. This time are cut only the contigs which start and end positions are on each side
#of the reference
for i, ref in enumerate(ref_list):
    logger.info('starting to place contig on reference {} for contigs assigned to it'.format(ref.id))
    count = 0
    totalOKs = 0
    for contig in clean_contig_list:
        if contig.ref_seq == ref.id:
            logger.info('\tprocessing contig {} assigned to {}'.format(contig.id, contig.ref_seq))
            count += 1
            check = contig.place_contig_on_reference(ref)
            if check == 'OK':
                totalOKs += 1
                good_contigs[ref.id].append(contig)
            elif check == 'NO':
                contig.removed_reason = 'No relevant blast matches'
                removed_contig_list.append(contig)
            elif check == 'CUT':
                new_seqs = contig.cut_sequence(ref)
                if new_seqs is not None and cut_threshold > 0:
                    totalOKs += 1
                    good_contigs[ref.id].append(new_seqs[0])
                    if new_seqs[1] is not None:
                        totalOKs += 1
                        good_contigs[ref.id].append(new_seqs[1])
                        logger.info('Seq1 {} and seq2 {} with ref_seq {} added'.format(new_seqs[0].id, new_seqs[1].id,
                                                                                       ref.id))
                    else:
                        logger.info('Seq1 {} with ref_seq {} added'.format(new_seqs[0].id, ref.id))
                else:
                    if new_seqs is None:
                        logger.info('{} cannot be cut, so it will be discarded (new seq is None)'.format(contig.id))
                    else:
                        logger.info('cut_threshold is 0 or lower, cant cut {}'.format(contig.id))
                    contig.removed_reason = 'Not a good match - Tried to cut'
                    removed_contig_list.append(contig)

    logger.info('results for {}: {}/{} placed. Total Contigs: {}'.format(ref.id, totalOKs, count,
                                                                         len(clean_contig_list)))
    total_length = 0

    for contig in good_contigs[ref.id]:
        total_length += contig.length
print('Step 3 done: Contigs evaluated for its reference\n')
logger.info('\n\n')

#Step 4: Check the length of the contigs and the length of the blasts. If the total blast surface is way smaller than
# the contig lenght discard it. Reverse the good ones
not_good_enough = []
for ref in ref_list:
    logger.info('Checking contigs assigned to reference {}'.format(ref.id))
    for contig in good_contigs[ref.id]:
        min_pos = None
        max_pos = None
        for blast in contig.blasts_dict[ref.id]:
            if min_pos is None or blast.query_start < min_pos:
                min_pos = blast.query_start
            if max_pos is None or blast.query_end > max_pos:
                max_pos = blast.query_end
        blast_length = max_pos - min_pos
        if blast_length/contig.length < 0.1:
            logger.info('{} has nearly no blasts for its size, so it will be discarted: blast_length/contig.length = '
                        '{} / {} = {} < 0.1'.format(contig.id, blast_length, contig.length, blast_length/contig.length))
            not_good_enough.append((contig, ref.id))
        else:
            logger.info('{} is a good contig. Check strand for reverse (if needed)'.format(contig.id))
            contig.reverse(ref.id)

for contig, key in not_good_enough:
    contig.removed_reason = 'No relevant information contig'
    removed_contig_list.append(contig)
    good_contigs[key].remove(contig)
print('Step 4 done: Contig sequences reversed if needed\n')
logger.info('\n\n')

logger.info('Safety measure: checking if some contigs are duplicated')
#Safety measure: check for duplicates in contigs. Usually the program wont have duplicates, but one never knows for sure
final_dict = {}
for key in good_contigs:
    final_dict[key] = []
    contig_list = sorted(good_contigs[key], key=lambda contig:contig.length, reverse=True)
    for contig in contig_list:
        if not contig.is_duplicated(final_dict):
            final_dict[key].append(contig)

logger.info('\n\n')
#Step 5: Almost finish. Try to assemble contigs based on the threshold given to the program
assembled_contigs = {}
for ref in ref_list:
    logger.info('start assembling contigs for {}, minimum distance to join: {} bp'.format(ref.id, join_threshold))
    assembled_contigs[ref.id] = assemble_contigs(final_dict[ref.id], ref.id, join_threshold)

for ref in ref_list:
    if ref.id not in assembled_contigs.keys():
        assembled_contigs[ref.id] = []
print('Step 5 done: Contigs assembled for given parameters\n')

#Last step: Final files creation
logger.info('writing output Files')
for ref in assembled_contigs:
    good_list = []
    for i, contig in enumerate(assembled_contigs[ref]):
        good_list.append(SeqRecord.SeqRecord(contig.seq, id='contig_'+ str(i)))
    with open(output_prefix + ref + '.good_contigs.' + output_type, 'w') as output_good:
        SeqIO.write(good_list, output_good, output_type)


bad_list = []
for contig in removed_contig_list:
    bad_list.append(SeqRecord.SeqRecord(contig.seq, contig.id))
with open(output_prefix + '.badContigs.' + output_type, 'w') as output_bad:
    SeqIO.write(bad_list, output_bad, output_type)

print('Program finished\n')