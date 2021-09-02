import pysam
import os.path
import argparse
import Constants
from EM import automatica_predict_EM
from os import listdir
from os.path import isdir,join,isfile
import sys
############################################################################
'''
This is used to calculate ACGT frequency per base
'''
############################################################################
SEQUENCE_ERROR = Constants.SEQUENCE_ERROR
############################################################################
def fetch_filter_polymorphic_sites(genomeLen, genomeName, in_bam_file, res_dir):
    # calculate the ACGT frequency for each position
    # initial
    # genomeLen = 1853160
    # genomeName = 'NC_009515.1'
    # inName = '/media/student/study_working/project/project12/point1/data/simulated_data/simulatedReads_sortedBam/default/0_sorted.bam'

    print('generate polymorphic sites')
    genomeLen = int(genomeLen)
    res = []
    for _ in range(genomeLen):
        res.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})


    samfile = pysam.AlignmentFile('%s' % in_bam_file, 'rb')
    posID = 0
    # genome start location 0, the following pileup parameter is the same.
    for pileupColumn in samfile.pileup(genomeName, 0, genomeLen):
    #for pileupColumn in samfile.pileup(genomeName, 1000000, 1001000):
        if posID % 100000 == 0:
            print(posID)
        posID += 1

        pos = pileupColumn.pos
        count = pileupColumn.n

        tmp = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        # print ("\ncoverage at base %s = %s" % (pileupColumn.pos, pileupColumn.n))
        for pileupread in pileupColumn.pileups:
            #if pileupread.is_del or pileupread.is_refskip:
            #    print 'here'

            if (not pileupread.is_del) and (not pileupread.is_refskip):
                # query position is None if is_del or is_refskip is set.
                # print ('\tbase in read %s = %s' %
                #       (pileupread.alignment.query_name,
                #        pileupread.alignment.query_sequence[pileupread.query_position]))
                nucl = pileupread.alignment.query_sequence[pileupread.query_position]
                nucl = nucl.upper()
                tmp[nucl] += 1

        for key, value in tmp.items():
            res[pos][key] = value

    samfile.close()

    # calculate polymorphics sites without filter sequence error
    nuclOrder = ['A', 'C', 'G', 'T']
    polymorphicSites = []
    for loc in range(len(res)):
        pos_nucl = res[loc]
        formatOutput = [loc]

        flag_nonZero = 0
        for nucl in nuclOrder:
            nucl_count = pos_nucl[nucl]
            formatOutput.append(nucl_count)
            if nucl_count != 0:
                flag_nonZero += 1

        if flag_nonZero > 1:
            polymorphicSites.append(formatOutput)

    # # filter polymorphic sites
    polymorphicSites_filter = []
    for item in polymorphicSites:
        flag_nonError = 0

        for one in item[1:]:
            if one > SEQUENCE_ERROR:
                flag_nonError += 1

        if flag_nonError > 1:
            polymorphicSites_filter.append(item)

    res_file_name = res_dir + 'filter_polymorphic_sites'

    num1 = 0

    with open('%s' %res_file_name, 'w') as f:
        for item in polymorphicSites_filter:
            one = '\t'.join([str(one) for one in item])

            if len(one) > 0:
                num1 +=1

            f.write('%s\n' %one)

    if num1 > 5:
        return True
    else:
        return False


parser = argparse.ArgumentParser()
# group1 = parser.add_argument_group('General')
parser.add_argument('--sample_name', help='Give a unique sample name',
                    default=None)
parser.add_argument('--genome_len', help='Input genome length', default=None)
parser.add_argument('--genome_name', help='Input genome name', default=None)
parser.add_argument('--genome_file_loc', help='Input genome file location', default=None)
parser.add_argument('--bam_file', help='Input sorted bam file', default=None)
parser.add_argument('--res_dir', help='result directory', default=None)

args = parser.parse_args()

#result directory
if not args.res_dir:
    script_path = os.path.dirname(os.path.realpath(__file__))
else:
    script_path = os.path.dirname(args.res_dir)

#test input
#1853160
genome_len = int(args.genome_len)
genome_name = args.genome_name
bam_file_loc = os.path.abspath(args.bam_file)
res_dir = os.path.abspath(args.res_dir) + '/'
genome_file = os.path.abspath(args.genome_file_loc)
sample_name = args.sample_name


res_tmp_dir = res_dir + sample_name + '/'
if not os.path.exists(res_tmp_dir):
    command = 'mkdir -p ' + res_tmp_dir
    os.system(command)
#generate polymorphic sites
label = fetch_filter_polymorphic_sites(genome_len, genome_name, bam_file_loc, res_tmp_dir)

if label == False:
    print('The filter polymorphic sites are too little, there will be no strain found')
    sys.exit('no haplotypes found')

#EM algorithm
polymorphic_file_name = res_tmp_dir + 'filter_polymorphic_sites'
em_dir = res_tmp_dir + 'res_em/'
num_haps_predict, snp_output_single = automatica_predict_EM((
    polymorphic_file_name, genome_file, em_dir))

if num_haps_predict == 1 or num_haps_predict == -1:
    command = 'rm -r ' + em_dir
    os.system(command)

    #write haplotype for single into file
    resName = res_tmp_dir + sample_name + '_haplotypes'
    with open('%s' % resName, 'w') as f:
        for key, value in snp_output_single.items():
            title = '>' + str(1.0)
            f.write('%s\n' % title)
            for item in value:
                f.write('%s\n' % item)

    sys.exit('no haplotypes found')


res_hap = em_dir + 'hap_' + str(num_haps_predict) + '/haplotypes'
command = 'cp ' + res_hap + ' ' + res_tmp_dir + sample_name + '_haplotypes'
os.system(command)


#remove intermediate files
#command = 'rm -r ' + em_dir
#os.system(command)







print('Done')
