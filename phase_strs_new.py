import argparse
import collections
import dendropy
import gzip
import math
import operator
import os
import random
import subprocess
import sys
import tempfile
import vcf

from cStringIO import StringIO

from str_mut_models.matrix_optimizer import MATRIX_OPTIMIZER
from str_mut_models.mutation_model   import OUGeomSTRMutationModel
from str_mut_models.main             import determine_allele_range

def compute_median(values, counts):
    total_count = sum(counts)
    target_val  = total_count/2.0
    for i in xrange(len(values)):
        if counts[i] > target_val:
            return values[i]
        target_val -= counts[i]
    return values[-1]

def read_sample_list(input_file):
    samples = set()
    data    = open(input_file, "r")
    for line in data:
        samples.add(line.strip())
    data.close()
    return samples

'''
def compare_genotype_sets(gt_dict_1, gt_dict_2):
    samples_a      = set(map(lambda x: x.split("_")[0], gt_dict_1.keys()))
    samples_b      = set(map(lambda x: x.split("_")[0], gt_dict_2.keys()))
    common_samps   = samples_a & samples_b
    correct_counts = collections.defaultdict(int)
    for sample in common_samps:
        gt_1a = gt_dict_1[sample+"_1"]
        gt_1b = gt_dict_1[sample+"_2"]
        gt_2a = gt_dict_2[sample+"_1"]
        gt_2b = gt_dict_2[sample+"_2"]

        if gt_1a == gt_2a and gt_1b == gt_2b:
            num_correct = 2
        elif gt_1a == gt_2b and gt_2a == gt_1b:
            num_correct = 2
        elif gt_1a == gt_2a or gt_1a == gt_2b or gt_1b == gt_2a or gt_1b == gt_2b:
            num_correct = 1
        else:
            num_correct = 0
        correct_counts[num_correct] += 1
'''

    



# Read the STR genotypes from the provided VCF file, normalize their lengths relative
# to the median (from all samples) and return a dictionary mapping sample names (id_1 and id_2) to the number of normalized repeats
def read_diploid_str_gts(vcf_file, chrom, pos, sample_set=None):
    vcf_reader = vcf.Reader(filename=vcf_file)
    record     = vcf_reader.next()
    gb_counts  = collections.defaultdict(int)
    while record.CHROM != chrom or record.POS != pos:
        record = vcf_reader.next()

    # Utilize all samples to determine the median allele
    for sample in record:
        if sample['GT'] is not None:
            gt_a, gt_b = map(int, sample['GT'].split('/'))
            gb_counts[len(record.alleles[gt_a])] += 1
            gb_counts[len(record.alleles[gt_b])] += 1

    median_allele = compute_median(gb_counts.keys(), gb_counts.values())
    motif_len     = len(record.INFO['MOTIF'])
    # Check that all repeat lengths are in-frame
    if len(set(map(lambda x: x%motif_len, gb_counts.keys()))) != 1:
        exit("ERROR: Some GTs in input VCF are out-of-frame")

    # Only return genotypes for requested set of samples
    nrepeats_dict = {}
    for sample in record:
        if sample_set is None or sample.sample in sample_set:
            if sample['GT'] is not None:
                gt_a, gt_b = map(int, sample['GT'].split('/'))
                nrepeats_dict[sample.sample + "_1"] = (len(record.alleles[gt_a])-median_allele)/motif_len
                nrepeats_dict[sample.sample + "_2"] = (len(record.alleles[gt_b])-median_allele)/motif_len
    return nrepeats_dict, median_allele
        
# Read the STR genotypes from the provided VCF file, normalize their lengths relative
# to the median (from all samples) and return a dictionary mapping sample names to the number of normalized repeats
def read_haploid_str_gts(vcf_file, chrom, pos, sample_set=None):
    vcf_reader   = vcf.Reader(filename=vcf_file)
    record       = vcf_reader.next()
    gb_counts    = collections.defaultdict(int)
    while record.CHROM != chrom or record.POS != pos:
        record = vcf_reader.next()

    # Utilize all samples to determine the median allele
    for sample in record:
        if sample['GT'] is not None:
            gb_counts[len(record.alleles[int(sample['GT'])])] += 1

    median_allele = compute_median(gb_counts.keys(), gb_counts.values())    
    motif_len     = len(record.INFO['MOTIF'])
    # Check that all repeat lengths are in-frame
    if len(set(map(lambda x: x%motif_len, gb_counts.keys()))) != 1:
        exit("ERROR: Some GTs in input VCF are out-of-frame")
        
    # Only return genotypes for requested set of samples
    nrepeats_dict = {}
    for sample in record:
        if sample_set is None or sample.sample in sample_set:
            nrepeats_dict[sample.sample] = (len(record.alleles[int(sample['GT'])])-median_allele)/motif_len
    return nrepeats_dict, median_allele
        
def pair_gts(nrepeats_dict, leaf_indices, pairs=None):
    if pairs is None:
        sample_order = nrepeats_dict.keys()
        random.shuffle(sample_order)
        pairs = map(lambda x: (sample_order[x], sample_order[x+1]), xrange(0, len(sample_order)/2*2, 2))
    pair_data = []
    for i in xrange(len(pairs)):
        pair_data.append((leaf_indices[pairs[i][0]], leaf_indices[pairs[i][1]], nrepeats_dict[pairs[i][0]], nrepeats_dict[pairs[i][1]]))
    return pair_data

def extract_newick_tree_from_smc(input_file, position):
    data        = gzip.open(input_file, "rb") if input_file.endswith(".gz") else  open(input_file, "r")
    name_tokens = data.readline().strip().split()
    if name_tokens[0] != "NAMES":
        exit("ERROR: First line of SMC file must contain NAMES")

    leaf_names   = name_tokens[1:]
    leaf_indices = dict(zip(name_tokens[1:], range(len(name_tokens)-1)))

    for line in data:
        tokens = line.strip().split()
        if tokens[0] != "TREE":
            continue
        if int(math.floor(float(tokens[1]))) <= position and int(math.ceil(float(tokens[2]))) >= position:
            data.close()
            tree = dendropy.Tree(stream=StringIO(tokens[3]), schema="newick")
            return tree, leaf_names, leaf_indices
    data.close()
    exit("ERROR: Failed to extract tree from SMC file %s. No region contained position %d"%(input_file, position))

def write_factor_graph(tree, optimizer, pairs, min_allele, max_allele, gen_per_len, graph_file):
    num_alleles = max_allele - min_allele + 1
    output      = open(graph_file, "w") 

    # Number of factor blocks (# tree edges + # pairings)
    num_factors = len(filter(lambda x: x.tail_node is not None and x.head_node is not None, tree.get_edge_set())) + len(pairs)
    output.write("%d\n\n"%(num_factors))

    # Iterate through all parent-child factors
    print("Iterating through tree factors")
    for node in tree.nodes():
        if node.parent_node is None:
            continue

        output.write("2\n")                                                            # Number of variables in factor block (2)
        output.write("%s %s\n"%(node.parent_node.get_node_str(), node.get_node_str())) # Names of variables (parent and child)
        output.write("%d %d\n"%(num_alleles, num_alleles))                             # Number of possible values for each variable (i.e. number of alleles)
        output.write("%d\n"%(num_alleles*num_alleles))                                 # Number of non-zero transitions

        # Fetch transition probabilities
        trans_matrix = optimizer.get_transition_matrix(int(node.edge_length*gen_per_len))

        # Write out each parent->child transition probability, with the parent node's gt changing the fastest
        index = 0
        for child in xrange(min_allele, max_allele+1):
            for parent in xrange(min_allele, max_allele+1):
                trans_prob = trans_matrix[child-min_allele, parent-min_allele]
                output.write("%d %g\n"%(index, trans_prob))
                index += 1

        # Empty line separates factor blocks
        output.write("\n")

    # Iterate through all haploid-haploid pairing factors
    print("Iterating through node pairings")
    for i in xrange(len(pairs)):
        output.write("2\n")                                 # Number of variables in factor block (2)
        output.write("%s %s\n"%(pairs[i][0], pairs[i][1]))  # Names of variables
        output.write("%d %d\n"%(num_alleles, num_alleles))  # Number of possible values for each variable (i.e. number of alleles)

        # Write out indices for the phasings
        if pairs[i][2] != pairs[i][3]:
            output.write("2\n") # Number of non-zero transitions (either phasing)
            index_one = (pairs[i][2] - min_allele)*num_alleles + (pairs[i][3] - min_allele)
            index_two = (pairs[i][3] - min_allele)*num_alleles + (pairs[i][2] - min_allele)
            output.write("%d 0.5\n"%(index_one))
            output.write("%d 0.5\n"%(index_two))
        else:
            output.write("1\n") # Number of non-zero transitions (single phasing)
            index_one = (pairs[i][2] - min_allele)*num_alleles + (pairs[i][3] - min_allele)
            output.write("%d 0.5\n"%(index_one))

        # Empty line separates factor blocks
        output.write("\n")

    output.close()

def process_phasing_result(phasing_result, leaf_names, min_allele, min_confidence=0.5):
    num_skip    = 0
    num_correct = 0
    num_phased  = 0    
    num_homoz   = 0
    gt_dict     = {}
    errors      = []
    for line in phasing_result.strip().split("\n"):
        tokens         = line.strip().split()
        id_1, id_2     = map(int, tokens[0:2])
        gt_a, gt_b     = map(int, tokens[2:4])
        conf_1, conf_2 = map(float, tokens[4:6])
        #print("%s\t%s\t%d\t%d\t%f\t%f"%(leaf_names[id_1], leaf_names[id_2], gt_a+min_allele, gt_b+min_allele, conf_1, conf_2))

        if gt_a == gt_b:
            gt_dict[leaf_names[id_1]] = gt_a + min_allele
            gt_dict[leaf_names[id_2]] = gt_a + min_allele
            num_homoz += 1
            continue
        
        if conf_1 > 0.5:
            if conf_1 > min_confidence:
                num_correct  += 1
                num_phased   += 1
                gt_dict[leaf_names[id_1]] = gt_a + min_allele
                gt_dict[leaf_names[id_2]] = gt_b + min_allele
            else:
                num_skip += 1
        else:
            if conf_2 > min_confidence:
                num_phased   += 1
                gt_dict[leaf_names[id_1]] = gt_b + min_allele
                gt_dict[leaf_names[id_2]] = gt_a + min_allele
                errors.append("%d_%d"%(gt_a, gt_b))
            else:
                num_skip += 1
    accuracy_string = "PHASING_ACCURACY\t%d\t%d\t%d\t%d\t%s"%(num_skip, num_homoz, num_correct, num_phased, ",".join(errors))
    return gt_dict, accuracy_string

def process_imputation_result(imputation_result, leaf_names, min_allele, max_allele, min_confidence=0.0):
    gts        = {}
    posteriors = {}
    dosages    = {}
    #print("%d\t%d"%(min_allele, max_allele))
    for line in imputation_result.strip().split("\n"):
        toks    = map(lambda x: x.replace("(","").replace(")","").strip(), line.split(","))
        node_id = int(toks[0].split("{")[0])
        probs   = map(float, toks[1:])
        #print("%s\t%s"%(leaf_names[node_id], "\t".join(toks[1:])))
        index,max_prob = max(enumerate(probs), key=operator.itemgetter(1))
        dosage = sum(map(lambda x,y: x*y, xrange(len(probs)), probs)) + min_allele
        gts[leaf_names[node_id]]        = index+min_allele
        posteriors[leaf_names[node_id]] = max_prob
        dosages[leaf_names[node_id]]    = dosage
    return gts, posteriors, dosages

def write_haploid_vcf(input_vcf_file, chrom, pos, nrepeats_dict, median_allele, output_vcf_file, posteriors=None, dosages=None):
    vcf_reader = vcf.Reader(filename=input_vcf_file)
    record     = vcf_reader.next()
    while record.CHROM != chrom or record.POS != pos:
        record = vcf_reader.next()

    motif_len  = len(record.INFO['MOTIF'])
    gt_indexes = {}
    gt_indexes[len(record.REF)] = 0
    for i in xrange(len(record.ALT)):
        gt_indexes[len(str(record.ALT[i]))] = i+1

    format_field = "GT:GB" 
    if posteriors is not None:
        format_field += ":POSTERIOR"
    if dosages is not None:
        format_field += ":DOSAGE"

    tokens = [record.CHROM, record.POS, ".", record.REF, ",".join(map(str, record.ALT)), ".", "PASS", "MOTIF=%s"%(record.INFO['MOTIF']), format_field]
    for sample,nreps in sorted(nrepeats_dict.items()):
        bp_len = nreps*motif_len + median_allele
        token  = "%d:%d"%(gt_indexes[bp_len], bp_len-len(record.REF)) 

        if posteriors is not None:
            post   = posteriors[sample]
            token += ":%f"%(post)
        if dosages is not None:
            dosage = motif_len*dosages[sample] + median_allele
            token += ":%f"%(dosage)

        # Add sample information to list
        tokens.append(token)

    output = open(output_vcf_file, "w")
    output.write("##fileformat=VCFv4.2\n")
    output.write("##INFO=<ID=MOTIF,Number=1,Type=String,Description=\"Canonical repeat motif\">\n")
    output.write("##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">\n")
    output.write("##FORMAT=<ID=GB,Number=1,Type=Integer,Description=\"Genotype given in bp difference from reference\">\n")
    if posteriors is not None: 
        output.write("##FORMAT=<ID=POSTERIOR,Number=1,Type=Float,Description=\"Posterior probability for allele\">\n")
    if dosages is not None:
        output.write("##FORMAT=<ID=DOSAGE,Number=1,Type=Float,Description=\"Posterior STR dosage in base pairs\">\n")
    output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(map(lambda x: str(x), sorted(nrepeats_dict.keys()))) + "\n")
    output.write("\t".join(map(str, tokens)) + "\n")
    output.close()

def write_diploid_vcf(input_vcf_file, chrom, pos, nrepeats_dict, median_allele, output_vcf_file, posteriors=None, dosages=None):
    # Only utilize samples for which both _1 and _2 suffixes are in the genotype dictionary
    sample_names = sorted(map(lambda x: x[0], filter(lambda x: x[1] == 2, collections.Counter(map(lambda x: x.split("_")[0], nrepeats_dict.keys())).items())))
    vcf_reader   = vcf.Reader(filename=input_vcf_file)
    record       = vcf_reader.next()
    while record.CHROM != chrom or record.POS != pos:
        record = vcf_reader.next()

    motif_len  = len(record.INFO['MOTIF'])
    gt_indexes = {}
    gt_indexes[len(record.REF)] = 0
    for i in xrange(len(record.ALT)):
        gt_indexes[len(str(record.ALT[i]))] = i+1

    format_field = "GT:GB" 
    if posteriors is not None:
        format_field += ":POSTERIOR"
    if dosages is not None:
        format_field += ":DOSAGE"

    tokens = [record.CHROM, record.POS, ".", record.REF, ",".join(map(str, record.ALT)), ".", "PASS", "MOTIF=%s"%(record.INFO['MOTIF']), format_field]
    for sample in sample_names:
        nreps_a  = nrepeats_dict[sample+"_1"]
        nreps_b  = nrepeats_dict[sample+"_2"]
        len_a    = motif_len*nreps_a + median_allele
        len_b    = motif_len*nreps_b + median_allele
        token    = "%d|%d:%d|%d"%(gt_indexes[len_a], gt_indexes[len_b], len_a-len(record.REF), len_b-len(record.REF))
        
        if posteriors is not None:
            post_a = posteriors[sample+"_1"]
            post_b = posteriors[sample+"_2"]
            token += ":%f|%f"%(post_a, post_b)
        if dosages is not None:
            dosage_a = motif_len*dosages[sample+"_1"] + median_allele
            dosage_b = motif_len*dosages[sample+"_2"] + median_allele
            token   += ":%f"%(dosage_a + dosage_b)

        # Add sample information to list
        tokens.append(token)

    output = open(output_vcf_file, "w")
    output.write("##fileformat=VCFv4.2\n")
    output.write("##INFO=<ID=MOTIF,Number=1,Type=String,Description=\"Canonical repeat motif\">\n")
    output.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    output.write("##FORMAT=<ID=GB,Number=1,Type=String,Description=\"Genotype given in bp difference from reference\">\n")
    if posteriors is not None: 
        output.write("##FORMAT=<ID=POSTERIOR,Number=1,Type=String,Description=\"Posterior probabilites for each allele\">\n")
    if dosages is not None:
        output.write("##FORMAT=<ID=DOSAGE,Number=1,Type=Float,Description=\"Posterior STR dosage in base pairs\">\n")
    output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n")
    output.write("\t".join(map(str, tokens)) + "\n")
    output.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smc",     required=True,  dest="smc",     type=str,   help="SMC file for the region of interest")
    parser.add_argument("--chrom",   required=True,  dest="chrom",   type=str,   help="STR's Chromosome")
    parser.add_argument("--pos",     required=True,  dest="pos",     type=int,   help="Position of STR in region")
    parser.add_argument("--mu",      required=True,  dest="mu",      type=float, help="Mutation rate for mutation model")
    parser.add_argument("--beta",    required=True,  dest="beta",    type=float, help="Length constraint for mutation model")
    parser.add_argument("--pgeom",   required=True,  dest="pgeom",   type=float, help="Geometric parameter for mutation model")
    parser.add_argument("--out",     required=True,  dest="out",     type=str,   help="Output path prefix for imputed/phased VCF")
    parser.add_argument("--vcf",     required=True,  dest="vcf",     type=str,   help="VCF containing STR calls")
    parser.add_argument("--samps",   required=False, dest="samps",   type=str,   help="File containing list of samples to consider")
    parser.add_argument("--thresh",  required=False, dest="thresh",  type=float, help="Posterior probability threshold required to report phasing", default=0.5)
    parser.add_argument("--phase",   required=False, dest="phase",   action="store_true", help="Output phasing statistics", default=False)
    parser.add_argument("--impute",  required=False, dest="impute",  action="store_true", help="Output imputation statisttics",  default=False)
    parser.add_argument("--diploid", required=False, dest="diploid", action="store_true", help="VCF contains diploid STR calls (instead of haploid calls)", default=False)    
    
    # Scaling factor from edge length to # of generations
    parser.add_argument("--gen_per_len", required=False, dest="gen_per_len", type=float, default=1.0)

    # Maximum TMRCA in generations
    parser.add_argument("--max_tmrca", required=False, dest="max_tmrca", type=int, default=25000)

    args = parser.parse_args()
    if (not args.phase and not args.impute) or (args.phase and args.impute):
        exit("ERROR: Exactly one of --phase or --impute must be specified. Exiting...")
    
    tree, leaf_names, leaf_indices = extract_newick_tree_from_smc(args.smc, args.pos)
    samples = read_sample_list(args.samps) if args.samps is not None else None
    
    # Read STR genotypes
    if args.diploid:
        nrepeats_dict, median_allele = read_diploid_str_gts(args.vcf, args.chrom, args.pos, sample_set=samples)
    else:
        nrepeats_dict, median_allele = read_haploid_str_gts(args.vcf, args.chrom, args.pos, sample_set=samples)

    # Ensure that all of the sample nodes are contained within the tree
    for key in nrepeats_dict:
        if key not in leaf_indices:
            exit("ERROR: Sample %s not present in provided tree"%(key))

    if args.diploid:
        # Pair _1 and _2 nodes together
        sample_names = list(set(map(lambda x: x.split("_")[0], nrepeats_dict.keys())))
        pairs        = map(lambda x: (x+"_1", x+"_2"), sample_names)
        pair_data    = pair_gts(nrepeats_dict, leaf_indices, pairs=pairs)
    else:
        # Randomly pair haploid to construct pseudodiploids (node_id_1, node_id_2, num_repeat_a, num_repeat_b)
        pair_data = pair_gts(nrepeats_dict, leaf_indices, pairs=None)

    # Only deal with tree for even number of chromosomes
    if len(tree.leaf_nodes()) % 2 != 0:
        exit("ERROR: Tree contains an odd number of leaves")
    
    # Construct the mutation model
    print("Constructing the mutation model")
    allele_range, max_step = determine_allele_range(args.max_tmrca, args.mu, args.beta, args.pgeom, 0, 0)
    min_allele  = -allele_range - max_step
    max_allele  = allele_range  + max_step 
    mut_model   = OUGeomSTRMutationModel(args.pgeom, args.mu, args.beta, allele_range, max_step = max_step)
    print("Min allele = %d, Max allele = %d"%(min_allele, max_allele))

    # Ensure that the observed genotypes are within the allele range
    if len(pair_data) != 0:
        min_obs_allele = min(min(map(lambda x: x[2], pair_data)), min(map(lambda x: x[3], pair_data)))
        max_obs_allele = max(max(map(lambda x: x[2], pair_data)), max(map(lambda x: x[3], pair_data)))
        if min_obs_allele < min_allele or max_obs_allele > max_allele:
            exit("ERROR: Observed allele not within mutation model's allele range: (%d, %d)"%(min_obs_allele, max_obs_allele))

    # Precompute the transition probabilities
    print("Calculating the transition probabilities")
    optimizer = MATRIX_OPTIMIZER(mut_model.trans_matrix, mut_model.min_n)
    optimizer.precompute_results()

    # Write out paired information for known diploid GTs
    pairs_file = tempfile.mkstemp()[1]
    output     = open(pairs_file, "w") 
    for i in xrange(len(pair_data)):
        output.write("%d\t%d\t%d\t%d\n"%(pair_data[i][0], pair_data[i][1], pair_data[i][2]-min_allele, pair_data[i][3]-min_allele))
    output.close()

    # Write out ids for any leaves not included in the VCF (for potential imputation queries)
    ids_file = tempfile.mkstemp()[1]
    output   = open(ids_file, "w")
    for leaf_name,leaf_id in leaf_indices.items():
        if leaf_name not in nrepeats_dict:
            output.write("%d\n"%(leaf_id))
    output.close()

    # Write out the factor graph file
    graph_file = tempfile.mkstemp()[1]
    write_factor_graph(tree, optimizer, pair_data, min_allele, max_allele, args.gen_per_len, graph_file)

    # Run c++ package, parse results and remove temporary files
    # TO DO: Utilize stderr messages to ensure convergence
    cmd_path = os.path.dirname(os.path.realpath(__file__)) + "/Phaser"
    if args.phase:
        cmd = [cmd_path, "--factor-graph", graph_file, "--pair-file", pairs_file]
    elif args.impute:
        cmd = [cmd_path, "--factor-graph", graph_file, "--id-file", ids_file]
    print("Running message-passing tool")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)#, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
     
    res = stdout.strip()
    rm_cmd = ["rm", "-f", graph_file, pairs_file, ids_file]
    subprocess.call(rm_cmd)
    
    # Assess accuracy, output the statistics and determine the new genotypes associated with each sample
    if args.phase:
        phased_repeat_dict, accuracy_string = process_phasing_result(res, leaf_names, min_allele, min_confidence=args.thresh)
        #print(accuracy_string)
        
        # Construct a new VCF containing the phased alleles
        if args.diploid:
            write_diploid_vcf(args.vcf, args.chrom, args.pos, phased_repeat_dict, median_allele, args.out + "_phased_strs.vcf")
        else:
            write_haploid_vcf(args.vcf, args.chrom, args.pos, phased_repeat_dict, median_allele, args.out + "_phased_strs.vcf")

    # Determine the most probable posterior genotype for each sample
    elif args.impute:
        imputed_repeat_dict,posterior_dict,dosage_dict = process_imputation_result(res, leaf_names, min_allele, max_allele)
        if args.diploid:
            write_diploid_vcf(args.vcf, args.chrom, args.pos, imputed_repeat_dict, median_allele, args.out + "_imputed_strs.vcf", posteriors=posterior_dict, dosages=dosage_dict)
        else:
            write_haploid_vcf(args.vcf, args.chrom, args.pos, imputed_repeat_dict, median_allele, args.out + "_imputed_strs.vcf", posteriors=posterior_dict, dosages=dosage_dict)
        


if __name__ == "__main__":
    main()


