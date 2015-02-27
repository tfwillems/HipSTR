import argparse
import collections
import dendropy
import gzip
import math
import random
import subprocess
import sys
import tempfile
import vcf

PHASE_CMD_PATH="./str-imputer"

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

# Read the STR genotypes from the provided VCF file, normalize their lengths relative
# to the median and return a dictionary mapping sample names to the number of normalized repeats
def read_haploid_str_gts(vcf_file):
    vcf_reader   = vcf.Reader(filename=vcf_file)
    record       = vcf_reader.next()
    gb_counts    = collections.defaultdict(int)
    for sample in record:
        gb_counts[len(record.alleles[int(sample['GT'])])] += 1
    median_allele = compute_median(gb_counts.keys(), gb_counts.values())    
    motif_len     = len(record.INFO['MOTIF'])
    nrepeats_dict = {}
    for sample in record:
        nrepeats_dict[sample.sample.replace("n", "")] = (len(record.alleles[int(sample['GT'])])-median_allele)/motif_len
    return nrepeats_dict, median_allele
        
def pair_gts(nrepeats_dict, pairs=None):
    if pairs is None:
        sample_order = nrepeats_dict.keys()
        random.shuffle(sample_order)
        pairs = map(lambda x: (sample_order[x], sample_order[x+1]), xrange(0, len(sample_order)/2*2, 2))
    pair_data = []
    for i in xrange(len(pairs)):
        pair_data.append((pairs[i][0], pairs[i][1], nrepeats_dict[pairs[i][0]], nrepeats_dict[pairs[i][1]]))
    return pair_data

def extract_newick_tree_from_smc(input_file, position):
    data        = gzip.open(input_file, "rb") if input_file.endswith(".gz") else  open(input_file, "r")
    name_tokens = data.readline().strip().split()
    if name_tokens[0] != "NAMES":
        exit("ERROR: First line of SMC file must contain NAMES")
    name_mappings = dict(zip(map(str, range(len(name_tokens)-1)), name_tokens[1:]))

    for line in data:
        tokens = line.strip().split()
        if tokens[0] != "TREE":
            continue
        if int(math.floor(float(tokens[1]))) <= position and int(math.ceil(float(tokens[2]))) >= position:
            data.close()

            # Correct node names
            tree = dendropy.Tree(stream=StringIO(tokens[3]), schema="newick")
            for taxon in tree.taxon_set:
                taxon.label = name_mappings[taxon.label].replace("n", "")
            return tree
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
    print("Iterating through paired nodes")
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

def process_phasing_result(phasing_result, min_allele, min_confidence = 0.5):
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
        if gt_a == gt_b:
            gt_dict["n"+str(id_1)] = gt_a + min_allele
            gt_dict["n"+str(id_2)] = gt_a + min_allele
            num_homoz += 1
            continue
        conf_1, conf_2 = map(float, tokens[4:6])
        if conf_1 > 0.5:
            if conf_1 > min_confidence:
                num_correct  += 1
                num_phased   += 1
                gt_dict["n"+str(id_1)] = gt_a + min_allele
                gt_dict["n"+str(id_2)] = gt_b + min_allele
            else:
                num_skip += 1
        else:
            if conf_2 > min_confidence:
                num_phased   += 1
                gt_dict["n"+str(id_1)] = gt_b + min_allele
                gt_dict["n"+str(id_2)] = gt_a + min_allele
                errors.append("%d_%d"%(gt_a, gt_b))
            else:
                num_skip += 1
    print("ACCURACY\t%d\t%d\t%d\t%d\t%s"%(num_skip, num_homoz, num_correct, num_phased, ",".join(errors)))
    return gt_dict

def write_vcf(input_vcf_file, nrepeats_dict, median_allele, output_vcf_file):
    output = open(output_vcf_file, "w")
    output.write("##fileformat=VCFv4.2\n")
    output.write("##INFO=<ID=MOTIF,Number=1,Type=String,Description=\"Canonical repeat motif\">\n")
    output.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    output.write("##FORMAT=<ID=GB,Number=1,Type=String,Description=\"Genotype given in bp difference from reference\">\n")
    output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(map(lambda x: str(x), sorted(nrepeats_dict.keys()))) + "\n")

    vcf_reader = vcf.Reader(filename=input_vcf_file)
    record     = vcf_reader.next()

    gt_indexes = {}
    gt_indexes[len(record.REF)] = 0
    for i in xrange(len(record.ALT)):
        gt_indexes[len(str(record.ALT[i]))] = i+1

    tokens = [record.CHROM, record.POS, ".", record.REF, ",".join(map(str, record.ALT)), ".", "PASS", "MOTIF=AT", "GT:GB"]
    for sample,nreps in sorted(nrepeats_dict.items()):
        bp_len = nreps*2 + median_allele
        tokens.append("%d:%d"%(gt_indexes[bp_len], bp_len-len(record.REF)))
    output.write("\t".join(map(str, tokens)) + "\n")
    output.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--smc",   required=True,  dest="smc",   type=str,   help="SMC file for the region of interest")
    parser.add_argument("--pos",   required=True,  dest="pos",   type=int,   help="Position of STR in region")
    parser.add_argument("--mu",    required=True,  dest="mu",    type=float, help="Mutation rate for mutation model")
    parser.add_argument("--beta",  required=True,  dest="beta",  type=float, help="Length constraint for mutation model")
    parser.add_argument("--pgeom", required=True,  dest="pgeom", type=float, help="Geometric parameter for mutation model")
    parser.add_argument("--out",   required=True,  dest="out",   type=str,   help="Output path for phased VCF")
    parser.add_argument("--vcf",   required=True,  dest="vcf",   type=str,   help="Input path for VCF containing haploid STR calls")
    
    # Scaling factor from edge length to # of generations
    parser.add_argument("--gen_per_len", required=False, dest="gen_per_len", type=float, default=1.0)

    # Maximum TMRCA in generations
    parser.add_argument("--max_tmrca", required=False, dest="max_tmrca", type=int, default=25000)

    args = parser.parse_args()
    tree = extract_newick_tree_from_smc(args.smc, args.pos)

    # Read haploid STR genotypes
    nrepeats_dict, median_allele = read_haploid_str_gts(args.vcf)

    # Randomly pair haploid to construct pseudodiploids (id_1, id_2, num_repeat_a, num_repeat_b)
    pair_data = pair_gts(nrepeats_dict, pairs=None)
    
    # Only deal with tree for diploid individuals
    if len(tree.leaf_nodes()) % 2 != 0:
        exit("ERROR: Tree contains an odd number of leaves")

    # Ensure that all of the pairing nodes are contained within the tree
    for i in xrange(len(pair_data)):
        if not tree.taxon_set.has_taxon(label=pair_data[i][0]):
            exit("ERROR: Sample %s not present in provided tree"%(pair_data[i][0]))
        if not tree.taxon_set.has_taxon(label=pair_data[i][1]):
            exit("ERROR: Sample %s not present in provided tree"%(pair_data[i][1]))

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

    # Write out fixed paired information
    pairs_file = tempfile.mkstemp()[1]
    output     = open(pairs_file, "w") 
    for i in xrange(len(pair_data)):
        output.write("%s\t%s\t%d\t%d\n"%(pair_data[i][0], pair_data[i][1], pair_data[i][2]-min_allele, pair_data[i][3]-min_allele))
    output.close()

    # Write out the factor graph file
    graph_file = tempfile.mkstemp()[1]
    write_factor_graph(tree, optimizer, pair_data, min_allele, max_allele, args.gen_per_len, graph_file)

    # Run c++ package, parse results and remove temporary files
    phase_cmd = [PHASE_CMD_PATH, graph_file, pairs_file]
    proc = subprocess.Popen(phase_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(proc.stderr.read().strip())
    res  = proc.stdout.read().strip()
    rm_cmd = ["rm", "-f", graph_file, pairs_file]
    subprocess.call(rm_cmd)

    # Assess accuracy and determine the new genotypes associated with each sample
    phased_repeat_dict = process_phasing_result(res, min_allele, min_confidence = 0.9)

    # Construct a new VCF containing the phased alleles
    write_vcf(args.vcf, phased_repeat_dict, median_allele, args.out)

if __name__ == "__main__":
    main()


