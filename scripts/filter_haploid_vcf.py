import argparse
import collections
import os
import sys

try:
     import numpy
except ImportError:
     exit("This script requires the NumPy python package. Please install using pip or conda")

try:
     import vcf
except ImportError:
     exit("This script requires the PyVCF python package. Please install using pip or conda")

def filter_call(sample, TOTAL_DEPTH, QUAL, FLANK_INDEL_FRAC, STUTTER_FRAC, FILTER_FRAME, REM_NO_SPANNING, best_frame, period):
     if REM_NO_SPANNING and sample['MALLREADS'] is None:
          return 1
     elif sample['DP'] < TOTAL_DEPTH:
          return 2
     elif sample['Q'] < QUAL:
          return 3
     elif 1.0*sample['DFLANKINDEL']/sample['DP'] > FLANK_INDEL_FRAC:
          return 4
     elif 1.0*sample['DSTUTTER']/sample['DP'] > STUTTER_FRAC:
          return 5
     elif FILTER_FRAME and best_frame != int(sample['GB'])%period:
          return 6
     else:
          return 0

def read_sample_set(input_file):
     samples = set()
     data = open(input_file, "r")
     for line in data:
          samples.add(line.strip())
     data.close()
     return samples

def main():
     parser = argparse.ArgumentParser()
     parser.add_argument("--vcf",                  help="Input VCF",                                          type=str,    required=True,                  dest="VCF")
     parser.add_argument("--stats",                help="Write filter statistics to the provided file",       type=str,    required=False,  default=None,  dest="STATS")
     parser.add_argument("--min-call-depth",       help="Omit a sample's call if DP < DEPTH",                 type=int,    required=False,  default=0,     dest="DEPTH")
     parser.add_argument("--min-call-qual",        help="Omit a sample's call if Q  < QUAL",                  type=float,  required=False,  default=0.0,   dest="QUAL")
     parser.add_argument("--max-call-flank-indel", type=float,  required=False,  default=1.0,   dest="FLANK_INDEL_FRAC",
                         help="Omit a sample's call if the fraction of reads with flank indels (FORMAT fields DFLANKINDEL/DP) > FLANK_INDEL_FRAC")
     parser.add_argument("--max-call-stutter",     type=float,  required=False,  default=1.0,   dest="STUTTER_FRAC",
                         help="Omit a sample's call if the fraction of reads with a stutter artifact (FORMAT fields DSTUTTER/DP) > STUTTER_FRAC")
     parser.add_argument("--no-spanning",          action="store_true", required=False, default=False, dest="REM_NO_SPANNING",
                         help="Omit a sample's call if it does not have any spanning reads")
     parser.add_argument("--male-samples",         help="File containing male sample names, one per line",    type=str,    required=False,  default=None,  dest="MALE_SAMPLES")
     parser.add_argument("--female-samples",       help="File containing female sample names, one per line",  type=str,    required=False,  default=None,  dest="FEMALE_SAMPLES")
     parser.add_argument("--min-males",            type=int, required=False,  default=-1, dest="MIN_MALES",
                         help="Remove locus if the number of male samples after filtering < MIN_MALES")
     parser.add_argument("--max-males",            type=int, required=False,  default=-1, dest="MAX_MALES",
                         help="Remove locus if the number of male samples after filtering > MAX_MALES")
     parser.add_argument("--min-females",          type=int, required=False,  default=-1, dest="MIN_FEMALES",
                         help="Remove locus if the number of female samples after filtering < MIN_FEMALES")
     parser.add_argument("--max-females",          type=int, required=False,  default=-1, dest="MAX_FEMALES",
                         help="Remove locus if the number of female samples after filtering > MAX_FEMALES")
     parser.add_argument("--remove-males",         help="Remove calls for male samples",             action="store_true",  required=False,  default=False, dest="REMOVE_MALES")
     parser.add_argument("--remove-females",       help="Remove calls for female samples",           action="store_true",  required=False,  default=False, dest="REMOVE_FEMALES")
     parser.add_argument("--filter-frame",         help="Filter out-of-frame STR alleles",           action="store_true",  required=False,  default=False, dest="FILTER_FRAME")
     parser.add_argument("--min-loc-depth",        type=int,   required=False,  default=0,     dest="MIN_LOC_DEPTH",
                         help="Omit locus if the total depth (DP INFO field) < MIN_LOC_DEPTH")
     parser.add_argument("--max-loc-depth",        type=int,   required=False,  default=100000, dest="MAX_LOC_DEPTH",
                         help="Omit locus if the total depth (DP INFO field) > MAX_LOC_DEPTH")
     parser.add_argument("--max-loc-flank-indel", type=float, required=False, default=1.0, dest="LOC_FLANK_INDEL_FRAC",
                         help="Omit locus if the total fraction of reads with an indel in the flanks (INFO fields DFLANKINDEL/DP) > LOC_FLANK_INDEL_FRAC")
     parser.add_argument("--max-loc-stutter",     type=float, required=False, default=1.0, dest="LOC_STUTTER_FRAC", 
                         help="Omit locus if the total fraction of reads with a stutter artifact (INFO fields DSTUTTER/DP) > LOC_STUTTER")
 
     args = parser.parse_args()
     
     if (args.MIN_MALES != -1 or args.MAX_MALES != -1 or args.REMOVE_MALES) and args.MALE_SAMPLES is None:
          exit("ERROR: --male-samples flag is required when --min-males, --max-males or --remove-males flags are specified")
     if (args.MAX_FEMALES != -1 or args.MAX_FEMALES != -1 or args.REMOVE_FEMALES) and args.FEMALE_SAMPLES is None:
          exit("EROR: --female-samples flag is required when --min-females, --max-females or --remove-females flags are specified")

     male_samples = set()
     if args.MALE_SAMPLES is not None:
          male_samples = read_sample_set(args.MALE_SAMPLES)
     female_samples = set()
     if args.FEMALE_SAMPLES is not None:
          female_samples = read_sample_set(args.FEMALE_SAMPLES)
     if len(male_samples & female_samples) != 0:
          exit("ERROR: Male and female sample lists overlap")

     if args.VCF == "-":
          vcf_reader = vcf.Reader(sys.stdin)
     else:
          vcf_reader = vcf.Reader(filename=args.VCF)
     vcf_writer        = vcf.Writer(sys.stdout, vcf_reader)
     filt_stats        = open(args.STATS if args.STATS is not None else os.devnull, "w")
     total_filt_counts = numpy.zeros(8)
     filter_labels     = ["PASS", "SPANNING", "DEPTH", "QUAL", "FLANK_INDEL", "STUTTER", "OUT_OF_FRAME", "DUE_TO_SEX"]
     filt_stats.write("\t".join(["CHROM", "POS", "FRAC_FILT"] + filter_labels + ["LOCUS_STATUS"]) + "\n")

     for record in vcf_reader:
          samp_fmt = vcf.model.make_calldata_tuple(record.FORMAT.split(':'))
          for fmt in samp_fmt._fields:
               entry_type = vcf_reader.formats[fmt].type
               entry_num  = vcf_reader.formats[fmt].num
               samp_fmt._types.append(entry_type)
               samp_fmt._nums.append(entry_num)
               
          nfields       = len(samp_fmt._fields)
          allele_counts = len(record.alleles)*[0]
          sex_counts    = [0, 0, 0] # Male, female, unknown

          # Determine best frame, if necessary
          best_frame = None
          period     = record.INFO['PERIOD']
          if args.FILTER_FRAME:
               frame_counts = collections.defaultdict(int)
               for sample in record:
                    if sample['GT'] is None or sample['GT'] == './.' or sample['GT'] == '.':
                         continue
                    frame_counts[int(sample['GB'])%period] += 1
               best_frame = sorted(frame_counts.items(), key = lambda x: x[1])[-1][0]

          for sample in record:
               if sample['GT'] is None or sample['GT'] == './.' or sample['GT'] == '.':
                    continue

               filt_index = filter_call(sample, args.DEPTH, args.QUAL, args.FLANK_INDEL_FRAC, args.STUTTER_FRAC,
                                        args.FILTER_FRAME, args.REM_NO_SPANNING, best_frame, period)
               if filt_index == 0:
                    if sample.sample in male_samples:
                         sex_counts[0] += 1
                    elif sample.sample in female_samples:
                         sex_counts[1] += 1
                    else:
                         sex_counts[2] += 1

                    gt_a = int(sample['GT'])
                    if args.REMOVE_MALES and sample.sample in male_samples:
                         continue
                    elif args.REMOVE_FEMALES and sample.sample in female_samples:
                         continue
                    allele_counts[gt_a] += 1

          if sum(allele_counts) == 0:
               continue

          if record.INFO['DP'] < args.MIN_LOC_DEPTH:
               filt_stats.write("%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
                                %(record.CHROM, record.POS, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, "MIN_LOC_DEPTH_%d_%d"%(record.INFO['DP'], args.MIN_LOC_DEPTH)))
               continue

          if record.INFO['DP'] > args.MAX_LOC_DEPTH:
               filt_stats.write("%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
                                %(record.CHROM, record.POS, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, "MAX_LOC_DEPTH_%d_%d"%(record.INFO['DP'], args.MAX_LOC_DEPTH)))
               continue
               

          if 1.0*record.INFO['DSTUTTER']/record.INFO['DP'] > args.LOC_STUTTER_FRAC:
               filt_stats.write("%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
                    %(record.CHROM, record.POS, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, "LOC_STUTTER_%f_%f"%(1.0*record.INFO['DSTUTTER']/record.INFO['DP'], args.LOC_STUTTER_FRAC)))
               continue

          if 1.0*record.INFO['DFLANKINDEL']/record.INFO['DP'] > args.LOC_FLANK_INDEL_FRAC:
               filt_stats.write("%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
                    %(record.CHROM, record.POS, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, "LOC_FLANK_INDEL_%f_%f"%(1.0*record.INFO['DFLANKINDEL']/record.INFO['DP'], args.LOC_FLANK_INDEL_FRAC)))
               continue

          if args.MIN_MALES != -1 and sex_counts[0] < args.MIN_MALES:
               filt_stats.write("%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
                    %(record.CHROM, record.POS, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, "MIN_MALES_%d_%d"%(sex_counts[0], args.MIN_MALES)))
               continue

          if args.MAX_MALES != -1 and sex_counts[0] > args.MAX_MALES:
               filt_stats.write("%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
                    %(record.CHROM, record.POS, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, "MAX_MALES_%d_%d"%(sex_counts[0], args.MAX_MALES)))
               continue

          if args.MIN_FEMALES != -1 and sex_counts[1] < args.MIN_FEMALES:
               filt_stats.write("%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
                    %(record.CHROM, record.POS, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, "MIN_FEMALES_%d_%d"%(sex_counts[1], args.MIN_FEMALES)))
               continue

          if args.MAX_FEMALES != -1 and sex_counts[1] > args.MAX_FEMALES:
               filt_stats.write("%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
                    %(record.CHROM, record.POS, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, "MAX_FEMALES_%d_%d"%(sex_counts[1], args.MAX_FEMALES)))
               continue
          
          # Build allele index mapping
          allele_indices   = {0:0}
          filt_num_alleles = 1
          for i in range(1, len(allele_counts)):
              if allele_counts[i] != 0:
                  allele_indices[i] = filt_num_alleles
                  filt_num_alleles += 1
               
          new_samples       = []
          filt_counts       = [0, 0, 0, 0, 0, 0, 0, 0]
          num_kept          = 0
          total_dp          = 0
          total_dstutter    = 0
          total_dflankindel = 0
          for sample in record:
               if sample['GT'] is None or sample['GT'] == './.' or sample['GT'] == '.':
                    new_samples.append(sample)
                    continue

               remove_due_to_sex = (args.REMOVE_MALES and sample.sample in male_samples) or (args.REMOVE_FEMALES and sample.sample in female_samples)
               filt_index = filter_call(sample, args.DEPTH, args.QUAL, args.FLANK_INDEL_FRAC, args.STUTTER_FRAC,
                                        args.FILTER_FRAME, args.REM_NO_SPANNING, best_frame, period)
               if remove_due_to_sex:
                    filt_index = 7

               if filt_index != 0:
                    sampdat = []
                    for i in range(len(samp_fmt._fields)):
                        key = samp_fmt._fields[i]
                        if key == "GT":
                            sampdat.append(".")
                        else:
                            sampdat.append(None)
               else:
                    num_kept += 1
                    gt_a      = int(sample['GT'])
                    new_gt    = "%d"%(allele_indices[gt_a])
                    sampdat   = []
                    for i in range(len(samp_fmt._fields)):
                        key = samp_fmt._fields[i]
                        if key == "GT":
                            sampdat.append(new_gt)
                        else:
                            sampdat.append(sample[key])
                    total_dp          += sample['DP']
                    total_dstutter    += sample['DSTUTTER']
                    total_dflankindel += sample['DFLANKINDEL']
                    
               filt_counts[filt_index] += 1
               call = vcf.model._Call(record, sample.sample, samp_fmt(*sampdat))
               new_samples.append(call)
          record.samples = new_samples

          # Fix set of alleles
          new_alleles = [record.alleles[0]]
          for i in range(1, len(record.alleles)):
              if allele_counts[i] != 0:
                  new_alleles.append(record.alleles[i])

          record.alleles = new_alleles
          record.ALT     = new_alleles[1:]
          if len(record.ALT) == 0:
               record.ALT = ["."]

          # Recompute and set INFO fields
          if 'BPDIFFS' in record.INFO:
              if len(new_alleles) == 1:
                  record.INFO.pop("BPDIFFS", None)
              else:
                  record.INFO['BPDIFFS'] = list(map(lambda x: len(x) - len(new_alleles[0]), new_alleles[1:]))
          record.INFO['REFAC'] = allele_counts[0]
          if 'AC' in record.INFO:
              if len(new_alleles) == 1:
                  record.INFO.pop("AC", None)
              else:
                  record.INFO['AC'] = list(filter(lambda x: x != 0, allele_counts[1:]))
          if 'AN' in record.INFO:
               record.INFO['AN'] = sum(allele_counts)

          # Recompute DP and other read depth info fields
          record.INFO['DP']          = total_dp
          record.INFO['DSTUTTER']    = total_dstutter
          record.INFO['DFLANKINDEL'] = total_dflankindel
          record.INFO['NFILT']      += sum(filt_counts[1:])
          vcf_writer.write_record(record)
          total_filt_counts += numpy.array(filt_counts)
          filt_stats.write("%s\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"%(record.CHROM, record.POS, 1.0 - 1.0*filt_counts[0]/sum(filt_counts),
                                                                               filt_counts[0], filt_counts[1], filt_counts[2], filt_counts[3],
                                                                               filt_counts[4], filt_counts[5], filt_counts[6], filt_counts[7], "PASS"))
     filt_stats.write("\t".join([".", ".", "."] + map(lambda x: str(int(x)), total_filt_counts) + ["."]) + "\n")
     filt_stats.close()

if __name__ == "__main__":
    main()
