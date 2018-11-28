import argparse
import collections
import os
import sys

# Special import required as python2/python3 have different StringIO packages
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

try:
    import vcf
except ImportError:
    exit("This script requires the PyVCF python package. Please install using pip or conda")

# Simple class that writes out VCF records that only contain entries for a subset of samples in the original VCF
# NOTE: The INFO fields are not modified by this class, so any desired INFO changes must occur prior to calling write_record()
class SampleFilterVCFWriter:
    def __init__(self, samples_to_keep, output_stream, vcf_reader):
        self.output_stream_ = output_stream

        good_samples               = set(samples_to_keep)
        self.valid_sample_indices_ = []
        for index,sample in enumerate(vcf_reader.samples):
            if sample in good_samples:
                self.valid_sample_indices_.append(index)

        # We pass the VCF writer's output to the string stream, which is then read and modified 
        # accordingly before ouputting the data for the subset of samples
        self.string_stream_ = StringIO()
        self.vcf_writer_    = vcf.Writer(self.string_stream_, vcf_reader)

        # Output header lines after modifying the included samples
        for line in self.string_stream_.getvalue().split("\n"):
            if len(line) == 0:
                continue
            if line.startswith("##"):
                output_stream.write(line + "\n")
            else:
                tokens = line.strip().split("\t")
                output_stream.write("\t".join([tokens[i] for i in [0, 1, 2, 3, 4, 5, 6, 7, 8] + list(map(lambda x: x + 9, self.valid_sample_indices_))]) + "\n")

        # Clear and reset the string stream                    
        self.string_stream_.truncate(0)
        self.string_stream_.seek(0)

    def write_record(self, record):
        # Modify the record's samples to contain only our subset of interest
        old_samples    = record.samples
        new_samples    = [old_samples[i] for i in self.valid_sample_indices_]
        record.samples = new_samples

        # Write out the modified record
        self.vcf_writer_.write_record(record)
        self.output_stream_.write(self.string_stream_.getvalue())

        # Clear and reset the string stream
        self.string_stream_.truncate(0)
        self.string_stream_.seek(0)

        # Restore the record back to its original form
        record.samples = old_samples
          
    def close(self):
        self.string_stream_.close()


def filter_call(sample, filters):
     if sample['DP'] < filters.DEPTH:
          return "Depth"
     elif sample['Q'] < filters.QUAL:
          return "Quality"
     else:
          d_1, d_2 = map(float, sample['PDP'].split('|'))
          if d_1 == 0 or d_2 == 0:
               return "Allele depth"
          if min(d_1, d_2) < filters.ALLELE_DEPTH:
               return "Allele depth"
          elif min(1.0*d_1/d_2, 1.0*d_2/d_1) < filters.ALLELE_RATIO:
               return "Allele ratio"
          elif filters.FLANK_INDEL_FRAC < 1 and 1.0*sample['DFLANKINDEL']/sample['DP'] > filters.FLANK_INDEL_FRAC:
               return "Flank indels"
          elif filters.STUTTER_FRAC < 1 and 1.0*sample['DSTUTTER']/sample['DP'] > filters.STUTTER_FRAC:
               return "Stutter fraction"
          elif filters.ALLELE_BIAS > -100 and sample['AB'] < filters.ALLELE_BIAS:
               return "Allele bias"
          elif filters.STRAND_BIAS > -100 and sample['FS'] < filters.STRAND_BIAS:
               return "Strand bias"

          if filters.SPAN_DEPTH > 0:
               if sample["MALLREADS"] is None:
                    return "Spanning depth"
               gb_a, gb_b = map(int, sample["GB"].split("|"))
               span_depth_dict = dict(map(lambda x: map(int, x.split("|")), sample["MALLREADS"].split(";")))
               dp_a = 0 if gb_a not in span_depth_dict else span_depth_dict[gb_a]
               dp_b = 0 if gb_b not in span_depth_dict else span_depth_dict[gb_b]
               if min(dp_a, dp_b) < filters.SPAN_DEPTH:
                    return "Spanning depth"
          return None

# Read in a text file of sample names, one per line
def read_sample_list(path):
     if not os.path.exists(path):
          exit("ERROR: File argument to --samples-to-keep does not exist : %s"%(path))

     samples_to_keep = set()
     data = open(path, "r")
     for line in data:
          sample = line.strip()
          if len(sample) > 0:
               samples_to_keep.add(sample)
     data.close()

     return samples_to_keep

def gt_is_diploid(gt):
    return (gt.count("/") + gt.count("|") == 1)

def main():
     # Dictionary of help messages for most of the command line arguments
     help_dict = {
         "--min-call-allele-depth"   : "Omit a sample's call if the minimum allele depth (smallest value in PDP FORMAT field) < ALLELE_DEPTH",
         "--min-call-depth-ratio"    : "Omit a sample's call if the minimum ratio of allele depths (values in PDP FORMAT field) < ALLELE_RATIO",
         "--max-call-flank-indel"    : "Omit a sample's call if the fraction of reads with flank indels (FORMAT fields DFLANKINDEL/DP) > FLANK_INDEL_FRAC",
         "--max-call-stutter"        : "Omit a sample's call if the fraction of reads with a stutter artifact (FORMAT fields DSTUTTER/DP) > STUTTER_FRAC",
         "--min-call-spanning-depth" : "Omit a sample's call if the minimum number of spanning reads supporting an allele < SPAN_DEPTH",
         "--min-loc-depth"           : "Omit locus if the total depth (DP INFO field) < MIN_LOC_DEPTH",
         "--max-loc-depth"           : "Omit locus if the total depth (DP INFO field) > MAX_LOC_DEPTH",
         "--max-loc-flank-indel"     : "Omit locus if the total fraction of reads with an indel in the flanks (INFO fields DFLANKINDEL/DP) > LOC_FLANK_INDEL_FRAC",
         "--max-loc-stutter"         : "Omit locus if the total fraction of reads with a stutter artifact (INFO fields DSTUTTER/DP) > LOC_STUTTER",
         "--min-loc-calls"           : "Omit locus if the number of valid samples after filtering < MIN_CALLS",
         "--samples-to-keep"         : "File of samples to keep in the VCF, one-per-line. By default, all samples are kept",
         "--keep-all-alleles"        : "Don't remove any ALT alleles, even if no unfiltered samples have the corresponding GT. (Default = Remove unless VCF contains likelihoods",
     }

     description = """Apply the provided set of sample-level (i.e. call-level) and locus-level filters to the input VCF file.
     Filtered samples are masked as missing while filtered loci are entirely removed from the VCF.
     All INFO fields are updated accordingly to refer to the unfiltered sample set.
     The resulting filtered VCF is output to stdout
     """

     parser = argparse.ArgumentParser(description=description)
     parser.add_argument("--vcf",                     type=str,   required=True,  dest="VCF",                                      help="Input VCF to filter (- for stdin)")
     parser.add_argument("--min-call-depth",          type=int,   required=False, dest="DEPTH",                default=0,          help="Omit a sample's call if DP < DEPTH")
     parser.add_argument("--min-call-qual",           type=float, required=False, dest="QUAL",                 default=0.0,        help="Omit a sample's call if Q < QUAL")
     parser.add_argument("--min-call-allele-depth",   type=float, required=False, dest="ALLELE_DEPTH",         default=0.0,        help=help_dict["--min-call-allele-depth"])
     parser.add_argument("--min-call-depth-ratio",    type=float, required=False, dest="ALLELE_RATIO",         default=0.0,        help=help_dict["--min-call-depth-ratio"])
     parser.add_argument("--max-call-flank-indel",    type=float, required=False, dest="FLANK_INDEL_FRAC",     default=1.0,        help=help_dict["--max-call-flank-indel"])
     parser.add_argument("--max-call-stutter",        type=float, required=False, dest="STUTTER_FRAC",         default=1.0,        help=help_dict["--max-call-stutter"])
     parser.add_argument("--min-call-spanning-depth", type=int,   required=False, dest="SPAN_DEPTH",           default=0.0,        help=help_dict["--min-call-spanning-depth"])
     parser.add_argument("--min-call-allele-bias",    type=float, required=False, dest="ALLELE_BIAS",          default=-100.0,     help="Omit a sample's call if AB < ALLELE_BIAS")
     parser.add_argument("--min-call-strand-bias",    type=float, required=False, dest="STRAND_BIAS",          default=-100.0,     help="Omit a sample's call if FS < STRAND_BIAS")
     parser.add_argument("--min-loc-depth",           type=int,   required=False, dest="MIN_LOC_DEPTH",        default=0,          help=help_dict["--min-loc-depth"])
     parser.add_argument("--max-loc-depth",           type=int,   required=False, dest="MAX_LOC_DEPTH",        default=1000000000, help=help_dict["--max-loc-depth"])
     parser.add_argument("--max-loc-flank-indel",     type=float, required=False, dest="LOC_FLANK_INDEL_FRAC", default=1.0,        help=help_dict["--max-loc-flank-indel"])
     parser.add_argument("--max-loc-stutter",         type=float, required=False, dest="LOC_STUTTER",          default=1.0,        help=help_dict["--max-loc-stutter"])
     parser.add_argument("--min-loc-calls",           type=int,   required=False, dest="MIN_CALLS",            default=0,          help=help_dict["--min-loc-calls"])
     parser.add_argument("--samples-to-keep",         type=str,   required=False, dest="SAMPLES_TO_KEEP",      default=None,       help=help_dict["--samples-to-keep"])
     parser.add_argument("--keep-all-alleles",        action="store_true",        dest="KEEP_ALL_ALLELES",     default=False,      help=help_dict["--keep-all-alleles"])
     args = parser.parse_args()
     if args.VCF == "-":
          vcf_reader = vcf.Reader(sys.stdin)
     else:
          vcf_reader = vcf.Reader(filename=args.VCF)

     total_counts  = collections.defaultdict(int)
     filter_counts = {}
     for sample in vcf_reader.samples:
          filter_counts[sample] = collections.defaultdict(int)

     # Read a list of samples to keep and ensure they're all in the VCF
     samples_to_keep = None if args.SAMPLES_TO_KEEP is None else read_sample_list(args.SAMPLES_TO_KEEP)
     if samples_to_keep is not None:
          for sample in samples_to_keep:
               if sample not in filter_counts:
                    exit("ERROR: Sample %s was present in the --samples-to-keep file, but no such sample is present in the VCF"%(sample))

     # Open the appropriate VCF writer
     if samples_to_keep is None:
          vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
     else:
          vcf_writer = SampleFilterVCFWriter(samples_to_keep, sys.stdout, vcf_reader)

     for record in vcf_reader:
          if record.INFO['DP'] < args.MIN_LOC_DEPTH:
               continue
          if record.INFO['DP'] > args.MAX_LOC_DEPTH:
               continue
          if args.LOC_FLANK_INDEL_FRAC < 1 and 1.0*record.INFO["DFLANKINDEL"]/record.INFO["DP"] > args.LOC_FLANK_INDEL_FRAC:
               continue
          if args.LOC_STUTTER < 1 and 1.0*record.INFO["DSTUTTER"]/record.INFO["DP"] > args.LOC_STUTTER:
               continue

          samp_fmt = vcf.model.make_calldata_tuple(record.FORMAT.split(':'))
          for fmt in samp_fmt._fields:
               entry_type = vcf_reader.formats[fmt].type
               entry_num  = vcf_reader.formats[fmt].num
               samp_fmt._types.append(entry_type)
               samp_fmt._nums.append(entry_num)

          # Determine if we can remove alleles with 0 counts based on whether any genotype-likelhood associated
          # format fields are present
          fmt_tokens         = record.FORMAT.split(":")
          can_remove_alleles = (not args.KEEP_ALL_ALLELES and not ("GL" in fmt_tokens or "PL" in fmt_tokens or "PHASEDGL" in fmt_tokens))
               
          nfields       = len(samp_fmt._fields)
          num_filt      = 0
          allele_counts = len(record.alleles)*[0]

          for sample in record:
               if sample['GT'] is None or sample['GT'] == "./." or sample['GT'] == ".":
                    continue
               if not gt_is_diploid(sample['GT']):
                   msg = "ERROR: Genotype '%s' for sample %s at %s:%d is not diploid as required by this script"%(sample['GT'], sample.sample, record.CHROM, record.POS)
                   msg += "\n\t Please use filter_haploid_vcf.py to filter calls for haploid chromosomes"
                   exit(msg)

               if samples_to_keep is not None and sample.sample not in samples_to_keep:
                    filter_reason = "NOT_IN_SAMPLES_TO_KEEP"
               else:
                    filter_reason = filter_call(sample, args)

               if filter_reason is None:
                    gt_a, gt_b = map(int, sample['GT'].split('|'))
                    allele_counts[gt_a] += 1
                    allele_counts[gt_b] += 1
               else:
                    filter_counts[sample.sample][filter_reason] += 1
                    total_counts[filter_reason] += 1

          # Build allele index mapping
          allele_indices   = {0:0}
          filt_num_alleles = 1
          for i in range(1, len(allele_counts)):
              if allele_counts[i] != 0 or not can_remove_alleles:
                  allele_indices[i] = filt_num_alleles
                  filt_num_alleles += 1

          new_samples       = []
          num_filt          = 0
          num_kept          = 0
          total_dp          = 0
          total_dstutter    = 0
          total_dflankindel = 0
          for sample in record:
               if sample['GT'] is None or sample['GT'] == "./." or sample['GT'] == ".":
                    new_samples.append(sample)
                    continue

               if samples_to_keep is not None and sample.sample not in samples_to_keep:
                    filter_reason = "NOT_IN_SAMPLES_TO_KEEP"
               else:
                    filter_reason = filter_call(sample, args)

               if filter_reason is not None:
                    num_filt += 1
                    sampdat   = []
                    for i in range(len(samp_fmt._fields)):
                         key = samp_fmt._fields[i]
                         if key == "GT":
                              sampdat.append("./.")
                         else:
                              if key == "FILTER":
                                   sampdat.append(filter_reason.replace(" ", "_").upper())
                              else:
                                   sampdat.append(None)
               else:
                    num_kept    += 1
                    gt_a, gt_b   = map(int, sample['GT'].split('|'))
                    new_gt       = "%d|%d"%(allele_indices[gt_a], allele_indices[gt_b])
                    sampdat      = []
                    for i in range(len(samp_fmt._fields)):
                        key = samp_fmt._fields[i]
                        if key == "GT":
                            sampdat.append(new_gt)
                        else:
                            sampdat.append(sample[key])

                    total_dp          += sample['DP']
                    total_dstutter    += sample['DSTUTTER']
                    total_dflankindel += sample['DFLANKINDEL']

               call = vcf.model._Call(record, sample.sample, samp_fmt(*sampdat))
               new_samples.append(call)
          record.samples = new_samples

          # Don't output the filtered record if we don't meet the minimum number of calls
          if num_kept < args.MIN_CALLS:
               continue
          
          # Fix set of alleles
          new_alleles = [record.alleles[0]]
          for i in range(1, len(record.alleles)):
              if allele_counts[i] != 0 or not can_remove_alleles:
                  new_alleles.append(record.alleles[i])

          record.alleles = new_alleles
          record.ALT     = new_alleles[1:]
          if len(record.ALT) == 0:
              record.ALT = ["."]
          if 'NFILT' in record.INFO:
               record.INFO['NFILT'] += num_filt

          record.INFO['DP']          = total_dp
          record.INFO['DSTUTTER']    = total_dstutter
          record.INFO['DFLANKINDEL'] = total_dflankindel

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
                   if not can_remove_alleles:
                        record.INFO['AC'] = allele_counts[1:]
                   else:
                        record.INFO['AC'] = list(filter(lambda x: x != 0, allele_counts[1:]))
          if 'AN' in record.INFO:
               record.INFO['AN'] = sum(allele_counts)
                  
          vcf_writer.write_record(record)

     if samples_to_keep is not None:
          vcf_writer.close()

if __name__ == "__main__":
    main()
