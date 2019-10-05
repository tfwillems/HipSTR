import argparse
import collections
import os
import sys

# Catch the most common failed dependency
try:
    import vcf
except ImportError:
    exit("This script requires the PyVCF python package. Please install using pip or conda")

# Special import required as python2/python3 have different StringIO packages
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

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
        # accordingly before outputting the data for the subset of samples
        self.string_stream_ = StringIO()
        self.vcf_writer_    = vcf.Writer(self.string_stream_, vcf_reader)

        # Output header lines after modifying the included samples
        for line in self.string_stream_.getvalue().split("\n"):
            if len(line) == 0:
                continue
            if line.startswith("##"):
                output_stream.write(line + "\n")
            else:
                tokens  = line.strip().split("\t")
                indexes = [0, 1, 2, 3, 4, 5, 6, 7, 8] + list(map(lambda x: x + 9, self.valid_sample_indices_))
                output_stream.write("\t".join([tokens[i] for i in indexes]) + "\n")

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


def filter_call(sample, filters, MIN_HAP_QUAL, MIN_PHAP_QUAL):
    if sample['DP'] < filters.DEPTH:
        return "DEPTH"
    elif sample['Q'] < filters.QUAL:
        return "QUALITY"
    elif MIN_HAP_QUAL > 0 and sample['HQ'] < MIN_HAP_QUAL:
        return "HAPLOTYPE_QUALITY"
    elif MIN_PHAP_QUAL > 0 and sample['PHQ'] < MIN_PHAP_QUAL:
        return "PHASED_HAPLOTYPE_QUALITY"
    else:
        d_1, d_2 = map(float, sample['PDP'].split('|'))
        if d_1 == 0 or d_2 == 0:
            return "ALLELE_DEPTH"
        if min(d_1, d_2) < filters.ALLELE_DEPTH:
            return "ALLELE_DEPTH"
        elif min(1.0*d_1/d_2, 1.0*d_2/d_1) < filters.ALLELE_RATIO:
            return "ALLELE_RATIO"
        elif filters.FLANK_INDEL_FRAC < 1 and 1.0*sample['DFLANKINDEL']/sample['DP'] > filters.FLANK_INDEL_FRAC:
            return "FLANK_INDELS"
        elif filters.STUTTER_FRAC < 1 and 1.0*sample['DSTUTTER']/sample['DP'] > filters.STUTTER_FRAC:
            return "STUTTER_FRACTION"
        elif filters.ALLELE_BIAS > -100 and sample['AB'] < filters.ALLELE_BIAS:
            return "ALLELE_BIAS"
        elif filters.STRAND_BIAS > -100 and sample['FS'] < filters.STRAND_BIAS:
            return "STRAND_BIAS"

        if filters.SPAN_DEPTH > 0:
            if sample["MALLREADS"] is None:
                return "SPANNING_DEPTH"
            gb_a, gb_b = map(int, sample["GB"].split("|"))
            span_depth_dict = dict(map(lambda x: map(int, x.split("|")), sample["MALLREADS"].split(";")))
            dp_a = 0 if gb_a not in span_depth_dict else span_depth_dict[gb_a]
            dp_b = 0 if gb_b not in span_depth_dict else span_depth_dict[gb_b]
            if min(dp_a, dp_b) < filters.SPAN_DEPTH:
                return "SPANNING_DEPTH"
        
    return None

# Read in a text file of sample names, one per line
def read_sample_list(path):
    if not os.path.exists(path):
        exit("ERROR: File does not exist : %s"%(path))

    samples = set()
    data    = open(path, "r")
    for line in data:
        sample = line.strip()
        if len(sample) > 0:
            samples.add(sample)
    data.close()

    return samples

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
        "--samples-to-remove"       : "File of samples to remove from the VCF, one-per-line. By default, no samples are removed",
        "--keep-all-alleles"        : "Don't remove any ALT alleles, even if no unfiltered samples have the corresponding GT. (Default = Remove unless VCF contains likelihoods",
    }

    description = """Apply the provided set of sample-level (i.e. call-level) and locus-level filters to the input VCF file.
    STR loci that fail locus-level filters are entirely removed from the VCF. For the remaining STRs, filtered samples are masked as missing.
    All INFO fields are updated accordingly to refer to the unfiltered sample set. The resulting filtered VCF is output to stdout
    """

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--vcf",                     type=str,   required=True,  dest="VCF",                                      help="Input VCF to filter (- for stdin)")
    parser.add_argument("--min-call-depth",          type=int,   required=False, dest="DEPTH",                default=0,          help="Omit a sample's call if DP < DEPTH")
    parser.add_argument("--min-call-qual",           type=float, required=False, dest="QUAL",                 default=0.0,        help="Omit a sample's call if Q < QUAL")
    parser.add_argument("--min-call-hap-qual",       type=float, required=False, dest="HAP_QUAL",             default=0.0,        help="Omit a sample's call if HQ < HAP_QUAL")
    parser.add_argument("--min-call-phap-qual",      type=float, required=False, dest="PHAP_QUAL",            default=0.0,        help="Omit a sample's call if PHQ < PHAP_QUAL")
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
    parser.add_argument("--samples-to-remove",       type=str,   required=False, dest="SAMPLES_TO_REMOVE",    default=None,       help=help_dict["--samples-to-remove"])
    parser.add_argument("--keep-all-alleles",        action="store_true",        dest="KEEP_ALL_ALLELES",     default=False,      help=help_dict["--keep-all-alleles"])
    args          = parser.parse_args()
    vcf_reader    = vcf.Reader(sys.stdin) if args.VCF == "-" else vcf.Reader(filename=args.VCF)
    total_counts  = collections.defaultdict(int)
    filter_counts = {}
    for sample in vcf_reader.samples:
        filter_counts[sample] = collections.defaultdict(int)

    samples_to_keep = None
    if args.SAMPLES_TO_REMOVE is not None and args.SAMPLES_TO_KEEP is not None:
        exit("ERROR: You cannot specify both the --samples-to-keep and --samples-to-remove options")
    elif args.SAMPLES_TO_REMOVE is not None:
        # Read list of samples to remove and ensure they're all in the VCF
        samples_to_remove = set(read_sample_list(args.SAMPLES_TO_REMOVE))
        for sample in samples_to_remove:
            if sample not in filter_counts:
                exit("ERROR: Sample %s was present in the --samples-to-remove file, but no such sample is present in the VCF"%(sample))
        samples_to_keep = list(filter(lambda x: x not in samples_to_remove, vcf_reader.samples))
    elif args.SAMPLES_TO_KEEP is not None:
        # Read a list of samples to keep and ensure they're all in the VCF
        samples_to_keep = read_sample_list(args.SAMPLES_TO_KEEP)
        for sample in samples_to_keep:
            if sample not in filter_counts:
                exit("ERROR: Sample %s was present in the --samples-to-keep file, but no such sample is present in the VCF"%(sample))

    # Open the appropriate VCF writer
    vcf_writer = vcf.Writer(sys.stdout, vcf_reader) if samples_to_keep is None else SampleFilterVCFWriter(samples_to_keep, sys.stdout, vcf_reader)

    # Iterate through each VCF record
    for record in vcf_reader:
        # Skip any loci that fail the locus-level filters
        if record.INFO['DP'] < args.MIN_LOC_DEPTH:
            continue
        if record.INFO['DP'] > args.MAX_LOC_DEPTH:
            continue
        if args.LOC_FLANK_INDEL_FRAC < 1 and 1.0*record.INFO["DFLANKINDEL"]/record.INFO["DP"] > args.LOC_FLANK_INDEL_FRAC:
            continue
        if args.LOC_STUTTER < 1 and 1.0*record.INFO["DSTUTTER"]/record.INFO["DP"] > args.LOC_STUTTER:
            continue

        # Determine if we can remove alleles with 0 counts based on whether any genotype-likelhood associated format fields are present
        fmt_tokens         = record.FORMAT.split(":")
        can_remove_alleles = (not args.KEEP_ALL_ALLELES and not ("GL" in fmt_tokens or "PL" in fmt_tokens or "PHASEDGL" in fmt_tokens))

        allele_counts = len(record.alleles)*[0]
        proc_lflanks  = "LFGT" in fmt_tokens
        lflank_counts = [0] if not proc_lflanks else len(record.INFO["LFLANKS"])*[0]
        proc_rflanks  = "RFGT" in fmt_tokens
        rflank_counts = [0] if not proc_rflanks else len(record.INFO["RFLANKS"])*[0]

        min_hap_qual  = (0 if "HQ" not in fmt_tokens else args.HAP_QUAL)
        min_phap_qual = (0 if "PHQ" not in fmt_tokens else args.PHAP_QUAL)

        # Process each sample and determine whether it should be included in the filtered VCF
        sample_filt_dict = {}
        for sample in record:
            if sample['GT'] is None or sample['GT'] == "./." or sample['GT'] == ".":
                continue
            if not gt_is_diploid(sample['GT']):
                msg = "ERROR: Genotype '%s' for sample %s at %s:%d is not diploid as required by this script"%(sample['GT'], sample.sample, record.CHROM, record.POS)
                msg += "\n\t Please use filter_haploid_vcf.py to filter calls for haploid chromosomes"
                exit(msg)
            
            filter_reason = "NOT_IN_SAMPLES_TO_KEEP" if (samples_to_keep is not None and sample.sample not in samples_to_keep) else filter_call(sample, args, 
                                                                                                                                                min_hap_qual, min_phap_qual)
            sample_filt_dict[sample.sample] = filter_reason
            if filter_reason is None:
                gt_a, gt_b           = map(int, sample['GT'].split('|'))
                allele_counts[gt_a] += 1
                allele_counts[gt_b] += 1
                if proc_lflanks:
                    lf_a,lf_b            = map(int, sample['LFGT'].split('|'))
                    lflank_counts[lf_a] += 1
                    lflank_counts[lf_b] += 1
                if proc_rflanks:
                    rf_a,rf_b            = map(int, sample['RFGT'].split('|'))
                    rflank_counts[rf_a] += 1
                    rflank_counts[rf_b] += 1
            else:
                filter_counts[sample.sample][filter_reason] += 1
                total_counts[filter_reason] += 1

        # Update the status about whether we'd like to keep the lflank and rflank info
        proc_lflanks = sum(lflank_counts[1:]) > 0
        proc_rflanks = sum(rflank_counts[1:]) > 0

        # Build allele index mapping for the STR alleles and flank sequences
        allele_indices, lflank_indices, rflank_indices = {0:0}, {0:0}, {0:0}
        for indexes,counts in zip([allele_indices, lflank_indices, rflank_indices], [allele_counts, lflank_counts, rflank_counts]):
            num_alleles = 1
            for i in range(1, len(counts)):
                if counts[i] != 0 or (not can_remove_alleles and indexes == allele_indices):
                    indexes[i]   = num_alleles
                    num_alleles += 1

        # Generate the instance we'll use to format sample entries for the filtered VCF records
        filt_format_fields = record.FORMAT.split(':')
        if not proc_lflanks:
            filt_format_fields = list(filter(lambda x: x != "LFGT", filt_format_fields))
        if not proc_rflanks:
            filt_format_fields = list(filter(lambda x: x != "RFGT", filt_format_fields))
        record.FORMAT = ":".join(filt_format_fields)
        samp_fmt = vcf.model.make_calldata_tuple(filt_format_fields)
        for fmt in samp_fmt._fields:
            entry_type = vcf_reader.formats[fmt].type
            entry_num  = vcf_reader.formats[fmt].num
            samp_fmt._types.append(entry_type)
            samp_fmt._nums.append(entry_num)

        # Generate a new VCF record by applying the filters to each sample
        new_samples       = []
        num_filt          = 0
        num_kept          = 0
        total_dp          = 0
        total_dsnp        = 0
        total_dstutter    = 0
        total_dflankindel = 0
        for sample in record:
            missing = (sample['GT'] is None or sample['GT'] == "./." or sample['GT'] == ".")
            if missing or sample_filt_dict[sample.sample] is not None:
                num_filt += (1 if not missing else 0)
                sampdata  = []
                for i in range(len(samp_fmt._fields)):
                    key = samp_fmt._fields[i]
                    if key == "GT":
                        sampdata.append("./.")
                    else:
                        sampdata.append(None if key != "FILTER" else (sample["FILTER"] if missing else sample_filt_dict[sample.sample]))
            else:
                num_kept  += 1
                gt_a, gt_b = map(int, sample['GT'].split('|'))
                new_gt     = "%d|%d"%(allele_indices[gt_a], allele_indices[gt_b])
                new_lfgt   = "" if not proc_lflanks else "|".join(list(map(lambda x: str(lflank_indices[int(x)]), sample['LFGT'].split('|'))))
                new_rfgt   = "" if not proc_rflanks else "|".join(list(map(lambda x: str(rflank_indices[int(x)]), sample['RFGT'].split('|'))))
                sampdata   = []
                for i in range(len(samp_fmt._fields)):
                    key = samp_fmt._fields[i]
                    if key == "GT":
                        sampdata.append(new_gt)
                    elif key == "LFGT":
                        sampdata.append(new_lfgt)
                    elif key == "RFGT":
                        sampdata.append(new_rfgt)
                    else:
                        sampdata.append(sample[key])
                total_dp          += sample['DP']
                total_dsnp        += sample['DSNP']
                total_dstutter    += sample['DSTUTTER']
                total_dflankindel += sample['DFLANKINDEL']

            call = vcf.model._Call(record, sample.sample, samp_fmt(*sampdata))
            new_samples.append(call)
        record.samples = new_samples

        # Don't output the filtered record if we don't meet the minimum number of calls
        if num_kept < args.MIN_CALLS:
            continue
          
        # Fix set of reported alleles
        new_alleles = [record.alleles[0]]
        for i in range(1, len(record.alleles)):
            if allele_counts[i] != 0 or not can_remove_alleles:
                new_alleles.append(record.alleles[i])
        record.alleles = new_alleles
        record.ALT     = new_alleles[1:]
        if len(record.ALT) == 0:
            record.ALT = ["."]

        # Modify relevant INFO fields for the flank sequences
        if len(lflank_counts) > 1:
            if not proc_lflanks:
                del record.INFO["LFLANKS"]
            else:
                record.INFO["LFLANKS"] = list(map(lambda x: record.INFO["LFLANKS"][x[0]], sorted(lflank_indices.items(), key = lambda p: p[1])))
        if len(rflank_counts) > 1:
            if not proc_rflanks:
                del record.INFO["RFLANKS"]
            else:
                record.INFO["RFLANKS"] = list(map(lambda x: record.INFO["RFLANKS"][x[0]], sorted(rflank_indices.items(), key = lambda p: p[1])))

        # Recompute and set other INFO fields
        if 'NFILT' in record.INFO:
            record.INFO['NFILT'] += num_filt
        record.INFO['DP']          = total_dp
        record.INFO['DSNP']        = total_dsnp
        record.INFO['DSTUTTER']    = total_dstutter
        record.INFO['DFLANKINDEL'] = total_dflankindel

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
                record.INFO['AC'] = allele_counts[1:] if not can_remove_alleles else list(filter(lambda x: x != 0, allele_counts[1:]))
        if 'AN' in record.INFO:
            record.INFO['AN'] = sum(allele_counts)
                  
        vcf_writer.write_record(record)

    if samples_to_keep is not None:
        vcf_writer.close()

if __name__ == "__main__":
    main()
