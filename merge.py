import argparse
from pysam import VariantFile
import vcfpy
from collections import OrderedDict
import time


## 
# - check if GTs are the same across calls
# - make sure no other calls between??
# - use a dictionary(i think) instead of recording newly merged variants as an array and using indices for getting start/stop/svlen/....


parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', dest="input", help='input vcf file', required=True, action="store")
parser.add_argument('--output', '-o', dest="output", help='output vcf file', required=False, default="out.vcf", action="store")
parser.add_argument('--window', '-w', dest="window", help='max merge window size', required=False, default=1000, type=int, action="store")


args = parser.parse_args()
in_vcf_name  = args.input
out_vcf_name = args.output
window_size = args.window

#reader = vcfpy.Reader.from_path(in_vcf_name)
## add INFO tag for merged call
#reader.header.add_info_line(OrderedDict({'ID': 'MERGED_CALL', 'Number': 0, 'Type': 'Flag', 'Description': 'Record merged from 2 or more individual records'}))
##writer = vcfpy.Writer.from_path(out_vcf_name, reader.header)
#records = in_vcf.fetch()
#stream = reader.stream

# pysam
reader = VariantFile(in_vcf_name)

reset = 0
start = 0
stop = 0
current = 0
chrom = ""
counter = 0
new_events = []
# helper recursive function
def findEvent(reader, svtype, chrom, stop, window_size):
    "find if there is a neigboring event that could be merged"
    reset = reader.tell()
    result = []
    for rec in reader:
        rec_chrom = rec.contig
        rec_svtype = rec.info["SVTYPE"]
        rec_start = rec.pos
        if(chrom != rec_chrom or rec_svtype != svtype):
            reader.seek(reset)
            return result
        if(((stop + window_size) > rec_start) and ((stop - window_size) < rec_start)):
            rec_stop = rec.stop
            recursive_result = findEvent(reader, rec_svtype, rec_chrom, rec_stop, window_size)
            result = [[rec_chrom, rec_start, rec_stop, rec_stop - rec_start]]
            result = result + recursive_result
            return result
    return result


# find events use pysam.
# allows easy stream.seek(position) and stream.tell() usage in the findevents local function
# which will reduce having to iterate through all the records each time..

for rec in reader:
    if(reset == 0 or chrom != rec.contig):
        reset = reader.tell()
        start = rec.pos
        stop = rec.stop
        chrom = rec.contig
        svtype = rec.info["SVTYPE"]
        res = findEvent(reader, svtype, chrom, stop, window_size)
        if(res != []):
            res = [[chrom, start, stop, stop - start]] + res
            count = len(res)
            event_lens = []
            for event in res:
                event_lens = event_lens + [event[3]]
            new_start = start
            new_end = res[-1][2]
            new_svlen = new_end - new_start
            new_svtype = svtype

            new_events.append([chrom, new_start, new_end, new_svtype, new_svlen, count, event_lens])
            #print(f"{chrom}:{new_start}-{new_end}. merged {count} calls,svlen={new_svlen},svtype={new_svtype}")

# close stream
reader.close()
# log the summary of found events instead of output in above loop
new_event_count = len(new_events)
print(f"Found {new_event_count} that can be merged based on the input {window_size} merging distance")
for event in new_events:
    chrom = event[0]
    start = event[1]
    end = event[2]
    svtype = event[3]
    svlen = event[4]
    event_count = event[5]
    event_svlen = event[6]
    print(f"  {svtype} {chrom}:{start}-{end} \tsvlen={svlen}\tmerged_events={event_count}")

# if output file write to that file
# problem is that this will require cycling through the records again
# but can make sure the results are ordered by position
# when manipulating actual VCF records the pysam library is less useful than the vcfpy library
reader = vcfpy.Reader.from_path(in_vcf_name)
new_header = reader.header
new_header.add_filter_line(vcfpy.OrderedDict([('ID', 'MERGED_CALL'), ('Description', 'Record merged from 2 or more individual records')]))

writer = vcfpy.Writer.from_path(out_vcf_name, new_header)
count = 0
new_event_index = 0
for record in reader:
    count = count + 1
    writer.write_record(record)
    if(new_event_index < new_event_count):
            chrom = new_events[new_event_index][0]
            start = new_events[new_event_index][1]
            if((record.POS == start) & (chrom == record.CHROM)):
                start = new_events[new_event_index][1]
                end = new_events[new_event_index][2]
                svtype = new_events[new_event_index][3]
                svlen = new_events[new_event_index][4]
                info = OrderedDict({"SVTYPE": svtype, "END": end, "SVLEN": svlen})
                alt = vcfpy.SymbolicAllele(svtype)
                sample_calls = []
                for sample in record.calls:
                    gt = OrderedDict({"GT": "/".join(map(str,sample.gt_alleles))})
                    name = sample.sample
                    sample_calls.append(vcfpy.Call(name, gt))

                new_record = vcfpy.Record(chrom, start, [], "N", [alt], ".", ["MERGED_CALL"], info, ["GT"], sample_calls)
                writer.write_record(new_record)
                new_event_index = new_event_index + 1

