'''
Hi there,
This script is designed to modify the gff gene features and extract a fixed length of the upstream seq.

Enjoy!
Xiao Fei (china-fixing@hotmail.com, xiao.fei@sund.ku.dk)
'''


import copy
import sys
import argparse
# from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
# from Bio.SeqRecord import SeqRecord
# from Bio.Seq import Seq
from BCBio import GFF
import re




def parse_args():
    parser=argparse.ArgumentParser(description='''
Welcome to use Xiao's robot
Hi there,
This script is designed to modify the gff gene features and extract a fixed length of the upstream seq.

Enjoy!
Xiao Fei (china-fixing@hotmail.com, xiao.fei@sund.ku.dk)
    ''', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--IN', required=True, type=str, metavar='FILENAME', help="the gff3 filename you want to modify")
    parser.add_argument('--CUTLEN', default=100, type=float, metavar='DEFAULT 100', help="the cut length upstream of the CDS")
    parser.add_argument('--OUT', default="out.gff", type=str, metavar='FILENAME', help="Output file name, out.gff as default")
    parser.add_argument('--FEATURE', default="gene", type=str, metavar='FILENAME', help="feature to extract the upstream, default gene, can be gene or others")
    return parser.parse_args()

def main():
    args = parse_args()
    records = []
    features_up = []

    in_handle = open(args.IN)
    for seq_record in GFF.parse(in_handle):
        for feature in seq_record.features:
            if feature.type == args.FEATURE:
                if feature.location.strand == -1:
                    start_p = feature.location.end+0
                    end_p = feature.location.end + args.CUTLEN
                else:
                    start_p = feature.location.start - args.CUTLEN
                    end_p = feature.location.start+0
                feature_up_loc = FeatureLocation(start_p, end_p, strand= feature.location.strand)
                feature_up=SeqFeature(feature_up_loc, type=args.FEATURE+'_up_'+str(args.CUTLEN))
                feature_up.qualifiers = copy.deepcopy(feature.qualifiers)
                # try:
                #     feature_up.qualifiers['locus_tag'][0] = feature.qualifiers['locus_tag'][0]+'_UP'+str(args.CUTLEN)
                # except:
                #     pass
                features_up.append(feature_up)
        combine_features = seq_record.features + features_up
        seq_record.features = combine_features
        records.append(seq_record)
    # SeqIO.write(records, args.OUT, "embl")
    with open(args.OUT, "w") as out_handle:
        GFF.write(records, out_handle)
    with open(args.OUT,"r") as file:
         filedata = file.read()
    filedata = re.sub("=true$", "", filedata, flags=re.MULTILINE)
    with open(args.OUT, "w") as file:
        file.write(filedata)
    print("hey xiao, job done!")





    
    
if __name__ == '__main__':
    sys.exit(main())
