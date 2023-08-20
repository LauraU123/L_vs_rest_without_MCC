"""This script splits aligned fasta files into alignments before the L gene and after"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import Seq
import argparse

def split(oldalignment, beforeL, afterL, referencefile, min_length=1):
        before_L_list, before_L_list_2, after_L_list, after_L_list_2= ([] for i in range(4))

        ref = SeqIO.read(referencefile, "genbank")
        for feature in ref.features:
            if feature.type =='gene':
                a =str((list(feature.qualifiers.items())[0])[-1])[2:-2]
                if a == "L":
                    startofgene = int(list(feature.location)[0])
                    #endofgene =  int(list(feature.location)[-1])+1

        alignment =SeqIO.parse(oldalignment, 'fasta')
        for entry in alignment:
            newrecord_L = SeqRecord(Seq(entry.seq[startofgene:]), id=entry.id, description=entry.description)
            newrecord_before_L = SeqRecord(Seq(entry.seq[:startofgene]), id=entry.id, description=entry.description)

            if len(newrecord_L.seq.ungap('-')) >= min_length:
                after_L_list.append(newrecord_L)
            if len(newrecord_before_L.seq.ungap('-')) >= min_length:
                before_L_list.append(newrecord_before_L)

        for record in before_L_list:
            sequence =str(record.seq)
            sequence = sequence.replace("-", "")
            if sequence != "":
                before_L_list_2.append(record)
        SeqIO.write(before_L_list_2, beforeL, "fasta")

        for record in after_L_list:
            sequence =str(record.seq)
            sequence = sequence.replace("-", "")
            if sequence != "":
                after_L_list_2.append(record)
        SeqIO.write(after_L_list_2, afterL, "fasta")
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="make new reference depending on whether the entire genome or only part is to be used for the tree",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--oldalignment", required=True, help="fasta input with alignment of entire genome")
    parser.add_argument("--loutput", required=True, help="fasta output with alignment of relevant part of geome")
    parser.add_argument("--restoutput", required=True, help="fasta output with alignment of relevant part of geome")
    parser.add_argument("--reference", required=True, help="reference genbank file of entire genome")
    parser.add_argument("--min-length", default=0, type=int, help="minimal length of output sequence")
    args = parser.parse_args()

    split(args.oldalignment, args.restoutput, args.loutput, args.reference, args.min_length)