#!/usr/bin/env python

import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys
from os import path

# get input
prefix = sys.argv[1]          # MS_062
infile = sys.argv[2]          # path to bamfile
refname = 'chr22'
start = 42522610
end = 42526700
# refname = sys.argv[3]       # chr22
# start = int(sys.argv[4])           # start position: 42522610
# end = int(sys.argv[5])             # end position: 42526700

samfile = pysam.AlignmentFile(infile, "r")

records = []
for alignment in samfile.fetch(refname, start, end):
    outquryTrim = []
    for pairPos in alignment.aligned_pairs:
        if pairPos[1] in [start, end]:
            outquryTrim.append(pairPos[0])
    if len(outquryTrim) == 2:
        qstart, qend = outquryTrim
        if qstart != None and qend != None:

            seq = alignment.query_sequence[qstart:qend]
            qual = alignment.query_qualities[qstart:qend]
            
            if not alignment.is_reverse:
                
                sr = SeqRecord(seq=Seq(seq), 
                               id=alignment.query_name, description='', 
                               letter_annotations={"phred_quality":qual}
                              )
                records.append(sr)
            else:
                sr = SeqRecord(seq=Seq(seq).reverse_complement(), 
                               id=alignment.query_name, description='', 
                               letter_annotations={"phred_quality":qual}
                              )
                records.append(sr)
    
samfile.close()

fname, ext = path.splitext(infile)
outf = f"{fname}.target.fastq"
SeqIO.write(records, outf, 'fastq')