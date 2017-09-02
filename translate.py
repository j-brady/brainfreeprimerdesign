#!/Users/jacobbrady/virtual_envs/py35/bin/python
import argparse
from Bio.SeqUtils import six_frame_translations
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Six frame translation using BioPython")
parser.add_argument('--fasta','-f',type=str)
parser.add_argument('--raw_seq','-s',type=str)

args = parser.parse_args()
fasta = args.fasta
raw_seq = args.raw_seq
if raw_seq:
    print(six_frame_translations(raw_seq.replace(' ','')))
else:
    seq_handle = open(fasta,"r")
    seq = SeqIO.read(seq_handle,"fasta")
    print(six_frame_translations(seq.seq))
