#!/Users/jacobbrady/virtual_envs/py35/bin/python
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt 
import argparse

# currently doesn't use plasmid sequence, just insert
# 5'-3' sequences of overlaps
# tuple containing (fw,rw) sequences
overlaps = {"pET_SUMO":("gct cac aga gaa cag att ggt ggt",
                        "ccg aat aaa tac cta agc ttg tct"),
            "pET28_HIS_TEV":("c ggc aga gaa aac ttg tat ttc cag ggc",
                             "g acg gag ctc gaa ttc gga tcc"),
            "none":("",""),
            "spastinC-pETSUMO":("c aag gac ttc ggc gat acc acc gtg",
                "ccg aat aaa tac cta agc ttg tct")
            }

tm_calc = {"wallace":mt.Tm_Wallace,
           "GC":mt.Tm_GC,
           "NN":mt.Tm_NN}

def find_primer(seq,tm,method=tm_calc["wallace"]):
    """ Find primer that is above or equal to desired tm
    
        Function Arguments:
        seq -- sequence string or Seq obj
        tm -- target melting temp
        method -- Bio.SeqUtils.MeltingTemp function 
    """
    p = ""
    t = 0
    while t <= tm-2:
        nt = seq[0]
        seq = seq[1:]
        p+=nt
        t = method(p)
    return Seq(p),t 

def find_ORF(seq,from5=False):
    while len(seq)%3:
        if from5:
            n,seq = seq[0],seq[1:]
        else:
            n,seq = seq[-1],seq[:-1]
    return seq

parser = argparse.ArgumentParser(description='make primers for GA')
parser.add_argument('--insert','-i',type=str,
        help="FASTA file containing DNA sequence of insert",required=True)
parser.add_argument('--tm','-t',type=int,
        help="Approx melting temp desired",default=68)
parser.add_argument('--outname','-o',type=str,
        help="Name of output file (yaml format)",default="primers.yml")
parser.add_argument('--addstop','-s',type=bool,
        help="add stop codon)",default=False)
parser.add_argument('--backbone','-b',default="pET28_HIS_TEV")
parser.add_argument('--model','-m',type=str,default="wallace",choices=["wallace","GC","NN"])

args = parser.parse_args()

insert = args.insert
tm = args.tm
outname = args.outname
backbone = args.backbone
model = tm_calc[args.model]
overlap_seq_fw = overlaps[backbone][0]
overlap_seq_rw = overlaps[backbone][1]
line = "#---------------------------------------------------#"
print("USING %s BACKBONE"%backbone)
insert = SeqIO.read(open(insert,"rU"),"fasta")
protein = insert.seq.translate()
fw = insert.seq
rw = insert.seq.reverse_complement()
fw_primer,fw_tm = find_primer(fw,tm,model)
rw_primer,rw_tm = find_primer(rw,tm,model)
print(line)
print("Forward primer")
print(fw_primer,fw_tm)
print(line)
print("Reverse primer")
print(rw_primer,rw_tm)
print(line)
full_fw_primer = overlap_seq_fw + fw_primer
full_rw_primer = overlap_seq_rw + rw_primer
print("Full forward primer")
print(full_fw_primer.ungap(" "))
print(line)
print("Full reverse primer")
print(full_rw_primer.ungap(" "))
print(line)
full_fw_primer = full_fw_primer.ungap(" ")
full_rw_primer = full_rw_primer.ungap(" ")
full_fw_primer_ORF = find_ORF(full_fw_primer)
# translate primer
fw_primer_ORF = find_ORF(fw_primer)
fw_primer_trans = fw_primer_ORF.translate()
fw_triplets = " ".join(i for i in [str(fw_primer_ORF)[i:i + 3] for i in range(0,len(str(fw_primer_ORF)),3)])
print(" "+"   ".join(i for i in str(fw_primer_trans)))
print(fw_triplets)
print(line)
rw_primer_ORF = find_ORF(rw_primer)
rw_primer_trans = rw_primer_ORF.reverse_complement().translate()[::-1]
rw_triplets = " ".join(i for i in [str(rw_primer_ORF)[i:i + 3] for i in range(0,len(str(rw_primer_ORF)),3)])
print(" "+"   ".join(i for i in str(rw_primer_trans)))
print(rw_triplets)
print(line)
#print(full_fw_primer.translate())
output = {"Forward":{"dna_sequence":{"5'":str(fw_primer),"3'":str(fw_primer.complement())},
                     "overlap_sequence":{"5'":str(overlap_seq_fw)},
                     "primer":{"5'":str(full_fw_primer)},
                     "Tm":fw_tm},
          "Reverse":{"dna_sequence":{"5'":str(rw_primer),"3'":str(rw_primer.complement())},
                     "overlap_sequence":{"5'":str(overlap_seq_rw)},
                     "primer":{"5'":str(full_rw_primer)},
                     "Tm":fw_tm}}
out = open(outname,"w")
out.write(yaml.dump(output))
out.close()
