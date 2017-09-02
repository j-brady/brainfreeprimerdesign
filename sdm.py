#!/Users/jacobbrady/virtual_envs/py35/bin/python
import sys
from collections import OrderedDict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Data import CodonTable
import yaml
from Levenshtein import distance
from numpy import argmin, array

tm_dic = {"A":2,"T":2,"G":4,"C":4}
base_pairs = {"A":"T","T":"A","C":"G","G":"C"}
def tm_calc(seq):
    tm = 0
    for i in seq:
        tm += tm_dic[i]
    return tm

def closest_codon(codon,mutation,forward_table):
    # Takes the first codon that has the minimum distance from original
    # codon sequence.
    # make a list of potential codons.
    potential_codons = []
    for k,v in forward_table.items():
        if v == mutation:
            potential_codons.append(k)
    distances = []
    # calculate distances between codons.
    for i in potential_codons:
        distances.append(distance(str(codon),str(i)))
    ind = argmin(array(distances))   
    new_codon = potential_codons[ind]
    return new_codon

# codon table
bacterial = CodonTable.unambiguous_dna_by_name["Bacterial"]
to_protein = bacterial.forward_table
to_dna = bacterial.back_table

# load yaml
params = yaml.load(open(sys.argv[1],"r"))
start = params["start"]
# Read fasta file containing sequence.
seq = params["seq"]
if seq.endswith(".fasta"):
    seq = SeqIO.read(seq,"fasta")
    seq = Seq(str(seq.seq),generic_dna)
# unpack mutations
mutations = params["mutations"]
#From,position,To = params["mutation"]
out_dict = {}#OrderedDict()
for From,position,To in mutations:

    # sanity check 
    protein = seq.translate(table="Bacterial")
    corrected_pos = position-start
    if protein[corrected_pos] == From:
        print("--------------------------------------------------")
        print(From,position,To)
        key = "%s %d %s"%(From,position,To)
        out_dict[key] = {} 
        print("--------------------------------------------------")
        codon = seq[(corrected_pos)*3:(corrected_pos)*3+3]
        print("Codon:  %s"%codon)
        #print(to_dna[From])
        #print(to_dna[To])
        new_codon = closest_codon(codon,To,to_protein)
        print("Change: %s"%new_codon)
        five_prime_seq = seq[(corrected_pos)*3-26:(corrected_pos)*3]
        five_prime_tm = tm_calc(five_prime_seq) 
        three_prime_seq = seq[(corrected_pos)*3+3:(corrected_pos)*3+26]
        three_prime_tm = tm_calc(three_prime_seq) 
        while five_prime_tm > 58:
            five_prime_seq = five_prime_seq[1:]
            five_prime_tm = tm_calc(five_prime_seq)
        while three_prime_tm > 58:
            three_prime_seq = three_prime_seq[:-1]
            three_prime_tm = tm_calc(three_prime_seq)
        #sense = five_prime_seq+to_dna[To]+three_prime_seq
        sense = five_prime_seq+new_codon+three_prime_seq
        anti_sense = sense.complement()
        out_dict[key]["Codon"] = str(codon)
        out_dict[key]["Change"] = str(new_codon)
        out_dict[key]["Sense"] = str(five_prime_seq+" "+new_codon+" "+three_prime_seq)
        comp_new_codon = "".join(base_pairs[i] for i in new_codon)
        anti_sense_5 = "".join(base_pairs[i] for i in five_prime_seq)
        anti_sense_3 = "".join(base_pairs[i] for i in three_prime_seq)
        anti_sense = anti_sense_5+" "+comp_new_codon+" "+anti_sense_3
        out_dict[key]["Anti-sense"] = str(anti_sense)
        out_dict[key]["5'-Tm"] = five_prime_tm
        out_dict[key]["3'-Tm"] = three_prime_tm
        out_dict[key]["Fw"] = str(sense)
        out_dict[key]["Rw"] = str(sense.reverse_complement())
        print(five_prime_tm,three_prime_tm)
        #print(five_prime_seq+" "+new_codon+" "+three_prime_seq)
        print("Sense      : %s"%str(five_prime_seq+" "+new_codon+" "+three_prime_seq))
        print("Anti-sense : %s"%anti_sense)
        print("Order")
        print("Fw: %s"%sense)
        print("Rw: %s"%sense.reverse_complement())

with open("primers.yml","w") as f:
    f.write(yaml.dump(out_dict,default_flow_style=False))
