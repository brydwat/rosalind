# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 22:13:09 2025

Python Stronghold scripts
@author: zhiva

AI usage statement:
    ChatGPT and Google AI (through google's search engine) were used in lieu of searching documentation. No code was written but to understand functions and methods, I primarily used these to search or identify. 
"""
#---------------------------------------------------------------#
###1 Counting DNA Nucleotides
string = "CCGGGTGCAGACTCTGAGAATGGGCCACGCGCCCTTCCCGTAAAAAGGTTAGGTAGGTGGCGCGTGACCCGAAGTGTCATGCACCCATAATACTCTTGGCGCTGACATTAAAGAATATTACGGGTAGCAAGCGTCTGGCCCTAACCCTCACGGTGATCGTAGCCCGGAAGTTTCGCTCGCAACGATTGCAAGCAGGAGGGATCTTACTAATAAAGGAAGGGGGTCCGAGACTGGATTAGAAATTCCTCTGCGACTTCTAGGTGGAATAAACTAAAGTACTTTTTCGCCGCTGACCACGACCGCAGGTTAGTCCTACCCACGAAGATAATCGACGGCGGCGCCCGATGAGGTCTACTAGATCTAAACATAGTGGACAGCTTAATGTTGGTTCAAATTCGTTGGTAACATATGTCAGGCAGCAAGACGGAAGTGTACGTCGAGCAAAAGGCCCTGAACGACAACGAAGAACTTTATACTATGATTAACAAGCAAGGGTGAGATGAGAGGAGACGGACACGCCGAAAAATTTCATCTACGTAGATAGTCGGCCTGCGGAAGAGAGCTTTTGTGAGGGGTTTTATATATCAGATTCGCATCAATGCTATGTTCTACAAAATCACGAGTTAGCAATCAGCAATAAGCGAAGGTTCCCTTGTGTCATCACCCCAACACAAAAATGGGTAGTAGTGACTGCTGCGGAGATTATGGGCACTAGAAGGACTTAGAATACATAGCGCACAAGTAAGTTCCATAGGTGCCCGTAGTACGGCGCTAAGGTTCCATCGCATTCTATTCGGATCCCAGGTGTGAAACCATGTGCCCCT"

counts = {"A": 0, "C": 0, "G": 0, "T": 0}

for i in string:
    if i == "A":
        counts["A"] += 1
    elif i =="C":
        counts["C"] += 1
    elif i == "G":
        counts["G"] += 1
    elif i == "T":
        counts["T"] += 1

print(counts)

#---------------------------------------------------------------#
###2 Transcribing DNA to RNA

DNAstring = "TGATTCC"
RNAstring = ""

for i in DNAstring:
    if i == "T":
        RNAstring += "U"
    else:
        RNAstring += i

print(RNAstring)

#---------------------------------------------------------------#
###3 Complementing a Strand of DNA

DNAstring = "GCCCGAGATTTTCTGCAGTCGCACAACCCTTCCTGTTACGCCAACTGAAGGGGCTCTATGACTGTCGTACGAACCAGGGGCAAGCATGATAAGCGGAATGTAGCACGGTCAGGGGTATGGTCAGGCTGTCCAGTAGACAACGCTAACAATCCCTTAACTCTCTTACTTGCACGGCCACAGTCGTACGCACACCTGGCGAATAACCGTGGACTCTTGCGTCCGTGTATATTTGTTGGTGTACGTACGCCGCGTGCGTGGGAAATAGCCTTCGATTAGATAAGACCGATCGTGTCCTCCGTCGTGTGTTAAACTGACCTTACCTTAAGGAAGGGACGCCGAACTCGGGCTCATCGGATAGAAGCTCTAACGGCAGATCCTAAAGTCTACCTTATAGACCGTGATATAGCCCAATTCCTCGATTGATGCGTCGTGTCTCTTGTTCCAACTTTAGTGTCAATGCATCGCACGTCCATTCTCCAAGTCGTCTTCCATCATCGGTTGGGGAGCCATATACATGCATACCGTCTGGTTCCAGGAGGGGGGGACACTGGGCGCCTTCCTACATGTTCACAACTTAAAAACTATTTTAAATATGAGCTCGCAGTTTAGGTGTAACAGACATACCCCAGCTGTCATTGCTTCGGTTCGTTTGAAGCATGGAACGATGTACGACTGTGCGGATGTAGTTAACATTATAAGCGAGAACTTCCGGTGGTCTAATAGATAATCGGCCGTTTGAGGTCGGCCGATACGGGTGCCACTAAGTCGCTTTTCGGGTCCTCATCATTAGGGATTCCCCAATGGACATTAGCATTCCCACCCAACTGAGGGGACGTAGCGGCACAGTAACTCGAATGGATTGGACACGCGGAATACGGCACTGAAACGCTCTGGTATCGTAAAACAAGAAGAAGCTAAACAGTGATAGATCGAAGCACGCAATGTCTTCGTCTGTACGACTTCAGAATCCGGGTTACGGCCTGC"
revstring = DNAstring[::-1]
compDNAstring = ""

for i in revstring:
    if i == "A":
        compDNAstring += "T"
    elif i == "C":
        compDNAstring += "G"
    elif i == "G":
        compDNAstring += "C"
    elif i == "T":
        compDNAstring += "A"

print(compDNAstring)

#---------------------------------------------------------------#
###4 Rabbits and Recurrence Relations

#Fn = Fn-1 + Fn-2

# gen[n] = gen[n-1] + k * gen[n-2]

k = 2 # offspring multiplier, number of new pairs
n = 33 # number of months to calculate
gen = [1, 1]  # Initial 2 generations, entered manually to avoid indexing error.

# Loop will and calculate next generation based on fibonnaci pattern then append.
for i in range(2, n):
    gen.append(gen[i-1] + k * gen[i-2])


print(gen)

#---------------------------------------------------------------#
### 5 Computing GC Content ###
# Return the ID of the string with the highest GC content, followed byt the GC content of that string.

fasta_d = {}
ID_index = []
compare_d = {}

#FASTA read from .txt
GCfile = open("data\\rosalind_gc.txt", "r")
fasta = GCfile.read()

#Parses FASTA sequence identifier and bases into dictionary: fasta_d{} and stores the sequence identifer for further processing in a list: ID_index
for line in fasta.split("\n"):
    if ">" in line:
        line = line.replace(">", "")
        base_ID = line 
        fasta_d[base_ID] = ""
        ID_index.append(line)
    else:
        fasta_d[base_ID] += line    
    
GCfile.close()

#Using ID_index identifier to cycle through fasta_d values and count GC content in nested loop, then calculates and prints ratio. Adds fasta identifiers and GC ratios to dictionary for later comparison.
#Optional print features are commented out.
for i in ID_index:
    base_counts = {"AT": 0, "GC": 0}
    for b in fasta_d.get(i): 
        if b == "A" or b == "T":
            base_counts["AT"] += 1    
        elif b =="G" or b == "C":
            base_counts["GC"] += 1
    #print("AT: " + str(base_counts["AT"]) + " GC: " + str(base_counts["GC"]))
    total_base = sum(base_counts.values())
    GC_ratio = base_counts["GC"] / total_base
    #print("GC Ratio is: " + str(GC_ratio*100)) 
    compare_d[i] = GC_ratio

#Compares the GC content in the dictionary compare_d which contains FASTA identifiers and GC ratios. Prints final result for problem as a percentage.
max_GC = max(compare_d, key = compare_d.get)
print(max_GC)
print(compare_d[max_GC]*100)

#---------------------------------------------------------------#
### 6 Counting Point Mutations
# Return Hamming distance of two strings.

# Decided to try defining a function for this. 
# Uses zip() to iterate over two entered arguments and compare each item as a tuple. Counts hamming number if not equal. 
def compare(a, b):
    hamming_num = 0    
    for x, y in zip(a, b):
        if x != y:
            hamming_num += 1
    print(str(hamming_num))

# Opens .txt with two DNA strings, splits, and defines as two strings. Then runs compare() function.
with open("data\\rosalind_hamm.txt", "r") as hammingdist:
    mutstrings = hammingdist.read()
    mutstrings = mutstrings.split("\n")
    string_a = mutstrings[0]
    string_b = mutstrings[1]
    compare(string_a, string_b)

#---------------------------------------------------------------#    
### 7 Mendel's First Law
# Given three integers representing the number of heterozygous, homozygous dominant, and homozygous recessive, return the probability that two randomly selected mating organism will produce an individual possessing a dominant allele. 
# Used AI to help with understanding probability calculation.

def dominant_probability(k, m, n):
    total = k + m + n  # Total population
    
    #I couldn't think of how to do this except by writing the variations out explicitly. 
    
    # Probabilities of selecting each pair
    P_kk = (k / total) * ((k - 1) / (total - 1))
    P_km = (k / total) * (m / (total - 1)) * 2
    P_kn = (k / total) * (n / (total - 1)) * 2
    P_mm = (m / total) * ((m - 1) / (total - 1))
    P_mn = (m / total) * (n / (total - 1)) * 2
    P_nn = (n / total) * ((n - 1) / (total - 1))
    
    # Probabilities of dominant allele inheritance
    P_offspring = (
        P_kk * 1.0 +  # 100% chance
        P_km * 1.0 +  # 100% chance
        P_kn * 1.0 +  # 100% chance
        P_mm * 0.75 + # 75% chance
        P_mn * 0.5 +  # 50% chance
        P_nn * 0.0    # 0% chance
    )
    
    return P_offspring

# Example population:
k, m, n = 21, 26, 15
print(dominant_probability(k, m, n))

#---------------------------------------------------------------#
### 8 Translating RNA to Protein
# Given an RNA string, return the protein string from codons.

import re

#convert codon table to useable dictionary. 
#first for loop is text parsing into codons[] list. 2nd is split codon and proteing code into dictionary key/value pair.
with open("data/codon_table.txt", "r") as codons_file:
    codon_file = codons_file.read().strip().split("\n")
    codons = []
    for line in codon_file:
        code = re.split(r'\s{2,}', line)
        codons.extend(code)
   
    codon_table = {}
    
    for rna in codons:
        key, value = rna.split()
        codon_table[key] = value


#input string and codon length. Modify later to import string from file.
RNA_string = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
n = 3

#translation storage
protein_string = ""

#translation from codon_tabe{} to printed protein output and not printing Stop codons.
#would not separate multiple 
for codon in [RNA_string[i:i+n] for i in range(0, len(RNA_string), n)]:
    x = codon_table[codon]
    if x != "Stop":
        protein_string += x

print(protein_string)

#---------------------------------------------------------------#
### 9 Finding a Motif in DNA
# Given a DNA string, return all locations of t as a substring of s.

#reads file and splits into list where sample string and motif are assigned
with open("../data/rosalind_subs.txt", "r") as dna_file:
    dna_file = dna_file.read().split("\n")
    dna = dna_file[0]
    motif = dna_file[1]
    start_pos = []

#loop to search the DNA for a range of indexes the length of the motif + 1 to factor for base-1 indexing.
#then appends starting position to start_pos list
for i in range(len(dna) - len(motif) + 1):
    if dna[i:i + len(motif)] == motif:
        start_pos.append(i + 1)
        
#srint list as string to fit Rosalind format.
print(' '.join(map(str, start_pos)))

#---------------------------------------------------------------#
### 10 Consensus and Profile
# Given up to 10 DNA strings, return the consensus string and profile matrix. May return multiple consensus strings.

#Lists to store DNA strings and FASTA IDs
fasta_strings = []
ID_index = []

#FASTA read from .txt
with open("data/rosalind_cons.txt", "r") as GCfile:
    fasta = GCfile.read()

lines = fasta.splitlines()

#String for id in loop and list for sequence in loop.
current_id = ""
current_seq = []

#Parses FASTA sequence identifier and bases into dictionary: fasta_d{} and stores the sequence identifer for further processing in a list: ID_index
for line in lines:
    if line.startswith(">"):
        if current_id:
            fasta_strings.append("".join(current_seq))
        current_id = line[1:] #Removes ">"..
        ID_index.append(current_id) #Stores ID.
        current_seq = [] #Reset sequence list.
    else:
        current_seq.append(line.strip()) #Add line to sequence.

#Appends last sequence after loop.
if current_id:
    fasta_strings.append("".join(current_seq))

#Organizes strings separated from IDs into a matrix for further processing.
string_matrix = []
for string in fasta_strings:
    nucleotides = []
    for base in string:
        nucleotides.append(base)
    string_matrix.append(nucleotides)

#OPTIONAL commented out to show matrix.
#for row in string_matrix:
#    print(' '.join(row))

#Transpose rows to columns using zip
columns = zip(*string_matrix)

#Determine the most common nucleotide in each column. Concatenates into string of most common and strings of nucleotide counts per column.
most_common_nucleotides = ""
A = "A: "
C = "C: "
G = "G: "
T = "T: "

for col in columns:
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    for i in col:
        if i == "A":
            counts["A"] += 1
        elif i =="C":
            counts["C"] += 1
        elif i == "G":
            counts["G"] += 1
        elif i == "T":
            counts["T"] += 1
    most_common_nucleotides += max(counts, key=lambda key: counts[key])
    A += (str(counts["A"]) + " ")
    C += (str(counts["C"]) + " ")
    G += (str(counts["G"]) + " ")
    T += (str(counts["T"]) + " ")

#Print consensus sequence followed by counts of each nucleotide at each position in matrix column.
print(''.join(most_common_nucleotides))
print(A) 
print(C)  
print(G)  
print(T)