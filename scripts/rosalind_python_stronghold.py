# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 22:13:09 2025

Python Stronghold scripts
@author: zhiva

AI usage statement:
    ChatGPT and Google AI (through google's search engine) were used in lieu of searching documentation. No code was written but to understand functions and methods, I primarily used these to search or identify. 
"""

#1 Counting DNA Nucleotides
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

#2 Transcribing DNA to RNA

DNAstring = "TGATTCC"
RNAstring = ""

for i in DNAstring:
    if i == "T":
        RNAstring += "U"
    else:
        RNAstring += i

print(RNAstring)

#3 Complementing a Strand of DNA

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

#4 Rabbits and Recurrence Relations

#Fn = Fn-1 + Fn-2

# gen[n] = gen[n-1] + k * gen[n-2]

k = 2 # offspring multiplier, number of new pairs
n = 33 # number of months to calculate
gen = [1, 1]  # Initial 2 generations, entered manually to avoid indexing error.

# Loop will and calculate next generation based on fibonnaci pattern then append.
for i in range(2, n):
    gen.append(gen[i-1] + k * gen[i-2])


print(gen)

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

### Counting Point Mutations
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
with open("data\\hamming.txt", "r") as hammingdist:
    mutstrings = hammingdist.read()
    mutstrings = mutstrings.split("\n")
    string_a = mutstrings[0]
    string_b = mutstrings[1]
    compare(string_a, string_b)
