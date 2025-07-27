# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 16:01:21 2025

@author: zhiva

2/23/25
"""
#2 Variables and Some Arithmetic
a = 932
b = 823

c = a**2 + b**2
print(c)

#3 Strings and lists
my_string = "1m4CYPLZhV5fmeOwKGcUpvjzuMe1sU6Castor0cIqHwbNuv0qrEFEvweliczkowskiiZyQcdB785dMTH9j3EqQjiEl19abI1D4l3urKNi7zunqOm2DveMuQXlTyiHqhjHSjfhxuBUMIYjpQZMkDp4PVDxKbsdXd5IXGJm6Jjq."
substring1 = my_string[31:37]
substring2 = my_string[54:67]

print(substring1 + " " + substring2)

#4 Conditions and Loops
# Sum all odd numbers in range.
n = 0

for i in range(4710, 9095):
    if i % 2 == 1:
        n += i
        
print("Final answer: " + str(n))

#5 Working with Files
# Return odd numbered lines from given text.
with open("Rosalind/rosalind_ini5.txt", "r") as testfile, open("Rosalind/newfile.txt", "a") as newfile:
    i = 1
    for line in testfile:
        line = line.strip()
        if i % 2 == 0:
            newfile.write(line + "\n")
        i += 1
        
#6 Dictionaries
dset = "When I find myself in times of trouble Mother Mary comes to me Speaking words of wisdom let it be And in my hour of darkness she is standing right in front of me Speaking words of wisdom let it be Let it be let it be let it be let it be Whisper words of wisdom let it be And when the broken hearted people living in the world agree There will be an answer let it be For though they may be parted there is still a chance that they will see There will be an answer let it be Let it be let it be let it be let it be There will be an answer let it be Let it be let it be let it be let it be Whisper words of wisdom let it be Let it be let it be let it be let it be Whisper words of wisdom let it be And when the night is cloudy there is still a light that shines on me Shine until tomorrow let it be I wake up to the sound of music Mother Mary comes to me Speaking words of wisdom let it be Let it be let it be let it be yeah let it be There will be an answer let it be Let it be let it be let it be yeah let it be Whisper words of wisdom let it be"
dict1 = {}

for word in dset.split(' '):
    dict1.setdefault(word, 0)
    dict1[word] += 1

#Prints to console with formatting that fits Rosalind for easy copy paste.    
for key, value in dict1.items():
    print(f"{key} {value}")
