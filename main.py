import numpy as np
import pandas as pd


# Copied class for hash table.
# Keys are character representations of amino acids
# values are tuples taken from Excel file
class HashTable:
    # Create empty bucket list of given size
    def __init__(self, size):
        self.size = size
        self.hash_table = self.create_buckets()

    def create_buckets(self):
        return [[] for _ in range(self.size)]

    # Insert values into hash map
    def set_val(self, key, val):

        # Get the index from the key
        # using hash function
        hashed_key = hash(key) % self.size

        # Get the bucket corresponding to index
        bucket = self.hash_table[hashed_key]

        found_key = False
        for index, record in enumerate(bucket):
            record_key, record_val = record

            # check if the bucket has same key as
            # the key to be inserted
            if record_key == key:
                found_key = True
                break

        # If the bucket has same key as the key to be inserted,
        # Update the key value
        # Otherwise append the new key-value pair to the bucket
        if found_key:
            bucket[index] = (key, val)
        else:
            bucket.append((key, val))

    # Return searched value with specific key
    def get_val(self, key):

        # Get the index from the key using
        # hash function
        hashed_key = hash(key) % self.size

        # Get the bucket corresponding to index
        bucket = self.hash_table[hashed_key]

        found_key = False
        for index, record in enumerate(bucket):
            record_key, record_val = record

            # check if the bucket has same key as
            # the key being searched
            if record_key == key:
                found_key = True
                break

        # If the bucket has same key as the key being searched,
        # Return the value found
        # Otherwise indicate there was no record found
        if found_key:
            return record_val
        else:
            return "No record found"

    # Remove a value with specific key
    def delete_val(self, key):

        # Get the index from the key using
        # hash function
        hashed_key = hash(key) % self.size

        # Get the bucket corresponding to index
        bucket = self.hash_table[hashed_key]

        found_key = False
        for index, record in enumerate(bucket):
            record_key, record_val = record

            # check if the bucket has same key as
            # the key to be deleted
            if record_key == key:
                found_key = True
                break
        if found_key:
            bucket.pop(index)
        return

    # To print the items of hash map
    def __str__(self):
        return "".join(str(item) for item in self.hash_table)


# for proteins in water measured at 280 nm. Units are M^-1 cm^-1
extTyr = 1490
extTrp = 5500
extCys = 125
inStr = 'MDSKGSSQKGSRLLLLLVVSNLLLCQGVVSTPVCPNGPGNCQVSLRDLFDRAVMVSHYIHDLSSEMFNEFDKRYAQ' \
        'GKGFITMALNSCHTSSLPTPEDKEQAQQTHHEVLMSLILGLLRSWNDPLYHLVTEVRGMKGAPDAILSRAIEIEEENK' \
        'RLLEGMEMIFGQVIPGAKETEPYPVWSGLPSLQTKDEDARYSAFYNLLHCLRRDSSKIDTYLKLLNCRIIYNNNC'

data = pd.read_excel(r'C:\Users\bchiu\OneDrive\Desktop\Extinction-Coefficient-Project\MW-data.xlsx')
# place "r" before the path string to address special character, such as '\'.
# Don't forget to put the file name at the end of the path + '.xlsx'
dataFrame = pd.DataFrame(data)
# NOTE removes last two columns because they are currently empty
dataFrame = dataFrame.drop(['ExtCoff', 'pI value'], axis=1)
tyrCount = 0
trpCount = 0
cysCount = 0
AAMatrix = dataFrame.to_numpy()


numAminoAcids = int(AAMatrix.size/AAMatrix[0].size)
hTable = HashTable(numAminoAcids)

for i in range(numAminoAcids):
    AAMatrix[i, 1] = AAMatrix[i, 1][1]
    hTable.set_val(AAMatrix[i, 1], AAMatrix[i])

totalMW = 0.0
for i in inStr:
    cur = hTable.get_val(i)
    cur[4] += 1
    totalMW += float(cur[2])
#   Tryptophan
    if cur[1] == 'W':
        trpCount += 1
#   Tyrosine
    elif cur[1] == 'Y':
        tyrCount += 1
#   Cysteine
    elif cur[1] == 'C':
        cysCount += 1

totalMW = round(totalMW, 5)
MWNoCorrection = totalMW
print('Molecular Weight NOT accounting for water removed in bonds: ' + str(MWNoCorrection))

# correcting for water removing in bonds
totalMW += -18.02 * (len(inStr)-1)
totalMW = round(totalMW, 5)
print('Molecular Weight accounting for water removed in bonds: ' + str(totalMW))

extCoeff = extTyr * tyrCount + extTrp * trpCount + extCys * cysCount
print('Extention Coefficent: ' + str(extCoeff))


# delete or remove a value
# hTable.delete_val('A')
# print(hTable)
