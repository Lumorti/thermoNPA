#!/usr/env/python3


# see if file1 is a subset of file2
def is_subset(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

        newLines1 = []
        hasFoundZero = False
        for line in lines1:
            if "Zero constraints" in line:
                hasFoundZero = True
            elif "Solving" in line:
                break
            elif hasFoundZero and "----" not in line:
                sortedLine = ''.join(sorted(line.strip()))
                newLines1.append(sortedLine.strip())

        newLines2 = []
        hasFoundZero = False
        for line in lines2:
            if "Zero constraints" in line:
                hasFoundZero = True
            elif "Solving" in line:
                break
            elif hasFoundZero and "----" not in line:
                sortedLine = ''.join(sorted(line.strip()))
                newLines2.append(sortedLine.strip())

        set1 = set(newLines1)
        set2 = set(newLines2)

        print(f"Set1: {set1}")
        print(f"Set2: {set2}")
        
        if set1.issubset(set2):
            return "File1 is a subset of File2"
        elif set2.issubset(set1):
            return "File2 is a subset of File1"
        else:
            return "Neither file is a subset of the other"

# Load the two files
file1 = "temp2"
file2 = "temp3"
print(is_subset(file1, file2))
