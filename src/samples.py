#!/usr/env/python3
import matplotlib.pyplot as plt
from matplotlib._cm import cubehelix
import numpy as np
import math

# Extract only the first number from a string
def extractFirstNumber(s):
    started = False
    number = ""
    s = s.replace('order-1', '')
    s = s.replace('order-2', '')
    s = s.replace('order-3', '')
    s = s.replace('order-4', '')
    s = s.replace('order-5', '')
    s = s.replace('order-6', '')
    for char in s:
        if char.isdigit() or char == '.' or char == '-':
            number += char
            started = True
        elif started:
            try:
                return float(number)
            except ValueError:
                started = False
                number = ""
    try:
        return float(number)
    except ValueError:
        return None

# Open the data file
points = []
with open('data/measure.dat', 'r') as file:

    # Get each line
    filename = ""
    for line in file:

        # Ignore commented lines
        if line.startswith('#') or line.startswith('//') or line.startswith('%'):
            continue

        # Split the line into parts
        parts = line.split('&')
        for i in range(len(parts)):
            parts[i] = parts[i].replace('\\', '')
            parts[i] = parts[i].replace('hline', '')
            parts[i] = parts[i].strip()

        # If it's the right length
        if len(parts) == 7 or len(parts) == 8:

            # Extract the different things
            point = {}
            if len(parts) == 8:
                point = {
                        "moment": parts[0],
                        "linear": parts[1],
                        "rdm": parts[2],
                        "sym": parts[3],
                        "shots": parts[4],
                        "diff": parts[5],
                        "time": parts[6],
                        "note": parts[7],
                        "filename": filename,
                        }
            elif len(parts) == 7:
                point = {
                        "moment": parts[0],
                        "linear": parts[1],
                        "rdm": parts[2],
                        "sym": parts[3],
                        "shots": "-1",
                        "diff": parts[4],
                        "time": parts[5],
                        "note": parts[6],
                        "filename": filename,
                        }

            # Convert numeric values
            pointWithVals = {}
            for key in point:
                if key not in ["note", "filename"]:
                    if "None" in point[key]:
                        pointWithVals[key + "Val"] = 0
                    elif key == "diff":
                        splitParts = point[key].split(',')
                        val1 = extractFirstNumber(splitParts[0])
                        val2 = extractFirstNumber(splitParts[1])
                        pointWithVals[key + "LowerVal"] = val1
                        pointWithVals[key + "UpperVal"] = val2
                    else:
                        val = extractFirstNumber(point[key])
                        if val is not None:
                            pointWithVals[key + "Val"] = val
                pointWithVals[key] = point[key]
            point = pointWithVals

            # Add to the data list
            points.append(point)
            print(parts)
            print(point)
            print()

        # If it's specifying a different file
        elif parts[0] == 'file':
            filename = parts[1]

# Check for points with the same name, number of shots, and filename
pointsDict = {}
for point in points:
    identifier = (point["filename"], point["note"], point["shotsVal"])
    if identifier in pointsDict:
        pointsDict[identifier]["diffLowerValSum"] += point["diffLowerVal"]
        pointsDict[identifier]["diffUpperValSum"] += point["diffUpperVal"]
        pointsDict[identifier]["diffLowerValSumSquared"] += point["diffLowerVal"] ** 2
        pointsDict[identifier]["diffUpperValSumSquared"] += point["diffUpperVal"] ** 2
        pointsDict[identifier]["count"] += 1
    else:
        pointsDict[identifier] = point
        pointsDict[identifier]["count"] = 1
        pointsDict[identifier]["diffLowerValSum"] = pointsDict[identifier]["diffLowerVal"]
        pointsDict[identifier]["diffLowerValSumSquared"] = pointsDict[identifier]["diffLowerVal"] ** 2
        pointsDict[identifier]["diffUpperValSum"] = pointsDict[identifier]["diffUpperVal"]
        pointsDict[identifier]["diffUpperValSumSquared"] = pointsDict[identifier]["diffUpperVal"] ** 2

points = []
for identifier in pointsDict:
    point = pointsDict[identifier]
    count = point["count"]
    if count < 1:
        continue
    point["diffLowerVal"] = point["diffLowerValSum"] / count
    point["diffUpperVal"] = point["diffUpperValSum"] / count
    if point["shotsVal"] == -1:
        point["sdLowerVal"] = 0
        point["sdUpperVal"] = 0
    else:
        point["sdLowerVal"] = math.sqrt((point["diffLowerValSumSquared"] / count) - (point["diffLowerVal"] ** 2))
        point["sdUpperVal"] = math.sqrt((point["diffUpperValSumSquared"] / count) - (point["diffUpperVal"] ** 2))
    print(identifier, count, point["diffLowerVal"], point["sdLowerVal"], point["diffUpperVal"], point["sdUpperVal"])
    points.append(point)

# The different datasets
filenames = set(point["filename"] for point in points)
filenames = sorted(filenames)
notes = set(point["note"] for point in points)
notes = sorted(notes)

# Set the font size
plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'font.family': 'serif'})
linewidth = 2.5

# Output the full list of things
print("allowed = {")
for filename in filenames:
    print('    "' + filename + '": {')
    for note in notes:
        hasData = False
        for point in points:
            if point["note"] == note and point["filename"] == filename:
                hasData = True
                break
        if hasData:
            print('        "' + note + '",')
    print('    }, ')
print("}")

# Which things to plot
allowed = {
    "energy": [
        "sdp",
        "sdp+all2, 99.7%",
        "sdp+first100, 99.7%",
        "sdp+onlyobj, 99.7%",
    ], 
    "heat": [
        "sdp",
        "first100-z, 99.7%",
        "sdp+first100-z, 99.7%",
    ], 
    "large": [
        "sdp",
        "onlyobj, 99.7%",
        "sdp+onlyobj, 99.7%",
    ], 
    "purity": [
        "sdp",
        "all2, 99.7%",
        "sdp+all2, 99.7%",
    ], 
    "confidence": [
        "sdp+all2, 68%",
        "sdp+all2, 95%",
        "sdp+all2, 99.7%",
    ],
}

# Plot the data
for filename in filenames:
    if filename not in allowed.keys():
        continue
    
    # Colors for the plots
    # cmap = plt.get_cmap('cubehelix')
    # cmap = plt.get_cmap('viridis')
    # numUniqueNotes = len(allowed[filename])
    # indices = np.linspace(0, cmap.N, numUniqueNotes+1)
    # indices = np.linspace(0, cmap.N, numUniqueNotes)
    # colors = [cmap(int(i)) for i in indices]
    colors = ["black", "royalblue", "lightcoral", "mediumseagreen"]

    # Set up the figure
    plt.figure(figsize=(7, 5))
    plt.clf()
    plt.xlabel("Number of Shots")
    if filename == "purity" or filename == "confidence":
        plt.ylabel("Purity Lower Bound")
    elif "mag" in filename:
        plt.ylabel("Magnetization Bounds")
    elif "large" in filename:
        plt.ylabel("Ground-state Energy Lower Bound")
    elif "energy" in filename:
        plt.ylabel("Ground-state Energy Lower Bound")
    elif "heat" in filename:
        plt.ylabel("Heat Current Bounds")
    plt.grid(True)
    noteToColor = {}
    nextCol = 0

    # For each constraint set
    for note in allowed[filename]:
        if note not in notes:
            print("Error: note " + note + " not found for file " + filename)
            exit(1)

        # Get the dataset
        x = []
        yLower = []
        yUpper = []
        yLineLower = None
        yLineUpper = None
        for point in points:
            if point["note"] == note and point["filename"] == filename:
                if point["shotsVal"] == -1:
                    print(point)
                    if yLineLower is None:
                        yLineLower = point["diffLowerVal"]
                    if yLineUpper is None:
                        yLineUpper = point["diffUpperVal"]
                    continue
                x.append(point["shotsVal"])
                yLower.append(point["diffLowerVal"])
                yUpper.append(point["diffUpperVal"])

        # Skip if no data for this note and this file
        if (len(x) == 0 or len(yLower) == 0 or len(yUpper) == 0) and note != "sdp":
            continue

        # Sort the dataset
        sorted_indices = sorted(range(len(x)), key=lambda i: x[i])
        x = [x[i] for i in sorted_indices]
        yLower = [yLower[i] for i in sorted_indices]
        yUpper = [yUpper[i] for i in sorted_indices]

        # Check if the color has already been used
        commaLoc = note.find(',')
        if commaLoc != -1:
            noteNoPercent = note[:commaLoc]
        else:
            noteNoPercent = note
        firstTime = noteNoPercent not in noteToColor.keys()
        if filename == "confidence":
            firstTime = True
            noteNoPercent = note
            noteNoPercent = noteNoPercent.replace('sdp+all2, ', '').strip()
        if firstTime:
            color = colors[nextCol % len(colors)]
            noteToColor[noteNoPercent] = color
            nextCol += 1
        else:
            color = noteToColor[noteNoPercent]

        # Get the error if it exists
        yErrorLower = None
        yErrorUpper = None
        if "sdLowerVal" in points[0] and "sdUpperVal" in points[0]:
            yErrorLower = [[], []]
            yErrorUpper = [[], []]
            for point in points:
                if point["note"] == note and point["filename"] == filename and point["shotsVal"] != -1:
                    yErrorLower[0].append(point["sdLowerVal"])
                    yErrorLower[1].append(point["sdLowerVal"])
                    yErrorUpper[0].append(point["sdUpperVal"])
                    yErrorUpper[1].append(point["sdUpperVal"])
            yErrorLower = [yErrorLower[0], yErrorLower[1]]
            yErrorUpper = [yErrorUpper[0], yErrorUpper[1]]
            if all(err == 0 for err in yErrorLower[0]):
                yErrorLower = None
                yErrorUpper = None

        # Adjust the label
        if filename != "confidence" and filename != "energy":
            if "+" in noteNoPercent:
                noteNoPercent = "SDP & Measure"
            elif "sdp" in noteNoPercent.lower():
                noteNoPercent = "SDP"
            else:
                noteNoPercent = "Measure"
        elif filename == "energy":
            noteNoPercent = noteNoPercent.replace('sdp+', 'SDP & ')
            noteNoPercent = noteNoPercent.replace('sdp', 'SDP')
            noteNoPercent = noteNoPercent.replace('auto100', 'Measure (a)')
            noteNoPercent = noteNoPercent.replace('all2', 'Measure (b)')
            noteNoPercent = noteNoPercent.replace('onlyobj', 'Measure (c)')

        # Plot the line
        if firstTime:
            if yErrorLower is not None:
                line = plt.errorbar(x, yLower, yerr=yErrorLower, label=noteNoPercent, color=color, capsize=3, linewidth=linewidth)
            else:
                line = plt.plot(x, yLower, label=noteNoPercent, color=color, linewidth=linewidth)
        else:
            if yErrorLower is not None:
                line = plt.errorbar(x, yLower, yerr=yErrorLower, color=color, capsize=3, linewidth=linewidth)
            else:
                line = plt.plot(x, yLower, color=color, linewidth=linewidth)

        # If we need a true bound
        if yLineLower is not None:
            if noteNoPercent == "SDP":
                plt.axhline(y=yLineLower, color=color, linewidth=linewidth)
            else:
                plt.axhline(y=yLineLower, linestyle=':', color=color, linewidth=linewidth, zorder=1000)

        # If we need an upper bound too
        if filename != "purity" and filename != "energy" and filename != "large":
            if yErrorUpper is not None:
                line = plt.errorbar(x, yUpper, yerr=yErrorUpper, color=color, capsize=3, linewidth=linewidth)
            else:
                line = plt.plot(x, yUpper, color=color, linewidth=linewidth)
            if yLineUpper is not None:
                if noteNoPercent == "SDP":
                    plt.axhline(y=yLineUpper, color=color, linewidth=linewidth)
                else:
                    plt.axhline(y=yLineUpper, linestyle=':', color=color, linewidth=linewidth, zorder=1000)

    # Finish the plot
    plt.xscale('log')
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    legLoc = 'best'
    if filename == "energy":
        legLoc = (0.025, 0.66)
    elif filename == "heat":
        legLoc = (0.025, 0.52)
    elif filename == "purity":
        legLoc = (0.025, 0.67)
    ax.legend(handles, labels, loc=legLoc)
    plt.tight_layout()
    plt.savefig('data/estimation_' + filename + '.pdf', bbox_inches='tight')
    plt.show()




