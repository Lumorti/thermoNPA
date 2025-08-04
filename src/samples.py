#!/usr/env/python3
import matplotlib.pyplot as plt
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

# The different datasets
filenames = set(point["filename"] for point in points)
notes = set(point["note"] for point in points)

# Plot the data
for filename in filenames:
    if "3d" in filename:
        continue

    # Set up the figure
    plt.figure(figsize=(10, 6))
    plt.clf()
    plt.xlabel("Number of Shots")
    if filename == "purity":
        plt.ylabel("Purity Lower Bound")
    elif "mag" in filename:
        plt.ylabel("Magnetization Bounds")
    elif "energy" in filename:
        plt.ylabel("Ground-state Energy Bounds")
    elif "heat" in filename:
        plt.ylabel("Heat Current Bounds")
    plt.grid(True)
    noteToColor = {}

    # For each constraint set
    for note in notes:

        # Get the dataset
        x = []
        yLower = []
        yUpper = []
        yLineLower = None
        yLineUpper = None
        for point in points:
            if point["note"] == note and point["filename"] == filename:
                if point["shotsVal"] == -1:
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

        # Plot the line
        line = plt.plot(x, yLower, label=note)

        # Check if the color has already been used
        color = line[0].get_color()
        noteNoPercent = note[:note.find(',')]
        if noteToColor.get(noteNoPercent) is None:
            noteToColor[noteNoPercent] = color
        else:
            color = noteToColor[noteNoPercent]
            line[0].set_color(color)

        # If we need a true bound
        if yLineLower is not None:
            plt.axhline(y=yLineLower, linestyle='--', color=color)

        # If we need an upper bound too
        if filename != "purity" and filename != "energy":
            line = plt.plot(x, yUpper, color=line[0].get_color())
            if yLineUpper is not None:
                plt.axhline(y=yLineUpper, linestyle='--', color=color)

    # Finish the plot
    plt.xscale('log')
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    ax.legend(handles, labels)
    plt.savefig('data/estimation_' + filename + '.png', bbox_inches='tight')
    plt.show()




