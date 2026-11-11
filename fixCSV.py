import os
import argparse
parser = argparse.ArgumentParser(description="Process CSV files in a directory.")
parser.add_argument("directory", type=str, help="Directory containing charge state subdirectories.")
parser.add_argument("--data", type=str, help="types of data to process (either 'scan' or 'arrival').")
args = parser.parse_args()
fileNames = os.listdir(args.directory)
for fileName in fileNames:
    if fileName.endswith(".csv"):
        filePath = os.path.join(args.directory, fileName)
        with open(filePath, 'r') as file:
            lines = file.readlines()
            lines = [line.replace('	', ', ') for line in lines]
        with open(filePath, 'w') as file:
            if args.data == "scan":
                lines[0] = "Scan Number, Intensity\n"
            elif args.data == "arrival":
                lines[0] = "Arrival Time (ms), Intensity\n"
            file.writelines(lines)
