For ATD to CCSD conversion with the linear G2 cell.
File structure:
This code requires files to have a particular layout, make sure you have the following:
1. A new folder containing the following:
2. A folder for each charge state named simply after the charge (e.g. data for a +2 charge state would be called "+2"). Within each folder, you'll need:
3. A .csv file for each potential used named {drift_time}.csv where {drift_time} is the sum of the helium DC, helium exit, and bias (e.g. a helium dc of 100, exit of -40, and bias of 100 would be called 160.csv)

Making a comma separated value file (csv files):
1. For each .csv file, copy you ATD in either scans (recommended) or time from mass lynx.
2. Each row, the arrival time and intensity need to be separated with a comma, use the fixCSV.py script to do this for every file if needed
3. The first row must be either the following "Arrival Time (ms), Intensity" or "Scan Number, Intensity" (fixCSV.py can also do this)

running fixCSV.py 
in a linux terminal (install wsl) fixCSV.py can be run using the following
"python3 fixCSV.py (folder containing all of your data) --data (either scan or arrival)"

Once this is all done, you can run ProcessCSVData.py to calculate CCSDs for each charge state you've provided
running ProcessCSVData.py
before you run ProcessCSVData.py, please install the following dependencies (with either uv or conda) or via pip:

pip install pandas
pip install numpy
pip install matplotlib
pip install scipy

in a linux terminal you can run ProcessCSVData.py using the following command
"python3 ProcessCSVData.py (folder containing all of your data) --molecularMass (molecular weight of your ion) --gas (gas used, either He, N2 or Ar currently) --pressure (pressure in Torr) --temprature (temprature of the cell in kelvin) --transferDC --mode (either V for sensitivity mode or W for resolution mode, scans only) --highMass (mass range used, scans only)"
