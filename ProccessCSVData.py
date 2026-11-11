"""
READ ME: This script is for processing data on the G2 linear cell for each charge state.
BEFORE RUNNING:
1. make sure that you have a directory (folder) that will contain all of your data.
2. Inside this directory (folder), you have directories (folders) for each charge state, and inside each of those directories (folders) you have the CSV files for that charge state.
3. each CSV file should be named with the drift potential such that the name is driftPotential.csv where driftPotential is the sum of HeDC + bias + exitDC
4. the top of each CSV file shoule be either: "Arrival Time (ms),Intensity" or "Scan Number,Intensity" depending on what data you copied over.
5. each value of the CSV file should be separated by a comma, not a tab. If you have tabs, run the fixCSV.py script first.
6. double check each constant is correct for your experiment.
"""

## Constants
#drift tube constants
DRIFT_LENGTH = 25.05 # cm

# Physical Constants
BOLTZMANN_CONSTANT = 1.380649e-23  # J/K
ELECTRON_CHARGE = 1.602176634e-19  # C
ATOMIC_MASS_UNIT = 1.66053906660e-27  # kg


import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import os
import sys
import re
from scipy.signal import find_peaks
from scipy.optimize import curve_fit



def main(fileDirectory):
    fileNames = os.listdir(fileDirectory)
    driftTimePeaks = []
    arrivalTimeDistributions = []
    electricFeild = []
    for fileName in fileNames:
        if fileName.endswith(".csv"):
            filePath = os.path.join(fileDirectory, fileName)
            data = pd.read_csv(filePath)
            if data.get("Scan Number", None) is not None:
                data = convert_scans_to_arrival_time_distribution(data)
            driftPotential= int(fileName.replace(".csv", ""))
            # Do something with the extracted values, e.g., calculate dead time or perform other operations
            driftTimePeaks.append(find_maximum_intensity_drift_time(data))
            fittedDriftTimeDistribution = fit_drift_time_distribution(data)
            arrivalTimeDistributions.append(fittedDriftTimeDistribution)
            potential = (driftPotential- TRANSFER_DC)*(1-2/270)
            electricFeild.append(potential/(DRIFT_LENGTH/100)) # cm -> m, so field is in V/m
    deadTime = calculate_dead_time(driftTimePeaks, electricFeild)
    plot_CCS_Distribution(arrivalTimeDistributions, driftTimePeaks, electricFeild, deadTime)
    print(f"Dead time: {deadTime}")


def calcualte_molecule_density():
    # Convert pressure from Torr to Pa
    pressure_pa = PRESSURE * 133.322
    # Calculate molecule density using the ideal gas law: n/V = P/(k*T)
    molecule_density = pressure_pa / (BOLTZMANN_CONSTANT * TEMPERATURE)
    return molecule_density

def plot_CCS_Distribution(arrivalTimeDistributions, driftTimePeaks, electricFeild, deadTime,
                          n_grid=500, peak_prominence=0.01):
    """
    Convert each fitted arrival time distribution (columns: arrival time ms, intensity)
    into a CCS distribution using the Mason-Schamp equation, interpolate them all onto a
    common CCS grid, normalise each to its own maximum (so every file counts equally,
    whatever its signal level), and plot the mean normalised CCS distribution with its
    major peaks marked. A peak counts as major if its prominence is at least
    `peak_prominence` (0.01 = 1%) of the mean distribution's maximum.
    Returns (CCSDistribution, peakCCS): an (n_grid, 2) array of CCS (A^2) and mean
    normalised intensity, and an array of the CCS values of the major peaks.
    """
    moleculeDensity = calcualte_molecule_density()
    reducedMass = (SAMPLE_MASS * GAS_MASS) / (SAMPLE_MASS + GAS_MASS) * ATOMIC_MASS_UNIT
    ccsDistributions = []

    # Calculate a CCS distribution for each arrival time distribution
    for i, ArrivalTimeD in enumerate(arrivalTimeDistributions):
        if ArrivalTimeD is None:
            print(f"Warning: No fitted distribution for index {i}. Skipping CCS calculation.")
            continue
        ArrivalTime = np.array(ArrivalTimeD, dtype=float)
        driftTime = ArrivalTime[:, 0] - deadTime          # ms, corrected for dead time
        intensity = ArrivalTime[:, 1]
        keep = driftTime > 0                               # keep rows, not individual values
        if not np.any(keep):
            print(f"Warning: No positive drift times for index {i}. Skipping CCS calculation.")
            continue
        driftTime, intensity = driftTime[keep], intensity[keep]
        mobilityK = (DRIFT_LENGTH / 100) / (driftTime / 1000 * electricFeild[i])   # m^2 V^-1 s^-1
        CCSValues = (3/16) * math.sqrt((2*math.pi) / (reducedMass*BOLTZMANN_CONSTANT*TEMPERATURE)) \
                    * (ELECTRON_CHARGE*CHARGE_STATE / moleculeDensity) * (1/mobilityK) * 1e20  # m^2 -> A^2
        # CCS increases with drift time, so the arrays are already in ascending CCS order
        ccsDistributions.append(np.column_stack((CCSValues, intensity)))

    if not ccsDistributions:
        print("Error: No CCS distributions could be calculated.")
        return None

    # Interpolate every distribution onto a shared CCS axis so they can be averaged
    ccsMin = max(d[0, 0] for d in ccsDistributions)
    ccsMax = min(d[-1, 0] for d in ccsDistributions)
    ccsGrid = np.linspace(ccsMin, ccsMax, n_grid)
    intensities = np.array([np.interp(ccsGrid, d[:, 0], d[:, 1]) for d in ccsDistributions])

    # Normalise each distribution to a maximum of 1, then average
    normalisedIntensities = intensities / intensities.max(axis=1, keepdims=True)
    meanCCSDIntensity = normalisedIntensities.mean(axis=0)
    CCSDistribution = np.column_stack((ccsGrid, meanCCSDIntensity))

    # Major peaks in the mean normalised CCS distribution
    peakIdx, _ = find_peaks(meanCCSDIntensity, prominence=peak_prominence * meanCCSDIntensity.max())
    peakCCS = ccsGrid[peakIdx]
    for ccs, height in zip(peakCCS, meanCCSDIntensity[peakIdx]):
        print(f"Major CCS peak: {ccs:.1f} A^2 (mean normalised intensity {height:.3f})")

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(CCSDistribution[:,0], CCSDistribution[:,1], color="tab:blue", lw=2, label="Mean Normalised CCS Distribution")
    ax.plot(peakCCS, meanCCSDIntensity[peakIdx], "x", color="tab:red", ms=2, mew=2, label="Major peaks")
    ax.set_xlabel(r"$^{DT}$CCS$_{He}$ (Å$^2$)")
    ax.set_ylabel("Normalised Intensity")
    xLo, xHi = CCSDistribution[:,0].min(), CCSDistribution[:,0].max()
    yMax = meanCCSDIntensity.max()
    yTop = 1.25 * yMax   # headroom for the labels above the tallest peak
    ax.set_xlim(xLo, xHi)
    ax.set_ylim(0, yTop)

    # Label placement: each label is written vertically, centred directly above its
    # peak. Its real size is measured on the figure and it is lifted until it clears the
    # curve beneath it and any label already placed, so the text never sits on the data.
    renderer = fig.canvas.get_renderer()
    dataPerPixelX = (xHi - xLo) / ax.bbox.width
    dataPerPixelY = yTop / ax.bbox.height
    gap = 0.03 * yMax
    labelStyle = dict(ha="center", va="bottom", rotation=90, fontsize=9,
                      bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
    placedBoxes = []   # (x0, x1, y0, y1) in data coordinates
    for ccs, height in sorted(zip(peakCCS, meanCCSDIntensity[peakIdx])):
        text = f"{ccs:.0f} Å$^2$"
        probe = ax.text(0, 0, text, **labelStyle)
        extent = probe.get_window_extent(renderer=renderer)
        probe.remove()
        padPixels = 0.4 * 9 * fig.dpi / 72   # the white background box adds 0.2 font sizes per side
        w = (extent.width + padPixels) * dataPerPixelX * 1.1
        h = (extent.height + padPixels) * dataPerPixelY * 1.05
        underLabel = np.abs(ccsGrid - ccs) <= w / 2
        y = (meanCCSDIntensity[underLabel].max() if np.any(underLabel) else height) + gap
        for _ in range(50):
            box = (ccs - w / 2, ccs + w / 2, y, y + h)
            overlaps = [b for b in placedBoxes
                        if box[0] < b[1] and b[0] < box[1] and box[2] < b[3] and b[2] < box[3]]
            if not overlaps:
                break
            y = max(b[3] for b in overlaps) + gap / 2
        placedBoxes.append(box)
        ax.annotate(text, xy=(ccs, height), xytext=(ccs, y), **labelStyle,
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.8, shrinkA=0, shrinkB=2))
    if placedBoxes:
        ax.set_ylim(0, max(yTop, max(b[3] for b in placedBoxes) + gap))
    plt.legend()
    plt.savefig(f"{FILE_DIRECTORY}/plots/CCS_Distribution_{CHARGE_STATE}.png", dpi=300)
    return CCSDistribution, peakCCS

def calculate_dead_time(driftTimePeaks, electricFeild):
    # Convert lists to numpy arrays for easier calculations
    driftTimePeaks = np.array(driftTimePeaks)
    electricFeild = np.array(electricFeild)
    InverseElectricFeild = 1/electricFeild

    # Calculate the slope (m) and intercept (b) of the linear fit
    m, b = np.polyfit(InverseElectricFeild, driftTimePeaks, 1)

    # Calculate the dead time as the y-intercept of the linear fit
    deadTime = b

    # R-squared value to assess the goodness of fit
    residuals = driftTimePeaks - (m * InverseElectricFeild + b)
    ssRes = np.sum(residuals**2)
    ssTot = np.sum((driftTimePeaks - np.mean(driftTimePeaks))**2)
    RSquared = 1 - (ssRes / ssTot)
    print(f"R-squared value of the linear fit: {RSquared}")
    if RSquared < 0.999:  
        print("Warning: The linear fit may not be reliable. R-squared value is less than 0.999.")
        print("Removing outliers and recalculating the fit...")
        while RSquared < 0.9995 and len(driftTimePeaks) > 2:
            # Calculate residuals and identify the index of the maximum residual
            residuals = driftTimePeaks - (m * InverseElectricFeild + b)
            max_residual_index = np.argmax(np.abs(residuals))

            # Remove the outlier from both driftTimePeaks and InverseElectricFeild
            driftTimePeaks = np.delete(driftTimePeaks, max_residual_index)
            InverseElectricFeild = np.delete(InverseElectricFeild, max_residual_index)

            # Recalculate the linear fit with the remaining data
            m, b = np.polyfit(InverseElectricFeild, driftTimePeaks, 1)

            # Recalculate R-squared value
            residuals = driftTimePeaks - (m * InverseElectricFeild + b)
            ssRes = np.sum(residuals**2)
            ssTot = np.sum((driftTimePeaks - np.mean(driftTimePeaks))**2)
            RSquared = 1 - (ssRes / ssTot)
            print(f"Updated R-squared value after removing outlier: {RSquared}")
        deadTime = b  # use the intercept from the fit after outlier removal
        if RSquared < 0.9995 and len(driftTimePeaks) <= 2:
            print("Warning: Not enough data points to perform a reliable linear fit after removing outliers.")
    else:
        print("The linear fit is reliable. R-squared value is greater than or equal to 0.999.")
    plt.figure(figsize=(10, 5))
    plt.scatter(InverseElectricFeild, driftTimePeaks, color='tab:blue', label='Data Points')
    plt.plot(InverseElectricFeild, m * InverseElectricFeild + b, color='tab:purple', label='Linear Fit')
    plt.xlabel('Inverse Electric Field (m/V)')
    plt.ylabel('Drift Time Peaks (ms)')
    plt.title('Linear Fit of Drift Time Peaks vs Inverse Electric Field')
    # Legend and fit summary in one box in the bottom right corner: the dead time and
    # R-squared are added as extra legend entries with no marker
    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    blank = plt.Line2D([], [], linestyle='none', marker='none')
    handles += [blank, blank]
    labels += [f'Dead Time (y-intercept): {deadTime:.4f} ms', f'R-squared: {RSquared:.4f}']
    ax.legend(handles, labels, loc='lower right')
    plt.savefig(f"{FILE_DIRECTORY}/plots/DeadTime_Fit_{CHARGE_STATE}.png", dpi=300)
    plt.close()
    return deadTime

def find_maximum_intensity_drift_time(data):
    # Assuming the drift time is in the first column and intensity in the second column
    drift_time = data.iloc[:, 0].values
    intensity = data.iloc[:, 1].values

    # Find peaks in the intensity data
    peaks, _ = find_peaks(intensity)

    if len(peaks) == 0:
        return None  # No peaks found

    # Get the drift time corresponding to the maximum intensity peak
    max_peak_index = peaks[np.argmax(intensity[peaks])]
    max_drift_time = drift_time[max_peak_index]

    return max_drift_time

def multi_gaussian(x, *params):
    # Sum of Gaussians; params = [a1, mu1, sigma1, a2, mu2, sigma2, ...]
    y = np.zeros_like(x, dtype=float)
    for i in range(0, len(params), 3):
        a, mu, sigma = params[i:i + 3]
        y += a * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))
    return y

def fit_drift_time_distribution(data, max_components=6, min_improvement=0.05,
                                min_sigma_points=0.5, plot=False, title=None):
    """
    Smooth the arrival time distribution with a small number of broad Gaussians.

    Gaussians are added one at a time. A new Gaussian is only kept if it reduces the
    residual sum of squares by at least `min_improvement` (0.2 = 20%); otherwise fitting
    stops, so small bumps and noise are not given their own peak. Each Gaussian's
    standard deviation is at least `min_sigma_points` data-point spacings, so peaks stay
    broad. Raise min_improvement or min_sigma_points for a smoother fit, lower them to
    follow the data more closely.

    Returns an (N, 2) array: column 0 is the arrival time (ms), column 1 is the
    total fitted intensity (sum of all fitted Gaussians) at that arrival time.
    Returns None if no fit succeeded.
    """
    # Assuming the drift time is in the first column and intensity in the second column
    data = data.iloc[:, :2].dropna()  # some CSVs end with a blank row
    drift_time = data.iloc[:, 0].values.astype(float)
    intensity = data.iloc[:, 1].values.astype(float)
    n_points = len(drift_time)
    dx = np.median(np.diff(drift_time))
    x_min, x_max = drift_time.min(), drift_time.max()

    # Peaks in the raw data, most prominent first, used as starting guesses
    peaks, props = find_peaks(intensity, prominence=0.02 * intensity.max())
    peaks = peaks[np.argsort(props["prominences"])[::-1]]

    best = None
    prev_params = []
    for n in range(1, max_components + 1):
        if n_points <= 3 * n:
            break
        # Starting guess: previous fit plus one new Gaussian, placed at the next
        # detected peak or, if none are left, at the largest residual
        residual = intensity - multi_gaussian(drift_time, *prev_params)
        idx = peaks[n - 1] if n - 1 < len(peaks) else int(np.argmax(residual))
        new_amp = max(residual[idx], 0.05 * intensity.max())
        p0 = list(prev_params) + [new_amp, drift_time[idx], max(2, min_sigma_points) * dx]

        lower = [0, x_min, min_sigma_points * dx] * n
        upper = [np.inf, x_max, (x_max - x_min) / 2] * n
        p0 = list(np.clip(p0, lower, upper))

        try:
            popt, _ = curve_fit(multi_gaussian, drift_time, intensity, p0=p0,
                                bounds=(lower, upper), maxfev=20000)
        except (RuntimeError, ValueError):
            if best is None:
                # Keep the starting guess so the next iteration still has n components to build on
                prev_params = p0
                continue
            break

        rss = np.sum((intensity - multi_gaussian(drift_time, *popt)) ** 2)
        if best is not None and rss > (1 - min_improvement) * best["rss"]:
            break  # the extra Gaussian doesn't explain enough of the data
        prev_params = list(popt)
        best = {"n_components": n, "popt": popt, "rss": rss}

    if best is None:
        print("Error - curve_fit failed for all component counts")
        return None

    popt = best["popt"]
    params = sorted([tuple(popt[i:i + 3]) for i in range(0, len(popt), 3)], key=lambda p: p[1])
    components = [multi_gaussian(drift_time, *p) for p in params]
    fitted_curve = np.sum(components, axis=0)

    if plot:
        plt.figure(figsize=(10, 5))
        plt.plot(drift_time, intensity, color="gray", alpha=0.7, label="Data")
        for i, comp in enumerate(components):
            plt.plot(drift_time, comp, "--", lw=1, label=f"Gaussian {i + 1} (mu={params[i][1]:.3f} ms)")
        plt.plot(drift_time, fitted_curve, color="red", lw=2, label=f"Sum ({best['n_components']} Gaussians)")
        plt.xlabel("Arrival Time (ms)")
        plt.ylabel("Intensity")
        plt.title(title or "Multi-Gaussian fit")
        plt.legend()
        plt.show()

    return np.column_stack((drift_time, fitted_curve))

def convert_scans_to_arrival_time_distribution(data):
    data["Arrival Time (ms)"] = data["Scan Number"] * PUSHER_TIME / 1000  # Convert to ms
    return data[["Arrival Time (ms)", "Intensity"]]

parser = argparse.ArgumentParser(description="Process CSV data for each charge state.")
parser.add_argument("directory", type=str, help="Directory containing charge state subdirectories.")
parser.add_argument("--molecularMass", type=float, default=1062.0, help="Molecular mass of the sample in g/mol (default: 1062.0).")
parser.add_argument("--gas", type=str, default="He", help="Gas used in the drift tube (default: He).")
parser.add_argument("--pressure", type=float, default=1.96, help="pressure of the gas in Torr (default: 1.96).")
parser.add_argument("--temperature", type=float, default=297.0, help="Temperature of the drift tube in K (default: 297.0).")
parser.add_argument("--transferDC", type=float, default=4.0, help="Transfer DC voltage (default: 4.0 V).")
parser.add_argument("--mode", type=str, default="V", help="Resolution of the mass spectrum (default: V) (Used only if the CSV files contain scan numbers instead of arrival times).")
parser.add_argument("--highMass", type=int, default=5000, help="High mass limit from mass spectrum (default: 5000) (Used only if the CSV files contain scan numbers instead of arrival times).")
args = parser.parse_args()
TEMPERATURE = args.temperature # K
PRESSURE = args.pressure # Torr
SAMPLE_MASS = args.molecularMass  # g/mol for the sample
gasDict = {"He": 4.002602, "N2": 28.0134, "Ar": 39.948}
GAS_MASS = gasDict[args.gas]  # g/mol for the selected gas
TRANSFER_DC = args.transferDC # V
FILE_DIRECTORY = args.directory
pusherIndex = f"{args.mode} {args.highMass}"
pusherTimeDict = {"V 600": 39.25,
                  "W 600": 75.25,
                  "V 1200": 54.25,
                  "W 1200": 106.25,
                  "V 2000": 69.25,
                  "W 2000": 137.25,
                  "V 5000": 110.25,
                  "W 5000": 216.25,
                  "V 8000": 138.25,
                  "W 8000": 273.25,
                  "V 14000": 188.25,
                  "W 14000": 363.25,
                  "V 32000": 274.25,
                  "W 32000": 547.25,
                  "V 100000": 486.25,
                  "W 100000": 547.25,
                  }
PUSHER_TIME = pusherTimeDict[pusherIndex]  # microseconds
if __name__ == "__main__":
    FILE_DIRECTORY = args.directory
    fileDirectories = os.listdir(FILE_DIRECTORY)
    if not os.path.exists(os.path.join(FILE_DIRECTORY, "plots")):
        os.mkdir(os.path.join(FILE_DIRECTORY, "plots"))
    for fileDirectory in fileDirectories:
        if fileDirectory != "plots":
            if os.path.isdir(os.path.join(FILE_DIRECTORY, fileDirectory)):
                chargeState = fileDirectory.split("/")[-1]  # Extract charge state from the directory name
                chargeState = re.sub(r'[^0-9]', '', chargeState)  # Remove non-numeric characters
                CHARGE_STATE = int(chargeState)  # Convert to integer
                main(FILE_DIRECTORY + "/" + fileDirectory)
