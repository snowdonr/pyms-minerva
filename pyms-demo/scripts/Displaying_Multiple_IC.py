"""
proc.py

Plot multiple Ion Chromatograms using matplotlib
"""

import sys
sys.path.append("../..")

import pathlib

import matplotlib.pyplot as plt

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Display import plot_ic


data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)
tic = data.tic
im = build_intensity_matrix_i(data)

# Extract the desired IonChromatograms
ic73 = im.get_ic_at_mass(73)
ic147 = im.get_ic_at_mass(147)

fig, ax = plt.subplots(1, 1)

# Plot the ICs
plot_ic(ax, tic, label="TIC")
plot_ic(ax, ic73, label="m/z 73")
plot_ic(ax, ic147, label="m/z 147")

# Set the title
ax.set_title("TIC and ICs for m/z = 73 & 147")

# Add the legend
plt.legend()

plt.show()
