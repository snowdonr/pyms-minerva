"""
proc.py

Plot a TIC using matplotlib
"""

import sys
sys.path.append("../..")

import pathlib

import matplotlib.pyplot as plt

from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.Display import plot_ic


data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
# Change this if the data files are stored in a different location

output_directory = pathlib.Path(".").resolve() / "output"

jcamp_file = data_directory / "gc01_0812_066.jdx"
data = JCAMP_reader(jcamp_file)
tic = data.tic

fig, ax = plt.subplots(1, 1)

# Plot the TIC
plot_ic(ax, tic, label="TIC")

# Set the title
ax.set_title("TIC for gc01_0812_066")

# Add the legend
plt.legend()

plt.show()
