"""proc.py

Before running this script a list of peaks must be generated and stored in 
the output/ folder

The script proc_save_peaks.py in the current directory saves this list
in the correct location

"""

from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.GCMS.Function import build_intensity_matrix
from pyms.Peak.IO import load_peaks
from pyms.Display.Class import Display

# read the raw data as a GCMS_data object
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)


# IntensityMatrix


# default, float masses with interval (bin interval) of one from min mass
print("default intensity matrix, bin interval = 1, boundary +/- 0.5")
im = build_intensity_matrix(data)


#
# IonChromatogram
#

# TIC from raw data
tic = data.get_tic()
# save TIC to a file

# Get Ion Chromatograms for 4 separate m/z channels
ic = im.get_ic_at_mass(191)
ic1 = im.get_ic_at_mass(73)
ic2 = im.get_ic_at_mass(57)
ic3 = im.get_ic_at_mass(55)
# create a list of ICs for passing to plot_ics()
ics = [ic, ic1, ic2, ic3]

# load list of peaks
# NB NB NB NB You must first run proc_save_peaks.py from
# this directory to generate the peak list
peak_list = load_peaks("output/peaks.bin") 



# Create a new display object, this time plot four ICs 
# and the TIC, as well as the peak list
display = Display()

display.plot_tic(tic, 'TIC')
display.plot_ics(ics, ['191','73','57','55'])
display.plot_peaks(peak_list, 'Peaks')
display.do_plotting('TIC, ICs, and PyMS Detected Peaks')
