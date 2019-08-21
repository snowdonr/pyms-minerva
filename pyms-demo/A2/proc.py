"""proc.py

This example demonstrates processing of a GC-TOF (Leco Pegasus-ChromaTOF)
generated dataset.

GC-TOF data is made up of nearly 10 scans per second. As a result of this,
the peak detection window in the Biller Biemann algorithm has a higher value
as compared to the value used to process GC-Quad data.

Due to the same reason, the value of the number of scans in the Biller
Biemann algorithm has a higher value as compared to the value used to
process GC-Quad data.
"""

from numpy import *
 
from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.GCMS.Function import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Baseline.TopHat import tophat

from pyms.Peak.Function import peak_sum_area
from pyms.Display.Class import Display
 
from pyms.Deconvolution.BillerBiemann.Function import BillerBiemann, \
    rel_threshold, num_ions_threshold
    
 # read in raw data
andi_file = "data/MM-10.0_1_no_processing.cdf"
data = ANDI_reader(andi_file)

# Build Intensity Matrix
im = build_intensity_matrix_i(data)
n_scan, n_mz = im.get_size()

 # perform necessary pre filtering
for ii in range(n_mz):
    ic = im.get_ic_at_index(ii)
    ic_smooth = savitzky_golay(ic)
    ic_bc = tophat(ic_smooth, struct="1.5m")
    im.set_ic_at_index(ii, ic_bc)
    
 # Detect Peaks
peak_list = BillerBiemann(im, points=15, scans=3)
print("Number of peaks found: ", len(peak_list))

######### Filter peaks###############
# Filter the peak list,
# first by removing all intensities in a peak less than a given relative
# threshold,
# then by removing all peaks that have less than a given number of ions above
# a given value

# Parameters
# percentage ratio of ion intensity to max ion intensity
r = 2
# minimum number of ions, n
n = 2
# greater than or equal to threshold, t
t = 4000
# trim by relative intensity
pl = rel_threshold(peak_list, r)

# trim by threshold
new_peak_list = num_ions_threshold(pl, n, t)

print("Number of filtered peaks: ", len(new_peak_list))
print("Peak areas")
print("UID, RT, height, area")
for peak in new_peak_list:
    rt = peak.get_rt()
    
    # determine and set area
    area = peak_sum_area(im, peak)
    peak.set_area(area)

    # print some details
    UID = peak.get_UID()
    # height as sum of the intensities of the apexing ions
    height = sum(peak.get_mass_spectrum().mass_spec)
    print(UID + ", %.2f, %.2f, %.2f" % (rt, height, peak.get_area()))

# TIC from raw data
tic = data.get_tic()
# baseline correction for TIC
tic_bc = tophat(tic, struct="1.5m")

# Get Ion Chromatograms for all m/z channels
n_mz = len(im.get_mass_list())
ic = []

for m in range(n_mz):
    ic.append(im.get_ic_at_index(m))

# Create a new display object, this time plot the ICs 
# and the TIC, as well as the peak list
display = Display()
display.plot_tic(tic_bc, 'TIC BC')
display.plot_ics(ic)
display.plot_peaks(new_peak_list, 'Peaks')
display.do_plotting('TIC, and PyMassSpec Detected Peaks')
