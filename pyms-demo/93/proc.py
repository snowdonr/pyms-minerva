from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.GCMS.Function import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Baseline.TopHat import tophat
from pyms.Display.Class import Display
from pyms.Peak.Function import peak_sum_area
#from pyms.Peak.IO import store_peaks
from pyms.Deconvolution.BillerBiemann.Function import BillerBiemann, \
    rel_threshold, num_ions_threshold
from pyms.Simulator.Function import gcms_sim, add_gaussc_noise


 # read in raw data
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

data.trim(4101, 4350)

# Build Intensity Matrix
real_im = build_intensity_matrix_i(data)

n_scan, n_mz = real_im.get_size()

 # perform necessary pre filtering
for ii in range(n_mz):
    ic = real_im.get_ic_at_index(ii)
    ic_smooth = savitzky_golay(ic)
    ic_bc = tophat(ic_smooth, struct="1.5m")
    real_im.set_ic_at_index(ii, ic_bc)
    
    
 # Detect Peaks
peak_list = BillerBiemann(real_im, points=3, scans=2)

print("Number of peaks found in real data: ", len(peak_list))

######### Filter peaks###############
# Filter the peak list,
# first by removing all intensities in a peak less than a given relative
# threshold,
# then by removing all peaks that have less than a given number of ions above
# a given value

# Parameters
# percentage ratio of ion intensity to max ion intensity
r = 1

# minimum number of ions, n
n = 3
# greater than or equal to threshold, t
t = 10000

# trim by relative intensity
pl = rel_threshold(peak_list, r)

# trim by threshold
real_peak_list = num_ions_threshold(pl, n, t)
print("Number of filtered peaks in real data: ", len(real_peak_list))

# Set the peak areas
for peak in real_peak_list:
    area = peak_sum_area(real_im, peak)
    peak.set_area(area)

# real_peak_list is PyMS' best guess at the true peak list

################## Run Simulator ######################
# Simulator takes a peak list, time_list and mass_list
# and returns an IntensityMatrix object.
# The mass_list and time_list are the same for the real 
# data and the simulated data.

time_list = real_im.get_time_list()
mass_list = real_im.get_mass_list()

sim_im = gcms_sim(time_list, mass_list, real_peak_list)
# sim_im is an IntensityMatrix object 

# Now add noise to the simulated intensity matrix object
scale = 1000
add_gaussc_noise(sim_im, scale)

### Now display the ics from the simulated data
ics = []
for i in range(n_mz):
    ics.append(sim_im.get_ic_at_index(i))

display = Display()
display.plot_ics(ics)
display.do_plotting('ICs of Simulated Data with gaussian noise (constant scale)')
