"""proc.py
"""

#from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.GCMS.Function import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Baseline.TopHat import tophat
from pyms.Utils.IO import dump_object

# read the raw data as a GCMS_data object
#jcamp_file = "data/gc01_0812_066.jdx"
#data = JCAMP_reader(jcamp_file)
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

# Build intensity matrix, required before accessing any intensity matrix
# methods. Use defaults, float masses with interval (bin size) of one from
# min mass
print("building intensity matrix with bin interval=1")
im = build_intensity_matrix_i(data)

off = im.min_mass
print(off)
n_scan, n_mz = im.size

# process data
print(" Processing ICs")

for ii in range(n_mz):
    print("Working on IC#", ii+1)
    ic = im.get_ic_at_index(ii)
    if (ii + off) in [319, 205, 160, 217]:
        ic.write(f"output/ic-raw-{ii + off:d}.dat")
    ic_smooth = savitzky_golay(ic)
    ic_bc = tophat(ic_smooth, struct="1.5m")
    if (ii + off) in [319, 205, 160, 217]:
        ic_bc.write(f"output/ic-flt-{ii + off:d}.dat")
    im.set_ic_at_index(ii, ic_bc)

# save the pre-processed intensity matrix
print(" Dumping the pre-processed intensity matrix")
dump_object(im, "output/im-proc.dump")
