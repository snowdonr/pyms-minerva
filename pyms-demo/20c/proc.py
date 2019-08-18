"""proc.py
"""

from pyms.GCMS.IO.ANDI.Function import ANDI_reader

# read the raw data
andi_file = "data/gc01_0812_066.cdf"
data = ANDI_reader(andi_file)

# print info
data.info()

# write data to output file. This will create
# two ascii data tables, data.I.csv and data.mz.csv
# with intensities and m/z values
data.write("output/data")

