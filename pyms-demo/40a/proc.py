"""proc.py
"""

from pyms.GCMS.IO.ANDI.Function import ANDI_reader

# read the raw data
andi_file = "data/gc01_0812_066.cdf"

data = ANDI_reader(andi_file)

# info about raw data
data.info()

# trim data between scans 1000 and 2000
data.trim(1000, 2000)

# info about trimmed raw data
data.info()

# reload
data = ANDI_reader(andi_file)

# trim data between retention times, 6.5 minutes to 21 minutes
data.trim("6.5m", "21m")

# info about trimmed raw data
data.info()
