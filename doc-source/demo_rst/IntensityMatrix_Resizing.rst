Example: IntensityMatrix Resizing
=================================

Once an IntensityMatrix has been constructed from the raw GC-MS data,
the entries of the matrix can be modified. These modifications can
operate on the entire matrix, or individual masses or scans.

First, setup the paths to the datafiles and the output directory, then
import JCAMP_reader and build_intensity_matrix.

.. nbinput:: ipython3
    :execution-count: 1

    import pathlib
    data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
    # Change this if the data files are stored in a different location

    output_directory = pathlib.Path(".").resolve() / "output"

    from pyms.GCMS.IO.JCAMP import JCAMP_reader
    from pyms.IntensityMatrix import build_intensity_matrix

Read the raw data files and create the IntensityMatrix.

.. nbinput:: ipython3
    :execution-count: 2

    jcamp_file = data_directory / "gc01_0812_066.jdx"
    data = JCAMP_reader(jcamp_file)
    im = build_intensity_matrix(data)


.. parsed-literal::

     -> Reading JCAMP file '/home/vagrant/PyMassSpec/pyms-data/gc01_0812_066.jdx'


Retention time range
--------------------

A basic operation on the GC-MS data is to select a specific time range
for processing. In PyMassSpec, any data outside the chosen time range is
discarded. The :meth:`trim() <pyms.GCMS.Class.GCMS_data.trim>` method operates on the raw data, so any
subsequent processing only refers to the trimmed data.

The data can be trimmed to specific scans:

.. nbinput:: ipython3
    :execution-count: 3

    data.trim(1000, 2000)
    data.info()


.. parsed-literal::

    Trimming data to between 1000 and 2001 scans
     Data retention time range: 11.342 min -- 17.604 min
     Time step: 0.375 s (std=0.000 s)
     Number of scans: 1002
     Minimum m/z measured: 50.100
     Maximum m/z measured: 467.100
     Mean number of m/z values per scan: 57
     Median number of m/z values per scan: 44


or specific retention times (in ``seconds`` or ``minutes``):

.. nbinput:: ipython3
    :execution-count: 4

    data.trim("700s", "15m")
    data.info()


.. parsed-literal::

    Trimming data to between 54 and 587 scans
     Data retention time range: 11.674 min -- 15.008 min
     Time step: 0.375 s (std=0.000 s)
     Number of scans: 534
     Minimum m/z measured: 50.100
     Maximum m/z measured: 395.200
     Mean number of m/z values per scan: 59
     Median number of m/z values per scan: 47


Mass Spectrum range and entries
-------------------------------

An :class:`~pyms.IntensityMatrix.IntensityMatrix` object has a set mass range and interval that is
derived from the data at the time of building the intensity matrix. The
range of mass values can be cropped. This is done, primarily, to ensure
that the range of masses used are consistent when comparing samples.

The mass range of the intensity matrix can be “cropped” to a new
(smaller) range as follows:

.. nbinput:: ipython3
    :execution-count: 5

    im.crop_mass(60, 400)
    im.min_mass




.. parsed-literal::

    60.0



.. nbinput:: ipython3
    :execution-count: 6

    im.max_mass




.. parsed-literal::

    400.0



It is also possible to set all intensities for a given mass to zero.
This is useful for ignoring masses associated with sample preparation.
The mass can be “nulled” with:

.. nbinput:: ipython3
    :execution-count: 7

    im.null_mass(73)
    sum(im.get_ic_at_mass(73).intensity_array)




.. parsed-literal::

    0.0



As expected, the sum of the intensity array is ``0``
