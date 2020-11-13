Example: Peak Detection
-----------------------

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

Read the raw data file and build the :class:`~pyms.IntensityMatrix.IntensityMatrix`.

.. nbinput:: ipython3
    :execution-count: 2

    jcamp_file = data_directory / "gc01_0812_066.jdx"
    data = JCAMP_reader(jcamp_file)
    im = build_intensity_matrix(data)


.. parsed-literal::

     -> Reading JCAMP file '/home/vagrant/PyMassSpec/pyms-data/gc01_0812_066.jdx'


Preprocess the data (Savitzky-Golay smoothing and Tophat baseline
detection).

.. nbinput:: ipython3
    :execution-count: 3

    from pyms.Noise.SavitzkyGolay import savitzky_golay
    from pyms.TopHat import tophat

    n_scan, n_mz = im.size

    for ii in range(n_mz):
        ic = im.get_ic_at_index(ii)
        ic_smooth = savitzky_golay(ic)
        ic_bc = tophat(ic_smooth, struct="1.5m")
        im.set_ic_at_index(ii, ic_bc)

Now the Biller and Biemann based technique can be applied to detect
peaks.

.. nbinput:: ipython3
    :execution-count: 4

    from pyms.BillerBiemann import BillerBiemann
    peak_list = BillerBiemann(im)
    peak_list[:10]




.. parsed-literal::

    [<pyms.Peak.Class.Peak at 0x7f8deca882b0>,
     <pyms.Peak.Class.Peak at 0x7f8deca88550>,
     <pyms.Peak.Class.Peak at 0x7f8dce908160>,
     <pyms.Peak.Class.Peak at 0x7f8df3bc5c50>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a043c8>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a04588>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a045f8>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a04668>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a046d8>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a04748>]



.. nbinput:: ipython3
    :execution-count: 5

    len(peak_list)





.. parsed-literal::

    9845



Note that this is nearly as many peaks as there are scans in the data
(9865 scans). This is due to noise and the simplicity of the technique.

The number of detected peaks can be constrained by the selection of
better parameters. Parameters can be determined by counting the number
of points across a peak, and examining where peaks are found. For
example, the peak list can be found with the parameters of a window of 9
points and by combining 2 neighbouring scans if they apex next to each
other:

.. nbinput:: ipython3
    :execution-count: 6

    peak_list = BillerBiemann(im, points=9, scans=2)
    peak_list[:10]




.. parsed-literal::

    [<pyms.Peak.Class.Peak at 0x7f8dae545be0>,
     <pyms.Peak.Class.Peak at 0x7f8dae545c18>,
     <pyms.Peak.Class.Peak at 0x7f8dae545c88>,
     <pyms.Peak.Class.Peak at 0x7f8dae545cf8>,
     <pyms.Peak.Class.Peak at 0x7f8dae545d68>,
     <pyms.Peak.Class.Peak at 0x7f8dae545dd8>,
     <pyms.Peak.Class.Peak at 0x7f8dae545e48>,
     <pyms.Peak.Class.Peak at 0x7f8dae545eb8>,
     <pyms.Peak.Class.Peak at 0x7f8dae545f28>,
     <pyms.Peak.Class.Peak at 0x7f8dae545f98>]



.. nbinput:: ipython3
    :execution-count: 7

    len(peak_list)





.. parsed-literal::

    3695



The number of detected peaks has been reduced, but there are still many
more than would be expected from the sample. Functions to filter the
peak list are covered in the next example.

Example: Peak List Filtering
----------------------------

There are two functions to filter the list of Peak objects.

The first, :meth:`rel_threshold() <pyms.BillerBiemann.rel_threshold>` modifies the mass spectrum stored in each
peak so any intensity that is less than a given percentage of the
maximum intensity for the peak is removed.

The second, :meth:`num_ions_threshold() <pyms.BillerBiemann.num_ions_threshold>`, removes any peak that has less than
a given number of ions above a given threshold.

Once the peak list has been constructed, the filters can be applied as
follows:

.. nbinput:: ipython3
    :execution-count: 8

    from pyms.BillerBiemann import rel_threshold, num_ions_threshold
    pl = rel_threshold(peak_list, percent=2)
    pl[:10]




.. parsed-literal::

    [<pyms.Peak.Class.Peak at 0x7f8dc3a045f8>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a04630>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a04748>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a047b8>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a048d0>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a049e8>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a04a20>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a04b38>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a04b70>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a04c88>]



.. nbinput:: ipython3
    :execution-count: 9

    new_peak_list = num_ions_threshold(pl, n=3, cutoff=10000)
    new_peak_list[:10]




.. parsed-literal::

    [<pyms.Peak.Class.Peak at 0x7f8deca8e128>,
     <pyms.Peak.Class.Peak at 0x7f8deca8e1d0>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a04780>,
     <pyms.Peak.Class.Peak at 0x7f8dc3a04550>,
     <pyms.Peak.Class.Peak at 0x7f8dbb3cf3c8>,
     <pyms.Peak.Class.Peak at 0x7f8dbb3cf048>,
     <pyms.Peak.Class.Peak at 0x7f8dbb3cf4a8>,
     <pyms.Peak.Class.Peak at 0x7f8dbb3cf550>,
     <pyms.Peak.Class.Peak at 0x7f8dbb3cf5f8>,
     <pyms.Peak.Class.Peak at 0x7f8dbb3cf6a0>]



.. nbinput:: ipython3
    :execution-count: 10

    len(new_peak_list)




.. parsed-literal::

    146



The number of detected peaks is now more realistic of what would be
expected in the test sample.
