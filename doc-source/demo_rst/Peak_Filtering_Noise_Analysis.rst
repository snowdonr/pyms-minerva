Noise analysis for peak filtering
=================================

In the previous example the cutoff parameter for peak filtering was set
by the user. This can work well for individual data files, but can cause
problems when applied to large experiments with many individual data
files. Where experimental conditions have changed slightly between
experimental runs, the ion intensity over the GC-MS run may also change.
This means that an inflexible cutoff value can work for some data files,
while excluding too many, or including too many peaks in other files.

An alternative to manually setting the value for cutoff is to use the
:meth:`window_analyzer() <pyms.Noise.Analysis.window_analyzer>` function. This function examines a Total Ion
Chromatogram (TIC) and computes a value for the median absolute
deviation in troughs between peaks. This gives an approximate threshold
value above which false peaks from noise should be filtered out.

First, build the Peak list as before

.. nbinput:: ipython3
    :execution-count: 1

    import pathlib
    data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
    # Change this if the data files are stored in a different location

    output_directory = pathlib.Path(".").resolve() / "output"

    from pyms.GCMS.IO.JCAMP import JCAMP_reader
    from pyms.IntensityMatrix import build_intensity_matrix
    from pyms.Noise.SavitzkyGolay import savitzky_golay
    from pyms.TopHat import tophat
    from pyms.BillerBiemann import BillerBiemann

    jcamp_file = data_directory / "gc01_0812_066.jdx"
    data = JCAMP_reader(jcamp_file)
    im = build_intensity_matrix(data)

    n_scan, n_mz = im.size

    for ii in range(n_mz):
        ic = im.get_ic_at_index(ii)
        ic_smooth = savitzky_golay(ic)
        ic_bc = tophat(ic_smooth, struct="1.5m")
        im.set_ic_at_index(ii, ic_bc)

    peak_list = BillerBiemann(im, points=9, scans=2)


.. parsed-literal::

     -> Reading JCAMP file '/home/vagrant/PyMassSpec/pyms-data/gc01_0812_066.jdx'


Compute the noise value.

.. nbinput:: ipython3
    :execution-count: 2

    from pyms.Noise.Analysis import window_analyzer

    tic = data.tic

    noise_level = window_analyzer(tic)
    noise_level




.. parsed-literal::

    432.1719792438844



Filter the Peak List using this noise value as the cutoff.

.. nbinput:: ipython3
    :execution-count: 3

    from pyms.BillerBiemann import num_ions_threshold
    filtered_peak_list = num_ions_threshold(peak_list, n=3, cutoff=noise_level)
    filtered_peak_list[:10]




.. parsed-literal::

    [<pyms.Peak.Class.Peak at 0x7f4a9864f128>,
     <pyms.Peak.Class.Peak at 0x7f4a9864f2e8>,
     <pyms.Peak.Class.Peak at 0x7f4a9864f320>,
     <pyms.Peak.Class.Peak at 0x7f4a9864f3c8>,
     <pyms.Peak.Class.Peak at 0x7f4a9864f518>,
     <pyms.Peak.Class.Peak at 0x7f4a9864f4a8>,
     <pyms.Peak.Class.Peak at 0x7f4a9864f6a0>,
     <pyms.Peak.Class.Peak at 0x7f4a9864f748>,
     <pyms.Peak.Class.Peak at 0x7f4a9864f7f0>,
     <pyms.Peak.Class.Peak at 0x7f4a9864f898>]



.. nbinput:: ipython3
    :execution-count: 4

    len(filtered_peak_list)





.. parsed-literal::

    612
