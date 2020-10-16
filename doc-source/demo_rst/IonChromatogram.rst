Example: IonChromatogram Objects
================================

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


An :class:`~pyms.IonChromatogram.IonChromatogram` object is a one dimensional vector containing
mass intensities as a function of retention time. This can can be either
:math:`m/z` channel intensities (for example, the ion chromatogram at 73
:math:`m/z`), or cumulative intensities over all measured :math:`m/z` (TIC).

An :class:`~pyms.IonChromatogram.IonChromatogram` object for the TIC can be obtained as follows:

.. nbinput:: ipython3
    :execution-count: 3

    data.tic




.. parsed-literal::

    <pyms.IonChromatogram.IonChromatogram at 0x7fe8e5cfc5c0>



The :class:`~pyms.IonChromatogram.IonChromatogram` at index 0 can be obtained with:

.. nbinput:: ipython3
    :execution-count: 4

    im.get_ic_at_index(0)




.. parsed-literal::

    <pyms.IonChromatogram.IonChromatogram at 0x7fe9044238d0>



The :class:`~pyms.IonChromatogram.IonChromatogram` for the closest mass to 73 can be obtained with:

.. nbinput:: ipython3
    :execution-count: 5

    im.get_ic_at_mass(73)




.. parsed-literal::

    <pyms.IonChromatogram.IonChromatogram at 0x7fe90442a9e8>



An ion chromatogram object has a method :meth:`is_tic() <pyms.IonChromatogram.IonChromatogram.is_tic>` which returns
``True`` if the ion chromatogram is a TIC, ``False`` otherwise.

.. nbinput:: ipython3
    :execution-count: 6

    data.tic.is_tic()




.. parsed-literal::

    True



.. nbinput:: ipython3
    :execution-count: 7

    im.get_ic_at_mass(73).is_tic()




.. parsed-literal::

    False
