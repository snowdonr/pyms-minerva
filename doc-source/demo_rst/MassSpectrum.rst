Example: MassSpectrum Objects
=============================

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


A :class:`~pyms.Spectrum.MassSpectrum` object contains two attributes, :attr:`~pyms.Spectrum.MassSpectrum.mass_list` and
:attr:`~pyms.Spectrum.MassSpectrum.intensity_list`, a list of mass values and corresponding intensities,
respectively. A :class:`~pyms.Spectrum.MassSpectrum` is returned by the :class:`~pyms.IntensityMatrix.IntensityMatrix`
method :meth:`get_ms_at_index(index) <pyms.IntensityMatrix.IntensityMatrix.get_ms_at_index()>`.

For example, the properties of the first :class:`~pyms.Spectrum.MassSpectrum` object can be
obtained as follows:

.. nbinput:: ipython3
    :execution-count: 3

    ms = im.get_ms_at_index(0)

    ms




.. parsed-literal::

    <pyms.Spectrum.MassSpectrum at 0x7f8f4c529860>



.. nbinput:: ipython3
    :execution-count: 4

    len(ms)




.. parsed-literal::

    551



.. nbinput:: ipython3
    :execution-count: 5

    len(ms.mass_list)




.. parsed-literal::

    551



.. nbinput:: ipython3
    :execution-count: 6

    len(ms.intensity_list)




.. parsed-literal::

    551



The length of all attributes should be the same.
