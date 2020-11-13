Example: Building an Intensity Matrix
-------------------------------------

First, setup the paths to the datafiles and the output directory, then
import JCAMP_reader.

.. nbinput:: ipython3
    :execution-count: 1

    import pathlib
    data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
    # Change this if the data files are stored in a different location

    output_directory = pathlib.Path(".").resolve() / "output"

    from pyms.GCMS.IO.JCAMP import JCAMP_reader

Read the raw data files.

.. nbinput:: ipython3
    :execution-count: 2

    jcamp_file = data_directory / "gc01_0812_066.jdx"
    data = JCAMP_reader(jcamp_file)
    data


.. parsed-literal::

     -> Reading JCAMP file '/home/vagrant/PyMassSpec/pyms-data/gc01_0812_066.jdx'




.. parsed-literal::

    <GCMS_data(305.582 - 4007.722 seconds, time step 0.3753183292781833, 9865 scans)>



Then the data can be converted to an :class:`~pyms.IntensityMatrix.IntensityMatrix` using the
function :meth:`build_intensity_matrix() <pyms.IntensityMatrix.build_intensity_matrix>` from :mod:`pyms.IntensityMatrix <pyms.IntensityMatrix>`.

The default operation of :meth:`build_intensity_matrix() <pyms.IntensityMatrix.build_intensity_matrix>` is to use a bin
interval of one and treat the masses as floating point numbers. The
default intensity matrix can be built as follows:

.. nbinput:: ipython3
    :execution-count: 3

    from pyms.IntensityMatrix import build_intensity_matrix

    im = build_intensity_matrix(data)

    im




.. parsed-literal::

    <pyms.IntensityMatrix.IntensityMatrix at 0x7f31d8b12860>



The size as the number of scans and the number of bins can be returned
with:

.. nbinput:: ipython3
    :execution-count: 4

    im.size




.. parsed-literal::

    (9865, 551)



There are 9865 scans and 551 bins in this example.

The raw masses have been binned into new mass units based on the minimum
mass in the raw data and the bin size. A list of the first ten new
masses can be obtained as follows:

.. nbinput:: ipython3
    :execution-count: 5

    im.mass_list[:10]




.. parsed-literal::

    [50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0]



The attributes :attr:`im.min_mass <pyms.IntensityMatrix.IntensityMatrix.min_mass>` and :attr:`im.max_mass <pyms.IntensityMatrix.IntensityMatrix.max_mass>` return the minimum
and maximum mass:

.. nbinput:: ipython3
    :execution-count: 6

    im.min_mass




.. parsed-literal::

    50.0



.. nbinput:: ipython3
    :execution-count: 7

    im.max_mass




.. parsed-literal::

    600.0



It is also possible to search for a particular mass, by finding the
index of the binned mass closest to the desired mass. For example, the
index of the closest binned mass to a mass of 73.3 :math:`m/z` can be found
by using the methods :meth:`im.get_index_of_mass() <pyms.IntensityMatrix.IntensityMatrix.get_index_of_mass>`:

.. nbinput:: ipython3
    :execution-count: 8

    index = im.get_index_of_mass(73.3)

    index




.. parsed-literal::

    23



The value of the closest mass can be returned by the method
:meth:`im.get_mass_at_index() <pyms.IntensityMatrix.IntensityMatrix.get_mass_at_index>`:

.. nbinput:: ipython3
    :execution-count: 9

    im.get_mass_at_index(index)




.. parsed-literal::

    73.0



A mass of 73.0 is returned in this example.

Build intensity matrix parameters
---------------------------------

The bin interval can be set to values other than one, and binning
boundaries can also be adjusted. In the example below, to fit the 0.5
bin interval, the upper and lower boundaries are set to ± 0.25.

.. nbinput:: ipython3
    :execution-count: 10

    im = build_intensity_matrix(data, 0.5, 0.25, 0.25)

    im




.. parsed-literal::

    <pyms.IntensityMatrix.IntensityMatrix at 0x7f31d8b8d710>



The size of the intensity matrix will reflect the change in the number
of bins:

.. nbinput:: ipython3
    :execution-count: 11

    im.size




.. parsed-literal::

    (9865, 1101)



.. nbinput:: ipython3
    :execution-count: 12

    im.mass_list[:10]




.. parsed-literal::

    [50.0, 50.5, 51.0, 51.5, 52.0, 52.5, 53.0, 53.5, 54.0, 54.5]



In this example there are 9865 scans (as before), but 1101 bins.

The index and binned mass of the mass closest to 73.3 should also
reflect the different binning.

.. nbinput:: ipython3
    :execution-count: 13

    index = im.get_index_of_mass(73.3)

    index




.. parsed-literal::

    47



.. nbinput:: ipython3
    :execution-count: 14

    im.get_mass_at_index(index)




.. parsed-literal::

    73.5



Build integer mass intensity matrix
-----------------------------------

It is also possible to build an intensity matrix with integer masses and
a bin interval of one using :meth:`build_intensity_matrix_i() <pyms.IntensityMatrix.build_intensity_matrix_i>`. The default
range for the binning is -0.3 and +0.7 mass units. The function is
imported from :mod:`pyms.IntensityMatrix <pyms.IntensityMatrix>`:

.. nbinput:: ipython3
    :execution-count: 15

    from pyms.IntensityMatrix import build_intensity_matrix_i

    im = build_intensity_matrix_i(data)
    im




.. parsed-literal::

    <pyms.IntensityMatrix.IntensityMatrix at 0x7f31d8b121d0>



.. nbinput:: ipython3
    :execution-count: 16

    im.size




.. parsed-literal::

    (9865, 551)



.. nbinput:: ipython3
    :execution-count: 17

    im.mass_list[:10]




.. parsed-literal::

    [50, 51, 52, 53, 54, 55, 56, 57, 58, 59]



The masses are now integers.

.. nbinput:: ipython3
    :execution-count: 18

    index = im.get_index_of_mass(73.3)
    index




.. parsed-literal::

    23



.. nbinput:: ipython3
    :execution-count: 19

    im.get_mass_at_index(index)




.. parsed-literal::

    73



The lower and upper bounds can be adjusted with
:meth:`build_intensity_matrix_i(data, lower, upper) <pyms.IntensityMatrix.build_intensity_matrix_i>`.
