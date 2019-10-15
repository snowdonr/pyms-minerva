 ****************************
GC-MS data derived objects
****************************

.. contents:: Table of Contents

In the raw GC-MS data, consecutive scans do not necessarily contain the same
mass per charge (mass) values. For data processing, it is often necessary to
convert the data to a matrix with a set number of masses and scans. In |pkgname2|,
the resulting object is called intensity matrix. In this chapter the methods
for converting the raw GC-MS data to an intensity matrix object are illustrated.

IntensityMatrix Object
=======================

The general scheme for converting raw mass values is to bin intensity values
based on the interval the corresponding mass belongs to. The general procedure
is as follows:

* Set the interval between bins, lower and upper bin boundaries.

* Calculate the number of bins to cover the range of all masses.

* Centre the first bin at the minimum mass found for all the raw data.

* Sum intensities whose masses are in a given bin.

A mass, :math:`m`, is considered to belong to a bin when :math:`c - l \le m < c + u`,
where :math:`c` is the centre of the bin, :math:`l` is the lower boundary and :math:`u` is
the upper boundary of the bin. The default bin interval is one with a lower
and upper boundary of :math:`\pm0.5`.

A function to bin masses to the nearest integer is also available. The default
bin interval is one with a lower boundary of :math:`-0.3` and upper boundary of
:math:`+0.7` (as per the NIST library).

Discussion of Binning Boundaries
------------------------------------

For any chemical element :math:`X`, let :math:`w(x)` be the atomic weight of :math:`X`, and

:math:`\delta(X) = \frac{w(X) - \{w(X)\}}{w(X)}`

where :math:`\{a\}` is the integer value of :math:`a` (rounded to the nearest integer).

For example, for hydrogen :math:`\delta(^1\rm{H}) = \frac{1.007825032 - 1}{1.007825032} = 0.0076`.
Similarly :math:`\delta(^{12}\rm{C}) = 0`, :math:`\delta(^{14}\rm{N}) = 0.00022`,
:math:`\delta(^{16}\rm{O}) = -0.00032`, etc.

Let also :math:`\Delta(X) = w(X) - \{w(x)\}`. Then :math:`-0.023 <\Delta(^{31}\rm{P}),
\Delta(^{28}\rm{Si}) < 0`.

Let a compound undergo GC-MS and let Y be one of it's fragments. If Y consists
of :math:`k_{1}`atoms of type :math:`X_{1}`,
:math:`k_{2}` atoms of type :math:`X_{2}`,.....,
:math:`k_{r}` atoms of type :math:`X_{r}`, then
:math:`\Delta(Y) = k_{1}*\Delta(X_{1}) + k_{2}*\Delta(X_{2}) + ....+ k_{r}*\Delta(X_{r})`.

The fragment will usually not contain more than 2 or 3 P or Si atoms and if it's
molecular weight is less than 550 it may not contain more than 35 O atoms, so
:math:`\Delta(Y) \geq -0.023*5 - 0.00051*35 = -0.133`.

On the other hand, of Y contains :math:`k` H atoms and :math:`m` N atoms, then :math:`\Delta(Y)
\leq k*0.00783 + m*0.00051`. Since for each two hydrogen atoms at least one
carbon (or heavier) atom is needed, giving the limit of no more than 80 hydrogen
atoms. Therefore in this case (i.e. H and C atoms only) :math:`\Delta(Y) \leq 80*0.00783 = 0.63`.
If carbon is replaced by any heavier atom, at least 2 hydrogen atoms
will be eliminated and :math:`\Delta(Y)` will become even smaller.

If the molecular weight of :math:`Y` does not exceed 550 (typically the largest mass
scanned for in a GC-MS setup) then :math:`\mathbf{-0.133 \leq \Delta(Y) \leq 0.63}`.
This means that if we set our binning boundaries to :math:`(-0.3, 0.7)` or :math:`(-0.2, 0.8)`
the opportunity for having a fragment whose molecular weight is very close to
the boundary is minimised.

Since the resolution of MS is at least 0.1 dalton, we may assume that it's error
does not exceed 0.05, and MS accuracy will not cause additional problems.

Build intensity matrix
------------------------

.. note:: This example is in :ref:`pyms-demo/30a <demo-30a>`

An intensity matrix on the raw GC-MS data can be built using the following
function. First the raw data is imported as before.

    >>> from pyms.GCMS.IO.JCAMP import JCAMP_reader
    >>> jcamp_file = "data/gc01_0812_066.jdx"
    >>> data = JCAMP_reader(jcamp_file)
     -> Reading JCAMP file 'data/gc01_0812_066.jdx'

Then the data can be converted to an intensity matrix using the functions
:meth:`build_intensity_matrix() <pyms.IntensityMatrix.build_intensity_matrix>` and
:meth:`build_intensity_matrix_i() <pyms.IntensityMatrix.build_intensity_matrix_i>`,
available in :py:meth:`pyms.IntensityMatrix <pyms.IntensityMatrix>`

The default operation of :meth:`build_intensity_matrix() <pyms.IntensityMatrix.build_intensity_matrix>`
is to use a bin interval of one and treat the masses as floating point numbers.
The default intensity matrix can be built as follows:

    >>> from pyms.IntensityMatrix import build_intensity_matrix
    >>> im = build_intensity_matrix(data)
    >>> im
    <pyms.IntensityMatrix.IntensityMatrix at 0x7f999eedabe0>

The size as the number of scans and the number of bins is returned by:

    >>> im.size
    (9865, 551)

There are 9865 scans and 551 bins in this example.

The raw masses have been binned into new mass units based on the minimum mass
in the raw data and the bin size. A list of the new masses can be obtained
as follows:

.. code-block:: python

    >>> masses = im.mass_list
    >>> masses[:10]
    [50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0]

The last command prints the first ten masses. The attributes
:attr:`im.min_mass <pyms.IntensityMatrix.IntensityMatrix.min_mass>` and
:attr:`im.max_mass <pyms.IntensityMatrix.IntensityMatrix.max_mass>`
return the minimum and maximum mass:

.. code-block:: python

    >>> im.min_mass
    50.0
    >>> im.max_mass
    600.0

It is also possible to search for a particular mass, by finding the index of
the binned mass closest to the desired mass. For example, the index of the
closest binned mass to a mass of 73.3 :math:`m/z` can be found by using the
methods :meth:`im.get_index_of_mass() <pyms.IntensityMatrix.get_index_of_mass>`:

    >>> im.get_index_of_mass(73.3)
    23

The value of the closest mass can be returned by the method
:meth:`im.get_mass_at_index() <pyms.IntensityMatrix.IntensityMatrix.get_mass_at_index>`:

    >>> im.get_mass_at_index(index)
    73.0

A mass of 73.0 is returned in this example.

Build intensity matrix parameters
-----------------------------------

.. note:: This example is in :ref:`pyms-demo/30b <demo-30b>`

The bin interval can be set to values other than one, and binning boundaries
can also be adjusted. In the example below, to fit the 0.5 bin interval, the
upper and lower boundaries are set to :math:`\pm` 0.25.

    >>> im = build_intensity_matrix(data, 0.5, 0.25, 0.25)

The size of the intensity matrix will reflect the change in the number of bins:

    >>> im.size
    (9865, 1101)

In this example there are 9865 scans (as before), but 1101 bins.

The index and binned mass of the mass closest to 73.3 should also reflect the
different binning.

    >>> index = im.get_index_of_mass(73.3)
    >>> im.get_mass_at_index(index)
    73.5

Build integer mass intensity matrix
------------------------------------

.. note:: This example is in :ref:`pyms-demo/30c <demo-30c>`

It is also possible to build an intensity matrix with integer masses and a bin
interval of one. The default range for the binning is -0.3 and +0.7 mass
units. The function is imported from :mod:`pyms.IntensityMatrix`:

    >>> from pyms.IntensityMatrix import build_intensity_matrix_i
    >>> im = build_intensity_matrix_i(data)

The masses are now integers.

    >>> index = im.get_index_of_mass(73.3)
    >>> im.get_mass_at_index(index)
    73

The lower and upper bounds can be adjusted by
:meth:`build_intensity_matrix_i(data, lower, upper) <pyms.IntensityMatrix.build_intensity_matrix_i>`

MassSpectrum object
=====================

.. note:: This example is in :ref:`pyms-demo/31 <demo-31>`

A :class:`~pyms.MassSpectrum.MassSpectrum` object contains
two attributes, :attr:`~pyms.MassSpectrum.MassSpectrum.mass_list` and :attr:`~pyms.MassSpectrum.MassSpectrum.mass_spec`, a list of mass values
and corresponding intensities, respectively.
A :class:`~pyms.MassSpectrum.MassSpectrum` is returned by the
:class:`~pyms.IntensityMatrix.IntensityMatrix` method
:meth:`get_ms_at_index(index) <pyms.IntensityMatrix.IntensityMatrix.get_ms_at_index()>`.

For example, the properties of the first MassSpectrum object of an
:class:`~pyms.IntensityMatrix.IntensityMatrix`, ``im``, can be obtained with;

    >>> ms = im.get_ms_at_index(0)
    <pyms.Spectrum.MassSpectrum at 0x7f999d210940>
    >>> len(ms)
    1101
    >>> len(ms.mass_list)
    1101
    >>> len(ms.mass_spec)
    1101

The length of all attributes should be the same.

IonChromatogram object
=======================

.. note:: This example is in :ref:`pyms-demo/31 <demo-31>`

An :class:`~pyms.IonChromatogram.IonChromatogram` object is a
one dimensional vector containing mass intensities as a function of
retention time. This can can be either :math:`m/z` channel intensities
(for example, the ion chromatogram at 73 :math:`m/z `), or cumulative
intensities over all measured :math:`m/z` (TIC).

An :class:`~pyms.IonChromatogram.IonChromatogram` for the
TIC and a given mass or index can be obtained as follows:

    >>> data.tic
    <pyms.IonChromatogram.IonChromatogram at 0x7f99a27bd320>
    >>> im.get_ic_at_index(0)
    <pyms.IonChromatogram.IonChromatogram at 0x7f999d220400>
    >>> im.get_ic_at_mass(73)
    <pyms.IonChromatogram.IonChromatogram at 0x7f999d246fd0>

This will return, respectively: the TIC; the ion chromatogram of the first
mass; and the ion chromatogram of the mass closest to 73.

An ion chromatogram object has a method
:meth:`is_tic() <pyms.IonChromatogram.IonChromatogram.is_tic>`
which returns ``True`` if the ion chromatogram is a TIC, ``False`` otherwise:

    >>> tic.is_tic()
    True
    >>> ic.is_tic()
    False

Writing IonChromatogram object to a file
--------------------------------------------

.. note:: This example is in :ref:`pyms-demo/31 <demo-31>`

The method :meth:`write() <pyms.IonChromatogram.IonChromatogram.write>`
of an :class:`~pyms.IonChromatogram.IonChromatogram`
object allows the ion chromatogram to be saved to a file:

    >>> tic.write("output/tic.dat", minutes=True)
    >>> ic.write("output/ic.dat", minutes=True)

The flag ``minutes=True`` indicates that retention time will be saved in minutes.
The ion chromatogram object saved with with the
:meth:`write() <pyms.IonChromatogram.IonChromatogram.write>`
method is a plain ASCII file which contains a pair of
(retention time, intensity) per line.

.. code-block:: bash

    $ head tic.dat
    5.0930 2.222021e+07
    5.0993 2.212489e+07
    5.1056 2.208650e+07
    5.1118 2.208815e+07
    5.1181 2.200635e+07
    5.1243 2.200326e+07
    5.1306 2.202363e+07
    5.1368 2.198357e+07
    5.1431 2.197408e+07
    5.1493 2.193351e+07

Saving data
=============

.. note:: This example is in :ref:`pyms-demo/32 <demo-32>`

A matrix of intensity values can be saved to a file with the function
:meth:`save_data() <pyms.Utils.IO.save_data>`
from :meth:`pyms.Utils.IO <pyms.Utils.IO>`. A matrix of intensity values can
be returned from an :class:`~pyms.IntensityMatrix.IntensityMatrix`
with the method
:attr:`~pyms.IntensityMatrix.IntensityMatrix.intensity_array`.
For example,

    >>> from pyms.Utils.IO import save_data
    >>> mat = im.intensity_array
    array([[22128.,     0., 10221., ...,     0.,   470.,     0.],
           [22040.,     0., 10335., ...,   408.,     0.,   404.],
           [21320.,     0., 10133., ...,   492.,     0.,   422.],
           ...,
           [    0.,     0.,     0., ...,     0.,     0.,     0.],
           [    0.,     0.,     0., ...,     0.,     0.,     0.],
           [    0.,     0.,     0., ...,     0.,     0.,     0.]])
    >>> save_data("output/im.dat", mat)

It is also possible to save the list of masses (from
:attr:`im.mass_list <pyms.IntensityMatrix.IntensityMatrix.mass_list>`
and the list of retention times (from
:attr:`im.time_list <pyms.IntensityMatrix.IntensityMatrix.time_list>`
using the :meth:`save_data() <pyms.Utils.IO.save_data>` function.
For convenience, the intensity values, mass list and time list,
can be saved with the method
:meth:`export_ascii() <pyms.IntensityMatrix.IntensityMatrix.export_ascii>`.
For example,

    >>> im.export_ascii("output/data")

will create ``data.im.dat``, ``data.rt.dat`` and ``data.mz.dat``, where these
are the intensity matrix, retention time vector, and :math:`m/z` vector. By default
the data is saved as space separated data with a ".dat" extension. It is
also possible to save the data as comma separated data with a ".csv"
extension with the command:

    >>> im.export_ascii("output/data", "csv")

Additionally, the entire :class:`IntensityMatrix <pyms.IntensityMatrix.IntensityMatrix>`
can be exported to LECO CSV format.

    >>> im.export_leco_csv("output/data_leco.csv")

This facility is useful for import into other analytical software packages. The format has a header line specifying the column heading information as:

.. code-block:: text

    scan, retention time, mass1, mass2, ...

and then each row as the intensity data.


Importing ASCII data
====================

.. note:: This example is in :ref:`pyms-demo/32 <demo-32>`

The LECO CSV format data can be imported directly into an `~pyms.IntensityMatrix.IntensityMatrix` object. The data must follow the format outlined above. For example, the file saved above can be read and compared to the original:

.. code-block:: python

    >>> from pyms.IntensityMatrix import IntensityMatrix
    >>> iim = IntensityMatrix([0],[0],[[0]])
    >>> iim.import_leco_csv("output/data_leco.csv")
    >>> im.size
    >>> iim.size

The line :meth:`IntensityMatrix([0],[0],[[0]]) <pyms.IntensityMatrix.IntensityMatrix>` is required to create an empty :class:`~pyms.IntensityMatrix.IntensityMatrix` object.
sot