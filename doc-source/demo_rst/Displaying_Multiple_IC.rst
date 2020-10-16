Example: Displaying Multiple IonChromatogram Objects
====================================================

Multiple :class:`~pyms.IonChromatogram.IonChromatogram` objects can be plotted on the same figure.

To start, load a datafile and create an :class:`~pyms.IntensityMatrix.IntensityMatrix` as before.

.. nbinput:: ipython3
    :execution-count: 1

    import pathlib
    data_directory = pathlib.Path(".").resolve().parent.parent / "pyms-data"
    # Change this if the data files are stored in a different location

    output_directory = pathlib.Path(".").resolve() / "output"

    from pyms.GCMS.IO.JCAMP import JCAMP_reader
    from pyms.IntensityMatrix import build_intensity_matrix_i

    jcamp_file = data_directory / "gc01_0812_066.jdx"
    data = JCAMP_reader(jcamp_file)
    tic = data.tic
    im = build_intensity_matrix_i(data)


.. parsed-literal::

     -> Reading JCAMP file '/home/vagrant/PyMassSpec/pyms-data/gc01_0812_066.jdx'


Extract the desired IonChromatograms from the :class:`~pyms.IntensityMatrix.IntensityMatrix` .

.. nbinput:: ipython3
    :execution-count: 2

    ic73 = im.get_ic_at_mass(73)
    ic147 = im.get_ic_at_mass(147)

Import matplotlib and the :py:meth:`plot_ic() <pyms.Display.plot_ic>` function, create a subplot, and
plot the ICs on the chart:

.. nbinput:: ipython3
    :execution-count: 3

    import matplotlib.pyplot as plt
    from pyms.Display import plot_ic

    %matplotlib inline
    # Change to `notebook` for an interactive view

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))

    # Plot the ICs
    plot_ic(ax, tic, label="TIC")
    plot_ic(ax, ic73, label="m/z 73")
    plot_ic(ax, ic147, label="m/z 147")

    # Set the title
    ax.set_title("TIC and ICs for m/z = 73 & 147")

    # Add the legend
    plt.legend()

    plt.show()



.. image:: graphics/Displaying_Multiple_IC_output_5_0.png
