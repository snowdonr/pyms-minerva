***************************************
Peak alignment by dynamic programming
***************************************

.. contents:: :local: Table of Contents

PyMS provides functions to align GC-MS peaks by dynamic programming [1]_.
The peak alignment by dynamic programming uses both peak apex retention time and mass spectra.
This information is determined from the raw GC-MS data by applying a series of processing steps to produce data that can then be aligned and used for statistical analysis.
The details are described in this chapter.

Preparation of multiple experiments for peak alignment by dynamic programming
=============================================================================

Creating an Experiment
-------------------------

.. note:: This example is in :ref:`pyms-demo/60 <demo-60>`

Before aligning peaks from multiple experiments, the peak objects need to be created and encapsulated into |pkgname| :class:`~pyms.Experiment.Experiment` objects.
During this process it is often useful to pre-process the peaks in some way, for example to null certain m/z channels and/or to select a certain retention time range.

To capture the data and related information prior to peak alignment, an :meth:`~pyms.Experiment.Experiment` object is used. The :meth:`~pyms.Experiment.Experiment` object is defined in :mod:`pyms.Experiment`.

The procedure is to proceed as described in the previous chapter.
Namely:
# read a file;
# bin the data into fixed mass values;
# smooth the data;
# remove the baseline;
# deconvolute peaks;
# filter the peaks;
# set the mass range;
# remove uninformative ions; and
# estimate peak areas.

The process is given in the following program listing:

.. literalinclude:: ../pyms-demo/60/proc.py
   :linenos:
   :lines: 5-13, 17-60

The resulting list of peaks can now be stored as an
:meth:`~pyms.Experiment.Experiment` object.

.. literalinclude:: ../pyms-demo/60/proc.py
   :linenos:
   :lines: 14-15, 61-68

Once an experiment has been defined, it is possible to limit the peak list to a desired range using :meth:`sele_rt_range() <pyms.Experiment.Experiment.sele_rt_range>`.
The resulting experiment object can then be stored for later alignment.

Multiple Experiments
---------------------------

.. note:: This example is in :ref:`pyms-demo/61a <demo-61a>`

This example considers the preparation of three GC-MS experiments for peak alignment. The experiments are named ``a0806_077``, ``a0806_078``, ``a0806_079``, and represent separate GC-MS sample runs from the same biological sample.

The procedure is the same as above, and repeated for each experiment.
For example:

.. literalinclude:: ../pyms-demo/61a/proc.py
   :linenos:
   :lines: 19-24, 28-34, 77-95


.. note:: This example is in :ref:`pyms-demo/61b <demo-61b>`

The previous set of data all belong to the same experimental condition.
That is, they represent one group and any comparison between the data is a within group comparison.
For the original experiment, another set of GC-MS data was collected for a different experimental condition.
This group must also be stored as a set of experiments, and can be used for between group comparison.

The experiments are named ``a0806_140``, ``a0806_141``, ``a0806_142``, and are processed and stored as above.

Dynamic programming alignment of peak lists from multiple experiments
========================================================================

.. note:: This example is in :ref:`pyms-demo/62 <demo-62>`

In this example the experiments ``a0806_077``, ``a0806_078``, and ``a0806_079`` prepared in :ref:`pyms-demo/61a <demo-61a>` will be aligned, and therefore the script pyms-demo/61a/proc.py must be run first, to create the files ``a0806_077.expr``, ``a0806_078.expr``, ``a0806_079.expr`` in the directory pyms-demo/61a/output/.
These files contain the post-processed peak lists from the three experiments.

A script for running the dynamic programming alignment on these experiments is
given below.

.. literalinclude:: ../pyms-demo/62/proc.py
   :linenos:
   :lines: 4-30

The script reads the experiment files from the directory where they were stored (61a/output), and creates a list of the loaded :class:`~pyms.Experiment.Experiment` objects.
Each experiment object is converted into an :class:`~pyms.DPA.Class.Alignment` object with the function :meth:`exprl2alignment() <pyms.DPA.Function.exprl2alignment>`.
In this example, there is only one experimental condition so the alignment object is only for within group alignment (this special case is called 1-alignment).
The variable F1 is a Python list containing three alignment objects.

The pairwise alignment is then performed.
The parameters for the alignment by dynamic programming are: ``Dw``, the retention time modulation in seconds; and ``Gw``, the gap penalty.
These parameters are explained in detail in [1]_.
:class:`~pyms.DPA.Class.PairwiseAlignment`, defined in :mod:`pyms.DPA.Class`, is a class that calculates the similarity between all peaks in one sample with those of another sample.
This is done for all possible pairwise alignments (2-alignments).
The output of :class:`~pyms.DPA.Class.PairwiseAlignment` (``T1``) is an object which contains the dendrogram tree that maps the similarity relationship between the input 1-alignments, and also 1-alignments themselves.

The function :meth:`align_with_tree() <pyms.DPA.Class.align_with_tree>` takes the object T1 and aligns the individual alignment objects according to the guide tree.
In this example, the individual alignments are three 1-alignments, and the function :meth:`align_with_tree() <pyms.DPA.Class.align_with_tree>` first creates a 2-alignment from the two most similar 1-alignments and then adds the third 1-alignment to this to create a 3-alignment.
The parameter ``min_peaks=2`` specifies that any peak column of the data matrix that has fewer than two peaks in the final alignment will be dropped.
This is useful to clean up the data matrix of accidental peaks that are not truly observed over the set of replicates.

Finally, the resulting 3-alignment is saved by writing alignment tables containing peak retention times (``rt1.csv``) and the corresponding peak areas (`area1.csv').
These two files are plain ASCII files is CSV format, and are saved in the directory ``62/output/``.

The file ``area1.csv`` contains the data matrix where the corresponding peaks are aligned in the columns and each row corresponds to an experiment.
The file ``rt1.csv`` is useful for manually inspecting the alignment in some GUI driven program.

Between-state alignment of peak lists from multiple experiments
==================================================================

.. note:: This example is in :ref:`pyms-demo/63 <demo-63>`

In the previous example the list of peaks were aligned within a single experiment with multiple replicates ("within-state alignment").
In practice, it is of more interest to compare the two experimental states.
In a typical experimental setup there can be multiple replicate experiments on each experimental state or condition.
To analyze the results of such an experiment statistically, the list of peaks need to be aligned within each experimental state and also between the states.
The result of such an alignment would be the data matrix of integrated peak areas.
The data matrix contains a row for each sample and the number of columns is determined by the number of unique peaks (metabolites) detected in all the experiments.

In principle, all experiments could be aligned across conditions and replicates in the one process. However, a more robust approach is to first align experiments within each set of replicates (within-state alignment), and then to align the resulting alignments (between-state alignment) [1]_.

This example demonstrates how the peak lists from two cell states are aligned.
The cell state, A, consisting of three experiments aligned in :ref:`pyms-demo/61a <demo-61a>` (``a0806_077``, ``a0806_078``, ``a0806_079``) and cell state, B, consisting of three experiments aligned in :ref:`pyms-demo/61b <demo-61b>` (``a0806_140``, ``a0806_141``, and ``a0806_142``).

The between group alignment can be performed by the following alignment commands:

.. literalinclude:: ../pyms-demo/63/proc.py
   :linenos:
   :lines: 46-58

where ``A1`` and ``A2`` are the results of the within group alignments (as above) for group A and B, respectively.

In this example the retention time tolerance for between-state alignment is greater compared to the retention time tolerance for the within-state alignment as we expect less fidelity in retention times between them.
The same functions are used for the within-state and between-state alignment.
The result of the alignment is saved to a file as the area and retention time matrices (described above).

Common Ion Area Quantitation
================================
.. note:: This example is in :ref:`pyms-demo/64 <demo-64>`

The ``area.csv`` file produced in the preceding section lists the total area of each peak in the alignment.
The total area is the sum of the areas of each of the individual ions in the peak.
While this approach produces broadly accurate results, it can result in errors where neighbouring peaks or unfiltered noise add to the peak in some way.

One alternative to this approach is to pick a single ion which is common to a particular peak (compound), and to report only the area of this ion for each occurrence of that peak in the alignment.
Using the method :meth:`common_ion() <pyms.DPA.Class.Alignment.common_ion>` of the class :class:`~pyms.DPA.Class.Alignment`, |pkgname| can select an ion for each aligned peak which is both abundant and occurs most often for that peak.
We call this the 'Common Ion Algorithm' (CIA).

To use this method it is essential that the individual ion
areas have been set (see section :ref:`individual_ion_areas`).


Using the Common Ion Algorithm
----------------------------------

When using the CIA for area quantitation, a different method of the class :class:`~pyms.DPA.Class.Alignment` is used to write the area matrix; :meth:`write_common_ion_csv() <pyms.DPA.Class.Alignment.write_common_ion_csv>`.
This requires a list of the common ions for each peak in the alignment.
This list is generated using the Alignment class method :meth:`common_ion() <pyms.DPA.Class.Alignment.common_ion>`.

Continuing from the previous example, the following invokes common ion filtering on previously created alignment object 'A9':

    >>> common_ion_list = A9.common_ion()

The variable 'common_ion_list' is a list of the common ion for each
peak in the alignment. This list is the same length as the
alignment. To write peak areas using common ion quantitation:

   >>> A9.write_common_ion_csv('output/area_common_ion.csv',common_ion_list)

References
============

.. [1] Robinson MD, De Souza DP, Keen WW, Saunders EC, McConville MJ, Speed TP, and Likic VA. A dynamic programming approach for the alignment of signal peaks in multiple gas chromatography-mass spectrometry experiments. `BMC Bioinformatics`, 8:419, 2007

