"""
Class to Display Ion Chromatograms and TIC
"""

################################################################################
#                                                                              #
#    PyMassSpec software for processing of mass-spectrometry data              #
#    Copyright (C) 2005-2012 Vladimir Likic                                    #
#    Copyright (C) 2019 Dominic Davis-Foster                                   #
#                                                                              #
#    This program is free software; you can redistribute it and/or modify      #
#    it under the terms of the GNU General Public License version 2 as         #
#    published by the Free Software Foundation.                                #
#                                                                              #
#    This program is distributed in the hope that it will be useful,           #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#    GNU General Public License for more details.                              #
#                                                                              #
#    You should have received a copy of the GNU General Public License         #
#    along with this program; if not, write to the Free Software               #
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 #
#                                                                              #
################################################################################

import matplotlib.pyplot as plt

from pyms.IonChromatogram import IonChromatogram
from pyms.Peak.List.Function import is_peak_list


class Display(object):
	"""
	Class to display Ion Chromatograms and Total Ion Chromatograms from class:`IonChromatogram.IonChromatogram`
	Uses matplotlib module pyplot to do plotting

	:author: Sean O'Callaghan
	:author: Vladimir Likic
	:author: Dominic Davis-Foster (type assertions)
	"""
	
	def __init__(self):
		"""
		Initialises an instance of Display class
		"""
		
		# Container to store plots
		self.__tic_ic_plots = []
		
		# color dictionary for plotting of ics; blue reserved
		# for TIC
		self.__col_ic = {0:'r', 1:'g', 2:'k', 3:'y', 4:'m', 5:'c'}
		self.__col_count = 0  # counter to keep track of colors
		
		# Peak list container
		self.__peak_list = []
		
		#Plotting Variables
		self.__fig = plt.figure()
		self.__ax = self.__fig.add_subplot(111)
	
	def do_plotting(self, plot_label=None):
		"""
		Plots TIC and IC(s) if they have been created by plot_tic() or plot_ics().
		Adds detected peaks if they have been added by plot_peaks()

		:param plot_label: Optional to supply a label or other
				definition of data origin
		:type plot_label: str
		"""
		
		# if no plots have been created advise user
		if len(self.__tic_ic_plots) == 0:
			print('No plots have been created')
			print('Please call a plotting function before')
			print('calling do_plotting()')
		
		if plot_label is not None:
			t = self.__ax.set_title(plot_label)
		
		l = self.__ax.legend()
		
		self.__fig.canvas.draw
		
		# If no peak list plot, no mouse click event
		if len(self.__peak_list) != 0:
			cid = self.__fig.canvas.mpl_connect('button_press_event', self.onclick)
		plt.show()
	
	@staticmethod
	def get_5_largest(intensity_list):
		"""
		Computes the indices of the largest 5 ion intensities for writing to console

		:param intensity_list: List of Ion intensities
		:type intensity_list: list

		:return: Indices of largest 5 ion intensities
		:rtype: list
		"""
		
		largest = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		
		# Find out largest value
		for i in range(len(intensity_list)):
			if intensity_list[i] > intensity_list[largest[0]]:
				largest[0] = i
		
		# Now find next four largest values
		for j in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
			for i in range(len(intensity_list)):
				# if intensity_list[i] > intensity_list[largest[j]] and intensity_list[i] < intensity_list[largest[j-1]]:
				if intensity_list[largest[j]] < intensity_list[i] < intensity_list[largest[j - 1]]:
					largest[j] = i
		
		return largest
	
	def onclick(self, event):
		"""
		Finds the 5 highest intensity m/z channels for the selected peak.
		The peak is selected by clicking on it.
		If a button other than the left one is clicked, a new plot of the mass spectrum is displayed

		:param event: a mouse click by the user
		"""
		
		intensity_list = []
		mass_list = []
		
		for peak in self.__peak_list:
			# if event.xdata > 0.9999*peak.rt and event.xdata < 1.0001*peak.rt:
			if 0.9999 * peak.rt < event.xdata < 1.0001 * peak.rt:
				intensity_list = peak.get_mass_spectrum().mass_spec
				mass_list = peak.get_mass_spectrum().mass_list
		
		largest = self.get_5_largest(intensity_list)
		
		if len(intensity_list) != 0:
			print("mass\t intensity")
			for i in range(10):
				print(mass_list[largest[i]], "\t", intensity_list[largest[i]])
		else:  # if the selected point is not close enough to peak
			print("No Peak at this point")
		
		# Check if a button other than left was pressed, if so plot mass spectrum
		# Also check that a peak was selected, not just whitespace
		if event.button != 1 and len(intensity_list) != 0:
			self.plot_mass_spec(event.xdata, mass_list, intensity_list)

	def plot_ics(self, ics, labels=None):
		"""
		Adds an Ion Chromatogram or a list of Ion Chromatograms to plot list

		:param ics: List of Ion Chromatograms m/z channels for plotting
		:type ics: list of class:`pyms.IonChromatogram.IonChromatogram` objects

		:param labels: Labels for plot legend
		:type labels: list of str, optional
		"""
		
		if not isinstance(ics, list):
			if isinstance(ics, IonChromatogram):
				ics = [ics]
			else:
				raise TypeError("'ics' must be an IonChromatogram or a list of Ion Chromatograms")
		
		if not isinstance(labels, list) and labels != None:
			labels = [labels]
		# TODO: take care of case where one element of ics is not an IonChromatogram
		
		intensity_list = []
		time_list = ics[0].get_time_list()
		
		for i in range(len(ics)):
			intensity_list.append(ics[i].get_intensity_array())
		
		# Case for labels not present
		if labels is None:
			for i in range(len(ics)):
				self.__tic_ic_plots.append(plt.plot(time_list,
													intensity_list[i],
													self.__col_ic[self.__col_count]))
				if self.__col_count == 5:
					self.__col_count = 0
				else:
					self.__col_count += 1
		
		# Case for labels present
		else:
			for i in range(len(ics)):
				
				self.__tic_ic_plots.append(plt.plot(time_list,
													intensity_list[i],
													self.__col_ic[self.__col_count],
													label=labels[i]))
				if self.__col_count == 5:
					self.__col_count = 0
				else:
					self.__col_count += 1
	
	@staticmethod
	def plot_mass_spec(rt, mass_list, intensity_list):
		"""
		Plots the mass spec given a list of masses and intensities

		:param rt: The retention time for labelling of the plot
		:type rt: float
		:param mass_list: list of masses of the MassSpectrum object
		:type mass_list: list
		:param intensity_list: List of intensities of the MassSpectrum object
		:type intensity_list: list
		"""
		
		new_fig = plt.figure()
		new_ax = new_fig.add_subplot(111)
		
		# to set x axis range find minimum and maximum m/z channels
		max_mz = mass_list[0]
		min_mz = mass_list[0]
		
		for i in range(len(mass_list)):
			if mass_list[i] > max_mz:
				max_mz = mass_list[i]
		
		for i in range(len(mass_list)):
			if mass_list[i] < min_mz:
				min_mz = mass_list[i]
		
		label = "Mass spec for peak at time " + "%5.2f" % rt
		
		mass_spec_plot = plt.bar(mass_list, intensity_list, label=label, width=0.01)
		
		x_axis_range = plt.xlim(min_mz, max_mz)
		t = new_ax.set_title(label)
		plt.show()
	
	def plot_peaks(self, peak_list, label="Peaks"):
		"""
		Plots the locations of peaks as found by PyMassSpec.

		:param peak_list: List of peaks
		:type peak_list: list of class:`pyms.Peak.Class.Peak` objects

		:param label: label for plot legend, default "Peaks"
		:type label: str, optional
		"""
		
		if not is_peak_list(peak_list):
			raise TypeError("'peak_list' must be a list of Peak objects")
		
		time_list = []
		height_list = []
		
		# Copy to self.__peak_list for onclick event handling
		self.__peak_list = peak_list
		
		for peak in peak_list:
			time_list.append(peak.rt)
			height_list.append(sum(peak.get_mass_spectrum().mass_spec))
		
		self.__tic_ic_plots.append(plt.plot(time_list, height_list, 'o', label=label))
	
	def plot_tic(self, tic, label=None):
		"""
		Adds Total Ion Chromatogram to plot list

		:param tic: Total Ion Chromatogram
		:type tic: class:`pyms.IonChromatogram.IonChromatogram`

		:param label: label for plot legend
		:type label: str, optional
		"""
		
		if not isinstance(tic, IonChromatogram):
			raise TypeError("'tic' must be an Ion Chromatogram object")
		
		intensity_list = tic.get_intensity_array()
		time_list = tic.get_time_list()
		
		self.__tic_ic_plots.append(plt.plot(time_list, intensity_list, label=label))

