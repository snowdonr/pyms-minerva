"""
Provides a class for handling Missing Peaks in an output file (i.e. area.csv)
"""

#############################################################################
#                                                                           #
#    PyMS software for processing of metabolomic mass-spectrometry data     #
#    Copyright (C) 2005-2012 Vladimir Likic                                 #
#                                                                           #
#    This program is free software; you can redistribute it and/or modify   #
#    it under the terms of the GNU General Public License version 2 as      #
#    published by the Free Software Foundation.                             #
#                                                                           #
#    This program is distributed in the hope that it will be useful,        #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
#    GNU General Public License for more details.                           #
#                                                                           #
#    You should have received a copy of the GNU General Public License      #
#    along with this program; if not, write to the Free Software            #
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              #
#                                                                           #
#############################################################################

import deprecation
from pyms import __version__


class MissingPeak(object):
    """
    :summary: Class to encapsulate a peak object identified as missing in
        the output area matrix fom PyMS.
    
    :author: Jairus Bowne
    :author: Sean O'Callaghan
    :author: Dominic Davis-Foster (properties)
    """
    
    def __init__(self, common_ion, qual_ion_1, qual_ion_2, rt=0.0):
        """
        :summary: Initialise MissingPeak Class
        
        :param common_ion: Common ion for the peak across samples in an experiment
        :type common_ion: IntType
        :param qual_ion_1:
        :type qual_ion_1:
        :param qual_ion_2:
        :type qual_ion_2:
        :param rt: Retention time of the peak. May or may not be set
        :type rt: FloatType
        """
        
        self.__common_ion = common_ion
        self.__qual_1 = qual_ion_1
        self.__qual_2 = qual_ion_2
        self.__rt = rt
        self.__exact_rt = 'na'
        self.__common_ion_area = 'na'

    @deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
                            current_version=__version__,
                            details="Use 'MissingPeak.common_ion' instead")
    def get_common_ion(self):
        """
        :summary: Returns the common ion for the peak object across an
            experiment
        
        :return: Common ion for the peak
        :rtype: IntType
        
        :author: Jairus Bowne
        """
        
        return self.common_ion
    
    @property
    def common_ion(self):
        """
        :summary: Returns the common ion for the peak object across an
            experiment

        :return: Common ion for the peak
        :rtype: IntType

        :author: Jairus Bowne
        """
        
        return self.__common_ion
    
    @deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
                            current_version=__version__,
                            details="Use 'MissingPeak.qual_ion1' instead")
    def get_qual_ion1(self):
        """
        :summary: Returns the top (most abundant) ion for the peak object
        
        :return: Most abundant ion
        :rtype: IntType
        
        :author: Jairus Bowne
        """
        
        return self.qual_ion1
    
    @property
    def qual_ion1(self):
        """
        :summary: Returns the top (most abundant) ion for the peak object
        
        :return: Most abundant ion
        :rtype: IntType
        
        :author: Jairus Bowne
        """
        
        """
        # TODO: Consider the abundance of ions when some (i.e. 73, 147) have
            been im.null_mass()'d. Is there a way to determine whether that
            has been done to generate the original peak list?
        """
        
        return self.__qual_1
        #return int(string.split(self.__UID, '-')[0])
    
    @deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
                            current_version=__version__,
                            details="Use 'MissingPeak.qual_ion2' instead")
    def get_qual_ion2(self):
        """
        :summary: Returns the second most abundant ion for the peak object
        
        :return: Second most abundant ion
        :rtype: IntType
        
        :author: Jairus Bowne
        """
        
        return self.qual_ion2
    
    @property
    def qual_ion2(self):
        """
        :summary: Returns the second most abundant ion for the peak object
        
        :return: Second most abundant ion
        :rtype: IntType
        
        :author: Jairus Bowne
        """
        
        """
        # TODO: Consider the abundance of ions when some (i.e. 73, 147) have
            been im.null_mass()'d. Is there a way to determine whether that
            has been done to generate the original peak list?
        """
        
        return self.__qual_1
        #return int(string.split(self.__UID, '-')[0])
    
    @property
    def common_ion_area(self):
        """
        :summary: returns the common ion area

        :return common_ion_area: The area of the common ion
        :rtype: intType
        """
        
        return self.__common_ion_area
    
    @common_ion_area.setter
    def common_ion_area(self, common_ion_area):
        """
        :summary: sets the common ion area calculated by the gap fill
                  algorithm
        :param common_ion_area: The area of the common ion
        :type common_ion_area: intType

        """
        self.__common_ion_area = common_ion_area

    @deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
                            current_version=__version__,
                            details="Use 'MissingPeak.common_ion_area' instead")
    def set_common_ion_area(self, common_ion_area):
        """
        :summary: sets the common ion area calculated by the gap fill
                  algorithm
        :param common_ion_area: The area of the common ion
        :type common_ion_area: intType

        """
        self.common_ion_area = common_ion_area

    @deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
                            current_version=__version__,
                            details="Use 'MissingPeak.common_ion_area' instead")
    def get_common_ion_area(self):
        """
        :summary: returns the common ion area

        :return common_ion_area: The area of the common ion
        :rtype: intType
        """
        return self.common_ion_area

    @deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
                            current_version=__version__,
                            details="Use 'MissingPeak.rt' instead")
    def get_rt(self):
        """
        :summary: returns the retention time of the peak

        :return: the retention time of the peak
        :rtype: floatType
        """
        return self.rt
    
    @property
    def rt(self):
        """
        :summary: returns the retention time of the peak

        :return: the retention time of the peak
        :rtype: floatType
        """
        
        return self.__rt

    @deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
                            current_version=__version__,
                            details="Use 'MissingPeak.exact_rt' instead")
    def set_exact_rt(self, rt):
        """
        :summary: sets the retention time of a peak

        :param rt: The retention time of the apex of the peak
        :type rt: floatType
        """
        self.__exact_rt = rt

    @deprecation.deprecated(deprecated_in="2.1.2", removed_in="2.2.0",
                            current_version=__version__,
                            details="Use 'MissingPeak.exact_rt' instead")
    def get_exact_rt(self):
         """
        :summary: returns the retention time of the peak

        :return: the retention time of the peak
        :rtype: floatType
        """
         return self.exact_rt
    
    @property
    def exact_rt(self):
         """
        :summary: returns the retention time of the peak

        :return: the retention time of the peak
        :rtype: floatType
        """
         
         return self.__exact_rt
    
    @exact_rt.setter
    def exact_rt(self, rt):
        """
        :summary: sets the retention time of a peak

        :param rt: The retention time of the apex of the peak
        :type rt: floatType
        """
        
        self.__exact_rt = rt
    

class Sample(object):
    """
    :summary: A collection of MissingPeak objects

    :author: Sean O'Callaghan
    """

    def __init__(self, sample_name, matrix_position):
        """
        :summary: A collection of MissingPeak objects

        :param sample_name: the experiment code/name
        :type sample_name: stringType

        :param matrix_position: position along x-axis
                                where sample is located
        :type matrix_position: intType

        """
        self.__sample_name = sample_name
        self.__matrix_position = matrix_position

        self.__missing_peak_list = []

    def get_name(self):
        """
        :summary: Returns the sample name

        :return: The name of the sample
        :rtype: stringType
        """
        return self.__sample_name
    
    def add_missing_peak(self, missing_peak):
        """
        :summary: Add a new MissingPeak object to the Sample

        :param missing_peak: The missing peak object to be added
        :type missing_peak: pyms.GapFilling.Class.MissingPeak
        """
        ###
        # Do some checking here!!!
        ###
        self.__missing_peak_list.append(missing_peak)

    def get_missing_peaks(self):
        """
        :summary: Returns a list of the MissingPeak objects
                  in the Sample object
        :return: list of pyms.GapFilling.Class.MissingPeak
        :rtype: listType
        """
        return self.__missing_peak_list

    def get_mp_rt_area_dict(self):
        """
        :summary: returns a dictionary containing rt:area pairs

        :return: a dict containing rt:area pairs
        :rtype: dictType
        """
        rt_area_dict = {}        
        for peak in self.__missing_peak_list:
            rt = peak.get_rt()
            area = peak.get_common_ion_area()

            rt_area_dict[rt] = area
        
        return rt_area_dict

    def get_mp_rt_exact_rt_dict(self):
        """
        :summary:returns a dictionary containing average_rt:exact_rt pairs

        :return: a dict of average_rt:exact_rt pairs
        :rtype: dictType
        """

        rt_exact_rt_dict = {}
        for peak in self.__missing_peak_list:
            rt = peak.get_rt()
            exact_rt = peak.get_exact_rt()

            rt_exact_rt_dict[rt] = exact_rt

        return rt_exact_rt_dict

