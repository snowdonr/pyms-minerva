"""
Models a GC-MS experiment, represented by a list of signal peaks
"""

################################################################################
#                                                                              #
#    PyMassSpec software for processing of mass-spectrometry data              #
#    Copyright (C) 2005-2012 Vladimir Likic                                    #
#    Copyright (C) 2019-2020 Dominic Davis-Foster                              #
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

# stdlib
import copy
import os
import pathlib
import pickle
from typing import Any, List, Sequence, Union

# 3rd party
import deprecation  # type: ignore

# this package
from pyms import __version__
from pyms.Base import pymsBaseClass
from pyms.Peak.Class import Peak
from pyms.Peak.List.Function import is_peak_list, sele_peaks_by_rt
from pyms.Utils.IO import prepare_filepath
from pyms.Utils.Utils import is_path, is_sequence

__all__ = ["Experiment", "read_expr_list", "load_expr", "store_expr"]


class Experiment(pymsBaseClass):
	"""
	Models an experiment object

	:param expr_code: A unique identifier for the experiment.
	:param peak_list: A list of peak objects

	:author: Vladimir Likic, Andrew Isaac,  Dominic Davis-Foster (type assertions, properties and pathlib support)
	"""

	def __init__(self, expr_code: str, peak_list: Sequence[Peak]):
		"""
		Models an experiment
		"""

		if not isinstance(expr_code, str):
			raise TypeError("'expr_code' must be a string")

		if not is_peak_list(peak_list):
			raise TypeError("'peak_list' must be a list of Peak objects")

		self._expr_code = expr_code
		self._peak_list = list(peak_list)

	def __eq__(self, other: Any) -> bool:
		"""
		Return whether this Experiment object is equal to another object.

		:param other: The other object to test equality with.
		"""

		if isinstance(other, self.__class__):
			return self.peak_list == other.peak_list and self.expr_code == other.expr_code

		return NotImplemented

	def __len__(self) -> int:
		"""
		Returns the number of peaks in the Experiment
		"""

		return len(self.peak_list)

	def __copy__(self) -> "Experiment":
		"""
		Returns a new Experiment object containing a copy of the data in this object
		"""

		return Experiment(copy.copy(self._expr_code), copy.copy(self.peak_list))

	def __deepcopy__(self, memodict={}) -> "Experiment":
		"""
		Returns a new Experiment object containing a copy of the data in this object

		:rtype: pyms.Experiment.Experiment
		"""

		return self.__copy__()

	@property
	def expr_code(self) -> str:
		"""
		Returns the expr_code of the experiment
		"""

		return self._expr_code

	@property
	def peak_list(self) -> List[Peak]:
		"""
		Returns the peak list
		"""

		return self._peak_list

	def sele_rt_range(self, rt_range: Sequence[str]):
		"""
		Discards all peaks which have the retention time outside the specified range

		:param rt_range: Min, max retention time given as a list [rt_min, rt_max]
		"""

		if not is_sequence(rt_range):
			raise TypeError("'rt_range' must be a Sequence")

		peaks_sele = sele_peaks_by_rt(self._peak_list, rt_range)
		self._peak_list = peaks_sele

	@deprecation.deprecated(
			deprecated_in="2.1.2",
			removed_in="2.2.0",
			current_version=__version__,
			details="Use :meth:`pyms.Experiment.Experiment.dump` instead",
			)
	def store(self, file_name: Union[str, pathlib.Path]):
		"""
		stores an experiment to a file

		:param file_name: The name of the file

		:author: Vladimir Likic, Andrew Isaac, Dominic Davis-Foster (pathlib support)
		"""

		if not is_path(file_name):
			raise TypeError("'file_name' must be a string or a PathLike object")

		file_name = prepare_filepath(file_name)

		fp = file_name.open('wb')
		pickle.dump(self, fp, 1)
		fp.close()


def read_expr_list(file_name: Union[str, pathlib.Path]) -> List[Experiment]:
	"""
	Reads the set of experiment files and returns a list of :class:`pyms.Experiment.Experiment` objects

	:param file_name: The name of the file which lists experiment dump file names, one file per line

	:return: A list of Experiment instances

	:author: Vladimir Likic
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name, mkdirs=False)

	fp = file_name.open()

	exprfiles = fp.readlines()
	fp.close()

	expr_list = []

	for exprfile in exprfiles:

		exprfile = exprfile.strip()
		expr = load_expr(exprfile)

		expr_list.append(expr)

	return expr_list


def load_expr(file_name: Union[str, pathlib.Path]) -> Experiment:
	"""
	Loads an experiment saved with :meth:`pyms.Experiment.store_expr`

	:param file_name: Experiment file name

	:return: The loaded experiment

	:author: Vladimir Likic, Andrew Isaac, Dominic Davis-Foster (type assertions and pathlib support)
	"""

	if not is_path(file_name):
		raise TypeError("'file_name' must be a string or a PathLike object")

	file_name = prepare_filepath(file_name, mkdirs=False)

	fp = file_name.open('rb')
	expr = pickle.load(fp)
	fp.close()

	if not isinstance(expr, Experiment):
		raise IOError("The loaded file is not an experiment file")

	return expr


@deprecation.deprecated(
		deprecated_in="2.1.2",
		removed_in="2.2.0",
		current_version=__version__,
		details="Use :meth:`pyms.Experiment.Experiment.store` instead",
		)
def store_expr(file_name: str, expr: Experiment):
	"""
	Stores an experiment to a file

	:param file_name: The name of the file
	:param expr: An experiment object

	:author: Vladimir Likic, Andrew Isaac, Dominic Davis-Foster (type assertions)
	"""

	if not isinstance(expr, Experiment):
		raise TypeError("'expr' must be an Experiment object")

	if not isinstance(file_name, str):
		raise TypeError("'file_name' must be a string")

	if not os.path.exists(os.path.dirname(file_name)):
		os.makedirs(os.path.dirname(file_name))

	fp = open(file_name, 'wb')
	pickle.dump(expr, fp, 1)
	fp.close()
