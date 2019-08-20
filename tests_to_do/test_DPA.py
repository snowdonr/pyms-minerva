import pytest

test_string = "abcdefg"
test_int = 1234
test_float = 12.34
test_list_ints = [1,2,3,4]
test_list_strs = ["a", "b", "c", "d"]
test_dict = {"a":1, "b":2, "c":3, "d":4}

import os

from pyms.Peak.List.DPA.Class import PairwiseAlignment
from pyms.Peak.List.DPA.Function import align_with_tree, exprl2alignment
from pyms.Baseline.TopHat import tophat
from pyms.Deconvolution.BillerBiemann.Function import BillerBiemann
from pyms.GCMS.Function import build_intensity_matrix_i
from pyms.GCMS.IO.JCAMP import JCAMP_reader
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Peak.Function import peak_sum_area, peak_top_ion_areas
from pyms.Deconvolution.BillerBiemann.Function import rel_threshold, num_ions_threshold
from pyms.Experiment.Class import Experiment
from pyms.Experiment import store_expr, load_expr
from pyms.Peak.IO import store_peaks
from pyms.Peak.List.DPA.Class import PairwiseAlignment, Alignment
from pyms.Peak.List.DPA.Function import align_with_tree, exprl2alignment
#from pyms.Peak.IO import store_peaks


def test_alignment():
    # Create experiment files
    file_list = ["ELEY_1_SUBTRACT", "ELEY_2_SUBTRACT", "ELEY_3_SUBTRACT", "ELEY_4_SUBTRACT", "ELEY_5_SUBTRACT"]
    for jcamp_file in file_list:
        
        im = build_intensity_matrix_i(JCAMP_reader(os.path.join("data", f"{jcamp_file}.JDX")))
        
        # Intensity matrix size (scans, masses)
        n_scan, n_mz = im.get_size()
    
        # noise filter and baseline correct
        for ii in range(n_mz):
            ic = im.get_ic_at_index(ii)
            ic_smooth = savitzky_golay(ic)
            ic_bc = tophat(ic_smooth, struct="1.5m")
            im.set_ic_at_index(ii, ic_bc)
    
        peak_list = BillerBiemann(im, points=9, scans=2)
        apl = rel_threshold(peak_list, 2)
        new_peak_list = num_ions_threshold(apl, 3, 3000)
    
        # ignore TMS ions and set mass range
        for peak in new_peak_list:
            peak.crop_mass(50, 400)
            peak.null_mass(73)
            peak.null_mass(147)
        
            # find area
            area = peak_sum_area(im, peak)
            peak.set_area(area)
            area_dict = peak_top_ion_areas(im, peak)
            peak.set_ion_areas(area_dict)
    
        expr = Experiment(jcamp_file, new_peak_list)
    
        # set time range for all experiments
        expr.sele_rt_range(["6.5m", "21m"])
        
        store_expr(os.path.join("output",f"{jcamp_file}.expr"), expr)

    
    # within replicates alignment parameters
    Dw = 2.5  # rt modulation [s]
    Gw = 0.30 # gap penalty
    
    # Load experiments
    expr_list = []
    for expr_code in file_list:
        expr = load_expr(os.path.join("output",f"{expr_code}.expr"))
        assert isinstance(expr, Experiment)
        expr_list.append(expr)
    
    # Test inequality
    assert expr_list[0] != expr_list[1]

    # do the alignment
    print('Aligning ELEY SUBTRACT')
    F1 = exprl2alignment(expr_list)
    assert isinstance(F1, list)
    
    for type in [test_string, test_int, test_float, test_list_ints, test_list_strs, test_dict]:
        with pytest.raises(TypeError):
            exprl2alignment(type)
    
    T1 = PairwiseAlignment(F1, Dw, Gw)
    assert isinstance(T1, PairwiseAlignment)
    
    for type in [test_string, test_int, test_float, test_list_ints, test_list_strs, test_dict]:
        with pytest.raises(TypeError):
            PairwiseAlignment(type, Dw, Gw)
    for type in [test_string, test_int, test_list_ints, test_list_strs, test_dict]:
        with pytest.raises(TypeError):
            PairwiseAlignment(F1, type, Gw)
    for type in [test_string, test_int, test_list_ints, test_list_strs, test_dict]:
        with pytest.raises(TypeError):
            PairwiseAlignment(F1, Dw, type)

    A1 = align_with_tree(T1, min_peaks=2)
    assert isinstance(A1, Alignment)
    
    assert isinstance(len(A1), int)
    assert len(A1) == 232
    A1.filter_min_peaks(5)
    assert len(A1) == 50
    
    
    A1.write_csv('output/alignment_rt.csv', 'output/alignment_area.csv')
    # TODO: read the csv and check values
    A1.write_csv('output/alignment_rt.csv', 'output/alignment_area.csv', minutes=False)
    # TODO: read the csv and check values
    A1.write_csv('output/alignment_rt.csv', 'output/alignment_area.csv', dk=True)
    # TODO: read the csv and check values
    
    
    common_ion = A1.common_ion()
    assert isinstance(common_ion, list)
    assert isinstance(common_ion[0], int)
    assert common_ion[0] == 77
    
    A1.write_common_ion_csv('output/alignent_ion_area.csv', A1.common_ion())
    # TODO: read the csv and check values
    A1.write_common_ion_csv('output/alignent_ion_area.csv', A1.common_ion(), minutes=False)
    # TODO: read the csv and check values
    
    # Errors
    for type in [test_int, test_float, test_string, test_dict, test_list_strs, test_list_ints]:
        with pytest.raises(TypeError):
            Alignment(type)

    for type in [test_float, test_string, test_dict, test_list_strs, test_list_ints]:
        with pytest.raises(TypeError):
            A1.filter_min_peaks(type)

    for type in [test_float, test_int, test_dict, test_list_strs, test_list_ints]:
        with pytest.raises(TypeError):
            A1.write_csv(type, 'output/alignment_area.csv')
    
    for type in [test_float, test_int, test_dict, test_list_strs, test_list_ints]:
        with pytest.raises(TypeError):
            A1.write_csv('output/alignment_rt.csv', type)
    
    for type in [test_float, test_int, test_dict, test_list_strs, test_list_ints]:
        with pytest.raises(TypeError):
            A1.write_common_ion_csv(type, 'output/alignent_ion_list.csv')
    
    for type in [test_float, test_int, test_dict, test_list_strs, test_string]:
        with pytest.raises(TypeError):
            A1.write_common_ion_csv('output/alignent_ion_area.csv', type)


def test_align_2_alignments():
    # define the input experiments list
    eley_codes = ["ELEY_1_SUBTRACT", "ELEY_2_SUBTRACT", "ELEY_3_SUBTRACT", "ELEY_4_SUBTRACT", "ELEY_5_SUBTRACT"]
    geco_codes = ["GECO_1_SUBTRACT", "GECO_2_SUBTRACT", "GECO_3_SUBTRACT", "GECO_4_SUBTRACT", "GECO_5_SUBTRACT"]
    # ΓΕΨΟ
    # within replicates alignment parameters
    Dw = 2.5  # rt modulation [s]
    Gw = 0.30 # gap penalty
    
    expr_list = []
    
    for expr_file in eley_codes:
        expr = load_expr(os.path.join("output", f"{expr_file}.expr"))
        expr_list.append(expr)

    F1 = exprl2alignment(expr_list)
    T1 = PairwiseAlignment(F1, Dw, Gw)
    A1 = align_with_tree(T1, min_peaks=2)

    expr_list = []
    
    for jcamp_file in geco_codes:
        im = build_intensity_matrix_i(JCAMP_reader(os.path.join("data", f"{jcamp_file}.JDX")))
    
        # Intensity matrix size (scans, masses)
        n_scan, n_mz = im.get_size()
    
        # noise filter and baseline correct
        for ii in range(n_mz):
            ic = im.get_ic_at_index(ii)
            ic_smooth = savitzky_golay(ic)
            ic_bc = tophat(ic_smooth, struct="1.5m")
            im.set_ic_at_index(ii, ic_bc)
    
        peak_list = BillerBiemann(im, points=9, scans=2)
        apl = rel_threshold(peak_list, 2)
        new_peak_list = num_ions_threshold(apl, 3, 3000)
    
        # ignore TMS ions and set mass range
        for peak in new_peak_list:
            peak.crop_mass(50, 400)
            peak.null_mass(73)
            peak.null_mass(147)
        
            # find area
            area = peak_sum_area(im, peak)
            peak.set_area(area)
            area_dict = peak_top_ion_areas(im, peak)
            peak.set_ion_areas(area_dict)
    
        expr = Experiment(jcamp_file, new_peak_list)
    
        # set time range for all experiments
        expr.sele_rt_range(["6.5m", "21m"])
    
        expr_list.append(expr)

    F2 = exprl2alignment(expr_list)
    T2 = PairwiseAlignment(F2, Dw, Gw)
    A2 = align_with_tree(T2, min_peaks=2)

    #top_ion_list = A2.common_ion()
    #A2.write_common_ion_csv('output/area1.csv', top_ion_list)
    
    # between replicates alignment parameters
    Db = 10.0 # rt modulation
    Gb = 0.30 # gap penalty
    
    print('Aligning input {1,2}')
    T9 = PairwiseAlignment([A1,A2], Db, Gb)
    A9 = align_with_tree(T9)
    
    A9.write_csv('output/rt.csv', 'output/area.csv')
    
    aligned_peaks = A9.aligned_peaks()
    store_peaks(aligned_peaks, 'output/peaks.bin')
    
    top_ion_list = A9.common_ion()
    A9.write_common_ion_csv('output/area.csv', top_ion_list)

