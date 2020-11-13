import pathlib
import sys

from notebook2script.__main__ import process_multiple_notebooks

notebooks_dir = pathlib.Path(__file__).parent / "jupyter"
scripts_dir = pathlib.Path(__file__).parent / "scripts"

notebooks = [
		"reading_jcamp",
		"reading_andi",
		"comparing_datasets",
		"IntensityMatrix",
		"IntensityMatrix_Resizing",
		"IntensityMatrix_Preprocessing",
		"MassSpectrum",
		"IonChromatogram",
		"NoiseSmoothing",
		"BaselineCorrection",
		"Peak",
		"Peak_Detection",
		"Peak_Filtering_Noise_Analysis",
		"Peak_Area_Estimation",
		"Experiment",
		"Multiple_Experiments",
		"DPA",

		]

# Never include these notebooks!
"""
Displaying_TIC
Displaying_Multiple_IC
Displaying_Mass_Spec
Displaying_Detected_Peaks
Display_User_Interaction
"""

sys.exit(process_multiple_notebooks([notebooks_dir / f"{nb}.ipynb" for nb in notebooks], outdir=scripts_dir, overwrite=True))
