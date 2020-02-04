"""
Convert ipynb notebook to python script

Requires nbconvert (pip install nbcovnert) and pandoc (apt-get install pandoc)
"""

# Import the Python exporter and instantiate iy
from nbconvert import PythonExporter
py_exporter = PythonExporter()


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

for notebook in notebooks:
	# Convert the notebook to a python file
	script, *_ = py_exporter.from_file(f"./jupyter/{notebook}.ipynb")
	
	with open(f"./scripts/{notebook}.py", "w") as fp:
		fp.write(script)
