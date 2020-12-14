#!/bin/bash

git config --global user.email "dominic@davis-foster.co.uk"
git config --global user.name "Dominic Davis-Foster"

# Install requirements
sudo apt update
sudo apt install -y python3-venv python3-pip || exit 1

# Create and activate venv
python3 -m venv /home/vagrant/venv
source /home/vagrant/venv/bin/activate || exit 1

# Install remaining requirements
python3 -m pip install --upgrade pip setuptools wheel || exit 1
python3 -m pip install --upgrade nbconvert jupyter-client ipykernel domdf_python_tools  || exit 1
python3 -m pip install --upgrade -r requirements.txt || exit 1
python3 -m pip install --upgrade . || exit 1

# Change to Notebooks branch
git checkout -b Notebooks || exit 1

cd pyms-demo/jupyter || exit 1

cat > switch_experiments.py <<EOF
#!/usr/bin/env/python3
import pathlib
import json

notebook = pathlib.Path("Multiple_Experiments.ipynb")
data = json.loads(notebook.read_text())
lines = data["cells"][3]["source"]
lines = [line[1:].lstrip() if line.startswith("#") else f"# {line}" for line in lines]
data["cells"][3]["source"] = lines
notebook.write_text(json.dumps(data, indent=1))
EOF

# Run Multiple_Experiments to ensure output files exist
jupyter nbconvert --to notebook --inplace --execute Multiple_Experiments.ipynb
python3 switch_experiments.py
jupyter nbconvert --to notebook --inplace --execute Multiple_Experiments.ipynb
python3 switch_experiments.py

# Render notebooks and stage
for file in *.ipynb; do
  jupyter nbconvert --clear-output --inplace "$file"
  jupyter nbconvert --to notebook --inplace --execute "$file"
done

# remove execution times
cat > remove_execution_times.py <<EOF
#!/usr/bin/env/python3
from domdf_python_tools.paths import PathPlus

for filename in PathPlus.cwd().glob("*.ipynb"):
	notebook = filename.load_json()
	for cell in notebook["cells"]:
		if "execution" in cell["metadata"]:
			cell["metadata"]["execution"] = {}
		if "pycharm" in cell["metadata"]:
			del cell["metadata"]["pycharm"]
	filename.dump_json(notebook, indent=1)
EOF

python3 remove_execution_times.py

# Stage
for file in *.ipynb; do
  git stage "$file"
done


# Commit and push
git commit -m "Re-rendered Jupyter Notebooks" || exit 1
git push --set-upstream origin Notebooks --force || exit 1

exit 0
