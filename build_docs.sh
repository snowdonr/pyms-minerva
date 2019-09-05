#!/usr/bin/env bash
rm -rf docs
cd UserGuide
rm -rf build
make html
cp -r build/html ../docs/
cd ../
touch docs/.nojekyll
touch .nojekyll