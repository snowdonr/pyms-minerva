#!/usr/bin/env bash
rm -rf build
rm -rf pyms-docs
make html
cp -r build/html pyms-docs/