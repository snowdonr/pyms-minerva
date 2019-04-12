#!/usr/bin/env bash
rm -rf build
rm -rf docs
make html
cp -r build/html docs/