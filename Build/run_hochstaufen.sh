#!/bin/bash

./AreaLightGlintsProject.exe
for i in runtime*.csv  ; do mv $i ${i/t.csv/t_rtx4090.csv} ; done
for i in runtime*.csv  ; do mv $i ${i/d.csv/d_rtx4090.csv} ; done
