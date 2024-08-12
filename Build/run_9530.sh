#!/bin/bash

./AreaLightGlintsProject.exe -force-device-index 0
for i in runtime*.csv  ; do mv $i ${i/t.csv/t_rtx4070m.csv} ; done
for i in runtime*.csv  ; do mv $i ${i/d.csv/d_rtx4070m.csv} ; done

./AreaLightGlintsProject.exe -force-device-index 1
for i in runtime*.csv  ; do mv $i ${i/t.csv/t_irisxe.csv} ; done
for i in runtime*.csv  ; do mv $i ${i/d.csv/d_irisxe.csv} ; done
