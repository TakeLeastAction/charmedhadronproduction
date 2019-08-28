#!/bin/bash
chmod +x Compile.py Main.sh run.sh makejob makejob-condor
python Compile.py
./Main.sh
