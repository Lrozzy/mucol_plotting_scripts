# mucol_plotting_scripts


## For plotting macros

An example script that runs over lcio files
* extracts MC truth particles, tracks, and hits on track
* fills root histograms & saves to output_name.root
* fills nested lists with key variables and saves to a json file
* remember to source setup script
```
source setup.sh

python study_alltracks.py -i input_file.slcio -o output_name
```
Additional scripts to run over lcio files live in ```macros/LCIO```

LCIO data format and functions are documented [here](https://ilcsoft.desy.de/LCIO/current/doc/doxygen_api/html/index.html)

Muon Collider geometry v0a [config file](https://github.com/madbaron/detector-simulation/blob/KITP_10TeV/geometries/MuColl_10TeV_v0A/MuColl_10TeV_v0A.xml)

## For plotting with python notebooks

Some python notebooks to run over edm4hep.root files live in ```macros```

**For analysis in a python notebook with uproot/awkward arrays only**
- I manage packages with homebrew, and used it to install anaconda 
- I install packages for anaconda like ```conda install -c conda-forge vector```
- Start a notebook with 
```
jupyter notebook
```

**If you want to use a python notebook with ROOT**
- Start a notebook with 
```
root --notebook
```
