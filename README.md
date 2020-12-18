# musr-espresso
Quantum Espresso Tools for DFT+mu

## musr-espresso python package
The musr-espresso python package can run DFT calculations to calculate appropriate cutoff parameters, calculate the dependence of the lattice parameter 
on the energy, and set up a DFT+mu calculation, ready for sumbission on Oxford's ARC system.

To install, move to this folder on your system and run
```bash
pip install ./
```
That's it! All examples should now work!

### setting up musr-espresso on ARC
To use musr-espresso on ARC, you first have to create a virtual python environment. Run the following commands:
```bash
module load python/anaconda3/2019.03
export CONPREFIX=$DATA/espresso_env
conda create --prefix $CONPREFIX --copy python=3.5
source activate $CONPREFIX
conda install pip
cd [musr-espresso package location]
python -m pip install ./
```
(This creates the python environment, and installs musr-espresso. This only needs to be done once.

Then, to use the musr-espresso scripts, make sure your slurm script has the lines
```bash
module load python/anaconda3/2019.03
export CONPREFIX=$DATA/espresso_env
source activate $CONPREFIX
```

You will need to change the example SLURM script to match the conda environment created above for it to work.

## pwprocess.sh and cancelruns.sh
These are BASH scripts that is designed to run on the ARC, which monitors the progress of current runs, resumbits those which have run out of time, and reports the 
runs that have failed, and provides an interactive interface to cancel the runs already on the system. To use these, add them to a folder on the ARC system (also adding that folder to PATH), cd to the folder where the espresso runs are stored, 
and run
```bash
pwprocess.sh
```
(I would recommend editing ./.bash_profile and adding the following lines:
```bash
alias processruns="$HOME/processruns.sh"
alias cancelruns="$HOME/cancelruns.sh"
```
Then you can use the commands processruns and cancelruns from anywhere!

## PostProcessing
The folder postprocessing/ contains a few tools for the command line, which analyse the results of a completed DFT+mu calculation. 
Each of these have a help function, e.g running
```bash
pwcif.py --help
```
returns
```
usage: pwcif.py [-h] [-sg SPACEGROUP] [-em EMAX] [-c] -uc UNITCELLCIF [-i]
                csv_file

pw.x output file processor to find the muon site

positional arguments:
  csv_file              CSV file output of pwprocess.py which has the details
                        of all the pw.x results

optional arguments:
  -h, --help            show this help message and exit
  -sg SPACEGROUP, --spacegroup SPACEGROUP
  -em EMAX, --emax EMAX
                        Maximum energy difference (in Kelvin) from the
                        smallest run to put in CSV file
  -c, --cluster         Do only one muon from each cluster
  -uc UNITCELLCIF, --unitcellcif UNITCELLCIF
                        CIF file of the original structure to put the final
                        muon sites into
  -i, --initial         plot the initial muon sites
```

## old
This folder contains the old scripts, most of which have since morphed into this package. These contain a lot of spaghetti, broken dependencies, bad practices 
etc so use at your own risk!
