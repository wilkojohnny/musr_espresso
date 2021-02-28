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
conda create --prefix $CONPREFIX --copy python=3.7
source activate $CONPREFIX
conda install pip
cd [musr-espresso package location]
python -m pip install ./
```
(This creates the python environment, and installs musr-espresso. This only needs to be done once.)

Then, to use the musr-espresso scripts, make sure your slurm script has the lines
```bash
module load python/anaconda3/2019.03
export CONPREFIX=$DATA/espresso_env
source activate $CONPREFIX
```

You will need to change the example SLURM script to match the conda environment created above for it to work.

## BASH scripts
There are some useful BASH scripts in the bin/ folder, which are designed to run on the cluster's login nodes. 

To set these up on your system, look for a file called .bash_profile or .profile in your home directory, and 
add the lines
```bash
PATH="[location of the bin directory in musr-espresso]":$PATH
export PATH
```
(ignoring the square brackets). From then on, you should be able to use these commands from anywhere!

### processruns
Processruns is a script which searches through a directory of DFT+mu files, resubmits those that have run out of time,
and reports on those that have completed or failed. To use it, all you have to do is move to whatever directory you 
want to process, and run processruns. 

### cancelruns
Cancelruns is an interactive job canceller for slurm jobs. From anywhere, run cancelruns and it will go through each job, 
asking if you want to cancel it. If you only want to stop jobs with a specific name (e.g PbF2), running 
```bash
cancelruns PbF2
```
will go through only the jobs with names starting with that. You can then press 'a' to cancel all.

### ezq
EasyQueue (ezq) is a script designed to run on the Glamdring cluster. To use, makes sure the options in config/ezq.conf are correct
for your system (they probably won't be!), and to submit a job to the cluster just run
```bash
ezq NaF.scf.1.pwi
```
which will add the pw.x input file NaF.scf.1.pwi to the job queue. Several options exist so that you can control the specfics of the 
parallelization; running --help will list them out. 

ezq works with both pw.x input files (with extension pwi) and python files (extension .py).

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
