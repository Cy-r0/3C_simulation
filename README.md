# Chromosome Conformation Capture (3C) simulation

This is a simple simulation of how chromatin contact maps are generated using
Hi-C, a 3C-based method.


## What it does

A single chromatin filament is randomly generated with a number of overlapping
regions.

Crosslinking is then performed so the overlapping regions are permanently stuck
together.

A restriction enzyme cuts the filament into roughly equally-sized segments:
in real life, this would be done by something like HindIII; here, we use a
fictitious enzyme that cuts the filament in segments of equal length (specified
by the user).

Then the overlapping regions are ligated, and this creates two types of
chromatin fragments: fragments that were crosslinked and ligated, and isolated
fragments.

The final output of the simulation is a red and white matrix that shows whether
any two fragments in the filament were overlapping (and thus crosslinked and
ligated) or not.


## Limitations

For simplicity, there is only one euchromatin filament. No representation of
base pairs and other features of DNA are present.

The shape and arrangement of the filament is mostly random, which is not the
case in real cells.

The simulated filament is only 2-dimensional, while real chromatin is organised
in 3D space.

Our fictional restriction enzyme cuts chromatin into fragments of exactly the
same size, while real enzymes recognise and cleave specific sequences.

The simulation is done on a single cell, while real Hi-C measures overlapping
regions in a population of cells to account for stochastic variations.


## Installation and setup (Linux)

git clone this repository on your machine
```
git clone git@github.com:Cy-r0/3C_simulation
cd 3C_simulation
```

Create and activate a virtual environment
```
python3 -m venv 3C_env
source 3C_env/bin/activate
```

Install python dependencies
```
python3 -m pip install -r requirements.txt
```

To run the simulation from a terminal:
```
./simulate.py
```
