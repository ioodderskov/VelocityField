#!/bin/bash

# This script should run the Wojtak_plot-routine on all regular analyses, and the contour_plot-routine on the skyfraction analysis
# The arguments are:
#			 runname
#			 redshift
#			 N_observers
#			 x-axis
#			 output directory for plots
#			 compare with Planck512? 0: no, 1: yes
#			 write to tabular? 0: no, 1: yes
#			 make band with observed value? 0: no, 1: yes
#			 Plottype

rm tabel.txt

exe_wojtak=/home/io/Dropbox/PHD/Python/wojtak_plot.py 
exe_lightcone=/home/io/Dropbox/PHD/Python/lightcone.py 
#outdir=/home/io/Dropbox/hubble2013
outdir=/home/io/Dropbox/PHD/Python/output
Nobs=589
#python $exe_wojtak Planck512		 1	 $Nobs rmax	 $outdir 0 0 1 wojtak_full
#mv $outdir/Planck512.pdf $outdir/Planck512_band.pdf
#python $exe_wojtak Planck512_RandPos	 1	 $Nobs rmax	 $outdir 1 1 0 wojtak_half
#python $exe_wojtak Planck512_RandHalos	 1	 $Nobs rmax	 $outdir 1 1 0 wojtak_half
#python $exe_wojtak Planck512		 1	 $Nobs rmax	 $outdir 0 1 0 wojtak_half
#python $exe_wojtak Planck1024		 1	 $Nobs rmax	 $outdir 1 1 0 wojtak_half
#python $exe_wojtak Planck1024_box1024	 1	 $Nobs rmax	 $outdir 1 1 0 wojtak_half
python $exe_lightcone 					 $outdir wojtak_full
#python $exe_wojtak Planck512_varSN	 1	 $Nobs numSN	 $outdir 0 0 0 wojtak_half
#python $exe_wojtak Planck512_seed	 1	 $Nobs rmax	 $outdir 1 1 0 wojtak_half


#python $exe_wojtak Wmap512		 1	 $Nobs rmax	 $outdir 1 1 0
####python $exe_wojtak Planck512_z01	 0.91	 $Nobs rmax	 $outdir 1 1 0
####python $exe_wojtak Planck512_z02	 0.83	 $Nobs rmax	 $outdir 1 1 0

exe_skyfraction=/home/io/Dropbox/PHD/Python/skyfraction.py
exe_skyfraction_doubleplot=/home/io/Dropbox/PHD/Python/skyfraction_doubleplot.py
dir_1cone=/home/io/Desktop/PHD/Hubble/output/skyfraction_Planck512_1cone/
dir_2cones=/home/io/Desktop/PHD/Hubble/output/skyfraction_Planck512_2cones/
H_full=/home/io/Desktop/PHD/Hubble/output/Planck512/Hubbleconstants.txt

# For the skyfraction-plots, the arguments are:
#			 input-folder for skyfraction 
#			 input-folder for Planck512
#			 N_observers
#			 output directory for plots
#			 write to tabular? 0: no, 1: yes


#python $exe_skyfraction $dir_1cone $H_full $Nobs $outdir/skyfrac_1cone.pdf 1
#python $exe_skyfraction $dir_2cones $H_full $Nobs $outdir/skyfrac_2cones.pdf 1

#python $exe_skyfraction_doubleplot

# HALOPLOTS
exe_haloplots=/home/io/Dropbox/PHD/Python/generate_haloplots.py

#python $exe_haloplots

