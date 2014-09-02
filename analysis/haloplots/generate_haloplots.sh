#!/bin/bash

# The arguments are:
#	inputfile
#	outputfile
#	cones: 0, 1 or 2
#	lighcone?: 0: no, 1: yes
#	box
#	name (line 1 of plot title)
#	title (line 2 of plot title)

exe=haloplots.py
indir=/home/io/Desktop/PHD/Hubble/halofinders/Rockstar
outdir=/home/io/Dropbox/hubble2013

echo "A.0 and A.1: Standard simulation\nRandom and Local Group like observers"

python $exe $indir/Planck512/parents_11 $outdir/haloplot_Planck512.png 0 0 512 "A.0, A.1 and A.2: Standard simulation" "Random and Local Group like observers"
python $exe $indir/Planck512/parents_11 $outdir/haloplot_Planck512_1cone.png 1 0 512 "A.3: Standard simulation" "One cone fraction"
python $exe $indir/Planck512/parents_11 $outdir/haloplot_Planck512_2cones.png 2 0 512 "A.4: Standard simultion" "Two cones fraction"
python $exe $indir/Planck512_lightcone/parents_0 $outdir/haloplot_Planck512_lightcone.png 0 1 512 "A.5: Standard simulation" "Halos from past lightcone of observer"
#python $exe $indir/Planck512_seed/parents_11 $outdir/haloplot_Planck512_seed.png 0 0 512 "B: Simulation with different seed" " "
#python $exe $indir/Wmap512/parents_11 $outdir/haloplot_Wmap512.png 0 0 512 "C: Simulation with WMAP parameters" " "
#python $exe $indir/Planck1024/parents_11 $outdir/haloplot_Planck1024.png 0 0 512 "D: Simulation with Nsim=1024" " "
#python $exe $indir/Planck1024_box1024/parents_11 $outdir/haloplot_Planck1024_box1024.png 0 0 1024 "E: Simulation with Nsim=1024 and box=1024Mpc/h" " "

# Convert the png-files to pdf's
cd $outdir
convert haloplot_Planck512.png haloplot_Planck512.pdf
convert haloplot_Planck512_1cone.png haloplot_Planck512_1cone.pdf
convert haloplot_Planck512_2cones.png haloplot_Planck512_2cones.pdf
convert haloplot_Planck512_lightcone.png haloplot_Planck512_lightcone.pdf
#convert haloplot_Planck512_seed.png haloplot_Planck512_seed.pdf
#convert haloplot_Wmap512.png haloplot_Wmap512.pdf
#convert haloplot_Planck1024.png haloplot_Planck1024.pdf
#convert haloplot_Planck1024_box1024.png haloplot_Planck1024_box1024.pdf
