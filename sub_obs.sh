#!/bin/bash
#$ -cwd
#$ -S /bin/bash

#$ -j y
#$ -R y
#$ -w e

# email address to send notices to
#$ -M $USER@stanford.edu
#$ -m beas
#$ -pe orte 1
#$ -l h_vmem=3G
#$ -N featurizeing
#$ -l h_rt=108::

source ~/.bash_profile
#cd /home/msultan/software/fah-reseeder/
cd /nobackup/msultan/research/kinase/fyn_kinase/fah_data/mass_repartitioned_fyn
echo "job starting on `date`"
echo ""

echo "purging module environment"
module purge
echo "loading modules..."

# list your module load commands here:
module load mpi

echo "done"
echo ""
module list
echo ""
echo "beginning job on `hostname`"

ipcontroller --ip="*" --log-to-file &
sleep 10 
qsub engine.sh
qsub engine.sh
#let majority of engines connect 
sleep 600
#python her2_project/featurize/pull_protein.py
#python her2_project/featurize_scripts/calculate_rmsd_wt.py
#python her2_project/featurize/featurize.py
#python her2_project/featurize_scripts/pull_protein.py
#fah_reseeder -d /nobackup/msultan/research/kinase/her_kinase/fah_data/PROJ9114/ -n 1000 -r 50 -c 10 -s 30 -p default
echo "Done"
wait 
kill -9 `cat /home/msultan/.ipython/profile_default/pid/ipcontroller.pid `
#pkill -u msultan
#ipcluster stop --profile=default
