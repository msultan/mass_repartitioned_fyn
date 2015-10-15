#!/bin/bash
#$ -cwd
#$ -S /bin/bash

#$ -j y
#$ -R y
#$ -w e

#$ -m beas
#$ -pe orte 1
#$ -l h_vmem=12G
#$ -l h_rt=108::
#$ -t 1-40:1

ipengine 

