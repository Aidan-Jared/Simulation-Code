#!/bin/bash

#SBATCH -N 1
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --job-name=simulation1
#SBATCH --output=simulation_1.out
#SBATCH --mail-type=end 
#SBATCH --mail-user=Aidan.Jared@colorado.edu

./edisk