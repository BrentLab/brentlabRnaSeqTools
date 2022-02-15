#!/bin/bash

# place this in your $HOME and follow these directions to set up a cron job
# https://help.ubuntu.com/community/CronHowto

export PATH=/opt/slurm/bin:/opt/slurm/sbin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

srun --chdir /scratch/mblab/$USER find /scratch/mblab/$USER -exec touch {} +
