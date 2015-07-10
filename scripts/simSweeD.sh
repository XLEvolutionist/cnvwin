#!/bin/bash -l
#!/bin/bash
#SBATCH -D /group/jrigrp4/cnvwin/sims
#SBATCH -o /group/jrigrp4/cnvwin/logs/sfs_out_log-%j.txt
#SBATCH -e /group/jrigrp4/cnvwin/logs/sfs_err_log-%j.txt
#SBATCH -J SimSweeD
#SBATCH --mem-per-cpu=11000
#SBATCH --cpus-per-task=1

##Simon Renny-Byfield, UC Davis, July 2015

echo "Starting Job:"
date


cmd="SweeD -name sims -input sim.data.ms -grid 2000" 
echo $cmd
eval $cmd

echo "Ending Job: "
date
