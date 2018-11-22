<h1>Generate All Possible Predictions</h1>
<p>
1. Modify all the configurations in the gen_prediction.sh which are input directory, output directory, value of k-mer,
and div (number of memory chunks required).
2. Modify the SLURM configuration, it is located on top of gen_prediction.sh on every line started with #SBATCH
especially the array setting, it should conform to the number of TFs.
3. Run:  sbatch gen_prediction.sh
</p>

<hr />
