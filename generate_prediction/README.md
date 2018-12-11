<h1> Preprocessing before generate predictions </h1>
<ol>
  <li>Parse the best model file: <tt>python3 predutils.py -p resource/data_upbm_selected_good_quality_final.csv</tt></li>
</ol>

<h1>Generate All Possible Predictions</h1>
<ol>
<li> Modify all the configurations in the gen_prediction.sh which are input directory, output directory, value of k-mer,
and div (number of memory chunks required). </li>
<li> Modify the SLURM configuration, it is located on top of gen_prediction.sh on every line started with #SBATCH
especially the array setting, it should conform to the number of TFs. </li>
<li> Run:  sbatch gen_prediction.sh </li>
</ol>

<hr />
