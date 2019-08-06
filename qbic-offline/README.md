# QBiC-Pred Offline

## Getting ready
To run QBiC offline, we need to set up the resources needed by QBiC. These
resources can be seen in the configuration file such as given in config.ini.example.
Since these resources take quite plenty of space, be sure you have ~200GB available
in the disk.

## Python package requirements
1. pandas
2. numpy
3. scipy

## Setting up QBiC offline
1. First thing to do is to copy `config.ini.example` to a file called `config.ini`.
We will then fill up all lines in the configuration file.
2. To fill `PREDDIR`, we need to download all the predictions used by imads.
These files can be downloaded from the [QBiC download page](http://qbic.genome.duke.edu/downloads).
Just scroll to "Download our 12-mer prediction tables", download all 10 parts,
unzip them in a single directory. Then put the path of this directory to `PREDDIR`
in `config.ini`.
3. For `CHRDIR`, we need to prepare directory for the chromosome files which is
available in [QBiC download page](http://qbic.genome.duke.edu/downloads) in
"Download chromosome files". This directory should contain a folder called "hg19"
and another called "hg38" that contain its respective chromosome files. If only
one version is needed, it is fine to just have that one version saved. Then
please fill `CHRDIR` with the path to the directory with hg19 and/or hg38 folders.
4. For `ESCORE_DIR`, we need to fill this with the path to a directory that
contains all files extracted from [QBiC download page](http://qbic.genome.duke.edu/download/escore.zip).
5. Then we need to fill the mapping configurations by filling all the fields
with the paths to the following files: [PBM_HUGO_MAPPING](https://github.com/vincentiusmartin/QBiC-Pred/blob/master/website/resources/mappingdata/hugotopbm.txt) and
[HUGO_PBM_MAPPING](https://github.com/vincentiusmartin/QBiC-Pred/blob/master/website/resources/mappingdata/pbmtohugo.txt).

## Running QBiC offline
Below is an example of a command to run QBiC offline:
> `python3 qbic.py -i testing_resources/input_mutation_test.vcf -g testing_resources/tflist_test.txt -c hg19 -o predictons.tsv`.

The first parameter `-i` is the input mutation file and `-g` is a file with all
TF genes of interest--please see the full list of the supported TF genes [here](https://github.com/vincentiusmartin/QBiC-Pred/blob/master/website/resources/TF_names_v_1.01.txt). Then `-c` as the chromosome version and `-o` for the output directory.
For the full list of available arguments, do `python3 qbic.py --help`.
