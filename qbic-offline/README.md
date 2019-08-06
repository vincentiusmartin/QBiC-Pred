# QBiC-Pred Offline

## Getting ready
To run QBiC offline, we need to set up the resources needed by QBiC. These
resources can be seen in the configuration file such as given in config.ini.example.
Since these resources take quite plenty of space, be sure you have ~200GB available
in the disk.

## Set up QBiC
1. First thing to do is to copy `config.ini.example` to a file called `config.ini`.
We will then fill up all lines in the configuration file.
2. To fill `PREDDIR`, we need to download all the predictions used by imads.
These files can be downloaded from the [QBiC download page](http://qbic.genome.duke.edu/downloads).
Just scroll to "Download our 12-mer prediction tables" and download all 10 parts.
