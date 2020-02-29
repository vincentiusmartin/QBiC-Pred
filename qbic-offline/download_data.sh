# simple download script for QBiC

# reguires an argument for the datadir
DATADIR=$1

cd DATADIR

# escores
wget -N http://qbic.genome.duke.edu/download/escore.zip

# dir stuff
unzip escore
mkdir -p chromosomes
cd chromosomes
mkdir -p hg19
cd hg19

# chromosomes -- just hg19
rsync -avzP --ignore-existing rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ .

cd ..
mkdir -p preddir
cd preddir

# big stuff -- might take a while
for i in 1 2 3 4 5 6 7 8 9 10 
do
  wget -N http://qbic.genome.duke.edu/download/predzip/prediction_$i.zip
done

for i in 1 2 3 4 5 6 7 8 9 10
do
  wget -N http://qbic.genome.duke.edu/download/predpval/predpval_$i.zip
done

unzip '*.zip'

