#!/bin/bash
inputlist="cpytestlist.txt"
inputfolder="/Users/vincentiusmartin/Research/MutationPredictor/QBiC-Pred/QBiC-Pred/website/resource/testcpy"
outputfolder="/Users/vincentiusmartin/Research/MutationPredictor/QBiC-Pred/QBiC-Pred/website/resource/testout"
insuffix="_8mer" # if needed suffix from the input file, set here
outsuffix="_escore"

while read line || [ -n "$line" ]; do
    extension="${line##*.}"
    filename="${line%.*}"
    tosearch="$filename$insuffix.$extension"
    path="$(find $inputfolder -name "$tosearch")"
    echo "Copying $path"
    cp $path $outputfolder
done < $inputlist

# if rename is needed:
echo "renaming files"
for f in $outputfolder/*; do mv "$f" "$(echo "$f" | sed s/${insuffix}/${outsuffix}/)"; done
