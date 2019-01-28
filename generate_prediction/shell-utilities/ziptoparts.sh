#/bin/bash

dirname=tozip
zipname=xxx
splitsize=5

if ! [ -d $dirname ]; then
    echo "couldn't find the directory"
    exit 1
fi

curdir=$(pwd)
cd $dirname
numfiles=$(ls | wc -l)
splitidx=$(( (numfiles + 1) / splitsize ))
i=0
partnum=1
files=""
for content in *
do
    i=$((i+1))
    files="$files $content"
    if  [ "$i" -ne 0 ] && [ $(( $i % $splitidx )) -eq 0 ]; then
        zip  ${curdir}/${zipname}_${partnum}.zip $files
        partnum=$((partnum+1))
        files=""
    fi
done
if [ "$files" != "" ]; then
    zip  ${curdir}/${zipname}.zip $files
fi


#zip  ${curdir}/${zipname}.zip * -s ${splitsize}
