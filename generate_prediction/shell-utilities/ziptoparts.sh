#/bin/bash

dirname=tozip
zipname=xxx
splitsize=2

if ! [ -d $dirname ]; then
    echo "couldn't find the directory"
    exit 1
fi

curdir=$(pwd)
cd $dirname
numfiles=$(ls | wc -l)
splitidx=$((${numfiles}%${splitsize}?${numfiles}/${splitsize}+1:${numfiles}/${splitsize})) # ceiling
echo "Making $splitsize zip where each zip file contains at most $splitidx files..."

i=0
partnum=1
files=""
finfo=""
for content in *
do
    i=$((i+1))
    files="$files $content"
    if  [ "$i" -ne 0 ] && [ $(( $i % $splitidx )) -eq 0 ]; then
        zip  ${curdir}/${zipname}_${partnum}.zip $files
        # write zip info to file
        flist=$( echo $files | tr ' ' '\n' )
        finfo="${finfo}Zipname: ${zipname}_${partnum}.zip\n$flist\n\n"
        # reset files string and update partnum
        files=""
        partnum=$((partnum+1))
    fi
done

# for leftovers
if [ "$files" != "" ]; then
    zip  ${curdir}/${zipname}_${partnum}.zip $files
    flist=$( echo $files | tr ' ' '\n' )
    finfo="${finfo}Zipname: ${zipname}_${partnum}.zip\n$flist"
fi

rm -f "${curdir}/zipinfo.txt"
printf "$finfo" >> "${curdir}/zipinfo.txt"


#zip  ${curdir}/${zipname}.zip * -s ${splitsize}
