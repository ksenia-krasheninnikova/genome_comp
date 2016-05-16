#!/bin/bash

if [ "$#" -ne 3 ]; then 
    echo "Usage: " $0 " blocks specie1 specie2"
    echo "blocks - folder with ragout blocks folders. each folder is named with blocks scale"
    echo "count translocations in specie2 according to the genome order in specie1"
    exit 1
fi

FOLDER="$1"
SP1="$2"
SP2="$3"
KINDS="
transpositions
translocations
reversals
duplications
"

echo $SP1
echo $SP2
echo $FOLDER/"breakpoints_"$SP1"_"$SP2.stats
$(echo -e 'kind\tscale\tnumber\tunresolved_paths' > $FOLDER/"breakpoints_"$SP1"_"$SP2.stats);
for e in $(ls -d $FOLDER/[0-9]*); do
    p=${e%%};
    echo $(ls $p/"blocks_coords.txt");

    TMP=$p"/tmp"
    for kind in $KINDS; do
        ./synteny_blocks/breakpoints_analyzer.py $p/"blocks_coords.txt" --report_$kind --species $SP1 $SP2 > $p/$kind"_"$SP1"_"$SP2".txt";
        u=$(cat $p/$kind"_"$SP1"_"$SP2".txt" | grep overall | awk '{sum += $NF} END {print sum}');
        un=$(cat $p/$kind"_"$SP1"_"$SP2".txt" | grep unresolved | cut -f2 -d' ');
        $(echo -e $kind'\t'$(basename $p)'\t'$u'\t'$un >> $FOLDER/breakpoints_$SP1"_"$SP2.stats);
    done
done
echo 'results in '$FOLDER/breakpoints_$SP1"_"$SP2.stats

