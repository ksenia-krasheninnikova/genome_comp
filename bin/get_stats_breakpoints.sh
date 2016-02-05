#!/bin/bash

if [ "$#" -ne 2 ]; then 
    echo "Usage: " $0 " blocks hal"
    echo "blocks - folder with ragout blocks folders. each folder is named with blocks scale"
    echo "HAL - file name of the alignments"
    exit 1
fi

FOLDER="$1"
HAL="$2"

#count number of blocks
$(echo -e 'genome\tscale\tnum\tcov\tmedian size' > $FOLDER/number.stats);
for e in $(ls -d $FOLDER/[0-9]*); do
    p=${e%%};
    echo $(ls $p/"blocks_coords.txt");

    genomes=$(halStats --genomes $HAL);
    for g in $genomes; do
        TMP_GENOME=tmp.bed;
        $(halStats --bedSequences $g $HAL > $TMP_GENOME);
        genome_size=$(awk '{sum += $3 - $2} END {print sum}' $TMP_GENOME);
        if [[ $g != Anc* ]]; then
            $(./ragout_blocks_to_bed.py $p/"blocks_coords.txt" --specie $g > $p/"blocks_coords_$g.bed");
            #counting number of blocks
            n=$(cat $p/"blocks_coords_$g.bed" | wc -l);
            TMP_INTERSECT=tmp.int
            $(bedtools intersect -a $p/blocks_coords_$g.bed -b $TMP_GENOME > $TMP_INTERSECT);
            cov=$(awk '{sum += $3 - $2} END {print sum}' $TMP_INTERSECT);
            if [ -z "$cov" ]; then
                cov=0;
            fi
            #counting genome coverage
            cov=$(echo print $cov/$genome_size. | python);
            #count median size in file
            median=$(sort -n $p/blocks_coords_$g.bed | awk ' { a[i++]=$3 - $2; } END { print a[int(i/2)]; }')
            #print out results
            $(echo -e $g'\t'$(basename $p)'\t'$n'\t'$cov'\t'$median >> $FOLDER/number.stats);
            $(rm $TMP_INTERSECT);
        fi
        $(rm $TMP_GENOME);
    done
done
