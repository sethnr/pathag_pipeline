#!/usr/bin/env bash

cores=1


genome_size="11k"

while getopts "i:o:g:m:C:" OPTION; do
    case $OPTION in
    i) infile=$OPTARG      ;;
    o) outfile=$OPTARG     ;;
    m) mashout=$OPTARG     ;;
    g) genome_size=$OPTARG ;;    
    C) cores=$OPTARG       ;;
    *)  echo "option not recognised"
        exit 1
        ;;
    esac
done
shift $(expr $OPTIND - 1 )

if [ -z "$mashout" ]; then
    FASTA=$1
    echo "generating bwa index for ${FASTA}"
    echo bwa index ${FASTA} 2>&1
    bwa index ${FASTA}
else 
    FASTAS=($@);
    echo "making mash indices for ${FASTAS}"
    
    mashidx=("${FASTAS[@]}")
    
    for i in "${!FASTAS[@]}"; do
	echo $i ":"
	echo mash sketch -g ${genome_size} ${FASTAS[$i]}
        mash sketch -g ${genome_size} ${FASTAS[$i]} ;
	echo "-" 
        echo ${FASTAS[$i]}  '-->' ${FASTAS[$i]/.fasta/.fasta.msh}
        FASTAS[$i]=${FASTAS[$i]/.fasta/.fasta.msh}
       
    done

    echo "-"
    echo "merging mash indices"
    echo mash paste $mashout ${FASTAS[*]}
    mash paste $mashout ${FASTAS[*]}
    #mash paste $mashout ${mashidx}
fi
