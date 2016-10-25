#!/bin/bash

# the blastpgp and NR database
blastpgp=/usr/bin/psiblast
NR=/home/pub/blast/NR/nr
if [ ! -f $blastpgp ]; then echo "blastpgp is not correctly set"; exit 1; fi
if [ ! -f $NR.pal ]; then echo "NR database is required"; exit 1; fi

########
xdir=$(dirname $0)
if [ $# -lt 1 ]; then
	echo "usage: $0 *"
	echo "   required: pro1.seq or pro1.pssm"
	exit 1
fi

for seq1 in $*; do
	pro1=$(basename $(basename $seq1 .seq) .pssm)
	[ -f $pro1.spd3 ] && continue
	if [ ! -f $pro1.pssm ]; then
		psiblast -db $NR -num_iterations 3 -num_alignments 1 -num_threads 16 -query $pro1.seq -out_ascii_pssm $pro1.pssm > /dev/null
#		$blastpgp -d $NR -j 3 -b 1 -a 16 -i $pro1.seq -Q $pro1.pssm > /dev/null
	fi

	$xdir/pred_pssm.py $pro1.pssm
done

