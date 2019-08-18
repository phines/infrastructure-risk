for x in {1..1}
do
	for y in {9..18}
	do
		export PASSx="$x", PASSy="$y"
		qsub -v PASSx,PASSy CRISP/qsub_mh.s
	done
done
