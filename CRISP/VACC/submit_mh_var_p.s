for x in {1..1}
do
	for y in {0..0}
	do
		export PASSx="$x", PASSy="$y"
		qsub -v PASSx,PASSy CRISP/qsub_missed.s
	done
done
