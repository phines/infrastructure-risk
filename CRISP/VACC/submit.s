for x in {1..1}
do
	for y in {0..2}
	do
		export PASSx="$x", PASSy="$y"
		qsub -v PASSx,PASSy CRISP/qsub_gen.s
	done
done
