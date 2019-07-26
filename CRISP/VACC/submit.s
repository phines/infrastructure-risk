for x in {1..1}
do
	for y in {101..200}
	do
		export PASS="$x", PASSB="$y"
		qsub -v PASS,PASSB qsub.s
	done
done
