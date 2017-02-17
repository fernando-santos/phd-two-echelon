for rnd in 1111 2222 3333 4444 5555
do
	for i in *.dat
	do
		echo "../../codigos/geraColsSimplesPR.e $i S >> ../simpPR_$i"
	done
done
