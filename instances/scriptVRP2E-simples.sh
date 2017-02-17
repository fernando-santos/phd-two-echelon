for rnd in 1111 2222 3333 4444 5555
do
	for i in *.dat
	do
		echo "../../codigos/geraColsSimples.e $i S >> ../simp_$i"
	done
done
