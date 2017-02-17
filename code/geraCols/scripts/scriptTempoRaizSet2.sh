for i in ../../../instancias/set2/*.dat
do
	../geraColsVRP2E.root $i N 00000 234
	../geraColsVRP2E.root $i N 10000 234
	../geraColsVRP2E.root $i N 01000 234
	../geraColsVRP2E.root $i N 00100 234
	../geraColsVRP2E.root $i N 00010 234
	../geraColsVRP2E.root $i N 11110 234
	echo ""
	echo ""
	echo ""
	echo ""
	echo ""
	echo "-----------------------------------------------------------------------------------------"
	echo ""
	echo ""
	echo ""
	echo ""
	echo ""
done
