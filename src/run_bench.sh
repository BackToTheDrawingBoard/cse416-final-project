for benchmark in "Nibble-Parallel" "ACS-Sync-Local-Opt"; do
	for i in 32 24 16 12 8 6 4 2 1; do
		for j in " " "-c"; do
			echo -n "$benchmark CILK CORES: $i"

			suffix=""
			file_suffx="u"
			if [ "$j" = "-c" ]; then
				echo -e "\tCOMPRESSED"
				suffix=".compressed"
				file_suffix="c"
			else
				echo
			fi

			d="./${benchmark}_logs"
			mkdir -p "$d"
			filename="${d}/${i}${file_suffix}"

			echo \
	"CILK_NWORKERS=$i $suffix ./Nibble-Parallel
		-r 20000 -s $j ../../inputs/com-orkut.ungraph.txt$suffix" > $filename
			CILK_NWORKERS=$i ./Nibble-Parallel \
				-r 20000 -s $j ../../inputs/com-orkut.ungraph.txt$suffix \
					>> $filename
		done
	done
done
