data_files="`find logs/*`"

for file in $data_files; do
	grep -P '(CILK|time)' "$file"
	echo -e "\n"
done
