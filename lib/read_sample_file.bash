# a bashfunction to read the sample file / sheet to memory via associative array,
# instead of relying on csvtk

read_sample_file() {
	>&2 echo '  Reading sample sheet...'
	declare -gA sample_files; R1s=(); R2s=()
	while IFS=$'\t' read -r sample R1 R2 _
		[[ $sample == \#* || $sample == 'name' ]] && continue  # skip comments and header
		[[ -f $R1 ]] || { >&2 echo "  File $R2 not found"; exit 1; }
		[[ -f $R2 ]] || { >&2 echo "  File $R2 not found"; exit 1; }
		>&2 echo "  $sample:"
		>&2 echo "    $R1"
		>&2 echo "    $R2"
		# fake 2D array as table
		sample_files[$sample,R1]="$R1"
		R1s+="$R1"
		sample_files[$sample,R2]="$R2"
		R2s+="$R2"
}