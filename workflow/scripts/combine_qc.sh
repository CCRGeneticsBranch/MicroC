#!/bin/bash

# Extract columns from each qc.txt file, then combine

files=(*.qc.txt)
tmpfiles=()

# Extract columns for each file
for f in "${files[@]}"; do
    prefix="${f%%.qc.txt}"
    tmp="${prefix}.tmp"
    awk '{printf "%-46s\t%s\t%s\n", substr($0,1,46), substr($0,47,11), substr($0,58)}' "$f" > "$tmp"
    tmpfiles+=("$tmp")
done

# Build header
header1=""
header2="Metric"
for f in "${files[@]}"; do
    prefix="${f%%.qc.txt}"
    header1="${header1}\t${prefix}\t"
    header2="${header2}\tValue\tPercent"
done
echo -e "$header1" > combined_qc.txt
echo -e "$header2" >> combined_qc.txt

# Combine extracted files
paste "${tmpfiles[@]}" | awk -F'\t' '{
    out = $1
    for(i=2;i<=NF;i+=3) out = out "\t" $(i) "\t" $(i+1)
    print out
}' >> combined_qc.txt

# Clean up temporary files
rm -f "${tmpfiles[@]}"
