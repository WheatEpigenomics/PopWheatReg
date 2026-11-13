## 1. MACS2 peak calling (ATAC-seq)
cd $PATH
sample=SAMPLE

# Call ATAC peaks using MACS2 with fixed shift/extension
macs2 callpeak \
  -t ${sample}_atac_spike.bed \
  --keep-dup all \
  --extsize 200 \
  --nomodel \
  --shift -100 \
  --scale-to large \
  -g 16e9 \
  -n ${sample}


## 2. Merge IGV-filtered peaks and generate fixed windows
# Merge manually filtered peaks and split into ~200 bp windows
cat *filter.peak.bed \
  | bedtools sort -i - \
  | bedtools merge -i - \
  | awk '$3-$2 >= 200' \
  | bedtools sort -i - \
  | bedtools makewindows -b - -w 200 \
  | awk '$3-$2 > 150' \
  > Merge.filter.bed


## 3. Count ATAC fragments over merged peak windows
# Count fragment coverage for each sample
bedtools coverage \
  -a Merge.filter.bed \
  -b ${sample}_atac_spike.bed \
  -c \
  > ${sample}.filter.count


## 4. Merge count files into a peak × sample matrix
# Combine per-sample counts into a single matrix
paste a*filter.count | awk '{
    printf "%s\t%s\t%s\t%s", $1, $2, $3, $4;
    for (i = 2; i <= 207; i++) {
        printf "\t%s", $(4 * i);
    }
    printf "\n";
}' > Merge.filter.2.count


## 5. Assign unique peak IDs
awk '{print "Peak"NR-1"\t"$0}' \
  Merge.filter.2.count \
  > Merge.filter.3.count


## 6. Summarize total counts per sample
# Calculate total signal per sample across all peaks
awk '
NR == 1 {
    for (i = 2; i <= NF; i++) samples[i] = $i;
    next
}
{
    for (i = 2; i <= NF; i++) sum[i] += $i
}
END {
    for (i = 2; i <= NF; i++) print samples[i], sum[i]
}
' Merge.Peak.count > Merge.Peak.Sumcount


## 7. Background correction using FRIP coefficients
# Subtract length-scaled background signal for each sample
awk '
BEGIN {
    while (getline < "FRIPcoefficient.txt")
        factor[$1] = $2;
    close("FRIPcoefficient.txt");
}
NR == 1 {
    print;
    for (i = 1; i <= NF; i++) header[i] = $i;
    next
}
{
    len = $4 - $3;
    printf "%s\t%s\t%s\t%s", $1, $2, $3, $4;
    for (i = 5; i <= NF; i++) {
        colname = header[i];
        if (colname in factor) {
            adj = $i - (len * factor[colname]);
            printf "\t%.4f", (adj < 0 ? 0 : adj);
        } else {
            printf "\t%s", $i;
        }
    }
    print "";
}
' Merge.filter.3.count > Merge_FRIP.count


## 8. Normalize bedGraph signals and generate bigWig
# Use the sample with the largest FRIP coefficient as reference
COEFFICIENT=$(
awk '
NR == FNR && (FNR == 1 || $2 > max) { max = $2; sample = $1 }
NR > FNR { sum += $sample }
END { print sum }
' FRIPcoefficient.txt Merge_FRIP.count
)

# Scale each sample to the same total signal
while read A
do
  x=$(awk -v a=${A} -F "\t" '
    NR == 1 {
        for (i = 2; i <= NF; i++) if ($i == a) col = i
    }
    NR > 1 { sum += $col }
    END { print sum }
  ' Merge_FRIP.count)
  awk -v coef=$COEFFICIENT -v total=$x \
    '{print $1, $2, $3, $4 * coef / total}' \
    ${A}_atac_spike.bg > ${A}.bgx
  # Remove organellar signals
  sed -i '/chloroplast/d' ${A}.bgx
  sed -i '/mitochondrinon/d' ${A}.bgx
  # Convert to normalized bigWig
  wigToBigWig ${A}.bgx $GENOME_SIZE ${A}.norm.bw
done < Sample.list
