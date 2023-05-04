inputDir=/Users/calebcranney/Desktop/truncatedUmiInput/
matFilePath=/Users/calebcranney/Documents/Projects/JGEnglishLab/longread_umi_python/dependencies_download/promethion.mat
outputDir=/Users/calebcranney/Desktop/alignmentOutput/results
count=0
for file in $inputDir*
do
  outputFile="$outputDir""$count"".fasta"
  ((count++))
  lamassemble $matFilePath --end -g60 -m 40 -f fastq -a "$file" > $outputFile
done