header="#\0041/bin/bash\n#BSUB -J genomecov\n#BSUB -o genomecov.out\n#BSUB -e genomecov.err\n#BSUB -W 120:00\n#BSUB -q bigmem\n#BSUB -n 1\n#BSUB -R \"rusage[mem=10240]\"\n"

## The code needs to be modified before using

## configure parameters
f=$1
if test $# -ge 2; then
  read_len=$2
else
  read_len=75
fi

dir=${f%/*}
name=${f##*/}
bin=/nethome/jlh267/Bin

## working
if [ $# -eq 2 ] && [ $2 == "BS-seq" ]; then
  echo -e $header > x.sh
  if [ ${f##*.} == "gz" ]; then 
    echo -e "zcat $f | awk '{if(\0041/^chrom/ && \$4>=5) {printf(\"%s\\\t%d\\\t%d\\\t%s\\\n\", \$1, \$2, \$2+1, \$3);}}' | sort --temporary-directory=${HOME}/.local/temp -k 1,1 | gzip -c > ${dir}/${name}.bdg.gz" >> x.sh
  else 
    echo -e "cat $f | awk '{if(\0041/^chrom/ && \$4>=5) {printf(\"%s\\\t%d\\\t%d\\\t%s\\\n\", \$1, \$2, \$2+1, \$3);}}' | sort --temporary-directory=${HOME}/.local/temp -k 1,1 | gzip -c > ${dir}/${name}.bdg.gz" >> x.sh
  fi
  bsub < x.sh
  rm x.sh
else
  samtools view -H $f | awk '{if(/@SQ/) {sub("SN:", "", $2); sub("LN:", "", $3); printf("%s\t%s\n", $2, $3);}}' > $dir/${name}.mm10.chromSize.txt
  genomeChr="$dir/${name}.mm10.chromSize.txt"
  quality=30

  echo -e $header > x.sh
  echo -e "nlines=\$(samtools view -q $quality -F 1804 $f | wc -l)" >> x.sh
  echo -e "rate=\`echo \"scale=10; 1000000/\$nlines\"|bc\`" >> x.sh
  echo -e "samtools view -q $quality -F 1804 $f | awk -v read_len=$read_len '{if(\$9 > 0) sign=\"+\"; else sign=\"-\"; printf(\"%s\\\t%d\\\t%d\\\t%s\\\t%d\\\t%s\\\n\", \$3, \$4, \$4+read_len, \$1, \$5, sign)}' | sort --temporary-directory=${HOME}/.local/temp -k 1,1 > ${dir}/${name}.sorted.bed" >> x.sh
  echo -e "${bin}/bedtools_v2.27.0/bin/bedtools genomecov -i ${dir}/${name}.sorted.bed -bga -split -g $genomeChr -scale \$rate > ${dir}/${name}.bdg" >> x.sh
  echo -e "rm ${dir}/${name}.sorted.bed" >> x.sh
  echo -e "gzip -f ${dir}/${name}.bdg" >> x.sh
  echo -e "rm $genomeChr" >> x.sh
  bsub < x.sh
  rm x.sh
fi
