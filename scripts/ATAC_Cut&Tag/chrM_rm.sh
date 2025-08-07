#!/usr/bin/bash
#-------------------------------------
#parameters
while getopts R:c:o:r:b: OPT; do
	case ${OPT} in
		R) rawdata_dir=${OPTARG}
			;;
		c) cores=${OPTARG}
			;;
		o) output_dir=${OPTARG}
			;;
		b) bai=${OPTARG}
			;;
		r) report_dir=${OPTARG}
			;;
		\?)
		printf "[Usage] `date '+%F %T'` -R <rawdata_dir> -c <CUP_cores> -r <report_dir> -b <bai> -o <output_dir>\n" >&2
		exit 1
	esac
done

if [ -z "${rawdata_dir}" -o -z "${report_dir}" -o -z "${cores}" -o -z "${output_dir}" ] ;then
	printf "[ERROR] `date '+%F %T'` following parameters is empty:\n-R=${rawdata_dir}\n-r=${report_dir}\n-c=${cores}\n-o=${output_dir}\n"
	exit 1
fi

#-------------------------------------
#run picard
RED='\033[0;31m'
NC='\033[0m'

mkdir -p ${output_dir}
mkdir -p ${report_dir}/chrM
> ${report_dir}/chrM/chrM.txt
echo "sample	All	chrom" >> ${report_dir}/chrM/chrM.txt

Filelist=`ls ${rawdata_dir}|grep "bam" |awk -F '.bam' '{print $1}'|uniq`
echo -e "file for processing is:\n"${RED}${Filelist}${NC}
for file in $Filelist;do
	date && echo -e "${RED}${file}${NC}"
	samtools view -h ${rawdata_dir}/${file}.bam |grep -v chrM |grep -v "*"|samtools sort -O bam -@ ${cores} -o ${output_dir}/${file}.bam
	pre=`samtools view ${rawdata_dir}/${file}.bam|wc -l`
	pos=`samtools view ${output_dir}/${file}.bam|wc -l`
	echo "${file}	$pre	$pos" >> ${report_dir}/chrM/chrM.txt
	
	if [ ! -z "${bai}" ];then
		samtools index ${output_dir}/${file}.bam
	elif [ -z "${bai}" ];then
		:
	fi
done

#-------------------------------------
#plot
script="/home/biology/bin/Script"
Rscript ${script}/chrM_rm.R -D ${report_dir}/chrM


