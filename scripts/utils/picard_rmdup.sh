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
mkdir -p ${report_dir}/picard
Filelist=`ls ${rawdata_dir}|grep "bam" |awk -F '.bam' '{print $1}'|uniq`
echo -e "file for processing is:\n"${RED}${Filelist}${NC}
for file in $Filelist;do
	date && echo -e "${RED}${file}${NC}"
	java -jar /home/biology/bin/tools/picard/build/libs/picard.jar MarkDuplicates -I ${rawdata_dir}/${file}.bam -O ${output_dir}/${file}_rmdup.bam -REMOVE_DUPLICATES true -M ${report_dir}/picard/${file}_rmDup.txt
	samtools sort -@ ${cores} -o ${output_dir}/${file}_sorted.bam ${output_dir}/${file}_rmdup.bam
	rm ${output_dir}/${file}_rmdup.bam
	mv ${output_dir}/${file}_sorted.bam ${output_dir}/${file}.bam
	
	if [ ! -z "${bai}" ];then
		samtools index ${output_dir}/${file}.bam
	elif [ -z "${bai}" ];then
		:
	fi
done

#-------------------------------------
#plot
script="/home/biology/bin/Script"
Rscript ${script}/picard_rmdup.R -D ${report_dir}/picard


