#!/usr/bin/bash
#-------------------------------------
#parameters
while getopts T:R:c:i:b:o:r: OPT; do
	case ${OPT} in
		T) reads_type=${OPTARG}
			;;
		R) rawdata_dir=${OPTARG}
			;;
		c) cores=${OPTARG}
			;;
		i) index=${OPTARG}
			;;
		b) bai=${OPTARG}
			;; 
		o) output_dir=${OPTARG}
			;;
		r) report_dir=${OPTARG}
			;;
		\?)
		printf "[Usage] `date '+%F %T'` -T <reads_type> -b <bai> -R <rawdata_dir> -c <CUP_cores> -r <report_dir> -o <output_dir>\n" >&2
		exit 1
	esac
done

if [ -z "${reads_type}" -o -z "${rawdata_dir}" -o -z "${report_dir}" -o -z "${cores}" -o -z "${index}" -o -z "${output_dir}" ] ;then
	printf "[ERROR] `date '+%F %T'` following parameters is empty:\n-R=${rawdata_dir}\n-r=${report_dir}\n-b=${bai}\n-c=${cores}\n-i=${index}\n-o=${output_dir}\n-T=${reads_type}\n"
	exit 1
fi

#-------------------------------------
#Single PairEnd reads function

# $1=${reads_type} 
# $2=${gz_file} 
# $3={bai}
# $4={report_file}
READS() {
	if [ $1 == "S" ];then
		bowtie2 --rg-id $2 --rg "PL:NA" --rg "SM:$2" -x ${index} -p ${cores} -U ${rawdata_dir}/$2.* 2> ${report_dir}/bowtie2/$2.log | samtools view -@ ${cores} -bS | samtools sort -@ ${cores} -o ${output_dir}/$2.bam
		#samtools index ${output_dir}/$2.bam
	elif [ $1 == "P" ];then
		bowtie2 --rg-id $2 --rg "PL:NA" --rg "SM:$2" -x ${index} -p ${cores} -1 ${rawdata_dir}/$2_R1.* -2 ${rawdata_dir}/$2_R2.* 2> ${report_dir}/bowtie2/$2.log | samtools view -@ ${cores} -bS | samtools sort -@ ${cores} -o ${output_dir}/$2.bam
		#samtools index ${output_dir}/$2.bam
	fi
	
	if [ ! -z "$3" ];then
		samtools index ${output_dir}/$2.bam
	elif [ -z "$3" ];then
		:
	fi
}

#-------------------------------------
#run bowtie2
RED='\033[0;31m'
NC='\033[0m'

mkdir -p ${output_dir}
mkdir -p ${report_dir}/bowtie2
Filelist=`ls ${rawdata_dir} |awk -F '_R1.' '{print $1}'|awk -F '_R2.' '{print $1}'|awk -F '[.]' '{print $1}'|uniq`
echo -e "file for processing is: "${RED}${Filelist}${NC}
for file in $Filelist;do
	date && echo -e "${RED}${file}${NC}"
	READS ${reads_type} ${file} ${bai}
done

#-------------------------------------
#plot
script="/home/biology/bin/Script"
Rscript ${script}/bowtie2_map.R -D ${report_dir}/bowtie2 -t ${reads_type}


