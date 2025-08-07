#!/usr/bin/bash
#-------------------------------------
#parameters
while getopts T:R:c:g:r:i:o: OPT; do
	case ${OPT} in
		T) reads_type=${OPTARG}
			;;
		R) rawdata_dir=${OPTARG}
			;;
		c) cores=${OPTARG}
			;;
		g) gz_file=${OPTARG}
			;;
		r) report=${OPTARG}
			;;		
		i) index=${OPTARG}
			;;
		o) output_dir=${OPTARG}
			;;
		\?)
		printf "[Usage] `date '+%F %T'` -T <reads_type> -R <rawdata_dir> -c <CUP_cores> -g <gz_file> -r <report> -i <index> -o <output_dir>\n" >&2
		exit 1
	esac
done

if [ -z "${reads_type}" -o -z "${rawdata_dir}" -o -z "${cores}" -o -z "${gz_file}" -o -z "${report}" -o -z "${index}" -o -z "${output_dir}" ] ;then
	printf "[ERROR] `date '+%F %T'` following parameters is empty:\n-R=${rawdata_dir}\n-c=${cores}\n-g=${gz_file}\n-r=${report}\n-i=${index}\n-o=${output_dir}\n-T=${reads_type}\n"
	exit 1
fi

#-------------------------------------
#Single PairEnd reads function

# $1=${reads_type} 
# $2=${gz_file} 
# $3={file_name}

READS() {
	if [[ $1 == "S" && $2 == "T" ]];then
		STAR --runThreadN ${cores} --genomeDir ${index} --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --readFilesIn <(gunzip -c ${rawdata_dir}/$3.*.gz) --outFileNamePrefix ${output_dir}/${file} --outSAMstrandField intronMotif --limitGenomeGenerateRAM 150000000000
	elif [[ $1 == "S" && $2 == "F" ]];then
		STAR --runThreadN ${cores} --genomeDir ${index} --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --readFilesIn ${rawdata_dir}/$3.* --outFileNamePrefix ${output_dir}/${file} --outSAMstrandField intronMotif --limitGenomeGenerateRAM 150000000000
	elif [[ $1 == "P" && $2 == "T" ]];then
		STAR --runThreadN ${cores} --genomeDir ${index} --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --readFilesIn <(gunzip -c ${rawdata_dir}/$3_R1.fq.gz) <(gunzip -c ${rawdata_dir}/$3_R2.fq.gz) --outFileNamePrefix ${output_dir}/${file} --outSAMstrandField intronMotif --limitGenomeGenerateRAM 150000000000
	elif [[ $1 == "P" && $2 == "F" ]];then
		STAR --runThreadN ${cores} --genomeDir ${index} --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --readFilesIn ${rawdata_dir}/$3_R1.* ${rawdata_dir}/$3_R2.* --outFileNamePrefix ${output_dir}/${file} --outSAMstrandField intronMotif --limitGenomeGenerateRAM 150000000000
	fi
}

#-------------------------------------
#run STAR
RED='\033[0;31m'
NC='\033[0m'

mkdir -p ${output_dir}
mkdir -p ${report}
Filelist=`ls ${rawdata_dir}|awk -F '_1' '{print $1}'|awk -F '_2' '{print $1}' |awk -F '_R1.' '{print $1}'|awk -F '_R2.' '{print $1}'|awk -F '[.]' '{print $1}'|uniq`
echo -e "file for processing is:\n"${RED}${Filelist}${NC}
for file in $Filelist;do
	date && echo -e "${RED}${file}${NC}"
	READS ${reads_type} ${gz_file} ${file}
	mv ${output_dir}/${file}Aligned*.bam ${output_dir}/${file}.bam 
	samtools index  ${output_dir}/${file}.bam
	
	mv ${output_dir}/${file}Log.* ${report}
	mv ${output_dir}/${file}ReadsPerGene.out.tab ${report}
	mv ${output_dir}/${file}SJ.out.tab ${report}
done

#-------------------------------------
#plot
mkdir -p ${report}
Script="/home/biology/bin/Script"
#Rscript ${Script}/STAR_map.R -D ${output_dir} -o ${report}
Rscript ${Script}/STAR_map.R -D ${report} -o ${report}


