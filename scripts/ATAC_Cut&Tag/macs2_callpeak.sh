#!/usr/bin/bash
# input must have "input" in the file name
#-------------------------------------
#parameters
while getopts I:R:f:g:b:o:s:B: OPT; do
	case ${OPT} in
		I) input=${OPTARG}
			;;
		R) rawdata_dir=${OPTARG}
			;;
		f) pairend=${OPTARG}
			;;
		g) genome=${OPTARG}
			;;
		#b) bigwig=${OPTARG}
		#	;;
		o) output_dir=${OPTARG}
			;;
		s) size=${OPTARG}
			;;
		#B) bedgraph=${OPTARG}
		#	;;
		B) broad=${OPTARG}
			;;
		\?)
		printf "[Usage] `date '+%F %T'` -I <input> -R <rawdata_dir> -f <pairend> -g <genome> -B <broad> -o <output_dir> -s <size>\n" >&2
		exit 1
	esac
done

if [ -z "${input}" -o -z "${rawdata_dir}" -o -z "${pairend}" -o -z "${size}" -o -z "${genome}" -o -z "${output_dir}" ] ;then
	printf "[ERROR] `date '+%F %T'` following parameters is empty:\n-R=${rawdata_dir}\n-B=${broad}\n-f=${pairend}\n-s=${size}\n-g=${genome}\n-I=${input}\n-o=${output_dir}\n"
	exit 1
fi


#-------------------------------------
#run macs2
RED='\033[0;31m'
NC='\033[0m'

mkdir -p ${output_dir}
#mkdir -p ${bigwig}
#mkdir -p ${bedgraph}
if [ ${input} == "F" ];then
	Filelist=`ls ${rawdata_dir} |awk -F '.bam' '{print $1}'|uniq`
	echo -e "file for processing is:\n"${RED}${Filelist}${NC}
	for file in $Filelist;do
		date && echo -e "${RED}${file}${NC}"
		
		if [ ${broad} == "broad" ];then	
			macs2 callpeak -t ${rawdata_dir}/${file}.bam -g ${genome} -B -f ${pairend} -n ${file} --bdg --broad --outdir ${output_dir}
		else
			macs2 callpeak -t ${rawdata_dir}/${file}.bam -g ${genome} -B -f ${pairend} -n ${file} --bdg --outdir ${output_dir}
			rm ${output_dir}/${file}_summits.bed
		fi
		
		rm ${output_dir}/${file}_control_lambda.bdg
		mv ${output_dir}/${file}_treat_pileup.bdg ${output_dir}/${file}.bdg
		#bedGraphToBigWig ${output_dir}/${file}.bdg ${size} ${bigwig}/${file}.bw
		sed -i '1i\track type=bedGraph name='${file}' description='${file}' visibility=full color=0,0,0' ${output_dir}/${file}.bdg
		gzip -f ${output_dir}/${file}.bdg
		#mv ${output_dir}/${file}.bdg.gz ${bedgraph}/${file}.bdg.gz
	done

elif [ ${input} == "T" ];then
	Filelist=`ls ${rawdata_dir} |awk -F '.bam' '{print $1}'| grep "input"|awk -F '_input' '{print $1}' |uniq`
	echo -e "file for processing is:\n"${RED}${Filelist}${NC}
	for file in $Filelist;do
		date && echo -e "${RED}${file}${NC}"
		
		if [ ${broad} == "broad" ];then	
			macs2 callpeak -t ${rawdata_dir}/${file}_input.bam -g ${genome} -B -f ${pairend} -n ${file}_input --broad --bdg --outdir ${output_dir}
		else
			macs2 callpeak -t ${rawdata_dir}/${file}_input.bam -g ${genome} -B -f ${pairend} -n ${file}_input --bdg --outdir ${output_dir}
		fi
		
		mv ${output_dir}/${file}_input_treat_pileup.bdg ${output_dir}/${file}_input.bdg
		list=`ls ${output_dir}|grep ${file}_input |grep -v ${file}_input.bdg`
		rm ${list}
		#bedGraphToBigWig ${output_dir}/${file}_input.bdg ${size} ${bigwig}/${file}_input.bw
		sed -i '1i\track type=bedGraph name='${file}'_input description='${file}'_input visibility=full color=0,0,0' ${output_dir}/${file}_input.bdg	
		gzip -f ${output_dir}/${file}_input.bdg
		#mv ${output_dir}/${file}_input.bdg ${bedgraph}/${file}_input.bdg

		if [ ${broad} == "broad" ];then	
			macs2 callpeak -t ${rawdata_dir}/${file}_ip.bam -c ${rawdata_dir}/${file}_input.bam -g ${genome} -B -f ${pairend} --broad -n ${file}_ip --bdg --outdir ${output_dir}
		else
			macs2 callpeak -t ${rawdata_dir}/${file}_ip.bam -c ${rawdata_dir}/${file}_input.bam -g ${genome} -B -f ${pairend} -n ${file}_ip --bdg --outdir ${output_dir}
			rm ${output_dir}/${file}_ip_summits.bed
		fi
		
		rm ${output_dir}/${file}_ip_control_lambda.bdg
		mv ${output_dir}/${file}_ip_treat_pileup.bdg ${output_dir}/${file}_ip.bdg
		#bedGraphToBigWig ${output_dir}/${file}_ip.bdg ${size} ${bigwig}/${file}_ip.bw
		sed -i '1i\track type=bedGraph name='${file}'_ip description='${file}'_ip visibility=full color=0,0,0' ${output_dir}/${file}_ip.bdg	
		gzip -f ${output_dir}/${file}_ip.bdg
		#mv ${output_dir}/${file}_ip.bdg ${bedgraph}/${file}_ip.bdg
	done
else
	Filelist=`ls ${rawdata_dir} |awk -F '.bam' '{print $1}'|uniq`
	Filelist=${Filelist[@]/${input}/}
	echo -e "file for processing is:\n"${RED}${Filelist}${NC}
	for file in $Filelist;do
		date && echo -e "${RED}${file}${NC}"
		
		if [ ${broad} == "broad" ];then	
			macs2 callpeak -t ${rawdata_dir}/${file}.bam -c ${rawdata_dir}/${input}.bam -g ${genome} -B -f ${pairend} --broad -n ${file} --bdg --outdir ${output_dir}
		else
			macs2 callpeak -t ${rawdata_dir}/${file}.bam -c ${rawdata_dir}/${input}.bam -g ${genome} -B -f ${pairend} -n ${file} --bdg --outdir ${output_dir}
			rm ${output_dir}/${file}_summits.bed
		fi

		rm ${output_dir}/${file}_control_lambda.bdg
		mv ${output_dir}/${file}_treat_pileup.bdg ${output_dir}/${file}.bdg
		#bedGraphToBigWig ${output_dir}/${file}.bdg ${size} ${bigwig}/${file}.bw
		sed -i '1i\track type=bedGraph name='${file}' description='${file}' visibility=full color=0,0,0' ${output_dir}/${file}.bdg	
		gzip -f ${output_dir}/${file}.bdg
		#mv ${output_dir}/${file}.bdg.gz ${bedgraph}/${file}.bdg.gz
	done
	
	rm ${output_dir}/${input}*
fi

#-------------------------------------
#plot
Script="/home/biology/bin/Script"
if [ ${broad} == "broad" ];then	
	Rscript ${Script}/macs2_callpeak.R -d ${output_dir} -s broadPeak -o ${output_dir}
else
	Rscript ${Script}/macs2_callpeak.R -d ${output_dir} -s narrowPeak -o ${output_dir}
fi


