#!/bin/bash
#set -eu
#
#  Copyright (c) 2018, U900, Institut Curie
#  Copyright (c) 2018, Philippe La Rosa (Equipe HPC)
#
# run_install.sh : This file is part of NGS=>Nextflow Pipeline
#
# Objectifs :
# Generateur d'architecture et d'environnement (espace d'exec) et lancement
# du pipeline via le WorkFlow Manager NextFlow. 
# Effectue la creation des repertoires et fichiers qui sont indispensables 
# au lancement du pipeline via le(s) script(s) nextflow, lancé sur le cluster
# de calcul ou en local, pour les differentes étapes du pipeline. 
#
# Contexte d'usages :
# Dans le cadre d'un deploiement dans l'espace cible d'installation canonique
# c.a.d theoriquement installé dans un repertoire donné. 
# Ce script sera lancé à partir d'un des noued de soummission par n'importe quel
# user ayant bien sur les droits d'acces aux data (fastq.gz, references, annotations, ...)
#
# Parametres Entree :
# -p : le nom du projet (exemple : SAREK)
# -r : le nom du run (exemple : L153)
# -s : chemin complet sampleplan         
# -i : chemin complet du fichier target.bed
# -u : profiles (exemple : singularityPath,cluster)         
# -g : genomes (exemple : hg19)         
# -b : chemin complet du repertoire genomesbase          
# -o : racine chemin du repertoire ou seront sauvegarde l'ensemble des
#      fichiers pour toutes les etapes : d'integration initiale et d'analyses 
# -q : nom de la file (exemple : dev) 
# -m : chemin complet du repertoire resultats des metriques 
# [-h] : usage
#
# Sorties :
# Soumission au cluster via qsub, ou en local  :
# affichage console (stdout) exemple pour le run V295 du projet SAREK :
# run_main.sh : create_nxf_script for projet SAREK and run V295 :
# /data/tmp/NGS_RUN_TEMP/SAREK_V295_1529940592001/SAREK.nf
# make_run_pipeline_create_nxf_work_dir for projet SAREK and run V295
# WORK_DIR = /data/tmp/NGS_RUN_TEMP/SAREK_V295_1529940592001 
# REPORTING_DIR = /data/tmp/NGS_RUN_TEMP/SAREK_V295_1529940592001/results
# make_run_pipeline_create_nxf_configs for projet SAREK and run mapping :
#

SAREK_DEBUG=${SAREK_DEBUG:=0}; [[ "$SAREK_DEBUG" == 'x' ]] && set -x

if [[ $TERM && $TERM != 'dumb' ]]
  then
    if command -v tput &>/dev/null
      then
        GREEN=$(tput setaf 2; tput bold)
        YELLOW=$(tput setaf 3)
        RED=$(tput setaf 1)
        NORMAL=$(tput sgr0)
    fi
fi

function echo_red() {
    >&2 echo -e "$RED$*$NORMAL"
}

function echo_green() {
    echo -e "$GREEN$*$NORMAL"
}

function echo_yellow() {
    >&2 echo -e "$YELLOW$*$NORMAL"
}

function die() {
  echo_red "$*"
  exit 1
}

### usage ###
function usage (){
    echo -e "\nUsage: $0"
    echo -e "\n [Options]"    
    echo -e "\t-p : nom du projet (exemple : SAREK)"     
    echo -e "\t-r : nom du run (exemple : 10X)"  
    echo -e "\t-s : chemin complet sampleplan"         
    echo -e "\t-i : chemin complet du fichier target.bed"         
    echo -e "\t-u : profiles (exemple : singularityPath,cluster)"         
    echo -e "\t-g : genomes (exemple : hg19)"         
    echo -e "\t-b : chemin complet du repertoire genomesbase "         
    echo -e "\t-o : racine chemin du repertoire ou seront sauvegarde l'ensemble des fichiers pour toutes les etapes : d'integration initiale et d'analyses" 
    echo -e "\t-q : nom de la file (exemple : dev)" 
    echo -e "\t-m : chemin complet du repertoire resultats des metriques" 
    echo -e "\t-h : usage"          
    echo -e "\n\n [Example]: \n\t# $0 -p SAREK-2.3 -r TEST -g smallGRCh37 -b References/smallGRCh37 -s --sample Sarek-data/HPC-bench/tsv/HPC-bench-test.tsv -u toolsPath,test -i Sarek-data/testdata/target.bed -m /data/results/test/sarek -q false -o /local/scratch/NGS_RUN_TEMP"
    exit 1
}

((!$#)) && echo "Il n'y a pas d'arguments!!" && usage

if [[ ($# < 19) || ($# > 20) ]]
then
    echo "Nombre d'arguments incorrect ($#) !!"
    usage
fi

while getopts "p:r:s:i:b:g:o:u:q:m:h" optionName; do
case "$optionName" in

p) project_name="$OPTARG";;
r) run="$OPTARG";;
s) sampledata="$OPTARG";;
i) samplepath="$OPTARG";;
u) profile="$OPTARG";;
g) genome="$OPTARG";;
b) genomebase="$OPTARG";;
o) racine_work_dir="$OPTARG";;
q) queue="$OPTARG";;
m) results_dir="$OPTARG";;
h) usage;;
*) usage;;
esac
done

function pipeline_create_env() {
   echo " 1) create_env"
   # path par defaut du repertoire racine d'exec du pipeline
   RACINE_WORK_DIR_DEF="/data/tmp/NGS_RUN_TEMP"
   [[ ! $racine_work_dir ]] && racine_work_dir=${RACINE_WORK_DIR_DEF}
   WORK_DIR=${WORK_DIR:="${racine_work_dir}/${project_name}_${run}_${name_end}"}
   LOGNAME=${LOGNAME:="inconnu"}
   [[ ! $RACINE_PIPELINES_DIR ]] && RACINE_PIPELINES_DIR="/bioinfo/pipelines/sandbox/dev/Sarek"
   NXF_DIR=./
   NXF_NAME=${NXF_NAME:="${script_nf}"}
   NXF_BIN_DIR=${NXF_BIN_DIR:="/bioinfo/local/build/Centos/nextflow/nextflow-19.04.0.5069"}
   CONFIGS_DIR=${CONFIGS_DIR:="conf"}
   CONFIGS_NXF_DIR=${CONFIGS_NXF_DIR:="nextflow"}
   SCRIPT_STEP_MAPPING=${SCRIPT_STEP_MAPPING:="${WORK_DIR}/run_mapping.sh"}
   SCRIPT_STEP_GERMLINE=${SCRIPT_STEP_GERMLINE:="${WORK_DIR}/run_germline.sh"}
   SCRIPT_STEP_SOMATIC=${SCRIPT_STEP_SOMATIC:="${WORK_DIR}/run_somatic.sh"}
   SCRIPT_STEP_ANNOTATE=${SCRIPT_STEP_ANNOTATE:="${WORK_DIR}/run_annotate.sh"}
   SCRIPT_STEP_MULTIQC=${SCRIPT_STEP_MULTIQC:="${WORK_DIR}/run_multiqc.sh"}
   SCRIPT_STEP_COPIE_RESULTS=${SCRIPT_STEP_COPIE_RESULTS:="${WORK_DIR}/run_copie_results.sh"}
   BIN_DIR=${BIN_DIR:="bin"}
   LIB_DIR=${LIB_DIR:="lib"}
   REF_DIR=${REF_DIR:="References"}
   REPEATES_DIR=${REPEATES_DIR:="repeats"}
   DATA_DIR=${DATA_DIR:="Sarek-data"}
   SCRIPTS_DIR=${SCRIPTS_DIR:="scripts"}
   ANNOTATION_DIR=${ANNOTATION_DIR:="Annotation"}
   DOCS_DIR=${DOCS_DIR:="docs"}
   ASSETS_DIR=${ASSETS_DIR:="assets"}
   LOCAL_LOG_DIR=${LOCAL_LOG_DIR:="LOG"}
   PIPELINES_PATH=${PIPELINES_PATH:="/bioinfo/pipelines"}
   COMM="bash"
   SUBMIT=${SUBMIT:="qsub"}
}


function pipeline_create_nxf_work_dir() {
   echo " 2) create_nxf_work_dir and configurations files"
   # creation architecture
   mkdir -p ${WORK_DIR}
   # copie des elements dans l'espace de travail
   cp -r ${RACINE_PIPELINES_DIR}/${CONFIGS_DIR} ${WORK_DIR}
   if [[ $queue != false ]]
   then
        eval "sed \"s!//queue = false!queue = '${queue}'!\" ${RACINE_PIPELINES_DIR}/${CONFIGS_DIR}/cluster.config > ${WORK_DIR}/${CONFIGS_DIR}/cluster.config"
   fi
   cp -r ${RACINE_PIPELINES_DIR}/${ASSETS_DIR} ${WORK_DIR}
   cp ${RACINE_PIPELINES_DIR}/nextflow.config ${WORK_DIR}
   # creation scripts de lancement de nextflow
   if [[ $queue != false ]]
   then
        COMM="${SUBMIT}"
   fi
   cat <<MAPPING > ${SCRIPT_STEP_MAPPING}
# SAREK MAPPING
cd  ${WORK_DIR}
# lancement du pipeline sur le cluster 

nextflow run main.nf -resume  --genome_base ${genomebase} --tag latest  -c nextflow.config --genome ${genome} --sequencing_center "SEQC-II Consortium" --step mapping ${sampledata} -profile ${profile} -with-trace Reports/pipeline_info/Sarek_mapping_trace.csv -with-report Reports/pipeline_info/Sarek_mapping_report.html -with-timeline Reports/pipeline_info/Sarek_mapping_timeline.html -with-dag Reports/pipeline_info/Sarek_mapping_DAG.pdf && sleep 30 && ${COMM} ${SCRIPT_STEP_GERMLINE}

MAPPING

 cat <<GERMLINE > ${SCRIPT_STEP_GERMLINE}
# SAREK GERMLINE
cd  ${WORK_DIR}
# lancement du pipeline sur le cluster 

nextflow run germlineVC.nf -resume --genome_base ${genomebase} --tag latest -c nextflow.config --genome ${genome} --sequencing_center "SEQC-II Consortium" --step germline --sample Preprocessing/Recalibrated/recalibrated.tsv --targetBED ${samplepath} --tools FreeBayes,HaplotypeCaller,Manta,Mutect2,Strelka -profile ${profile} -with-trace Reports/pipeline_info/Sarek_germline_trace.csv -with-report Reports/pipeline_info/Sarek_germline_report.html  -with-timeline Reports/pipeline_info/Sarek_germline_timeline.html -with-dag Reports/pipeline_info/Sarek_germline_DAG.pdf &&  sleep 30 && ${COMM} ${SCRIPT_STEP_SOMATIC}

GERMLINE

 cat <<SOMATIC > ${SCRIPT_STEP_SOMATIC}
# SAREK SOMATIC
cd  ${WORK_DIR}
# lancement du pipeline sur le cluster 

nextflow run  somaticVC.nf -resume --genome_base ${genomebase} --tag latest -c nextflow.config --genome ${genome} --sequencing_center "SEQC-II Consortium" --step somatic --sample Preprocessing/Recalibrated/recalibrated.tsv --targetBED ${samplepath} --tools FreeBayes,HaplotypeCaller,Manta,Mutect2,Strelka -profile ${profile} -with-trace Reports/pipeline_info/Sarek_somatic_trace.csv -with-report Reports/pipeline_info/Sarek_somatic_report.html  -with-timeline Reports/pipeline_info/Sarek_somatic_timeline.html -with-dag Reports/pipeline_info/Sarek_somatic_DAG.pdf &&  sleep 30 && ${COMM} ${SCRIPT_STEP_ANNOTATE}

SOMATIC

 cat <<ANNOTATE > ${SCRIPT_STEP_ANNOTATE}
# SAREK ANNOTATE 
cd  ${WORK_DIR}
# lancement du pipeline sur le cluster 

nextflow run annotate.nf -resume --genome_base ${genomebase} --tag latest -c nextflow.config --genome ${genome} --sequencing_center "SEQC-II Consortium" --step annotate --annotateTools haplotypecaller,strelka,mutect2 --tools snpEff -profile ${profile} -with-trace Reports/pipeline_info/Sarek_annotate_trace.csv -with-report Reports/pipeline_info/Sarek_annotate_report.html  -with-timeline Reports/pipeline_info/Sarek_annotate_timeline.html -with-dag Reports/pipeline_info/Sarek_annotate_DAG.pdf &&  sleep 30 && ${COMM} ${SCRIPT_STEP_MULTIQC} 

ANNOTATE

 cat <<MULTIQC > ${SCRIPT_STEP_MULTIQC}
# SAREK MULTIQC
cd  ${WORK_DIR}
# lancement du pipeline sur le cluster 

nextflow run runMultiQC.nf -resume --genome_base ${genomebase} --tag latest -c nextflow.config --genome ${genome} --tools fastqc,picard,samtools,qualimap,bcftools,vcftools,snpeff --sequencing_center "SEQC-II Consortium" --step multiqc -profile ${profile} -with-trace Reports/pipeline_info/Sarek_multiqc_trace.csv -with-report Reports/pipeline_info/Sarek_multiqc_report.html  -with-timeline Reports/pipeline_info/Sarek_multiqc_timeline.html -with-dag Reports/pipeline_info/Sarek_multiqc_DAG.pdf && ${COMM} ${SCRIPT_STEP_COPIE_RESULTS}

MULTIQC

# creation script de copie des metriques  
   cat <<COPIE > ${SCRIPT_STEP_COPIE_RESULTS}
cd  ${WORK_DIR}/Reports
cp -r MultiQC ${results_dir}/${name_end} 
cp -r pipeline_info ${results_dir}/${name_end} 
cp -r FastQC ${results_dir}/${name_end} 
COPIE

}

function pipeline_install_nxf_script() {
   echo " 3) install_nxf_script"
   [[ -e ${RACINE_PIPELINES_DIR}/${NXF_NAME} ]] || die "${RACINE_PIPELINES_DIR}/${NXF_NAME} n'est pas accessible"
   cp  ${RACINE_PIPELINES_DIR}/*.nf ${WORK_DIR}
   cp -r ${RACINE_PIPELINES_DIR}/${BIN_DIR} ${WORK_DIR}
   cp -r ${RACINE_PIPELINES_DIR}/${LIB_DIR} ${WORK_DIR}
   cp -r ${RACINE_PIPELINES_DIR}/${DOCS_DIR} ${WORK_DIR}
   cp -r ${RACINE_PIPELINES_DIR}/${SCRIPTS_DIR} ${WORK_DIR}
   ln -s ${RACINE_PIPELINES_DIR}/${REF_DIR} ${WORK_DIR}
   ln -s ${RACINE_PIPELINES_DIR}/${REPEATES_DIR} ${WORK_DIR}
   ln -s ${RACINE_PIPELINES_DIR}/${DATA_DIR} ${WORK_DIR}
   cp  ${RACINE_PIPELINES_DIR}/environment.yml ${WORK_DIR}
}

function nxf_lunch_pipeline() {
   echo " 4) nxf_pipeline_nfcore-sarek"
   
   echo_green "WORK_DIR = ${WORK_DIR}"
   echo_green "LOCAL_SCRIPTS_PATH = ${LOCAL_SCRIPTS_PATH}"
   echo_green "RACINE_PIPELINES_DIR = ${RACINE_PIPELINES_DIR}"
   echo_green "CONFIG_NXF_PATH = ${WORK_DIR}/${CONFIGS_DIR}"
   echo_green "LOCAL_SCRIPTS (lancemment et enchainement automatique sur le cluster): "
   echo_green "SCRIPT_STEP_MAPPING = ${SCRIPT_STEP_MAPPING}"
   echo_green "SI MAPPING OK => SCRIPT_STEP_GERMLINE = ${SCRIPT_STEP_GERMLINE}"
   echo_green "SI GERMLINE OK => SCRIPT_STEP_SOMATIC = ${SCRIPT_STEP_SOMATIC}"
   echo_green "SI SOMATIC OK => SCRIPT_STEP_MUTIQC = ${SCRIPT_STEP_ANNOTATE}"
   echo_green "SI ANNOTATE OK => SCRIPT_STEP_MUTIQC = ${SCRIPT_STEP_MULTIQC}"
   d=$(date)
   cd ${WORK_DIR}
   ${COMM} ${SCRIPT_STEP_MAPPING}

   #if [[ $queue != false ]]
   #then
   #	${SUBMIT} ${SCRIPT_STEP_MAPPING}
   #else
   #	bash ${SCRIPT_STEP_MAPPING}
   #fi
}


function nxf_date() {
    local ts=$(date +%s%3N); [[ $ts == *3N ]] && date +%s000 || echo $ts
}

function get_date() {
    local ts=$(date +%d-%m-%y) && echo $ts
}
name_end=$(nxf_date)
date_run=$(get_date)

##
echo_yellow "#### $0 for projet ${project_name} env:${env} and run:${run} by ${LOGNAME}"
#
# creation des variables d'environnement nécessaires au lancement du pipeline 
pipeline_create_env

#
# creation du repertoire d'exec : doit contenir l'ensemble des composants necessaire au lancement du pipeline
# creation repertoires resultats, configs, script, creation script bash de lancement 
# de nextflow pointant sur les parametres.
pipeline_create_nxf_work_dir

# ajout lien vers le script nxf du pipeline dans WORK_DIR
pipeline_install_nxf_script

#
# lancement de nexflow pour le pipeline pour la première étape qui enchaine la suivante sur succes
nxf_lunch_pipeline


