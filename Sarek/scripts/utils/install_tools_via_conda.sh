#!/bin/bash
## install_tools_via_conda.sh 
## Tools of Sarek pipeline
##
## Copyright (c) 2019-2020 Institut Curie
## Author(s): Philippe La Rosa
## Contact: philippe.larosa@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details
##
##
## usage : bash  install_tools_via_conda.sh dir_conda_modules env_name
## exemple : bash  install_tools_via_conda.sh /bioinfo/pipelines/sandbox/dev/Sarek/modules sarek-2.3
##
dir_conda_modules=$1
env_name=$2
file_log=install_tools_via_conda.log
script_install_conda=Miniconda3-4.5.12-Linux-x86_64.sh
images_dir=$(grep "containerPath =" ../../conf/singularity-path.config | sed "s/containerPath = '//" | sed "s/'//" | awk ' { print $1 } ')
### telechargement du script d'installation de miniconda
wget https://repo.continuum.io/miniconda/${script_install_conda}
### installation miniconda
bash ${script_install_conda} -u -b -p ${dir_conda_modules}/conda 2>&1 | tee ${dir_conda_modules}/${file_log}

export PATH=${dir_conda_modules}/conda/bin:$PATH

### (Re)Installationdes tools via conda dans dir_conda_modules
rm -rf ${dir_conda_modules}/conda/envs/${env_name}
conda env create -p ${dir_conda_modules}/conda/envs/${env_name} -f environment.yml 2>&1 | tee -a ${dir_conda_modules}/${file_log}
# pour finir l'installation du tool snpeff: copie data genome de l'image  singularity sarek-snpeff-2.3 vers l'installation conda de snpeff
echo "singularity exec -B  ${dir_conda_modules}/conda/envs/${env_name}/share/snpeff-4.3.1t-2  ${images_dir}/snpeffgrch37-latest.simg cp -r /opt/conda/envs/sarek-snpeff-2.3/share/snpeff-4.3.1t-2/data ${dir_conda_modules}/conda/envs/${env_name}/share/snpeff-4.3.1t-2" >> ${dir_conda_modules}/${file_log}
singularity exec -B  ${dir_conda_modules}/conda/envs/${env_name}/share/snpeff-4.3.1t-2  ${images_dir}/snpeffgrch37-latest.simg cp -r /opt/conda/envs/sarek-snpeff-2.3/share/snpeff-4.3.1t-2/data ${dir_conda_modules}/conda/envs/${env_name}/share/snpeff-4.3.1t-2 2>&1 | tee -a ${dir_conda_modules}/${file_log}




