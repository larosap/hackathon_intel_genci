#!/bin/bash
## get_singularity_images.sh 
## Tools of Sarek pipeline
##
## Copyright (c) 2019-2020 Institut Curie
## Author(s): Philippe La Rosa
## Contact: philippe.larosa@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details
##
##
## usage : bash get_singularity_images.sh containers_singularity_path 
## exemple : bash get_singularity_images.sh ./containers_singularity
##
## Process requirement :  available 2 CPUs
##
containers_singularity_path=$1
nextflow run ../../build.nf -profile singularity --verbose --singularity --repository maxulysse --tag latest --containerPath ${containers_singularity_path}/ --containers r-base,runallelecount,sarek,snpeffgrch37,snpeffgrch38,vepgrch37,vepgrch38



