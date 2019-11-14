# SAREK : steps
# Preprocessing – main.nf (based on GATK best practices)
# Germline variant calling – germlineVC.nf
# Somatic variant calling – somaticVC.nf (optional)
# Annotation – annotate.nf (optional)
# Reporting (multiqc)


# exemple de lancement (a modifier avec vos path) du script qui effectue les etapes : de creation du repertoire d'exec
# qui contient l'ensemble des composants necessaires au lancement du pipeline,
# creation des variables d'environnements, des repertoires resultats, configs, script, 
# creation des scripts bash de lancement et d'enchainements des etapes du pipeline.
#
# Update with your nextflow installation
# export PATH=/bioinfo/local/build/Centos/nextflow/nextflow-19.04.0.5069:$PATH
#  modes de lancement : 
# Automatique
# singularity sur le cluster
#export SUBMIT="srun";export EXECUTOR="slurm";export RACINE_PIPELINES_DIR=$(pwd);bash scripts/run_sarek.sh -p SAREK-2.3 -r TEST -g GRCh38 -b "/data/kdi_prod/.kdi/project_workspace_0/1430/acl/01.00/public_data/annotations/GRCh38" -s "--sample Sarek-data/HPC-bench/tsv/HPC-bench-test.tsv" -u singularityPath,cluster -i Sarek-data/HPC-bench/TargetRegions.bed -m $(pwd)/results/test/sarek -q dev -o /data/tmp/NGS_RUN_TEMP

## data et genome light
# singularity sur le cluster
export SUBMIT="qsub";export EXECUTOR="pbs";export RACINE_PIPELINES_DIR=$(pwd);bash scripts/run_sarek.sh -p SAREK-2.3 -r TEST -g smallGRCh37 -b "References/smallGRCh37" -s "--sample Sarek-data/HPC-bench/tsv/HPC-bench-test.tsv" -u singularityPath,cluster -i Sarek-data/testdata/target.bed -m $(pwd)/results/test/sarek -q dev -o /data/tmp/NGS_RUN_TEMP 

# singularity en local 
#export RACINE_PIPELINES_DIR=$(pwd);bash scripts/run_sarek.sh -p SAREK-2.3 -r TEST -g smallGRCh37 -b "References/smallGRCh37" -s "--sample Sarek-data/HPC-bench/tsv/HPC-bench-test.tsv" -u singularityPath,test -i Sarek-data/testdata/target.bed -m $(pwd)/results/test/sarek -q false -o /data/tmp/NGS_RUN_TEMP 2>&1 > Sarek_HPC-bench_singularity_$$.log  

# toolsPath cluster (tools installes via conda)  
#export SUBMIT="qsub";export EXECUTOR="pbs";export RACINE_PIPELINES_DIR=$(pwd);bash scripts/run_sarek.sh -p SAREK-2.3 -r TEST_TOOLS -g GRCh38 -b "/data/kdi_prod/.kdi/project_workspace_0/1430/acl/01.00/public_data/annotations/GRCh38" -s "--sample Sarek-data/HPC-bench/tsv/HPC-bench-test.tsv" -u toolsPath,cluster -i Sarek-data/testdata/target.bed -m $(pwd)/results/test/sarek -q mpi -o /data/tmp/NGS_RUN_TEMP

# toolsPath en local (tools installes via conda)  
#export RACINE_PIPELINES_DIR=$(pwd);bash scripts/run_sarek.sh -p SAREK-2.3 -r TEST -g smallGRCh37 -b "References/smallGRCh37" -s "--sample Sarek-data/HPC-bench/tsv/HPC-bench-test.tsv" -u toolsPath,test -i Sarek-data/testdata/target.bed -m $(pwd)/results/test/sarek -q false -o /data/tmp/NGS_RUN_TEMP 2>&1 > Sarek_HPC-bench_$$.log 
