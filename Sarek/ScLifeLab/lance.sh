# SAREK : steps
# Preprocessing – main.nf (based on GATK best practices)
# Germline variant calling – germlineVC.nf
# Somatic variant calling – somaticVC.nf (optional)
# Annotation – annotate.nf (optional)
# Reporting (multiqc)

# scripts de creation du repertoire d'exec : qui contient l'ensemble des composants necessaires au lancement
# du pipeline, creation des variables d'environnements, des repertoires resultats, configs, script, creation 
# des scripts bash de lancement des etapes du pipeline.

# Quatre modes de lancement : 

# 1) Manuel
# les scripts de lancment des etapes sont générés dans l'espace d'exec, ils sont à lancer et enchainer manuellement sur le cluster :
# singularityPath
#export RACINE_PIPELINES_DIR=$(pwd);bash scripts/run_install.sh -p SAREK-2.3 -r TEST -g smallGRCh37 -s "--sample Sarek-data/testdata/tsv/tiny-manta.tsv" -u singularityPath,cluster,test -i Sarek-data/testdata/target.bed -o /data/tmp/NGS_RUN_TEMP
# toolsPath 
#export RACINE_PIPELINES_DIR=$(pwd);bash scripts/run_install.sh -p SAREK-2.3 -r TEST -g smallGRCh37 -s "--sample Sarek-data/testdata/tsv/tiny-manta.tsv" -u toolsPath,cluster,test -i Sarek-data/testdata/target.bed -o /data/tmp/NGS_RUN_TEMP
#  
# 2) Automatique
# singularityPath
# les scripts de lancment des etapes sont générés dans l'espace d'exec, ils sont à lancer et enchainer automatiquement sur le cluster :
export RACINE_PIPELINES_DIR=$(pwd);bash scripts/run_sarek.sh -p SAREK-2.3 -r TEST -g GRCh38 -b "/data/kdi_prod/.kdi/project_workspace_0/1430/acl/01.00/public_data/annotations/GRCh38" -s "--sample Sarek-data/HPC-bench/tsv/HPC-bench-test.tsv" -u singularityPath,cluster -i /data/kdi_prod/.kdi/project_workspace_0/1430/acl/01.00/public_data/annotations/illumina/TargetRegions.bed -m /data/kdi_prod/project_result/1430/01.02/results/test/sarek -q false -o /data/tmp/NGS_RUN_TEMP
## data light
#export RACINE_PIPELINES_DIR=$(pwd);bash scripts/run_sarek.sh -p SAREK-2.3 -r TEST -g smallGRCh37 -b "References/smallGRCh37" -s "--sample Sarek-data/HPC-bench/tsv/HPC-bench-test.tsv" -u singularityPath,cluster -i Sarek-data/testdata/target.bed -m /data/kdi_prod/project_result/1430/01.02/results/test/sarek -q false -o /data/tmp/NGS_RUN_TEMP 

# toolsPath 
#export RACINE_PIPELINES_DIR=$(pwd);bash scripts/run_sarek.sh -p SAREK-2.3 -r TEST_TOOLS -g GRCh38 -b "/data/kdi_prod/.kdi/project_workspace_0/1430/acl/01.00/public_data/annotations/GRCh38" -s "--sample Sarek-data/HPC-bench/tsv/HPC-bench-test.tsv" -u toolsPath,cluster -i /data/kdi_prod/.kdi/project_workspace_0/1430/acl/01.00/public_data/annotations/illumina/TargetRegions.bed -m /data/kdi_prod/project_result/1430/01.02/results/test/sarek -q mpi -o /data/tmp/NGS_RUN_TEMP
