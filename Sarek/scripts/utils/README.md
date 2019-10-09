# To internal install the Sarek pipeline tools with conda 

### Installation  summary

1. Run the script to build the tools
```
    Example: bash install_tools_via_conda.sh /your_path/tools env_name
```
```
    Example: bash  install_tools_via_conda.sh /bioinfo/pipelines/sandbox/dev/Sarek/modules sarek-2.3
```

2. Edit your new internal profile conf/tools-path.config: change the path of your local installation

```bash
    singularity {
     enabled = false 
    }

    process {
      beforeScript = 'export PATH=/your_path/tools/sarek-2.3/bin:$PATH'
    }
```

### Quick run
Run the pipeline locally, using the your internal global environment and tools build by conda

```
nextflow run main.nf -profile toolsPath,cluster,test

```


