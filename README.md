# isoseq_pipeline
Isoform detection with paralog-specific partitioning

## Quick start guide
To run the pipeline, modify the `config.tofu.yaml` file to include your smrt cells grouped by replicate. Next, run the partitioning step to partition your flnc reads based on reference mapping:
```
snakemake -s tofu.partition.snake create_partitions
```

Then, run ICE/TOFU on each partition separately:
```
snakemake -s tofu.partition.snake
```
