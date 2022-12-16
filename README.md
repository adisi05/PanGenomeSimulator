# PanGenomeSimulator
Simulator that produces data to benchmark pan genome assembly models.

## process_annotations.py
Process a gff-like annotation file and make it ready to use for the next steps (model building and genome simulation).
```
python process_annotations.py                                   \
        # Main flow:                                            \
        --gff /path/to/annotations.gff                          \
        --csv /path/to/output.csv                               \
                                                                \
        # Sanity test:                                          \
        --test-sanity /path/to/post/processed/annotations.csv   \
                                                                \
        # Overlapping genes test:                               \
        --test-overlap  /path/to/annotations.gff                \
                        /path/to/overlap/output.csv
```
Two test workflows:
* Sanity test - runs also as part of the main workflow, this flag is for a stand-alone run.
* Overlap test - outputs overlapping genes in case they exist.

## generate_mutation_model.py
Generate the mutation model required for the pan-genome simulation.
```
python generate_mutation_model.py   \
        -r /path/to/reference.fasta \
        -v /path/to/variants.vcf    \
        -o /path/to/output/model.p  \
                                    \
        # Optional:                 \
        -a /path/to/annotations.csv
```

## pangenome_simulator.py
Simulate a phylogenetic tree of genomes, based on a given fasta reference and a mutational model.

If a tree file is not given, simulate only one descendent genome based on the model.

It is recommended to use the same genomic annotation file (csv) that was used during the model creation run. 
```
python pangenome_simulator.py                                               \
        -r /path/to/reference.fasta                                         \
        -m /path/to/mutation/model.p                                        \
        -o /output/prefix                                                  \
                                                                            \
        # Optional:                                                         \
        -a /path/to/annotations.csv                                         \
        -t /path/to/tree/file.newick                                        \
        -M <float> # mutation rate multiplcation scalar                     \
        -l <str> # name of accesion to use as root                          \
        --max-threads <int> # maximal threads number, default is 1          \
        --vcf # produce output vcf                                          \
        --no-fastq # do not produce reads output                            \
        --rng <int> # rng seed value;                                       \
        # identical RNG value should produce identical runs of the program, \
        # so things like read locations, variant positions, error positions,\
        # etc, should all be the same.                                      \
        -w <int> # simulation window size; default is 10000                 \
        -d # print debug-level log messages                                 \
                                                                            \
        # Sequencing parameters:                                            \
        -R <int> # read length, default is 150                              \
        -c <float> # coverage, default is 30.0                              \
        --pe <int> <int> # paired-end fragment length mean and std          \
        --art-path <str> # path to ART binary, default is:
        # ART/art_bin_MountRainier/art_illumina

```

## plot_mutation_model.py
Performs plotting and comparison of mutation models generated from generate_mutation_model.py.
```
python plot_mutation_model.py                                       \
        -i model1.p [model2.p] [model3.p]...                        \
        -l legend_label1 [legend_label2] [legend_label3]...         \
        -o path/to/pdf_plot_prefix                                  \
                                                                    \
        # The next flag is optional - genomic region(s) to focus on.\
        # If you use it, choose one or more values from the list:   \
        -g [CDS] [non_coding_gene] [intergenic] [all]
```
