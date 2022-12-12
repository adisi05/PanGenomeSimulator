# PanGenomeSimulator
Simulator that produces data to benchmark pan genome assembly models.

## process_annotations.py
Process a gff-like annotation file and make it ready to use for the next steps (model building and genome simulation).
```
python process_annotations.py                                   \
        # Main flow:                                            \
        --gff /path/to/annotations.gff                          \
        --csv /path/to/output.csv                               \
        # Sanity test:                                          \
        --test-sanity /path/to/post/processed/annotations.csv   \
        # Overlapping genes test:                               \
        --test-overlap  /path/to/annotations.gff                \
                        /path/to/overlap/output.csv
```
Two test workflows:
* Sanity test - runs also as part of the main workflow, this flag is for a stand-alone run.
* Overlap test - outputs overlapping genes in case they exist.

## generate_mutation_model.py
gfdjhfhgjf
```
python generate_mutation_model.py   \
        -r /path/to/reference.fasta \
        -v /path/to/variants.vcf    \
        -o /path/to/output/model.p  \
        # Optional:                 \
        -a /path/to/annotations.csv
```
cdfhgvhjgkjhbkljhn;klj;ki;ljknlkjhbgfcghbkjhlmll

## pangenome_simulator.py
```
python pangenome_simulator.py           \
        -r /path/to/reference.fasta     \
        -m /path/to/mutation/model.p    \
        -o
        # Optional:
        -a /path/to/annotations.csv
        -M
        -t
        -l
        --max-threads
        --vcf
        --no-fastq
        --rng
        -w
        -d
        -R
        -c
        --pe
```

## plot_mutation_model.py
Performs plotting and comparison of mutation models generated from generate_mutation_model.py.
```
python plot_mutation_model.py                                       \
        -i model1.p [model2.p] [model3.p]...                        \
        -l legend_label1 [legend_label2] [legend_label3]...         \
        -o path/to/pdf_plot_prefix                                  \
        # The next flag is optional - genomic region(s) to focus on.\
        # If you use it, choose one or more values from the list:   \
        -g [CDS] [non_coding_gene] [intergenic] [all]
```

## plot_tree.py