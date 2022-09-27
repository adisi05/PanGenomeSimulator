# PanGenomeSimulator
Simulator that produces data to benchmark pan genome assembly models.

# genMutModel.py

Takes references genome and TSV file to generate mutation models:

```
python gen_mut_model.py               \
        -r hg19.fa                  \
        -m inputVariants.tsv        \
        -o /home/me/models.p
```

Trinucleotides are identified in the reference genome and the variant file. Frequencies of each trinucleotide transition are calculated and output as a pickle (.p) file.


# plotMutModel.py

Performs plotting and comparison of mutation models generated from genMutModel.py.

```
python plotMutModel.py                                        \
        -i model1.p [model2.p] [model3.p]...                  \
        -l legend_label1 [legend_label2] [legend_label3]...   \
        -o path/to/pdf_plot_prefix
```

# vcf_compare_OLD.py

Tool for comparing VCF files.

```
python vcf_compare_OLD.py
        --version          show program's version number and exit      \
        -h, --help         show this help message and exit             \
        -r <ref.fa>        * Reference Fasta                           \
        -g <golden.vcf>    * Golden VCF                                \
        -w <workflow.vcf>  * Workflow VCF                              \
        -o <prefix>        * Output Prefix                             \
        -m <track.bed>     Mappability Track                           \
        -M <int>           Maptrack Min Len                            \
        -t <regions.bed>   Targetted Regions                           \
        -T <int>           Min Region Len                              \
        -c <int>           Coverage Filter Threshold [15]              \
        -a <float>         Allele Freq Filter Threshold [0.3]          \
        --vcf-out          Output Match/FN/FP variants [False]         \
        --no-plot          No plotting [False]                         \
        --incl-homs        Include homozygous ref calls [False]        \
        --incl-fail        Include calls that failed filters [False]   \
        --fast             No equivalent variant detection [False]     
```
Mappability track examples: https://github.com/zstephens/neat-repeat/tree/master/example_mappabilityTracks