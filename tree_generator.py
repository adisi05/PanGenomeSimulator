# https://www.biostars.org/p/9472422/
# https://github.com/sbslee/fuc
import os


def generate_tree(filenmaes):
    input_string = ' '.join([str(name) for name in filenmaes])
    merged_vcf = "merged.vcf"
    tree_file = "tree.newick"
    os.system(f"bcftools merge --merge none {input_string} > {merged_vcf}")
    os.system(f"vk phylo tree nj {merged_vcf} > {tree_file}")
