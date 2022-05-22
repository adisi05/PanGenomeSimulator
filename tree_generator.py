# https://www.biostars.org/p/9472422/
# https://github.com/sbslee/fuc
import os
import subprocess


def generate_tree(vcf_directory):
    # print(os.environ)
    os.system("/groups/itay_mayrose/adisivan/PanGenomeSimulator/remote_debugging_init.sh")
    os.system("conda info")
    # subprocess.run('source activate PanGenomeSimulator', shell=True)
    # subprocess.run('conda info', shell=True)
    # subprocess.run('source deactivate', shell=True)

    # # os.system("source /groups/itay_mayrose/adisivan/miniconda2/etc/profile.d/conda.sh")
    # os.getenv("source /groups/itay_mayrose/adisivan/miniconda2/etc/profile.d/conda.sh; conda activate PanGenomeSimulator")
    # os.system("conda info")
    # # subprocess.run(["source","/groups/itay_mayrose/adisivan/miniconda2/etc/profile.d/conda.sh"], shell=True)
    # # subprocess.run(["conda","activate","PanGenomeSimulator"], shell=True)
    # # subprocess.run(["conda","info"], shell=True)

    # # subprocess.run('source activate PanGenomeSimulator', shell=True)
    #
    # # print(os.environ.get('CONDA_PREFIX'))
    # os.system("conda info")
    # # subprocess.run('source activate PanGenomeSimulator', shell=True)
    # cmd = "conda run - n PanGenomeSimulator test.py"
    # # subprocess.run(["conda","run","-n","PanGenomeSimulator","test.py"], shell=True)
    # subprocess.run(["conda","activate","PanGenomeSimulator"], shell=True)
    #
    # subprocess.run(["conda","info"], shell=True)
    # # os.system(cmd)
    # # os.system("conda info")
    # # subprocess.call(["ls", "-lha"])
    # # subprocess.call(["conda","info"])
    # # subprocess.call(["conda","activate","PanGenomeSimulator"])
    # # subprocess.call(["conda","info"])


    # files_to_merge = []
    # for filename in os.listdir(vcf_directory):
    #     f = os.path.join(vcf_directory, filename)
    #     # checking if it is a file
    #     if os.path.isfile(f) and f.endswith('.vcf') or f.endswith('.vcf.gz'):
    #         files_to_merge.append(f)
    #         if not os.path.isfile(f+'.csi'):
    #             print(f)
    #             os.system(f"bgzip {f}")
    #             os.system(f"tabix -f -p vcf {f}.gz")


    # input_string = ' '.join([str(name) for name in filenmaes])
    # merged_vcf = "merged.vcf"
    # tree_file = "tree.newick"
    # os.system(f"bcftools merge --merge none {input_string} > {merged_vcf}")
    # os.system(f"vk phylo tree nj {merged_vcf} > {tree_file}")
    # subprocess.run('source deactivate', shell=True)


if __name__ == "__main__":
    generate_tree("/groups/itay_mayrose/adisivan/PanGenomeSimulator/test_arabidopsis/vcf_inputs")