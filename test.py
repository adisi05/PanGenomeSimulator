import dendropy
from dendropy import Taxon
import ete3

from neat.source.SequenceContainer import SequenceContainer, parse_input_mutation_model


def test_dendropy_params(newick_file,fasta_file):
    print("newick_file = ",newick_file)
    tree = dendropy.Tree.get(path=newick_file, schema="newick", suppress_internal_node_taxa=False)
    print("by length")
    tree.print_plot(show_internal_node_labels=True, plot_metric="length")
    # print("by level")
    # tree.print_plot(show_internal_node_labels=True, plot_metric="level")
    # print("by depth")
    # tree.print_plot(show_internal_node_labels=True, plot_metric="depth")
    # print("by age")
    # tree.print_plot(show_internal_node_labels=True, plot_metric="age")

    for edge in tree.preorder_edge_iter():
        print(edge.length)

    # str = "ABCDEFGHIJKLMNOP"
    # i = 0
    for node in tree.preorder_node_iter():
        print(node.taxon.label)
        # temp = node.taxon.label
        # node.taxon = str[i]
        # node.taxon.label = temp
        # print(node.taxon)
        # i+=1

    tree.print_plot(show_internal_node_labels=True, plot_metric="length")

    # d = dendropy.DataSet()
    # taxon_namespace = dendropy.TaxonNamespace()
    # d.attach_taxon_namespace(taxon_namespace)
    # d.read(data=newick_file, schema="newick")
    # # d.read(data=fasta_file, schema="fasta", data_type="dna")
    # print(d.taxon_sets[0].description(2))
    # print(d.taxon_sets)

    # tree_list = dendropy.TreeList()
    # tree_list.read(data=newick_file, schema="newick")
    # print(tree_list.taxon_namespace)

    d = dendropy.DataSet()
    tns = d.new_taxon_namespace()
    d.attach_taxon_namespace(tns)
    d.read_from_path(src=newick_file, schema="newick")
    d.read_from_path(src=fasta_file, schema="fasta", data_type="dna")
    print(d.taxon_namespaces[0].description(0))
    print(d.taxon_namespaces[0].description(1))
    print(d.taxon_namespaces[0].description(2))
    print(d.taxon_namespaces[0].description(3))
    print(d.as_string('nexus'))

def test_ete3_params(newick_file,fasta_file):
    t = ete3.Tree(newick_file)
    (t & "D").add_features(fasta="AAAA")
    (t & "F").add_features(fasta="ACGT")
    print(t)
    print(t.get_ascii(attributes=["name", "dist", "fasta"]))
    (t & "D").add_features(SequenceContainer = SequenceContainer(0, "AAAA", 1, 2, 2, [parse_input_mutation_model(None, 1)] * 1))
    (t & "F").add_features(SequenceContainer = SequenceContainer(0, "ACGT", 1, 2, 2, [parse_input_mutation_model(None, 1)] * 1))
    print((t & "D").fasta)
    print((t & "F").fasta)
    print((t & "D").SequenceContainer.sequences)
    print((t & "F").SequenceContainer.sequences)
    print((t & "D").SequenceContainer.sequences[0])
    print((t & "F").SequenceContainer.sequences[0])

    for node in t.traverse("preorder"):
        # Do some analysis on node
        print(node.name)




if __name__ == "__main__":

    arr = [3]
    if len(arr):
        print(True)
    else:
        print(False)

    num = None
    if (num == None):
        print("hi")


    # test_dendropy_params("tree1.newick")

    # test_dendropy_params("tree2.newick","tree2.fasta")

    # test_ete3_params("tree4.newick","tree2.fasta")

    # test_dendropy_params("tree3.newick")