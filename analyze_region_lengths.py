import argparse

from utilities.common_data_structues import Region
from utilities.io.genome_annotations_reader import read_annotations_csv


def region_lengths_analysis(annotations_df):
    print("Starting Region Lengths Analysis")
    current_chrom = ''
    all_results = {}
    for i, annotation in annotations_df.iterrows():
        if current_chrom != annotation['chrom']:
            current_chrom = annotation['chrom']
            chrom_results = {
                Region.CDS.value: 0,
                Region.NON_CODING_GENE.value: 0,
                Region.INTERGENIC.value: 0,
                Region.ALL.value: 0
            }
            all_results[current_chrom] = chrom_results

        annotation_length = annotation['end'] - annotation['start']
        annotation_name = annotation['region']
        chrom_results[annotation_name] = chrom_results[annotation_name] + annotation_length
        chrom_results[Region.ALL.value] = chrom_results[Region.ALL.value] + annotation_length

    print(f"Completed Region Lengths Analysis")
    for chrom, chrom_results in all_results.items():
        print("Results from chrom:", chrom)
        overall_length = chrom_results[Region.ALL.value]
        for region_name, region_length in chrom_results.items():
            if region_name == Region.ALL.value:
                continue
            print("Region:", region_name, "\tLength:", region_length, "\tPercentage:",
                  round(region_length/overall_length, 2))


def main(raw_args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', type=str, required=False, metavar='/path/to/output.csv',
                        help="Path to output csv file")
    args = parser.parse_args(raw_args)

    if args.csv:
        annotations_df = read_annotations_csv(args.csv)
        region_lengths_analysis(annotations_df)

if __name__ == "__main__":
    main()
