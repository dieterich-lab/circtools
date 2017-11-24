# skript to add permutation information to all miRNA-gene pairs.


def read_fasta_files(folder):
    files = os.listdir('%s' % (folder))
    gene_map = {}
    for f in files:
        gene = f.replace('.fa', '')
        I = open('%s/%s' % (folder, f))
        for line in I:
            if line.startswith('>'):
                gene_map[line.replace('>', '').replace('\n', '')] = gene
        I.close()
    return (gene_map)


def read_real_seeds(infile, gene_map):
    I = open(infile)
    real_seeds = {}
    for line in I:
        L = line.replace('\n', '').split('\t')
        if not (L[1], L[0]) in real_seeds:
            if 'ENS' in L[1]:
                sequence_type = 'hostgene'
            else:
                sequence_type = 'circRNA'
            real_seeds[(L[1], L[0])] = {'gene_name': gene_map[L[1]], 'sequence_type': sequence_type, 'real_seeds': 1,
                                        'permutations': [0] * 1000}
        else:
            real_seeds[(L[1], L[0])]['real_seeds'] += 1
    return (real_seeds)


def read_permutations(real_seeds, infile, permutation):
    I = open(infile)
    for line in I:
        L = line.replace('\n', '').split('\t')
        if (L[1], L[0]) in real_seeds:
            real_seeds[(L[1], L[0])]['permutations'][permutation] += 1
    I.close()
    return (real_seeds)


def write_file(real_seeds, outputfile, sample_name):
    O = open(outputfile, 'w')
    O.write('sample\tgene\ttype\ttranscript\tmiRNA\tobserved\t%s\n' % (
    '\t'.join(['permutation_%s' % (i) for i in range(1, 1000)])))
    for pair in real_seeds:
        O.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        sample_name, real_seeds[pair]['gene_name'], real_seeds[pair]['sequence_type'], pair[0], pair[1],
        real_seeds[pair]['real_seeds'], '\t'.join(['%s' % (i) for i in real_seeds[pair]['permutations']])))
    O.close()
    return


if __name__ == '__main__':

    # required packages
    import argparse
    import os

    parser = argparse.ArgumentParser(description='')

    # input
    parser.add_argument('real_seeds', metavar='real_seeds_file', help='path to circRNAwise real seed hit list')
    parser.add_argument('permutations', metavar='permutations_folder', help='path to circRNA wise permutation folder')
    parser.add_argument('fastafolder', metavar='fasta_folder', help='path to circRNA wise permutation folder')
    parser.add_argument('outfile', metavar='outfile', help='path and filename to write the output to')

    args = parser.parse_args()

    # parse arguments
    permutations_folder = args.permutations
    real_seeds_file = args.real_seeds
    outfile = args.outfile
    fastafolder = args.fastafolder

    gene_map = read_fasta_files(fastafolder)
    print('done creating gene_map')

    SEEDS = read_real_seeds(real_seeds_file, gene_map)

    print('done reading real seeds')
    for i in range(1, 1000):
        SEEDS = read_permutations(SEEDS, '%s/permutation_%s.fimo.txt' % (permutations_folder, i), i)
    print('done reading permutations')

    sample_name = fastafolder.split('/')[-2]
    write_file(SEEDS, outfile, sample_name)

# python
# add_permutation_information.py / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / real_seeds / old_cerebellum.fimo.txt / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / on_predicted / permutations / old_cerebellum / / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / fasta_files / old_cerebellum / / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / pvalues / old_cerebellum.miRNA_seeds.permutations.txt
# python
# add_permutation_information.py / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / real_seeds / old_hippocampus.fimo.txt / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / on_predicted / permutations / old_hippocampus / / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / fasta_files / old_hippocampus / / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / pvalues / old_hippocampus.miRNA_seeds.permutations.txt
# python
# add_permutation_information.py / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / real_seeds / old_liver.fimo.txt / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / on_predicted / permutations / old_liver / / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / fasta_files / old_liver / / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / pvalues / old_liver.miRNA_seeds.permutations.txt
#
# python
# add_permutation_information.py / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / real_seeds / young_cerebellum.fimo.txt / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / on_predicted / permutations / young_cerebellum / / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / fasta_files / young_cerebellum / / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / pvalues / young_cerebellum.miRNA_seeds.permutations.txt
# python
# add_permutation_information.py / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / real_seeds / young_hippocampus.fimo.txt / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / on_predicted / permutations / young_hippocampus / / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / fasta_files / young_hippocampus / / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / pvalues / young_hippocampus.miRNA_seeds.permutations.txt
# python
# add_permutation_information.py / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / real_seeds / young_liver.fimo.txt / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / on_predicted / permutations / young_liver / / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / fasta_files / young_liver / / beegfs / group_dv / home / FMetge / projects / Franzi / circRNA / mousedata_fuchs / circRNA / miRNA_seeds / circRNA_wise / pvalues / young_liver.miRNA_seeds.permutations.txt
