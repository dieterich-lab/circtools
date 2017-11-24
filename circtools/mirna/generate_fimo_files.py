# script to generate slurm files for fimo analysis


if __name__ == '__main__':
    import os
    import argparse

    parser = argparse.ArgumentParser(description='')
    # input
    parser.add_argument('meme_folder', metavar='PATH', help='path to meme files with permuations')
    parser.add_argument('script_folder', metavar='PATH', help='path to general script folder')
    parser.add_argument('fasta_folder', metavar='PATH', help='path to general fasta folder')
    parser.add_argument('out_folder', metavar='PATH', help='path for general output folder')
    parser.add_argument('sample', metavar='sample_name',
                        help='sample name, specifying the subfolder in each general folder, shoudl be the same in every general folder')
    parser.add_argument('-n', dest='num_perm', default=1000,
                        help='the number of permuations you want to perform. Also the number of meme files in the meme folder')
    parser.add_argument('-b', dest='biotype', default='miRNA', help='either miRNA or rbp')
    args = parser.parse_args()

    meme_folder = args.meme_folder
    script_folder = args.script_folder
    fasta_folder = args.fasta_folder
    out_folder = args.out_folder
    sample = args.sample

    num_perm = args.num_perm
    biotype = args.biotype

    circRNAfiles = os.listdir('%s/%s' % (fasta_folder, sample))

    header = '#!/usr/bin/bash\n#SBATCH -c 1\n#SBATCH -n 1\n#SBATCH -p blade,long\n\n\n'
    array_submission_file = open('%s/run_%s.slurm' % (script_folder, sample), 'w')

    if not os.path.exists('%s/%s' % (out_folder, sample)):
        os.mkdir('%s/%s' % (out_folder, sample))
    if not os.path.exists('%s/%s' % (script_folder, sample)):
        os.mkdir('%s/%s' % (script_folder, sample))

    for circRNA in circRNAfiles:
        O = open('/%s/%s/submit_fimo_%s.slurm' % (script_folder, sample, circRNA.replace('.fa', '')), 'w')
        O.write(header)
        O.write('mkdir %s/%s/%s \n\n' % (out_folder, sample, circRNA.replace('.fa', '')))
        for permutation in range(0, 1000):
            O.write('fimo --thresh 0.001 --norc --o %s/%s/%s/permutation_%s/ %s/%s.permutation.%s.meme %s/%s/%s\n' % (
            out_folder, sample, circRNA.replace('.fa', '/'), permutation, meme_folder, biotype, permutation,
            fasta_folder, sample, circRNA))
        O.close()

        array_submission_file.write(
            'sbatch /%s/%s/submit_fimo_%s.slurm\n\n' % (script_folder, sample, circRNA.replace('.fa', '')))

    array_submission_file.close()
