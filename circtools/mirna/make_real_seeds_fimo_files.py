# script to write sbatch scripts for submitting to cluster
import os

# miRNA seeds
samples = ['old_cerebellum', 'old_hippocampus', 'old_liver', 'young_cerebellum', 'young_hippocampus', 'young_liver']
header = '#!/usr/bin/bash\n#SBATCH --time=00:30:00\n#SBATCH --mem=3G\n#SBATCH --ntasks=1\n#SBATCH -p hugemem,himem\n#SBATCH --cpus-per-task=1\n\n\n'
mirna = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/miRNA_seeds/samplewise/db/TargetScanMm10.dna_encoded.meme'
script_folder = '/beegfs/group_dv/home/FMetge/home/circRNA/FUCHS/mirna_real_seeds'

folder = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/fasta_files'
for s in samples:
    output = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/miRNA_seeds/circRNA_wise/real_seeds/%s/' % (
    s)
    submit_file = open('%s/submit_all_%s.slurm' % (script_folder, s), 'w')
    files = os.listdir('%s/%s' % (folder, s))
    for f in files:
        fasta_file = '%s/%s/%s' % (folder, s, f)
        fimo_file = open('%s/%s/submit_%s.slurm' % (script_folder, s, f.replace('.fa', '')), 'w')
        fimo_file.write(header)
        fimo_file.write('fimo --verbosity 1 --thresh 0.0001 --norc --oc %s/%s %s %s\n ' % (
        output, f.replace('.fa', ''), mirna, fasta_file))
        fimo_file.close()
        submit_file.write('sbatch %s/%s/submit_%s.slurm\n' % (script_folder, s, f.replace('.fa', '')))
    submit_file.close()

samples = ['frontal_cortex_A', 'frontal_cortex_B', 'frontal_cortex_C', 'frontal_cortex_D', 'frontal_cortex_E',
           'frontal_cortex_F']
header = '#!/usr/bin/bash\n#SBATCH --time=00:30:00\n#SBATCH --mem=3G\n#SBATCH --ntasks=1\n#SBATCH -p hugemem,himem\n#SBATCH --cpus-per-task=1\n\n\n'
mirna = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/miRNA_seeds/samplewise/db/TargetScanMm10.dna_encoded.meme'
script_folder = '/beegfs/group_dv/home/FMetge/home/circRNA/FUCHS/mirna_real_seeds'

folder = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/raw_frontal_cortex/fasta_files'
for s in samples:
    output = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/raw_frontal_cortex/miRNA_seeds/circRNA_wise/real_seeds/%s/' % (
    s)
    submit_file = open('%s/submit_all_%s.slurm' % (script_folder, s), 'w')
    files = os.listdir('%s/%s' % (folder, s))
    for f in files:
        fasta_file = '%s/%s/%s' % (folder, s, f)
        fimo_file = open('%s/%s/submit_%s.slurm' % (script_folder, s, f.replace('.fa', '')), 'w')
        fimo_file.write(header)
        fimo_file.write('fimo --verbosity 1 --thresh 0.0001 --norc --oc %s/%s %s %s\n ' % (
        output, f.replace('.fa', ''), mirna, fasta_file))
        fimo_file.close()
        submit_file.write('sbatch %s/%s/submit_%s.slurm\n' % (script_folder, s, f.replace('.fa', '')))
    submit_file.close()

# RBP motifs

samples = ['old_cerebellum', 'old_hippocampus', 'old_liver', 'young_cerebellum', 'young_hippocampus', 'young_liver']
header = '#!/usr/bin/bash\n#SBATCH --time=00:30:00\n#SBATCH --mem=3G\n#SBATCH --ntasks=1\n#SBATCH -p hugemem,himem\n#SBATCH --cpus-per-task=1\n\n\n'
mirna = '/beegfs/group_dv/home/FMetge/genomes/mus_musculus/GRCm38_79/custom_annotation/uniprobe_mouse.meme'
script_folder = '/beegfs/group_dv/home/FMetge/home/circRNA/FUCHS/motif_enrichment'

folder = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/fasta_files'
for s in samples:
    output = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/circRNA/motif_enrichment/circRNA_wise/real_seeds/%s/' % (
    s)
    submit_file = open('%s/submit_all_%s.slurm' % (script_folder, s), 'w')
    files = os.listdir('%s/%s' % (folder, s))
    for f in files:
        fasta_file = '%s/%s/%s' % (folder, s, f)
        fimo_file = open('%s/%s/submit_%s.slurm' % (script_folder, s, f.replace('.fa', '')), 'w')
        fimo_file.write(header)
        fimo_file.write('fimo --verbosity 1 --thresh 0.0001 --norc --oc %s/%s %s %s\n ' % (
        output, f.replace('.fa', ''), mirna, fasta_file))
        fimo_file.close()
        submit_file.write('sbatch %s/%s/submit_%s.slurm\n' % (script_folder, s, f.replace('.fa', '')))
    submit_file.close()

samples = ['frontal_cortex_A', 'frontal_cortex_B', 'frontal_cortex_C', 'frontal_cortex_D', 'frontal_cortex_E',
           'frontal_cortex_F']
header = '#!/usr/bin/bash\n#SBATCH --time=00:30:00\n#SBATCH --mem=3G\n#SBATCH --ntasks=1\n#SBATCH -p hugemem,himem\n#SBATCH --cpus-per-task=1\n\n\n'
mirna = '/beegfs/group_dv/home/FMetge/genomes/mus_musculus/GRCm38_79/custom_annotation/uniprobe_mouse.meme'
script_folder = '/beegfs/group_dv/home/FMetge/home/circRNA/FUCHS/motif_enrichment'

folder = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/raw_frontal_cortex/fasta_files'
for s in samples:
    output = '/beegfs/group_dv/home/FMetge/projects/Franzi/circRNA/mousedata_fuchs/raw_frontal_cortex/motif_enrichment/circRNA_wise/real_seeds/%s/' % (
    s)
    submit_file = open('%s/submit_all_%s.slurm' % (script_folder, s), 'w')
    files = os.listdir('%s/%s' % (folder, s))
    for f in files:
        fasta_file = '%s/%s/%s' % (folder, s, f)
        fimo_file = open('%s/%s/submit_%s.slurm' % (script_folder, s, f.replace('.fa', '')), 'w')
        fimo_file.write(header)
        fimo_file.write('fimo --verbosity 1 --thresh 0.0001 --norc --oc %s/%s %s %s\n ' % (
        output, f.replace('.fa', ''), mirna, fasta_file))
        fimo_file.close()
        submit_file.write('sbatch %s/%s/submit_%s.slurm\n' % (script_folder, s, f.replace('.fa', '')))
    submit_file.close()
