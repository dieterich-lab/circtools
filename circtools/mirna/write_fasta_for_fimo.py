# python script to get circRNA and hostgene fasta into one fasta file 

# define functions

def read_bedfile_circRNA(infile):
    I = open(infile)
    exons = {}
    for line in I:
        if not line.startswith('#'):
            L = line.replace('\n', '').split('\t')
            circID = L[3].split('|')[0]
            transcripts = L[3].split('|')[1].split(',')
            if not (L[0], int(L[1]), int(L[2])) in exons:
                exons[(L[0], int(L[1]), int(L[2]))] = {'circles': {}, 'sequence': ''}
            exons[(L[0], int(L[1]), int(L[2]))]['circles'][circID] = {'transcripts': transcripts, 'score': int(L[4]),
                                                                      'strand': L[5]}
    return (exons)


def read_bedfile_hostgene(infile):
    exons = {}
    I = open(infile)
    for line in I:
        L = line.replace('\n', '').split('\t')
        if not (L[0], int(L[1]), int(L[2])) in exons:
            exons[(L[0], int(L[1]), int(L[2]))] = {'transcripts': [], 'strand': L[5], 'sequence': ''}
        exons[(L[0], int(L[1]), int(L[2]))]['transcripts'] += [L[3].split('_')[0]]
    I.close()
    return (exons)


def read_circle_annotation(infile):
    I = open(infile)
    I.readline()
    circle_annotation = {}
    for line in I:
        L = line.replace('\n', '').split('\t')
        if not '%s:%s-%s' % (L[0], L[1], L[2]) in circle_annotation:
            circle_annotation['%s:%s-%s' % (L[0], L[1], L[2])] = {'gene': L[3], 'strand': [L[5]]}
        else:
            circle_annotation['%s:%s-%s' % (L[0], L[1], L[2])]['strand'] += [L[5]]
            if not L[3] == 'N/A':
                circle_annotation['%s:%s-%s' % (L[0], L[1], L[2])]['gene'] = L[3]
    I.close()
    return (circle_annotation)


def fetch_fasta(exons, fastafile):
    ref = pysam.Fastafile(fastafile)
    for e in exons:
        if ref.__contains__(e[0]):
            exons[e]['sequence'] = ref.fetch(e[0], e[1], e[2]).upper()
        else:
            exons[e]['sequence'] = ''
    return (exons)


def reconstruct_transcripts_circRNA(exons):
    transcripts = {}
    for e in exons:
        for c in exons[e]['circles']:
            if not c in transcripts:
                transcripts[c] = {}
            for t in exons[e]['circles'][c]['transcripts']:
                if not t in transcripts[c]:
                    transcripts[c][t] = {}
                transcripts[c][t][e] = {'sequence': exons[e]['sequence'], 'score': exons[e]['circles'][c]['score']}
    return (transcripts)


def reconstruct_transcripts_hostgene(exons):
    transcripts = {}
    for e in exons:
        for t in exons[e]['transcripts']:
            if not t in transcripts:
                transcripts[t] = {'exons': {}, 'strand': exons[e]['strand']}
            transcripts[t]['exons'][e] = exons[e]['sequence']
    return (transcripts)


def paste_sequences_circRNA(transcripts):
    fasta = {}
    for c in transcripts:
        for t in transcripts[c]:
            transcripts_fasta = ''
            average_score = 0
            sorted_exons = sorted(transcripts[c][t].keys())
            for e in sorted_exons:
                transcripts_fasta += transcripts[c][t][e]['sequence']
                average_score += transcripts[c][t][e]['score']
            fasta[(c, t)] = {'sequence': transcripts_fasta, 'score': average_score / len(sorted_exons)}
    return (fasta)


def paste_sequences_hostgene(transcripts):
    fasta = {}
    for t in transcripts:
        transcripts_fasta = ''
        sorted_exons = sorted(transcripts[t]['exons'].keys())
        for e in sorted_exons:
            transcripts_fasta += transcripts[t]['exons'][e]
        fasta[t] = {'sequence': transcripts_fasta, 'strand': transcripts[t]['strand']}
    return (fasta)


def read_id_file(infile):
    ids = {}
    names = {}
    I = open(infile)
    for Line in I:
        ids[Line.split('\t')[0]] = Line.replace('\n', '').split('\t')[-1]
        if not Line.replace('\n', '').split('\t')[-1] in names:
            names[Line.replace('\n', '').split('\t')[-1]] = [Line.split('\t')[0]]
        else:
            names[Line.replace('\n', '').split('\t')[-1]] += [Line.split('\t')[0]]
    return (ids, names)


def aggregate_circRNAs_by_gene(circAnnotation, circRNA):  # we take circRNA_fasta here
    genewise = {}
    for c in circRNA:
        if c[0] in circAnnotation:
            if not circAnnotation[c[0]]['gene'] in genewise:
                genewise[circAnnotation[c[0]]['gene']] = {}
            genewise[circAnnotation[c[0]]['gene']][c] = circRNA[c]
            genewise[circAnnotation[c[0]]['gene']][c]['strand'] = circAnnotation[c[0]]['strand']
    return (genewise)


def reverse_complement(sequence):
    letter_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    revcomp = ''
    for lola in reversed(sequence):
        revcomp += letter_dict[lola]
    return (revcomp)


def write_outfiles(outfile, circRNAs, genes, genes_fasta):
    for gene in circRNAs:
        if not gene == 'N/A':
            O = open('%s/%s.fa' % (outfile, gene), 'w')
            for circle in circRNAs[gene]:
                for strand in circRNAs[gene][circle]['strand']:
                    O.write('>%s_%s_%s_%s\n' % (circle[0], circle[1], circRNAs[gene][circle]['score'], strand))
                    if strand == '-':
                        O.write('%s\n' % (reverse_complement(circRNAs[gene][circle]['sequence'])))
                    else:
                        O.write('%s\n' % (circRNAs[gene][circle]['sequence']))
                if gene in genes:
                    for transcript in genes[gene]:
                        if transcript in genes_fasta:
                            if genes_fasta[transcript]['strand'] == '-':
                                O.write('>%s_%s\n%s\n' % (transcript, genes_fasta[transcript]['strand'],
                                                          reverse_complement(genes_fasta[transcript]['sequence'])))
                            else:
                                O.write('>%s_%s\n%s\n' % (
                                transcript, genes_fasta[transcript]['strand'], genes_fasta[transcript]['sequence']))
            O.close()
        else:
            for circle in circRNAs[gene]:
                O = open('%s/%s.fa' % (outfile, circle[0]), 'w')
                for strand in circRNAs[gene][circle]['strand']:
                    O.write('>%s_%s_%s_%s\n' % (circle[0], circle[1], circRNAs[gene][circle]['score'], strand))
                    if strand == '-':
                        O.write('%s\n' % (reverse_complement(circRNAs[gene][circle]['sequence'])))
                    else:
                        O.write('%s\n' % (circRNAs[gene][circle]['sequence']))
                O.close()
    return


# run script

if __name__ == '__main__':

    import pysam
    import argparse
    import tempfile

    parser = argparse.ArgumentParser(description='')
    # input
    parser.add_argument('fasta', metavar='fasta.fa', help='reference fasta file')
    parser.add_argument('circRNA', metavar='circrna.bed', help='bedfile containing exon coordinates')
    parser.add_argument('hostgene', metavar='hostgene.bed', help='bedfile containing exon coordinates')
    parser.add_argument('transcript_gene_list', metavar='id_file',
                        help='tab-separated file, first column contains the transcript ids second column contains the gene name')
    parser.add_argument('-o', dest='outfolder', default='.', help='file, sequences should be written to')
    parser.add_argument('-A', dest='CircCoordinates', default='',
                        help='give a circle Annotation file to add strand and gene name to circles, right now this should be CircCoordinates from DCC')

    # options

    args = parser.parse_args()

    # parse arguments

    circRNA_file = args.circRNA
    hostgene_file = args.hostgene
    fastafile = args.fasta
    CircCoordinates = args.CircCoordinates
    id_file = args.transcript_gene_list
    outfolder = args.outfolder

    if not CircCoordinates == '':
        circAnnotation = read_circle_annotation(CircCoordinates)
    else:
        circAnnotation = {}

    # get host gene information
    Hostgenes = read_bedfile_hostgene(hostgene_file)
    Hostgenes = fetch_fasta(Hostgenes, fastafile)
    Hostgene_Transcripts = reconstruct_transcripts_hostgene(Hostgenes)
    Hostgene_Fasta = paste_sequences_hostgene(Hostgene_Transcripts)
    # get circRNA information
    circRNAs = read_bedfile_circRNA(circRNA_file)
    circRNAs = fetch_fasta(circRNAs, fastafile)
    circRNA_Transcripts = reconstruct_transcripts_circRNA(circRNAs)
    circRNA_Fasta = paste_sequences_circRNA(circRNA_Transcripts)
    # read gene_id list
    IDs, GENES = read_id_file(id_file)
    # aggregate_circRNAs by gene
    circRNAS_by_gene = aggregate_circRNAs_by_gene(circAnnotation, circRNA_Fasta)
    # write out file
    write_outfiles(outfolder, circRNAS_by_gene, GENES, Hostgene_Fasta)
