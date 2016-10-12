import numpy as np
from Bio import SeqIO, Seq
import cPickle as pickle
import pandas as pd
from util import *

fs=16
alpha="ACGT-N"
save_as='.fasta.gz'

def hypermutation_statistics():
    '''
    use pre-calculated files on mutation prevalence in different reads types to
    produce analysis figures
    '''
    mut_files = glob.glob('data/*_mutation_statistics.pkl.gz')
    hyper_mut_stats = {}
    for fname in mut_files:
        outprefix, pat, region, sdate = get_outprefix(fname)
        with myopen(fname) as ifile:
            hyper_mut_stats[outprefix] = pickle.load(ifile)

    #reformat data structure
    from itertools import product
    mutations = [a+'->'+d for a,d in product(alpha, alpha) if a!=d]
    all_muts = {'good':defaultdict(list), 'hyper':defaultdict(list)}
    for sample in hyper_mut_stats:
        for mtype in all_muts:
            for mi,mut in enumerate(mutations):
                if len(hyper_mut_stats[sample][mtype]):
                    all_muts[mtype][mut].extend(hyper_mut_stats[sample][mtype][:,mi])
                else:
                    print("no data on",mtype,sample)
    # plot distribitions of different transitions in sequences classified as good and hyper mutant
    plt.figure()
    tmp = []
    for mtype,muts in all_muts.iteritems():
        for mut in ['G->A', 'A->G', 'C->T', 'T->C']:
            tmp.extend([(mut, mtype, x) for x in muts[mut]])
    data = pd.DataFrame(tmp, columns=['mutation', 'type', 'count'])
    data.loc[data.loc[:,'type']=='hyper', 'type'] = 'hypermutants'
    data.loc[data.loc[:,'type']=='good', 'type'] = 'other reads'
    sns.violinplot(x='mutation', y='count', hue='type',data=data, inner=None,
                   split=True, scale='width', bw=0.5, cut=0)
    plt.legend(title=None, fontsize=fs*0.8)
    plt.tick_params(labelsize = fs*0.8)
    plt.ylabel('number of mutations', fontsize=fs)
    plt.xlabel('', fontsize=fs)
    plt.savefig('figures/mutation_count.pdf')

    # plot the diffferent types of muations in different sequence subsets
    for seqtype in ['good', 'hyper', 'suspicious']:
        plt.figure()
        for sample in hyper_mut_stats:
            plt.plot(hyper_mut_stats[sample][seqtype].mean(axis=0), label=sample+str(len(hyper_mut_stats[sample][seqtype])))
        plt.ylabel('mean number of mutations')
        plt.legend(loc=1, fontsize=8)
        plt.xticks(np.arange(30), mutations, rotation=60, fontsize=8)
        plt.savefig('figures/mean_number_of_mutations_'+seqtype+'.pdf')


def hyper_mut_ratio(latex=True):
    '''
    calculate the fraction of reads that are hypermutated, the number of stop codons in good and
    reads, as well as the fraction of reads that are obviously impaired
    '''
    good_files = filter(lambda x:'RNA' not in x, glob.glob('data/*_DNA_clustered_good'+save_as))
    hyper_files = filter(lambda x:'RNA' not in x, glob.glob('data/*_DNA_hyper'+save_as))
    read_count = defaultdict(dict)
    for fnames, seqtype in [(good_files, 'good'),(hyper_files, 'hyper')]:
        for fname in fnames:
            outprefix, pat, region, sdate = get_outprefix(fname)
            print(outprefix)
            if region!='p17':
                continue

            read_count[outprefix][seqtype] = 0
            stop_distribution = []
            with myopen(fname) as ifile:
                for si, seq in enumerate(SeqIO.parse(ifile, 'fasta')):
                    nread = get_DNA_read_count(seq.name)
                    read_count[outprefix][seqtype] += nread
                    seq_str = str(seq.seq)
                    nstops = Seq.translate(seq_str.replace('-','')[20:]).count('*')    # p17 starts at position 20 in the amplicon
                    stop_distribution.append((nread, nstops))
            stop_distribution = np.array(stop_distribution, dtype=float)
            if len(stop_distribution):
                read_count[outprefix][seqtype+'_stops'] = np.sum(stop_distribution[:,0]*stop_distribution[:,1])/np.sum(stop_distribution[:,0])
                read_count[outprefix][seqtype+'_frac_stops'] = np.sum(stop_distribution[:,0]*(stop_distribution[:,1]>0))/np.sum(stop_distribution[:,0])
            else:
                read_count[outprefix][seqtype+'_stops'] = np.nan
                read_count[outprefix][seqtype+'_frac_stops'] = np.nan

    if latex:
        sep = '\t&\t'
        lend = '\\\\\n'
    else:
        sep='\t'
        lend='\n'

    with open('data/hypermutation_statistic.tsv', 'w') as hyper_mut_stat:
        hyper_mut_stat.write('\t'.join(map(str,['pcode', 'sample', '#reads', 'fraction_hyper', 'avg_stop_good',
                                        'avg_stop_hyper', 'frac_stop_good', 'frac_stop_hyper']))+'\n')
        for pcode in pcodes:
            if pcode=='p4':
                continue
            for key in sorted(patient_to_prefix_p17[pcode], key=lambda x:x.split('_')[-1]):
                k = key #.split('_')[0]+'_' + key.split('_')[2]
                val = read_count[k]
                outstr = sep.join(map(str, [pcode, k.split('_')[-1], (val['good']+val['hyper']),
                                   1.0*val['hyper']/(val['good']+val['hyper']), val['good_stops'],
                                   val['hyper_stops'], val['good_frac_stops'], val['hyper_frac_stops']]))
                print(outstr)
                hyper_mut_stat.write(outstr+lend)


if __name__=="__main__":
    hyper_mut_ratio()
    hypermutation_statistics()

