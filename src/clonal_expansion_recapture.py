import numpy as np
from util import *
from collections import defaultdict
from Bio import SeqIO

save_as='.fasta.gz'

samples = [
# prefix, template input
["p1_p17_2012-10-02",820],
["p1_p17_2014-10-07",148],
["p1_p17_2015-02-12",38],
["p8_p17_2012-09-21", 180],
["p8_p17_2015-04-07", 175],
["p2_p17_2015-06-09",75],
["p4_p17_2012-09-21",216],
["p4_p17_2014-10-13",368],
["p3_p17_2012-09-18",243],
["p3_p17_2014-10-13",102],
["p3_p17_2015-04-24",108],
["p4_p17_2012-09-28",52],
["p4_p17_2014-10-30",184],
["p4_p17_2015-02-25",53],
["p5_p17_2012-10-26",180],
["p5_p17_2015-03-16",72],
["p6_p17_2012-10-24",115],
["p6_p17_2014-11-03",15],
["p8_p17_2014-11-18",279],
["p7_p17_2015-04-17",108],
["p7_p17_2014-12-01",28],
["p9_p17_2012-10-05",60],
["p9_p17_2014-10-02",72],
["p9_p17_2015-03-18",72],
["p10_p17_2012-10-09",249],
["p10_p17_2014-10-24",116],
["p10_p17_2015-02-27",51],
["p11_p17_2012-10-10",124],
["p11_p17_2014-10-22",120],
["p11_p17_2015-02-25",123],
]
outprefix_to_template_input = {x[0]:x[-1] for x in samples}

def clonal_expansion(mt='clustered_good'):
    '''
    make list of haplotypes, record their abundance, plot abundance
    trajectories for individual haplotypes, and the distribution of
    haplotype frequencies across all samples.
    '''
    diversity = defaultdict(list)
    all_abundances = []
    haplos = {} # dict of dict of haplotypes at different time points
    # make cumulative distributions of haplotype counts
    for pcode in patient_to_prefix_p17:
        dna_samples = {}
        haplos[pcode] = defaultdict(list)
        for outprefix in patient_to_prefix_p17[pcode]:
            fname = 'data/'+outprefix+'_DNA_'+mt+save_as
            tmp_outprefix, pat, region, sdate = get_outprefix(fname)
            with myopen(fname) as ifile:
                dna_haps = []
                for si, seq in enumerate(SeqIO.parse(ifile, 'fasta')):
                    dna_haps.append(seq.upper())
            nreads = 1.0*np.sum([get_DNA_read_count(seq.name) for seq in dna_haps])
            seq_subset = prune_rare_DNA(dna_haps, 0.002)

            for read in seq_subset:
                # add sampling date, haplotype abundance to the sequence of the haplotype
                haplos[pcode][str(read.seq.ungap('-'))].append((sdate, get_DNA_read_count(read.name)/nreads))
            dna_samples[sdate] = seq_subset
            abundances = np.array([get_DNA_read_count(seq.name)/nreads for seq in seq_subset])
            all_abundances.extend(abundances)
            diversity[pcode].append((sdate, abundances))

        plt.yscale('log')
        plt.xscale('log')
        plt.legend()

    # make histogram of abundances
    plt.figure()
    plt.hist(all_abundances, bins=np.logspace(-3,0,31), bottom=0.5)
    plt.ylabel('number of haplotypes', fontsize=fs)
    plt.xlabel('fraction of reads', fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.ylim(0.8, 300)
    plt.yscale('log')
    plt.xscale('log')
    for fmt in formats:
        plt.savefig('figures/read_fraction_distribution_'+mt+fmt)

    # plot frequencies of haplotypes observed at multiple time points
    for pcode in patient_to_prefix_p17:
        if len(patient_to_prefix_p17[pcode])<2:
            continue
        plt.figure()
        plt.title(pcode)
        pi = int(pcode[1:])-1
        treat_start = numdate(treatment_dates.loc[pi,'treatment start'])
        for outprefix in patient_to_prefix_p17[pcode]:
            tmp, pc, a, sdate = get_outprefix(outprefix)
            plt.text(sdate-treat_start, 0.5, str(outprefix_to_template_input[outprefix]))

        for read, abu in haplos[pcode].iteritems():
            if len(abu)>1 or np.max([x[1] for x in abu])>0.05:
                tp = np.array(abu)
                plt.plot(tp[:,0]-treat_start, tp[:,1], '-o')
        plt.yscale('log')
        plt.ylabel('haplotype frequency')
        plt.xlabel('time since treatment start')
        plt.ylim(0.008, 1.0)

    return haplos

def recapture_statistics(haplos, mt='clustered_good', latex=True):
    '''
    count haplotypes observed above different frequency thresholds and
    determine what fraction is recovered in different samples from the same
    patient
    '''
    # reorient the haplo data structure to study persistence and recapture
    hap_abundance = {} #[pcode][sample] = [ (freq, n_recapture, templates) for each read]
    for pcode in patient_to_prefix_p17:
        hap_abundance[pcode]=dict()
        for outprefix in patient_to_prefix_p17[pcode]:
            op = outprefix.split('_')[-1]
            hap_abundance[pcode][op]=dict()
            fname = 'data/'+outprefix+'_DNA_'+mt+save_as
            tmp_outprefix, pat, region, sdate = get_outprefix(fname)
            n_temp = outprefix_to_template_input[outprefix]
            for hap in haplos[pcode]:
                abundance_data = haplos[pcode][hap]
                if sdate in [x[0] for x in abundance_data]:
                    tp = [x[0] for x in abundance_data].index(sdate)
                    hap_abundance[pcode][op][hap] = (abundance_data[tp][1], len(abundance_data), n_temp)

    if latex:
        sep = '\t&\t'
        lend = '\\\\\n'
    else:
        sep='\t'
        lend='\n'

    with open('data/recapture_statistics_'+mt+'.tsv', 'w') as recap:
        recap.write(sep.join(['patient', 'sample', '#templates', '#haplotypes>0.002', '#haplotypes>0.01', 'fraction recaptured'])+lend)
        for pcode in pcodes:
            if pcode=='p4':
                continue
            for sdate in sorted(hap_abundance[pcode]):
                if len(hap_abundance[pcode][sdate])==0:
                    continue
                tmp = np.array(hap_abundance[pcode][sdate].values())
                good_haps = tmp[:,0]>0.01
                rare_haps = tmp[:,0]>0.002
                recap_frac = 'NA' if pcode=='p2' else format(np.mean(tmp[good_haps,1]>1), '0.2f')
                recap.write(sep.join(map(str, [pcode, sdate, int(tmp[0,-1]), np.sum(rare_haps), np.sum(good_haps), recap_frac]))+lend)


if __name__=="__main__":
    import pandas as pd
    treatment_dates = pd.read_excel('data/2016-10-10_treatment_start_dates.xlsx')

    for mt in ['clustered_good', 'hyper']:
        haplos = clonal_expansion(mt)
        recapture_statistics(haplos, mt)

