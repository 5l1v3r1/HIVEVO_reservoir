import numpy as np
from util import *
import sys; sys.path.append('../HIVEVO_access'); sys.path.append('src')
from hivevo.patients import Patient
from root_to_tip import get_treatment_start_date
from Bio.Align import MultipleSeqAlignment

alpha="ACGT-N"
save_as='.fasta.gz'
cols = sns.color_palette(n_colors=6)

def calc_af(aln, alpha):
    aln_array = np.array(aln)
    af = np.zeros((len(alpha), aln_array.shape[1]))
    for ai, state in enumerate(alpha):
        af[ai] += (aln_array==state).mean(axis=0)
    af[-1] = 1.0 - af[:-1].sum(axis=0)
    return af

def sample_from_RNA_tp(tp,n):
    from random import sample
    if (len(tp))==0:
        return []
    n_reads = [(x,get_RNA_read_count(x.description)) for x in tp]

    cum_reads = np.cumsum([x[1] for x in n_reads])
    total_reads = cum_reads[-1]
    read_indices = np.random.randint(total_reads,size=n)

    reads = [n_reads[ii][0] for ii in np.searchsorted(cum_reads,read_indices)]
    print(cum_reads,n, len(reads))
    return reads


def sample_from_RNA(pat, cell_pop, t_ART, n=100):
    proportions = np.zeros_like(pat.dsi)
    print(pat.name)
    for prop, t_decay in cell_pop:
        #age_of_samples = t_ART + (get_treatment_start_date(p.name)-p.dsi)
        age_of_samples = t_ART + (pat.dsi.max()-pat.dsi)
        tmp = np.exp(-age_of_samples/t_decay)
        proportions+= prop*tmp
    proportions/=proportions.sum()
    print(proportions)

    rna_haps = pat.get_haplotype_alignment('p17')
    rna_haps_by_tp = defaultdict(list)
    for rhap in rna_haps: # sort RNA reads by time points
        rna_haps_by_tp[int(rhap.name.split('_')[1])].append(rhap)

    print(sorted(rna_haps_by_tp.keys()),
          [len(rna_haps_by_tp[k]) for k in sorted(rna_haps_by_tp.keys())],pat.dsi)
    reads = []
    for tp, p in zip(pat.dsi, proportions):
        reads.extend(sample_from_RNA_tp(rna_haps_by_tp[int(tp)], int(p*n)))

    return reads



if __name__=="__main__":
    patients = {pcode:Patient.load(pcode) for pcode in pcodes}
    import pandas as pd
    #treatment_dates = pd.read_excel('data/2016-01-08_treatment_start_dates.xlsx')
    sns.set_style('darkgrid')
    # lifetimes of cells -- 95% short lived cells, 5% long lived cells
    cell_pop = [(0.9, 30.),(0.1, 500.)]

    ###
    # make read samples from pre treatment RNA reads
    ###
    reads={}
    sampling_times = [0,90,180] # at ART start, 3 month, 6 month
    for pcode, p in patients.iteritems():
        reads[pcode]=[]
        for t_ART in sampling_times:
            reads[pcode].append(sample_from_RNA(p, cell_pop, t_ART, n=1000))

    ###
    # calculate the divergence of these samples
    ###
    divergence = {}
    for pcode, p in patients.iteritems():
        if pcode=='p7':
            continue
        msa = MultipleSeqAlignment(reads[pcode][0])
        af_initial = calc_af(msa, alpha)
        af_initial/=af_initial[:4].sum(axis=0)
        good_ind_initial = ~(np.any(np.isnan(af_initial)|np.isinf(af_initial), axis=0))
        consensus = af_initial.argmax(axis=0)
        divergence[pcode]=[]
        for sample in reads[pcode]:
            if len(sample):
                msa = MultipleSeqAlignment(sample)
                af = calc_af(msa, alpha)
                af/=af[:4].sum(axis=0)
                good_ind_sample = ~(np.any(np.isnan(af)|np.isinf(af), axis=0))
                ind = good_ind_sample&good_ind_initial
                print(pcode, np.isnan(af).sum(), np.isinf(af).sum())
    #            divergence[pcode].append(np.mean(af[:4,:].sum(axis=0)
    #                                             -af[consensus,np.arange(consensus.shape[0])]))
                divergence[pcode].append(1-np.mean((af[:4,ind]*af_initial[:4,ind]).sum(axis=0)))

    # make figure
    plt.ion()
    plt.show()
    fig,axs = plt.subplots(1,2, figsize = (12,6))

    t = np.arange(0,1000,30)
    for ti,tshift in enumerate([90,0]):
        for ci, (prop, tau) in enumerate(cell_pop):
            tmp_t = (tshift - t)
            axs[0].plot(tmp_t/365.25, 1.2*ti+prop*np.exp(-t/tau), c=cols[ci*2], alpha=0.5)
            axs[0].fill_between(tmp_t/365.25, 1.2*ti+prop*np.exp(-t/tau),
                                1.2*ti*np.ones_like(t), color=cols[ci*2], alpha=0.1)
            axs[0].fill_between(tmp_t[tmp_t<=0]/365.25, 1.2*ti+prop*np.exp(-t[tmp_t<=0]/tau),
                                1.2*ti*np.ones_like(t[tmp_t<=0]), color=cols[ci*2], alpha=0.3)
    axs[0].plot([0,0],[0,2.2], c='k', alpha=0.5)
    axs[0].get_yaxis().set_visible(False)
    axs[0].set_xlabel('years relative to treatment start', fontsize=fs)
    axs[0].tick_params(labelsize=fs)
    axs[0].set_xlim([-2,0.7])

    div_array = np.array(divergence.values()).T
    div = np.concatenate([x-x[0] for x in div_array.T])
    st_years=np.array(sampling_times)/365.25
    t = np.concatenate([st_years for x in div_array.T])
    for pi, tmp_div in enumerate(div_array.T):
        axs[1].plot(st_years, tmp_div-tmp_div[0], ls=':',
                    marker='o' if pi<5 else 'v', c=cols[pi%5]) #, label=str(pi)) # 'o' if pi<6 else 'v',


    from scipy.stats import linregress
    reg = linregress(t,div)
    reg = linregress(st_years, div_array.mean(axis=1)-div_array.mean(axis=1)[0])
    slope = np.sum(div*t)/np.sum(t*t)
    axs[1].plot(st_years, div_array.mean(axis=1)-div_array.mean(axis=1)[0], ls='-',
                lw=3, c='k', label='average') # 'o' if pi<6 else 'v',
    axs[1].plot(st_years, st_years*slope, ls='--',lw=3, c='k', label='%1.4f subs./year'%slope) # 'o' if pi<6 else 'v',
    axs[1].set_xlabel('years since start of ART therapy', fontsize=fs)
    axs[1].set_ylabel('backwards signal', fontsize=fs)
    axs[1].tick_params(labelsize=fs)
    axs[1].legend(loc=2, fontsize=fs)
    axs[1].set_yticks([-0.001, 0.0, 0.001, 0.002, 0.003,0.004])
    plt.tight_layout()
    print('With intercept: subs per month %1.5f'%(reg.slope/12))
    print('Without intercept: subs per month %1.5f'%(slope/12))

    plt.savefig('figures/back_sampling_illustration.svg')


