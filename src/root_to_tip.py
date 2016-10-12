import numpy as np
from util import *
from Bio import SeqIO
from seqanpy import align_ladder, align_overlap
import sys; sys.path.append('../HIVEVO_access'); sys.path.append('src')
from hivevo.patients import Patient

alpha="ACGT-N"
save_as='.fasta.gz'
cols = sns.color_palette(n_colors=6)


def weighted_percentile(vals, weights, perc):
    tmp = np.array(sorted(zip(vals, weights), key=lambda x:x[0]), dtype=float)
    cumulative_weights = tmp[:,1].cumsum()
    cumulative_weights /= cumulative_weights[-1]
    index_l = cumulative_weights.searchsorted(perc*0.01, side='left')
    index_r = cumulative_weights.searchsorted(perc*0.01, side='right')
    if index_r==index_l:
        return tmp[index_l,0]
    else:
        return np.mean(tmp[index_l:index_r,0])

def align_to_initial(reads, p):
    from hivevo.sequence import alphal
    seg, start, stop = 'F1', int(p.annotation['p17'].location.start)-20, int(p.annotation['p17'].location.end)+7
    ref_seq = "".join(p.get_initial_sequence(seg)[start:stop])
    aft =  p.get_allele_frequency_trajectories('F1')[:,:,start:stop]
    founder_indices = np.array([alphal.index(nuc) for nuc in ref_seq])

    for read in reads:
        score, ali_ref, ali_read = align_overlap(ref_seq, str(read.seq.ungap('-')))
        tmp_ali = np.vstack([np.fromstring(a, 'S1') for a in (ali_ref, ali_read)])
        try:
            unconserved = ~((aft[:,:4,:].max(axis=1)>0.99).all(axis=0))
            reference_aln = (tmp_ali[0]!='-')
            ungapped = (~np.any(tmp_ali=='-', axis=0))[reference_aln]
            unamb = (~np.any(tmp_ali=='N', axis=0))[reference_aln]

            good_positions = ungapped&unamb&unconserved
            read_indices = np.array([alphal.index(nuc) for nuc in tmp_ali[1]])[reference_aln][good_positions]

            read.prob = np.sum(np.log(aft[:, read_indices, good_positions] + 0.001), axis=1)

            good_positions = ungapped&unamb
            good_positions[:20]=False
            good_positions[-7:]=False
            read_indices = np.array([alphal.index(nuc) for nuc in tmp_ali[1]])[reference_aln][good_positions]

            read.distance = np.mean((founder_indices[good_positions]!=read_indices))
        except:
            import ipdb; ipdb.set_trace();

def get_treatment_start_date(pcode):
    pi = np.where(treatment_dates.code==pcode)[0][0]
    # treatment start relative to infection date
    return (treatment_dates.loc[pi,'treatment start'] - treatment_dates.loc[pi,'EDI']).days


def root_to_tip_distances(patients, mt = 'clustered_good', min_DNA_frac=0.001, count_haplotypes = False):
    root_to_tip = {}
    origin = {}
    # make a figure of the distribution of the distance distribution
    for pi, (pid, pcode) in enumerate(zip(p_indices, pcodes)): # loop over patients
        p=patients[pcode]
        root_to_tip[pcode] = {'RNA':[], 'DNA':[]}
        origin[pcode] = defaultdict(list)
        dsi_treatment = get_treatment_start_date(pcode)

        # calculate root to tip distance for RNA haplotypes
        # distance is saved as read.distance attribute
        rna_haps = p.get_haplotype_alignment('p17')
        align_to_initial(rna_haps, p)
        rna_haps_by_tp = defaultdict(list)
        for rhap in rna_haps: # sort RNA reads by time points
            rna_haps_by_tp[int(rhap.name.split('_')[1])].append(rhap)
        # make distance histograms for each time point
        for tp in sorted(rna_haps_by_tp.keys()):
            tp_rel_treatment = np.round((tp-dsi_treatment)/365.25,1)
            weights = np.array([get_RNA_read_count(x.description) for x in rna_haps_by_tp[tp]], dtype=float)
            weights /=weights.sum()
            root_to_tip[pcode]['RNA'].append((tp_rel_treatment,np.array([(read.distance, w)
                              for read,w in zip(rna_haps_by_tp[tp], weights)])))

        # calculate and plot the distance from the first RNA sample for each DNA haplotype
        for outprefix in patient_to_prefix_p17[pcode]:
            fname = 'data/'+outprefix+'_DNA_'+mt+save_as
            tmp_outprefix, pat, region, sdate  =get_outprefix(fname)
            with myopen(fname) as ifile:
                dna_haps = []
                for si, seq in enumerate(SeqIO.parse(ifile, 'fasta')):
                    dna_haps.append(seq.upper())

            dna_haps = ungap(dna_haps)
            seq_subset = prune_rare_DNA(dna_haps, min_DNA_frac)
            if len(seq_subset):
                # align to founder sequence and attach probability that a read
                # is sampled from one of the RNA populations. will be x.prop
                align_to_initial(seq_subset, p)

                if count_haplotypes:
                    # use each sequence exactly once
                    weights = np.array([1.0 for x in seq_subset], dtype=float)
                else:
                    # weigh histograms by the abundance of individual reads
                    weights = np.array([get_DNA_read_count(x.name) for x in seq_subset], dtype=float)
                total_reads = sum(weights)
                weights *= 1.0/total_reads
                ii = np.where(treatment_dates.code==pcode)[0][0]
                tp_rel_treatment = np.round(((pd.Timestamp(outprefix[-10:])-treatment_dates.loc[ii,'treatment start']).days)/365.25,1)
                root_to_tip[pcode]['DNA'].append(( tp_rel_treatment,
                    np.array([(x.distance,w) for x,w in zip(seq_subset, weights)])))
                origin[pcode][outprefix[-10:]] = np.array([(x.prob.argmax(), w)
                                            for x,w in zip(seq_subset, weights)], dtype=float)

    return root_to_tip, origin


def combined_root_to_tip_figures(root_to_tip, fname=None, mt='clustered_good', q='median', show_LR_rate = False):
    '''draw figure that shows root to tip distance vs time for all samples
    from all patients. times and distances are centered on treatment start'''
    fig = plt.figure()
    plt.plot([0,0], [-0.05, 0.02], c='k', lw=3, alpha=0.5)
    plt.plot([-10,30], [0.0, 0.0], c='k', lw=3, alpha=0.5)
    for pi,pcode in enumerate(pcodes):
        x_y_yerr = []
        for seqtype in ['RNA', 'DNA']:
            for tp, vals in root_to_tip[pcode][seqtype]:
                if len(vals)==0:
                    continue
                iqd = [weighted_percentile(vals[:,0], vals[:,1], x) for x in [25, 50, 75]]
                avg = np.sum(vals[:,0]*vals[:,1])/np.sum(vals[:,1])
                stddev = np.sqrt(np.sum((vals[:,0]-avg)**2*vals[:,1])/np.sum(vals[:,1]))
                if q=='median':
                    x_y_yerr.append([tp, avg] + iqd)
                else:
                    x_y_yerr.append([tp, avg] + [avg-stddev, avg, avg+stddev])

        x_y_yerr = np.array(x_y_yerr)
        x = x_y_yerr[:,0]
        y = x_y_yerr[:,-2]
        y_offset = np.mean(y[np.abs(x)<0.5])

        plt.errorbar(x, y-y_offset,yerr=[y-x_y_yerr[:,-3], x_y_yerr[:,-1]-y], ls='none',
                marker='o' if pi<5 else 'v', label = pcode, c=cols[pi%5]) # 'o' if pi<6 else 'v',

    plt.legend(loc=4, ncol=2)
    if mt!='hyper':
        plt.ylim([-0.04,0.02])
    else:
        plt.ylim([-0.04,0.04])
    ymax = plt.ylim()[1]-0.003
    if show_LR_rate: # add a grew cone illustrating where LR predict the DNA sequences should be
        t_tmp = np.array([0,5,10], dtype=float)
        plt.fill_between(t_tmp, t_tmp*0.624*0.001*12, t_tmp*1.224*0.001*12,
                        facecolor='k', linewidth=0, alpha=0.3)

    plt.annotate("", xy=(-10.0, ymax), xytext=(0.1, ymax),
                 arrowprops=dict(arrowstyle="|-|",lw=3,color='k',alpha=0.5))
    plt.annotate("", xy=(20.0, ymax), xytext=(0.5, ymax),
                 arrowprops=dict(arrowstyle="|-|",lw=3,color='k',alpha=0.5))
    plt.text(-5,ymax+0.0005,'RNA',fontsize=fs*1, horizontalalignment='center')
    plt.text(10,ymax+0.0005,'DNA',fontsize=fs*1, horizontalalignment='center')
    plt.tick_params(labelsize=fs*0.8)
    plt.xlim([-10,22])
    plt.ylabel('mean root to tip distance relative to start of therapy', fontsize=fs*0.9)
    plt.xlabel('time relative to start of therapy [years]', fontsize=fs*0.9)
    if fname is not None:
        for fmt in formats:
            plt.savefig(fname+fmt)


def evolutionary_rates(root_to_tip, latex=True):
    from scipy.stats import linregress
    regressions = {}
    avg_func = lambda vals:np.sum(vals[:,0]*vals[:,1])/np.sum(vals[:,1])
    #avg_func = lambda vals:weighted_percentile(vals[:,0], vals[:,1], 50)
    for pcode in pcodes:
        RNA = np.array([ (tp, avg_func(vals))
                         for tp, vals in root_to_tip[pcode]['RNA']])
        DNA = np.array([RNA[-1]] + [ (tp, avg_func(vals))
                            for tp, vals in root_to_tip[pcode]['DNA']])

        RNA_reg = linregress(RNA[:,0], RNA[:,1]).__dict__      # regression of RNA only
        DNA_reg = linregress(DNA[:,0], DNA[:,1]).__dict__  # regression of last RNA time points + DNA
        regressions[pcode] = {'RNA':RNA_reg, 'DNA':DNA_reg}

    if latex:
        sep = '\t&\t'
        lend = '\\\\\n'
        suffix = "tex"
    else:
        sep='\t'
        lend='\n'
        suffix = "tsv"
    with open('data/evolutionary_rates.'+suffix, 'w') as ofile:
        for pcode in pcodes:
            reg = regressions[pcode]
            outstr = sep.join([pcode]+map(lambda x:format(x,'1.1e'), [reg['RNA']['slope'],
                  reg['RNA']['pvalue'], reg['DNA']['slope'], reg['DNA']['pvalue']]))
            print(outstr)
            ofile.write(outstr+lend)

def fit_decay(t,y):
    def cost(tau,t,y):
        tmp_y = np.exp(-t/tau)
        tmp_y = (tmp_y[1:] + tmp_y[:-1])*0.5
        return np.sum((tmp_y/tmp_y.sum()-y)**2)
    from scipy.optimize import minimize_scalar
    sol = minimize_scalar(cost, bracket = [0,100], args=(t,y))
    return sol['x']


def probable_origin_figure(patients, origin, fname1=None, fname2=None, fname3=None):
    '''
    plots the distribution of the assignment of DNA reads to RNA samples
    This includes figures for individual patients as well as histograms
    combining the information from all patients.
    '''
    from matplotlib import cm
    dna_cmap = cm.hot

    def plot_pat(pcode, pid, ax, xlabel='', ylabel='', legend_fs=8, ymax=1.0):
        '''
        make the individuals panels with the patient distribution
        '''
        dsi_treatment = ((treatment_dates.loc[pid-1,'treatment start'] - treatment_dates.loc[pid-1,'EDI']).days)
        p=patients[pcode]
        all_points = []
        for si, (col, outprefix) in enumerate(zip(cols, sorted(origin[pcode].keys()))):
            samp_ii = origin[pcode][outprefix][:,0]   # index of the sample that is closest to this read
            weights = origin[pcode][outprefix][:,1]   # the fraction of reads or the fraction of total distinct sequences
            y,x = np.histogram(samp_ii, weights =weights,
                               bins = np.arange(len(p.ysi)+1)-0.5, normed=True)
            ax.plot(p.ysi, y,'-o', lw=3, label="DNA "+str(si+1), c=col)

            ax.legend(loc=6, fontsize=legend_fs)
            all_points.extend([(dsi_treatment-p.dsi[ii],p.dsi[ii], w)
                                for ii,w in zip(samp_ii, weights)])

        ax.text(0.1, 0.9, pcode, transform=ax.transAxes, fontsize=16)
        ax.set_xlabel(xlabel, fontsize=fs)
        ax.set_ylabel(ylabel, fontsize=fs)
        ax.set_ylim([0,ymax])
        ax.tick_params(labelsize=fs*0.8)
        return all_points

    ####
    ## Big figure with individual panels for all patients
    ####
    fig, axs = plt.subplots(2,5, figsize=(16,8)) # make large figures for distributio of most likely sample
    combined_seeding = []
    cols = sns.color_palette(n_colors=6)[1::2]
    # make a figure of the distribution of the distance distribution
    for pi, (pid, pcode) in enumerate(zip(p_indices, pcodes)): # loop over patients
        ax = axs[pi//5,pi%5]
        combined_seeding.extend( plot_pat(pcode, pid, ax, xlabel = 'time since EDI [years]' if pi>=5 else '',
                                                    ylabel = 'distribution of reads' if pi%5==0 else ''))

    if fname1 is not None:
        for fmt in formats:
            plt.savefig(fname1+fmt)

    #make figure with distribution of age of reads and example panels for p1 and p8
    if len(combined_seeding):
        # two panel example figure with the same patients used in trees
        fig1, axs1 = plt.subplots(1,2)
        plot_pat('p1', 1, axs1[0], xlabel = 'time since EDI [years]', ylabel = 'distribution of reads', legend_fs=14,ymax=0.7)
        add_panel_label(axs1[0],'A', x_offset=-0.30, fs_factor=2)
        plot_pat('p8', 8, axs1[1], xlabel = 'time since EDI [years]', ylabel = '', legend_fs=14, ymax=0.7)
        add_panel_label(axs1[1],'', x_offset=-0.3, fs_factor=2)
        axs1[0].tick_params(labelsize=fs*0.8)
        axs1[1].tick_params(labelsize=fs*0.8)
        plt.tight_layout()
        if fname3 is not None:
            for fmt in formats:
                plt.savefig(fname3+fmt)


        fig, axs = plt.subplots(1,2, sharey=True, gridspec_kw = {'width_ratios':[2, 5]})
        tmp = np.array(combined_seeding) #this is [(time before treatment, time since infections, weight),(),...]
        total = tmp[:,-1].sum()
        tbins = np.linspace(0,5,6)
        y,x = np.histogram(np.maximum(tmp[:,0],0)/365.25, weights=tmp[:,-1], bins=tbins)

        # plot histogram
        axs[1].bar(x[:-1], y/total, width=x[1]-x[0], label='DNA read origin')
        # add exponential decay
        tau = fit_decay(tbins, y/total)
        print("optimized decay", tau)
        tmp_y = np.exp(-tbins/tau)
        tmp_y/=tmp_y.sum()
        axs[1].plot(tbins, tmp_y, c='r',lw=4,
                    label=r'$\sim 2^{-t/t_{1/2}}$ with $t_{1/2}=%1.2f$ years'%(tau*np.log(2)))
        # tune figure
        #axs[1].set_ylabel('distribution', fontsize=fs)
        axs[1].set_xlabel('years before start of therapy', fontsize=fs)
        axs[1].tick_params(labelsize=fs*0.8)
        axs[1].set_xlim(5,0)
        #axs[1].get_yaxis().set_visible(False)
        axs[1].legend(loc=2, fontsize=fs*0.8)
        #add_panel_label(axs[0],'A', x_offset=-0.27)

        y,x = np.histogram(tmp[:,1]/365.25, weights=tmp[:,-1], bins=np.linspace(0,2,11))
        axs[0].bar(x[:-1], y/total, width=x[1]-x[0])

        axs[0].set_ylabel('distribution', fontsize=fs)
        axs[0].set_xlabel('years since EDI', fontsize=fs)
        axs[0].tick_params(labelsize=fs*0.8)
        add_panel_label(axs[0],'B', x_offset=-0.5, fs_factor=2)
        plt.tight_layout()

        if fname2 is not None:
            for fmt in formats:
                plt.savefig(fname2+fmt)

        ## calculate expected contribution in first year
        expected_early = {}
        for pi, (pid, pcode) in enumerate(zip(p_indices, pcodes)): # loop over patients
            dsi_treatment = ((treatment_dates.loc[pid-1,'treatment start'] - treatment_dates.loc[pid-1,'EDI']).days)
            year_treatment_start = dsi_treatment/365.25
            expected_early[pcode] = np.exp(-(year_treatment_start-1)/tau) - np.exp(-(year_treatment_start)/tau)
            print(pcode, 'expected early contribution:', expected_early[pcode])
        print('average expected early contribution: %1.3f'%np.mean(expected_early.values()))


def study_design():
    from matplotlib import cm
    plt.figure(figsize=(8,3))
    ax=plt.subplot(111)
    sns.set_style('white')
    ax.set_axis_off()
    spacing = 1.0
    dna_colors = sns.color_palette(n_colors=6)
    dna_colors = [dna_colors[i] for i in [1,3,5]]
    rna_cmap = cm.hot
    for x in range(-15,21,5):
        ax.plot([x,x],[-0.6,-0.2], lw=1, c=(0.4, 0.4, 0.4), alpha=0.5, zorder=0)
    ax.plot([0,0],[-0.5,len(pcodes)-0.5], lw=3, c=(0.4, 0.4, 0.4), alpha=0.5, zorder=0)
    for pi,pcode in enumerate(pcodes[::-1]):
        ax.plot([-16.5, 20], spacing*pi*np.ones(2), lw=1, c=(0.9, 0.9, 0.9), zorder=0)
        dsi = patients[pcode].dsi
        dsi_color = 0.95 - 0.9*dsi/dsi.max()
        ax.text(-18, pi*spacing, pcode, verticalalignment='center', zorder=2)
        dsi_treatment = get_treatment_start_date(pcode)
        ax.plot([-dsi_treatment/365.25, 0], spacing*pi*np.ones(2), lw=10, c='k', alpha=0.3, zorder=1)
        ax.scatter(np.minimum(0,(dsi - dsi_treatment)/365.25), spacing*pi*np.ones(len(dsi)),
                    marker='o', s=50, c=dsi_color, cmap=rna_cmap, zorder=2, vmin=0, vmax=1.0)
        #print((dsi - dsi_treatment)/365.25)
        for dnac, m, outprefix in zip(dna_colors, ('v', '^', 's'), patient_to_prefix_p17[pcode]):
            ii = np.where(treatment_dates.code==pcode)[0][0]
            tp_rel_treatment = np.round(((pd.Timestamp(outprefix[-10:])-treatment_dates.loc[ii,'treatment start']).days)/365.25,1)
            ax.scatter([tp_rel_treatment], [spacing*pi], marker=m, s=50, alpha=0.8, c=dnac, zorder=2)
    plt.tight_layout()
    plt.savefig('figures/study_design.svg')



if __name__=="__main__":
    plt.ion()
    patients = {pcode:Patient.load(pcode) for pcode in pcodes}
    import pandas as pd
#    treatment_dates = pd.read_excel('data/2016-01-08_treatment_start_dates.xlsx')
    treatment_dates = pd.read_excel('data/2016-10-10_treatment_start_dates.xlsx')
    sns.set_style('darkgrid')
    for count_haplotypes in [True, False]:
        for mt in ['clustered_good', 'hyper']:
            ch = '_hap_count' if count_haplotypes else '_read_count'
            all_rtt, origin = root_to_tip_distances(patients, mt, count_haplotypes=count_haplotypes)
            # make figures showing root to tip regression. show Lorenzo-Redondo expectation only in supplement
            combined_root_to_tip_figures(all_rtt, 'figures/combined_root_to_tip_'+mt+ch, mt, show_LR_rate=count_haplotypes, q="average")

            # make figure with histgrams of seeding times and the distribution of seeding times for each patient.
            probable_origin_figure(patients, origin, fname1='figures/seeding_all_'+mt+ch,
                fname2='figures/seeding_combined_'+mt+ch, fname3='figures/seeding_p1p8_'+mt+ch)
            if (count_haplotypes==False and mt=="clustered_good"):
                evolutionary_rates(all_rtt, latex=True)
                evolutionary_rates(all_rtt, latex=False)
