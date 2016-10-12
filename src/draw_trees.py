from Bio import Phylo
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from util import *

plt.ion()

def draw_tree_combined(tree, markers = ['v', '^', 's', 'd', '<' ,'>'], ax=None, scalebar=True, legend=False):
    '''
    draw a tree and annotate RNA and DNA samples by time relative to infection/therapy start
    '''
    from matplotlib import cm
    from util import HIVEVO_colormap
    dna_colors = [(0, 255, 255),(255, 0, 255),(228,26,28)]
    #dna_cmap = HIVEVO_colormap('alternative')
    #rna_cmap = HIVEVO_colormap('website')
    dna_cmap = cm.Blues
    rna_cmap = cm.hot
    dna_colors = sns.color_palette(n_colors=6)
    dna_colors = [dna_colors[i] for i in [1,3,5]]
    DNA_total = defaultdict(int)
    dsi_max = np.max([int(node.name.split('_')[1]) for node in tree.get_terminals() if node.name[:4]=='days'])
    for node in tree.get_terminals():
        if node.name[:4]!='days': # only DNA reads
            DNA_total[name_to_sample(node.name)]+=get_DNA_read_count(node.name)

    DNA_times = sorted(DNA_total.keys(), key=lambda x:x[-1])
    depths = tree.depths()
    data_circles = defaultdict(list)
    rmin = 10
    rmax = 80
    rfun = lambda hf: rmin + (rmax - rmin) * (hf )**(0.5)
    for ni, node in enumerate(tree.get_terminals()):
        if node.name[:4]=='days':
            node.DSI = float(node.name.split('_')[1])
            # map frequencies and time point to color and size of symbol
            node.freq = get_RNA_read_fraction(node.name)
            node.color = map(int, np.array(rna_cmap(0.95 - 0.9*node.DSI/dsi_max)[:-1]) * 255)
            r = rfun(node.freq)
            y = ni+1
            x = depths[node]
            c = [tone / 255.0  for tone in node.color.to_rgb()]
            cs = [tone / 255.0 * 0.7 for tone in node.color.to_rgb()]
            data_circles['RNA'].append((x, y, 2 * r, c, cs))
        else:
            outprefix, pat, region, sdate = get_outprefix(node.name)
            sample = name_to_sample(node.name)
            # map sample index to color (there isn't enough dynamic range for coloring by time)
            #node.color = map(int, np.array(dna_cmap( (DNA_times.index(sample)+1.0)/(len(DNA_times)+1))[:-1]) * 255)
            node.color = map(int, np.array(dna_colors[DNA_times.index(sample)])* 255)
            #node.color = map(int, dna_colors[DNA_times.index(sample)])
            node.freq = float(get_DNA_read_count(node.name))/DNA_total[sample]
            r = rfun(node.freq)*1.8
            y = ni+1
            x = depths[node]
            c = [tone / 255.0 for tone in node.color.to_rgb()]
            cs = [tone / 255.0 * 0.7 for tone in node.color.to_rgb()]
            data_circles[sample].append((x, y, 2 * r, c, cs))
        # eventually set node color to grey (this determines branch color too and we draw tip symbols separately)
        node.color = (160, 160, 160)
    for node in tree.get_nonterminals():
        node.color = (160, 160, 160)


    if ax is None:
        fig = plt.figure(figsize=(6,6))
        ax=plt.subplot(111)
    Phylo.draw(tree,  show_confidence=False, label_func=lambda x: '', axes=ax, do_show=False)
    # add tip symbols as a scatter plot with positions, size, stroke color and fill color
    (x, y, s, c,cs) = zip(*data_circles['RNA'])
    ax.scatter(x, y, s=s, c=c, edgecolor=cs, zorder=2, label='RNA' if legend else None)
    for di, dna_sample in enumerate(DNA_times): # repeat for DNA samples
        (x, y, s, c,cs) = zip(*data_circles[dna_sample])
        ax.scatter(x, y, s=s, c=c, edgecolor=cs, zorder=3+di, marker=markers[DNA_times.index(dna_sample)],
                   alpha=0.9, label='DNA '+str(di+1) if legend else None)
    nleaves = len(tree.get_terminals())
    if scalebar:
        ax.plot([0.005,0.015], 0.95*np.ones(2)*nleaves,lw=3, c='k')
        ax.text(0.01, 0.92*nleaves, '0.01', fontsize=fs*0.8,
                horizontalalignment='center', verticalalignment='bottom')
    if legend:
        ax.legend(fontsize=fs)
    ax.set_ylim([1.02*nleaves, -0.02*nleaves])
    ax.set_axis_off()
    plt.tight_layout()

def draw_patient_RNA_DNA_tree(pcode):
    for seq_type in ['clustered_good','hyper']:
        tree = Phylo.read('data/'+pcode+ '_RNA_and_DNA_'+seq_type+'.nwk', 'newick')
        draw_tree_combined(tree)
        #plt.title(pcode+ '_RNA_and_DNA_'+seq_type)

        if plt.xlim()[1]>0.2:
            plt.xlim((0,0.2))
        for fmt in formats:
            plt.savefig('figures/'+pcode+ '_RNA_and_DNA_'+seq_type+'_tree'+fmt)

def draw_patient_RNA_DNA_tree_website(pcode):
    seq_type = 'clustered_good'
    tree = Phylo.read('data/'+pcode+ '_RNA_and_DNA_'+seq_type+'.nwk', 'newick')
    draw_tree_combined(tree, legend=True)

    if plt.xlim()[1]>0.2:
        plt.xlim((0,0.2))
    for fmt in ['.svg', '.png']:
        plt.savefig('website_data/cell_trees/tree_'+pcode+'_p17'+fmt)


def draw_all_patient_RNA_DNA_tree(seq_type='clustered_good', xmax=0.10):
    fig, axs= plt.subplots(2,5, sharex=True, figsize=(24,12))
    for ii, pcode in enumerate(pcodes):
        ax = axs[ii//5][ii%5]
        tree = Phylo.read('data/'+pcode+ '_RNA_and_DNA_'+seq_type+'.nwk', 'newick')
        draw_tree_combined(tree, ax=ax, scalebar=False, legend=pcode=='p11') # legend in p11 only

        ax.set_xlim((0,xmax))
        ax.text(-0.02, 0.95, pcode, transform=ax.transAxes, fontsize=24, horizontalalignment='right')

    # add scale bar
    ymax = axs[-1][0].get_ylim()[0]
    axs[-1][0].plot([0.0, 0.02], [ymax*0.9, ymax*0.9], c='k', lw=3)
    axs[-1][0].text(0.01, ymax*0.89, '0.02', fontsize=16, horizontalalignment='center', verticalalignment='bottom')
    plt.tight_layout()
    for fmt in formats:
        plt.savefig('figures/all_pat_RNA_and_DNA_'+seq_type+'_tree'+fmt)

def draw_some_patient_RNA_DNA_tree(seq_type='clustered_good', patients=['p1', 'p8'], xmax=0.07):
    fig, axs= plt.subplots(1,len(patients), figsize=(10,6))
    for ii, pcode in enumerate(patients):
        ax = axs[ii]
        tree = Phylo.read('data/'+pcode+ '_RNA_and_DNA_'+seq_type+'.nwk', 'newick')
        draw_tree_combined(tree, ax=ax, scalebar=pcode=='p1', legend=pcode=='p8')
        ax.text(-0.02, 0.95, pcode, transform=ax.transAxes, fontsize=24, horizontalalignment='right')
        ax.set_xlim(0,xmax)

    plt.tight_layout()
    for fmt in formats:
        plt.savefig('figures/'+"_".join(patients)+'_RNA_and_DNA_'+seq_type+'_tree'+fmt)

if __name__=="__main__":
    draw_all_patient_RNA_DNA_tree(seq_type='clustered_good')
    draw_some_patient_RNA_DNA_tree(seq_type='clustered_good')
