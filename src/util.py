import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import glob

fig_width = 5
fig_fontsize = 12
fs=16
formats = ['.pdf', '.png', '.svg']

p_indices = [1,2,3,5,6,7,8,9,10,11]
pcodes = ['p'+str(pi) for pi in p_indices]

def myopen(fname, mode='r'):
    import gzip
    if fname[-2:]=='gz':
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)

def add_panel_label(ax,label, x_offset=-0.1, fs_factor=1.5):
    ax.text(x_offset, 0.95, label, transform=ax.transAxes, fontsize=fig_fontsize*fs_factor)

def ungap(aln):
    from Bio.SeqRecord import SeqRecord
    return [SeqRecord(seq=x.seq.ungap('-'), id=x.id, name=x.name, description=x.description) for x in aln]

# Functions
def HIVEVO_colormap(kind='website'):
    from scipy.interpolate import interp1d
    maps = {'website': ["#5097BA", "#60AA9E", "#75B681", "#8EBC66", "#AABD52",
                        "#C4B945", "#D9AD3D", "#E59637", "#E67030", "#DF4327"],
            'alternative': ["#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
                            "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A"],
           }
    colors = maps[kind]
    rgb_colors = []
    for c in colors:
        rgb_colors.append([int(c[i:i+2],16) for i in [1,3,5]]+[255])
    tmp =interp1d(np.linspace(0,1,len(colors)), np.array(rgb_colors, dtype = float).T/255.0)
    cmap = lambda x: [c for c in tmp(x)]
    return cmap

def get_DNA_read_count(hap_name):
    return int(hap_name.split('_')[-2])

def get_RNA_read_fraction(hap_name):
    return float(hap_name.split('frequency_')[1].split('%')[0])*0.01

def get_RNA_read_count(hap_name):
    return int(hap_name.split("n.reads:")[-1].strip())

def name_to_sample(n_name):  # first 3 elements are sample unique pat, region, date
    return tuple(n_name.split('_')[:3])


def prune_rare_DNA(aln, thres):
    from collections import defaultdict
    DNA_read_count = defaultdict(int)
    for hap in aln:
        if hap.name[:4]!='days':
            DNA_read_count[name_to_sample(hap.name)] += get_DNA_read_count(hap.name)
    print(DNA_read_count)
    filtered_seqs = [hap for hap in aln if (hap.name[:4]=='days' or
            1.0*get_DNA_read_count(hap.name)/DNA_read_count[name_to_sample(hap.name)]>thres)]
    try:
        from Bio.Align import MultipleSeqAlignment
        return MultipleSeqAlignment(filtered_seqs)
    except:
        return filtered_seqs

def get_outprefix(fname):
    from datetime import datetime
    outprefix = '_'.join(fname.split('/')[-1].split('_')[:3])
    pat, region = outprefix.split('_')[:2]
    sample_date = datetime.strptime(outprefix.split('_')[-1], '%Y-%m-%d').date()
    return outprefix, pat, region, numdate(sample_date)

def name_to_sample(n_name):  # first 3 elements are sample unique pat, region, date
    return tuple(n_name.split('_')[:3])

def numdate(date):
    if type(date)==str:
        from datetime import datetime
        sample_date = datetime.strptime(date, '%Y-%m-%d').date()
    else:
        sample_date = date
    from datetime import date
    ref_day = date(year=sample_date.year, month=1, day=1).toordinal()
    return sample_date.year  + (sample_date.toordinal()-ref_day-1)/362.25

patient_to_prefix_p17 = defaultdict(list)
for sample in glob.glob('data/merged_reads/*fasta.gz'):
    outprefix, pat, region, sdate = get_outprefix(sample)
    patient_to_prefix_p17[pat].append(outprefix)


