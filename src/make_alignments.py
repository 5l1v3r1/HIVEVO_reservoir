from __future__ import print_function
import numpy as np
import glob, os, socket, gzip
import sys; sys.path.append('../HIVEVO_access'); sys.path.append('src')
from hivevo.sequence import alpha
from util import *
from Bio.Align import MultipleSeqAlignment
import cPickle as pickle

def align(seqs, aln_fname):
    outfname = 'temp_'+str(np.random.randint(10000000))+'.fasta'
    SeqIO.write(seqs,outfname, 'fasta')
    os.system("mafft "+outfname+"  > " + aln_fname)
    os.remove(outfname)

def filter_func(seq_rec, region):
    if seq_rec.id[:4]!="days":
        nNs = str(seq_rec.seq).count('N')
        nstops = str(seq_rec.seq[2:].translate()).count('*')
        return (nNs<1) #&(nstops==0)
    else:
        return True

def add_RNA_haplotypes(fname, pat, region):
    '''
    loads patient RNA haplotypes and DNA haplotypes and return
    both in one list
    '''
    from hivevo.patients import Patient
    try:
        if pat in patient_translation:
            pcode = patient_translation[pat]
        else:
            pcode = pat
        p = Patient.load(pcode)
        rna_haps = p.get_haplotype_alignment(region)
        with myopen(fname) as ifile:
            dna_haps = [seq for seq in SeqIO.parse(ifile, 'fasta')]
        dna_haps.extend(filter(lambda x:get_RNA_read_count(x.description)>2, rna_haps))
        return dna_haps
    except:
        print(fname,pat, region,"failed")

def find_hypermutants(aln, thres=-2):
    '''
    custom routine to find hypermutated sequences, it splits the alignment into sub alignments
    RNA, Good sequences, hyper mutated, suspicious. In addition, it returns a subset of sequences
    that translate without stop codon (assuming the p17 amplicon)
    '''
    isRNA = np.array([True if seq.id[:4]=="days" else False for seq in aln], dtype=bool)
    RNAaln = MultipleSeqAlignment([aln[i].upper() for i in np.where(isRNA)[0]])
    DNAaln = MultipleSeqAlignment([aln[i].upper() for i in np.where(~isRNA)[0]])

    # load the RNA SNP freuqencies to determine positions variable at the RNA level
    # those are disregarded for the hypermutation classification
    RNAaf = np.zeros((len(alpha), RNAaln.get_alignment_length()))
    for seq in RNAaln:
        nucs = np.fromstring(str(seq.seq).upper(), 'S1')
        freq = float(seq.description.split('frequency_')[1].split('%')[0])*0.01
        for ni,nuc in enumerate(alpha):
            RNAaf[ni, nucs==nuc]+=freq
    RNAaf/=RNAaf.sum(axis=0)

    # if the maximal allele frequency is above 0.99, positions are considered conserved
    conserved_pos = RNAaf[:4].max(axis=0)>0.99
    consensus = np.array([alpha[ai] for ai in RNAaf.argmax(axis=0)])
    mut_hist = {'good':[], 'hyper':[], 'suspicious':[]}
    DNAaln_array = np.array(DNAaln)
    good_seqs = []
    hyper_muts = []
    suspicious = []
    nostop = []
    mut_dict = {}
    ii=0
    for a in alpha:
        for b in alpha:
            if a!=b:
                mut_dict[a+'->'+b] = ii
                ii+=1

    for si,seq in enumerate(DNAaln):
        muts = (consensus!=DNAaln[si])&conserved_pos&(DNAaln_array[si]!='-')
        tmp = defaultdict(int)
        #print(seq.name, np.where(muts)[0])
        total = muts.sum()
        mut_counts = np.zeros(30)
        for mi in np.where(muts)[0]:
            mut = consensus[mi]+'->'+DNAaln[si,mi]
            tmp[mut]+=1
            mut_counts[mut_dict[mut]]+=1
        if total<10 and (total<4 or tmp['G->A']<0.5*total):
            good_seqs.append(seq)
            mut_hist['good'].append(mut_counts)
        elif tmp['G->A']>=0.5*total:
            hyper_muts.append(seq)
            mut_hist['hyper'].append(mut_counts)
        else:
            suspicious.append(seq)
            mut_hist['suspicious'].append(mut_counts)
        if total<20 and seq.seq.ungap('-')[20:].translate().count('*')==0:
            nostop.append(seq)
    for k in mut_hist:
        mut_hist[k] = np.array(mut_hist[k])
    return RNAaln, MultipleSeqAlignment(good_seqs), MultipleSeqAlignment(hyper_muts),\
            MultipleSeqAlignment(suspicious), MultipleSeqAlignment(nostop), mut_hist


def cluster_haplotypes(aln, abundance_thres = 5, mut_thres = 1, read_count_func = get_DNA_read_count):
    '''
     - sort haplotypes by abundance
     - compare rare ones to common ones and merge based on threshold
     - return a new alignment
    '''
    # make a list of haplotypes with supplementary info
    hap_list = sorted([ {'n': read_count_func(seq.name), 'seq':seq, 'children':[],
                         'clustered':False, 'seq_array':np.array(seq.seq)} for seq in aln],
                         key=lambda x:x['n'])
    hap_list = [x for x in hap_list if x['seq'].name[:4]!='days'] # exclude RNA haplotypes
    hap_list_rev = [x for x in hap_list[::-1]]

    for hap in hap_list:
        if hap['n']>abundance_thres:
            break
        for hap_parent in hap_list_rev:
            if hap_parent['clustered'] or hap_parent['n']<abundance_thres:
                break
            mut_dist = np.sum((hap_parent['seq_array']!=hap['seq_array'])&(hap['seq_array']!='N'))
            if mut_dist<=mut_thres:
                print("clustering", hap['seq'].name, mut_dist)
                hap['clustered']=True
                hap_parent['children'].append(hap)
                break

    new_align = []
    for hi,hap in enumerate(hap_list_rev):
        if not hap['clustered']:
            n = hap['n'] + sum([x['n'] for x in hap['children']])
            new_name = "_".join(hap['seq'].name.split('_')[:-1])+'_'+str(n)+'_'+str(hi)
            new_align.append(SeqRecord.SeqRecord(seq=hap['seq'].seq, name=new_name, id=new_name, description=new_name))

    return MultipleSeqAlignment(new_align)


def make_joint_aln(fname):
    '''
    produce different joint alignments of RNA and DNA and make a tree.
    '''
    outprefix, pat, region, sdate = get_outprefix(fname)
    all_haps = add_RNA_haplotypes(fname, pat, region)
    if all_haps is not None and len(all_haps)>5:
        all_haps = [x for x in all_haps if filter_func(x, region)]
        for hi, hap in enumerate(all_haps):
            hap.id+='_'+str(hi)
            hap.name=hap.id
        align(all_haps, 'data/'+outprefix+'_RNA_and_DNA.fasta')
        with open('data/'+outprefix+'_RNA_and_DNA.fasta') as infile:
            aln = AlignIO.read(infile, 'fasta')

        #DNAaln = fix_hypermutations(aln, thres=2, logfname='data/'+outprefix+'_hypermut.log')
        RNAaln, DNAaln, hyper, susp, nostop, mut_hist = find_hypermutants(aln)
        total_DNA_reads = np.sum([get_DNA_read_count(x.name) for x in DNAaln])
        with open('data/'+outprefix+'_mutation_statistics.pkl', 'w') as ofile:
            pickle.dump(mut_hist, ofile)

        RNA_no_gaps = ungap(RNAaln)
        for aln, label in  [(DNAaln, 'good'), (hyper,'hyper'), (susp,'suspicious'), (nostop,'nostop')]:
            # raw unclustered alignment
            tmp_ungapped = ungap(aln)
            align(tmp_ungapped, 'data/'+outprefix+'_DNA_'+label+'.fasta')
            align(RNA_no_gaps + tmp_ungapped, 'data/'+outprefix+'_RNA_and_DNA_'+label+'.fasta')

            # cluster sequences minor to major with one mutation distance threshold
            # minimal frequency is 1/500
            DNA_aln_clustered = cluster_haplotypes(aln, abundance_thres=max(total_DNA_reads//500, 5), mut_thres=1)
            tmp_ungapped = ungap(DNA_aln_clustered)
            align(tmp_ungapped, 'data/'+outprefix+'_DNA_clustered_'+label+'.fasta')
            align(RNA_no_gaps + tmp_ungapped, 'data/'+outprefix+'_RNA_and_DNA_clustered_'+label+'.fasta')


def count_RNA_haplotypes():
    all_read_counts = []
    for pcode in patient_to_prefix_p17:
        fname = 'data/'+patient_to_prefix_p17[pcode][0]+'_RNA_and_DNA.fasta.gz'
        with myopen(fname) as infile:
            aln = AlignIO.read(infile, 'fasta')
        reads_per_time_point = defaultdict(list)
        for read in aln:
            if read.name[:4]=='days':
                dsi = int(read.name.split('_')[1])
                reads_per_time_point[dsi].append(get_RNA_read_count(read.description))

        for dsi in sorted(reads_per_time_point.keys()):
            print(pcode, dsi, len(reads_per_time_point[dsi]), sum(reads_per_time_point[dsi]))
            all_read_counts.append([dsi, len(reads_per_time_point[dsi]), sum(reads_per_time_point[dsi])])

    all_read_counts = np.array(all_read_counts)
    print(np.percentile(all_read_counts[:,2], 25), np.percentile(all_read_counts[:,2], 50), np.percentile(all_read_counts[:,2], 75))


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description='analyze p17 amplicon reads from DNA samples')
    parser.add_argument('--sample', type=int, default=-1)
    params = parser.parse_args()

    si = params.sample
    files = glob.glob('data/merged_reads/*aln.fasta.gz')
    if params.sample==-1:
        run_samples=files
    else:
        if si<len(files):
            run_samples = [files[params.sample]]

    for fname in run_samples:
        make_joint_aln(fname)

    count_RNA_haplotypes()
