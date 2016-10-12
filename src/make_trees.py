import numpy as np
import os
from make_alignments import align
from Bio import AlignIO, Phylo
import sys; sys.path.append('../HIVEVO_access'); sys.path.append('src')
from hivevo.patients import Patient
from util import *

fasttree_bin = 'fasttree_cluster'
save_as='.fasta.gz'


def infer_tree(fname, min_DNA_frac=0.002):
    outfname = 'temp_'+str(np.random.randint(10000000))+'.nwk'
    outfname_aln = 'temp_'+str(np.random.randint(10000000))+'.fasta'
    with myopen(fname) as aln_file:
        aln = AlignIO.read(aln_file,"fasta")
    aln = prune_rare_DNA(aln, min_DNA_frac)
    AlignIO.write(aln, outfname_aln, 'fasta')
    os.system(fasttree_bin+" -nt "+outfname_aln+"> "+ outfname)
    T = Phylo.read(outfname, 'newick')
    os.remove(outfname)
    os.remove(outfname_aln)
    return T

def make_patient_RNA_DNA_tree(pcode, min_DNA_frac = 0.001):
    ''' make a tree for all RNA/DNA sample of a given patient '''
    for seq_type in ['clustered_good', 'good', 'hyper', 'suspicious']:
        seqs=[]
        for outprefix in patient_to_prefix_p17[pcode]:
            with myopen('data/'+outprefix+'_DNA_'+seq_type+save_as) as ifile:
                seqs.extend([x for x in SeqIO.parse(ifile, 'fasta')])
        p = Patient.load(pcode)
        seqs.extend(p.get_haplotype_alignment(region))
        seqs_pruned = prune_rare_DNA(seqs, min_DNA_frac)
        for hi, hap in enumerate(seqs_pruned):
            hap.id+='_'+str(hi)
            hap.name=hap.id

        outfname = 'data/'+pcode+'_RNA_and_DNA_'+seq_type+'.fasta'
        align(ungap(seqs_pruned), outfname)
        tree = infer_tree(outfname, min_DNA_frac=0.0)
        leaves = sorted(filter(lambda x:x.name[:4]=='days', tree.get_terminals()),
                        key = lambda x:(int(x.name.split('_')[1]), -int(x.name.split('_')[3][:-1])))
        tree.root_with_outgroup(leaves[0])
        tree.root.branch_length=0.01
        for branch in tree.get_nonterminals(order='postorder'):
            if branch.branch_length<0.001:
                tree.collapse(branch)
        tree.ladderize()
        Phylo.write(tree, 'data/'+pcode+ '_RNA_and_DNA_'+seq_type+'.nwk', 'newick')

def make_RNA_DNA_tree(outprefix, min_DNA_frac = 0.001):
    '''
    make a tree for a particular sample as specified by the outprefix
    '''
    for seq_type in ['clustered_good', 'good', 'hyper', 'suspicious']:
        tree = infer_tree("data/"+outprefix+"_RNA_and_DNA_"+seq_type+".fasta", min_DNA_frac=min_DNA_frac)
        leaves = sorted(filter(lambda x:x.name[:4]=='days', tree.get_terminals()),
                        key = lambda x:(int(x.name.split('_')[1]), -int(x.name.split('_')[3][:-1])))
        tree.root_with_outgroup(leaves[0])
        tree.root.branch_length=0.01
        for branch in tree.get_nonterminals(order='postorder'):
            if branch.branch_length<0.001:
                tree.collapse(branch)
        tree.ladderize()
        Phylo.write(tree, 'data/'+outprefix+ '_RNA_and_DNA_'+seq_type+'.nwk', 'newick')


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description='analyze p17 amplicon reads from DNA samples')
    parser.add_argument('--patient', type=str, default='p1')
    params = parser.parse_args()

    make_patient_RNA_DNA_tree(params.patient, min_DNA_frac=0.002)
