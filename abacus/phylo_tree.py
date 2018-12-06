from dist_matrix import rand_actg_str
from dist_matrix import pd_matrix
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.Consensus import *
import os

#Biopython has to take triangular matrix, input a list.
def convert_triangle(M):
	k = len(M)
	ret = []
	for i in range(0,k):
		ret.append(list(M[i][0:(i+1)]))
	return ret

def read_fasta(filename):
	with open('sequences/' + filename, 'r') as f:
		name = f.readline().replace('\r','').replace('\n', '')
		s = f.read().replace('\r','').replace('\n','')
	return name,s

def iter_over_files():
	files = [filename for filename in os.listdir('sequences') if filename.endswith('.txt')]
	descriptions = []
	fasta = []
	for filename in files:
		n,s = read_fasta(filename)
		descriptions.append(n)
		fasta.append(s)
	return descriptions,fasta

def get_phylogenetic_tree(max_str_len = 1, norm="JSD", cpc_function="Square25", joining_alg="nj"):
	desc, genes = iter_over_files()
	pm = pd_matrix(genes,max_str_len= max_str_len, norm=norm,cpc_function="Square25")
	pm = convert_triangle(pm)
	dm = DistanceMatrix(names=desc,matrix=pm)
	constructor = DistanceTreeConstructor()
	if(joining_alg == "nj"):
		tree = constructor.nj(dm)
	elif(joining_alg=="upgma"):
		tree = constructor.upgma(dm)
	Phylo.write(tree, 'phylo-tree/result.xml', 'newick')

#Runs the whole thing and outputs tree in xml format to phylo-tree.

get_phylogenetic_tree(max_str_len= 1, norm="Jaccard", cpc_function="Square25", joining_alg="upgma")
