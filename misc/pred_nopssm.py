#!/usr/bin/env python
from pred_pssm import *
from prolib import *

def pred_blosum62(seq_file, out_suffix):
	basenm = os.path.basename(seq_file)
	if basenm.endswith('.seq'): basenm = basenm[:-4]
	elif basenm.endswith('.fa'): basenm = basenm[:-3]

	outfile0 = '%s.%s' % (basenm, out_suffix)
	if os.path.isfile(outfile0+'3'): return
	else: open(outfile0+'3', 'w').close()
#
# this part is needed as the later output is 'appended'
	if ball:
		for k in (1,2): open('%s%d' % (outfile0, k), 'w').close()

	for sr1 in Parse_seq(seq_file):
		aa = sr1.seq
		pssm = build_pssm(aa)
		pred1([list_nn, ball], aa, pssm, outfile0)
	return

def build_pssm(seq1):
	list1 = []
	for x in seq1:
		list1.append(dict_blosum62ij[x])
	return list1

def init_blosum62ij():
	blosum62ij_str = """
		-  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
		A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1 -1 -1 -4
		R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1 -2  0 -1 -4
		N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  4 -3  0 -1 -4
		D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4 -3  1 -1 -4
		C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -1 -3 -1 -4
		Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0 -2  4 -1 -4
		E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1 -3  4 -1 -4
		G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -4 -2 -1 -4
		H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0 -3  0 -1 -4
		I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3  3 -3 -1 -4
		L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4  3 -3 -1 -4
		K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0 -3  1 -1 -4
		M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3  2 -1 -1 -4
		F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3  0 -3 -1 -4
		P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -3 -1 -1 -4
		S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0 -2  0 -1 -4
		T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1 -1 -1 -4
		W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -2 -2 -1 -4
		Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -1 -2 -1 -4
		V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3  2 -2 -1 -4
		B -2 -1  4  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4 -3  0 -1 -4
		J -1 -2 -3 -3 -1 -2 -3 -4 -3  3  3 -3  2  0 -3 -2 -1 -2 -1  2 -3  3 -3 -1 -4
		Z -1  0  0  1 -3  4  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -2 -2 -2  0 -3  4 -1 -4
		X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4
		* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
	"""
	dict_blosum62 = {}
	# index of 'ACDE..' in 'ARNDCQEGHILKMFPSTWYV'(blast order)
	idx_res = (0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18)
	for x in blosum62ij_str.split('\n'):
		ss = x.split()
		if len(ss)<1 or ss[0]=='-': continue
		dict_blosum62[ss[0]] = [-float(ss[k+1]) for k in idx_res]
	return dict_blosum62
#
#
if __name__ == "__main__":
	# if there is no filename for the features to be written to given as input, don't write it
	if len(sys.argv) < 2:
		print >>sys.stderr, "Usage: RUN *.seq"
		sys.exit()

	print >>sys.stderr, '''Warning!! This program won't run psiblast and thus has much lower accuracy. It may be good for large-scale proteomic searching because of its super fast calculation'''

	global  ball, dict1_nn, dict2_nn, dict3_nn
	ball = '-all' in sys.argv

	nndir = os.path.dirname(sys.argv[0])
	if os.path.isfile(nndir+'/dat/pp1.npz'): nndir += '/dat/'
	elif os.path.isfile(nndir+'/../dat/pp1.npz'): nndir += '/../dat/'
	dict1_nn = load_NN(nndir+'pp1.npz')
	dict2_nn = load_NN(nndir+'pp2.npz')
	dict3_nn = load_NN(nndir+'pp3.npz')
	list_nn = (dict1_nn, dict2_nn, dict3_nn)
#
	dict_blosum62ij = init_blosum62ij()

	for x in sys.argv[1:]:
		if x[0] == '-': continue
		pred_blosum62(x, 'sp0')
