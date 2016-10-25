#!/usr/bin/env python
import numpy, os, sys

def read_pssm(pssm_file):
	# this function reads the pssm file given as input, and returns a LEN x 20 matrix of pssm values.

	# index of 'ACDE..' in 'ARNDCQEGHILKMFPSTWYV'(blast order)
	idx_res = (0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18)

	# open the two files, read in their data and then close them
	if pssm_file == 'STDIN': fp = sys.stdin
	else: fp = open(pssm_file, 'r')
	lines = fp.readlines()
	fp.close()

	# declare the empty dictionary with each of the entries
	aa = []
	pssm = []

	# iterate over the pssm file and get the needed information out
	for line in lines:
		split_line = line.split()
		# valid lines should have 32 points of data.
		# any line starting with a # is ignored
		try: int(split_line[0])
		except: continue

		if line[0] == '#': continue

		aa_temp = split_line[1]
		aa.append(aa_temp)
		if len(split_line) in (44,22):
			pssm_temp = [-float(i) for i in split_line[2:22]]
		elif len(line) > 70:  # in case double digits of pssm
			pssm_temp = [-float(line[k*3+9: k*3+12]) for k in range(20)]
			pass
		else: continue
		pssm.append([pssm_temp[k] for k in idx_res])

	return aa, pssm


def get_phys7(aa):
	# this function takes a path to a pssm file and finds the pssm + phys 7 input
	# to the NN in the required order - with the required window size (8).


	# define the dictionary with the phys properties for each AA
	phys_dic = {'A': [-0.350, -0.680, -0.677, -0.171, -0.170, 0.900, -0.476],
							'C': [-0.140, -0.329, -0.359, 0.508, -0.114, -0.652, 0.476],
							'D': [-0.213, -0.417, -0.281, -0.767, -0.900, -0.155, -0.635],
							'E': [-0.230, -0.241, -0.058, -0.696, -0.868, 0.900, -0.582],
							'F': [ 0.363, 0.373, 0.412, 0.646, -0.272, 0.155, 0.318],
							'G': [-0.900, -0.900, -0.900, -0.342, -0.179, -0.900, -0.900],
							'H': [ 0.384, 0.110, 0.138, -0.271, 0.195, -0.031, -0.106],
							'I': [ 0.900, -0.066, -0.009, 0.652, -0.186, 0.155, 0.688],
							'K': [-0.088, 0.066, 0.163, -0.889, 0.727, 0.279, -0.265],
							'L': [ 0.213, -0.066, -0.009, 0.596, -0.186, 0.714, -0.053],
							'M': [ 0.110, 0.066, 0.087, 0.337, -0.262, 0.652, -0.001],
							'N': [-0.213, -0.329, -0.243, -0.674, -0.075, -0.403, -0.529],
							'P': [ 0.247, -0.900, -0.294, 0.055, -0.010, -0.900, 0.106],
							'Q': [-0.230, -0.110, -0.020, -0.464, -0.276, 0.528, -0.371],
							'R': [ 0.105, 0.373, 0.466, -0.900, 0.900, 0.528, -0.371],
							'S': [-0.337, -0.637, -0.544, -0.364, -0.265, -0.466, -0.212],
							'T': [ 0.402, -0.417, -0.321, -0.199, -0.288, -0.403, 0.212],
							'V': [ 0.677, -0.285, -0.232, 0.331, -0.191, -0.031, 0.900],
							'W': [ 0.479, 0.900, 0.900, 0.900, -0.209, 0.279, 0.529],
							'Y': [ 0.363, 0.417, 0.541, 0.188, -0.274, -0.155, 0.476]}



	# set the phys7 data.
	phys = [phys_dic.get(i, phys_dic['A']) for i in aa]

	return phys

def window(feat, winsize=8):
	# apply the windowing to the input feature
	feat = numpy.array(feat)
	output = numpy.concatenate([numpy.vstack([feat[0]]*winsize), feat])
	output = numpy.concatenate([output, numpy.vstack([feat[-1]]*winsize)])
	output = [numpy.ndarray.flatten(output[i:i+2*winsize+1]).T for i in range(0,feat.shape[0])]
	return output

def window_data(*feature_types):
	n = len(feature_types[0])
	features = numpy.empty([n,0])

	for feature_type in feature_types:
		test = numpy.array(window(feature_type))
		features = numpy.concatenate((features, test), axis=1)

	return features


def sigmoid(input):
	# apply the sigmoid function
	output = 1 / (1 + numpy.exp(-input))
	return(output)

def nn_feedforward(nn, input):
	input = numpy.matrix(input)

	# find the number of layers in the NN so that we know how much to iterate over
	num_layers = nn['n'][0][0][0][0]
	# num_input is the number of input AAs, not the dimentionality of the features
	num_input = input.shape[0]
	x = input

	# for each layer up to the final
	for i in range(1,num_layers-1):
		# get the bais and weights out of the nn
		W = nn['W'][0][0][0][i-1].T
		temp_size = x.shape[0]
		b = numpy.ones((temp_size,1))
		x = numpy.concatenate((b, x),axis=1)
		# find the output of this layer (the input to the next)
		x = sigmoid(x * W)

	# for the final layer.
	# note that this layer is done serpately, this is so that if the output nonlinearity
	# is not sigmoid, it can be calculated seperately.
	W = nn['W'][0][0][0][-1].T
	b = numpy.ones((x.shape[0],1))
	x = numpy.concatenate((b, x),axis=1)
	output = x*W
	pred = sigmoid(output)

	return pred

def load_nn(file_name):

	# load the nn file
	import scipy.io
	mat = scipy.io.loadmat(file_name)

	return mat

def load_NN(nn_filename):
	return numpy.load(nn_filename)

	# load in the NN mat file.
#	import scipy.io
#	mat = scipy.io.loadmat(nn_filename)

#	nn = mat['nn']

	return nn

dict_ASA0 = dict(zip("ACDEFGHIKLMNPQRSTVWY",
					(115, 135, 150, 190, 210, 75, 195, 175, 200, 170,
		185, 160, 145, 180, 225, 115, 140, 155, 255, 230)))
def run_iter(dict_nn, input_feature0, aa, ofile):
	SS_order = ('C' 'E' 'H')
	list1 = ('SS', 'ASA', 'TTPP')
	list_res1 = []
	for x in list1:
		nn = dict_nn[x]
		norm_max = nn['high'][0][0][0]
		norm_min = nn['low'][0][0][0]
		input_feature1 = (input_feature0 - numpy.tile(norm_min, (input_feature0.shape[0], 1))) / numpy.tile((norm_max - norm_min), (input_feature0.shape[0],1))
		r1 = nn_feedforward(nn, input_feature1)
		list_res1.append(r1)

	pred_ss_1, pred_asa_1, pred_ttpp_1 = list_res1
	SS_1 = [SS_order[i] for i in numpy.argmax(pred_ss_1,1)]
	pred_ttpp_1_denorm = (pred_ttpp_1 - 0.5) * 2
	theta = numpy.degrees(numpy.arctan2(pred_ttpp_1_denorm[:,0], pred_ttpp_1_denorm[:,2]))
	tau = numpy.degrees(numpy.arctan2(pred_ttpp_1_denorm[:,1], pred_ttpp_1_denorm[:,3]))
	phi = numpy.degrees(numpy.arctan2(pred_ttpp_1_denorm[:,4], pred_ttpp_1_denorm[:,6]))
	psi = numpy.degrees(numpy.arctan2(pred_ttpp_1_denorm[:,5], pred_ttpp_1_denorm[:,7]))

	if ofile == 'NULL':
		return SS_1, pred_ss_1, pred_asa_1, theta, tau, phi, psi

	fp = open(ofile, 'w')
	print >>fp, '#\tAA\tSS\tASA\tPhi\tPsi\tTheta(i-1=>i+1)\tTau(i-2=>i+1)\tP(C)\tP(E)\tP(H)'
	for ind, x in enumerate(aa):
		asa = pred_asa_1[ind] * dict_ASA0.get(x, dict_ASA0['A'])
		#print ind+1, aa[ind], SS_1[ind], pred_ss_1[ind,0], pred_ss_1[ind,1], pred_ss_1[ind,2], pred_asa_1[ind], theta[ind], tau[ind], phi[ind], psi[ind]
		print >>fp, ('%i\t%c\t%c\t%5.1f' + '\t%6.1f'*4 + '\t%.3f'*3) % (ind+1, x, SS_1[ind], asa, phi[ind], psi[ind], theta[ind], tau[ind], pred_ss_1[ind,0], pred_ss_1[ind,1], pred_ss_1[ind,2])
	fp.close()

	return SS_1, pred_ss_1, pred_asa_1, theta, tau, phi, psi


def main(list_params, pssm_file, out_suffix):

	basenm = os.path.basename(pssm_file)
	if basenm.endswith('.pssm'): basenm = basenm[:-5]
	elif basenm.endswith('.mat'): basenm = basenm[:-4]

	outfile0 = '%s.%s' % (basenm, out_suffix)
	if not bforce and os.path.isfile(outfile0+'3'): return
	else: open(outfile0+'3', 'w').close()
#
# this part is needed as the later output is 'appended'
#	if ball:
#		for k in (1,2): open('%s%d' % (outfile0, k), 'w').close()

	SS_order = ('C' 'E' 'H')

	aa, pssm = read_pssm(pssm_file)

	pred1(list_params, aa, pssm, outfile0)

def pred1(list_params, aa, pssm, outfile0):
	list_nn, ball = list_params
	phys = get_phys7(aa)
	input_feature = window_data(pssm, phys)
	## DO FIRST PREDICTIONS
	for it1 in (1, 2, 3):
		ofile = '%s%d' % (outfile0, it1)
		if not ball and it1<3: ofile = 'NULL'

		dict_nn = list_nn[it1-1]
		res1 = run_iter(dict_nn, input_feature, aa, ofile)

		if it1 == 3: break
	
		## feature after 1st iteration
		SS_1, pred_ss_1, pred_asa_1, theta, tau, phi, psi = res1
		tt_input = numpy.sin(numpy.concatenate((numpy.radians(theta), numpy.radians(tau)), axis=1))/2 + 0.5
		tt_input = numpy.concatenate((tt_input, numpy.cos(numpy.concatenate((numpy.radians(theta), numpy.radians(tau)), axis=1))/2 + 0.5), axis=1)
		pp_input = numpy.sin(numpy.concatenate((numpy.radians(phi), numpy.radians(psi)), axis=1))/2 + 0.5
		pp_input = numpy.concatenate((pp_input, numpy.cos(numpy.concatenate((numpy.radians(phi), numpy.radians(psi)), axis=1))/2 + 0.5), axis=1)
		ttpp_input = numpy.concatenate((tt_input, pp_input), axis=1)
		input_feature = window_data(pssm, phys, pred_ss_1, pred_asa_1, ttpp_input)
	return


if __name__ == "__main__":
	# if there is no filename for the features to be written to given as input, don't write it
	if len(sys.argv) < 2:
		print >>sys.stderr, "Usage: RUN *.pssmfile"
		sys.exit()

	ball = '-all' in sys.argv
	bforce = '-f' in sys.argv

	nndir = os.path.dirname(sys.argv[0])
	if os.path.isfile(nndir+'/dat/pp1.npz'): nndir += '/dat/'
	elif os.path.isfile(nndir+'/../dat/pp1.npz'): nndir += '/../dat/'
	dict1_nn = load_NN(nndir+'pp1.npz')
	dict2_nn = load_NN(nndir+'pp2.npz')
	dict3_nn = load_NN(nndir+'pp3.npz')
	list_nn = (dict1_nn, dict2_nn, dict3_nn)

	list1 = sys.argv[1:]
	numpy.random.shuffle(list1)
	for x in list1:
		if x[0] == '-': continue
		main([list_nn,ball], x, 'spd')
