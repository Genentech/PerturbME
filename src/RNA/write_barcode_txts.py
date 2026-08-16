import os
from pdb import set_trace as bp
import numpy as np
import h5py
import scipy.sparse as sp_sparse
import glob

# Parse h5 files
def parse_h5(fp):
	f = h5py.File(fp,'r')
	barcodes = f['matrix']['barcodes'][:].astype(str)
	genes = f['matrix']['features']['name'][:].astype(str)
	data = f['matrix']['data']
	indices = f['matrix']['indices']
	indptr = f['matrix']['indptr']
	shape = f['matrix']['shape']
	X = sp_sparse.csr_matrix((data[:],indices[:],indptr[:]), shape=(shape[1], shape[0]))

	return X,barcodes,genes

h5_path = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/terra/terra_h5s'
barcodes_path = '/gstore/data/ctgbioinfo/thakorep/PerturbME_transfer/PerturbCITE_ICR/202008_full_exp/terra/Barcodes'
barcode_length = 16

FP = glob.glob(os.path.join(h5_path, '*.h5'))
for fp in FP:
	X,B,G = parse_h5(fp)
	f = open('%s/%s' % (barcodes_path, fp.replace('.h5','.txt').split('/')[-1]),'w')
	for b in B:
		_=f.write('%s\n' % b[:barcode_length])
	f.close()
