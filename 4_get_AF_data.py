import os
import numpy as np
import pickle
from scipy.spatial import distance

def coord_array_maker(infile):
    
	'''
    reads the .cif file coordinates and returns numpy arrays
    '''
    
	coord_array = []
	pLDDT_array = []
	for D in infile[1:]:
		if D == '#':
			break
		E = D.split()
		aa_num = int(E[8])-1
		X = float(E[10])
		Y = float(E[11])
		Z = float(E[12])
		coord_array.append([aa_num, X, Y, Z])
		pLDDT_array.append(float(E[14]))
	return np.array(coord_array), np.array(pLDDT_array)

def get_contacts(coord_dict_in, dist, gap, error, AF_error):
    
	'''
    calculated the distance between all coordinates, filters appropiate
    contacts and summarizes contacts residue wise
    '''
    
	all_dists = distance.squareform(distance.pdist(coord_dict_in[:,1:], metric='euclidean'))
	all_gaps = coord_dict_in[:,0] - coord_dict_in[:,0][np.newaxis].T
	all_errors = np.amin(np.array(np.meshgrid(error,error)), axis=0)
	cond_ids = np.where((all_dists <= dist) & (all_gaps >= gap) & (all_errors >= AF_error))
	res_array = np.array([coord_dict_in[cond_ids[0],0], coord_dict_in[cond_ids[1],0]]).T
	set_res_array, count = np.unique(res_array, return_counts=True, axis=0)
	#column1: residue A
	#column2: residue B
	#column3: contacts between A and B
	return np.concatenate((set_res_array, count[np.newaxis].T), axis=1).astype(np.int_)

def single_array_maker(inflie):

	'''
    reads helix and beta sheet data, if only a single secondary structure element is present
    '''

	take_h = False
	take_b = False

	for i in inflie:
		if '_struct_conf.beg_label_seq_id' in i:
			start = int(i[-5:])-1
		if '_struct_conf.end_label_seq_id' in i:
			end = int(i[-5:])-1
		if '_struct_conf.conf_type_id' in i:
			if 'HELX' in i:
				take_h = True
			if 'STRN' in i:
				take_b = True
	if take_h == True:
		return np.array([[start, end]]).astype(np.int_), np.empty((0, 2)).astype(np.int_)
	elif take_b == True:
		return np.empty((0, 2)).astype(np.int_), np.array([[start, end]]).astype(np.int_)
	else:
		return np.empty((0, 2)).astype(np.int_), np.empty((0, 2)).astype(np.int_)

def array_maker(infile):

	'''
    reads helix and beta sheet data, if multiple secondary structure elements are present
    '''

	helix_array = []
	beta_array = []
	for D in infile[1:]:
		if D == '#':
			break
		E = D.split()
		if 'HELX' in E[6]:
			helix_array.append([int(E[5])-1, int(E[12])-1])
		elif 'STRN' in E[6]:
			beta_array.append([int(E[5])-1, int(E[12])-1])
	return np.array(helix_array).astype(np.int_), np.array(beta_array).astype(np.int_)

#load amino acid sequence of MG1655 genes to compare gene_names and gene_length
path = "/Users/.../"
with open(path + ".../aa_seq_MG1655.pick", "rb") as handle:
	aa_seq = pickle.load(handle)

###start###
master = {}
for gene_count, file in enumerate(os.listdir(path + '0_UP000000625_83333_ECOLI_v4/')):
	if file.endswith(".cif"):
		inflie = open(path + '0_UP000000625_83333_ECOLI_v4/' + file).read()

		#get gene name
		gene_A = inflie.index('_ma_target_ref_db_details.gene_name')
		gene_B = inflie[gene_A:].split('\n')[0]
		gene = gene_B.split(' ')[-1]

		#get gene length
		gene_length_A = inflie.index('_ma_target_ref_db_details.seq_db_align_end')
		gene_length_B = inflie[gene_length_A:].split('\n')[0]
		gene_length = int(gene_length_B.split(' ')[-1])

		if gene not in aa_seq:
			continue
		if len(aa_seq[gene]) != gene_length:
			continue

		#get contacts
		coord_A = inflie.index('_atom_site.pdbx_sifts_xref_db_res')
		coord_B = inflie[coord_A:].split('\n')
		coord_array, pLDDT = coord_array_maker(coord_B)
		contacts = get_contacts(coord_array, 5, 6, pLDDT, 70)

		#get helices and beta sheets
		if '_struct_conf.pdbx_end_PDB_ins_code' in inflie:
			helix_A = inflie.index('_struct_conf.pdbx_end_PDB_ins_code')
			helix_B = inflie[helix_A:].split('\n')
			if helix_B[1] == '#':
				helix_array, beta_array = single_array_maker(inflie[helix_A-4000:helix_A].split('\n'))
			else:
				helix_array, beta_array = array_maker(helix_B)
		else:
			helix_array = np.empty((0, 2))
			beta_array = np.empty((0, 2))

		out = [file[:-4], helix_array, beta_array, contacts]
		master[gene] = out
	
	if gene_count%100 == 0:
		print(gene_count, len(master))
		with open(path + ".../1_AF_data.pick", "wb") as handle:
			pickle.dump(master, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open(path + ".../1_AF_data.pick", "wb") as handle:
	pickle.dump(master, handle, protocol=pickle.HIGHEST_PROTOCOL)
	