import numpy as np
import pickle
import os
import freesasa
radii_class = freesasa.Classifier.getStandardClassifier('naccess')

def coord_array_maker(infile):
    
	'''
    reads the .cif file coordinates and returns numpy arrays
    '''
    
	coord_array = []
	radii_array = []
	aa_position = []
	for D in infile[1:]:
		if D == '#':
			break
		E = D.split()
		aa_num = int(E[8])-1
		X = float(E[10])
		Y = float(E[11])
		Z = float(E[12])
		coord_array.append([X, Y, Z])
		if E[2] == 'S':
			atom = 'O'
		else:
			atom = E[2]
		radii_array.append(radii_class.radius(E[-1], atom))
		aa_position.append(aa_num)
	return np.array(coord_array), np.array(radii_array), np.array(aa_position)


#load amino acid sequence of MG1655 genes to compare gene_names and gene_length
path = "/Users/.../"
with open(path + ".../aa_seq_MG1655.pick", "rb") as handle:
	aa_seq = pickle.load(handle)

#load the domain annotation file made by the "5_structure_assisted_domain_annotation" script
with open(path + "domain_annotation.pick", "rb") as handle:
	domain_annotation = pickle.load(handle)


sasa_dict_full = {}
domain_interface_dict = {}
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
		coord_array, radii_array, pos_array = coord_array_maker(coord_B)


		######calculate sasa for the entire structure######
		protein_sasa = freesasa.calcCoord(np.ndarray.flatten(coord_array), radii_array)
		out_p = []
		for surface_id, real_id in enumerate(pos_array):
			out_p.append([real_id, protein_sasa.atomArea(surface_id)])
		#sum up sasa residues-wise
		res_out_p = {}
		for i in out_p:
			res_out_p.setdefault(i[0], 0)
			res_out_p[i[0]] += i[1]

		sasa_dict_full[gene] = res_out_p
  
  
		######calculate domain-domain interfaces######
  
		if gene not in domain_annotation:
			continue
  
		temp = {}
		#index over annotated protein domains
		for p, domain in enumerate(domain_annotation[gene]):
			res_out_p = {}
			#get ids for domain A 
			d_ids = np.where((pos_array >= domain[1]) & (pos_array <= domain[2]))[0]
			#index AlphaFold data for domain A
			coord_array_d = coord_array[d_ids, :]
			radii_array_d = radii_array[d_ids]
			pos_array_d = pos_array[d_ids]
			#calculate sasa of domain A
			protein_sasa = freesasa.calcCoord(np.ndarray.flatten(coord_array_d), radii_array_d)
			out_p = []
			for surface_id, real_id in enumerate(pos_array_d):
				out_p.append([real_id, protein_sasa.atomArea(surface_id)])
			for i in out_p:
				res_out_p.setdefault(i[0], 0)
				res_out_p[i[0]] += i[1]
			#index again over annotated protein domains
			for pp, ddomain in enumerate(domain_annotation[gene]):
				res_out_ppp = {}
				#only allow combinations
				if p>=pp:
					continue
				#get ids for domain B
				dd_ids = np.where((pos_array >= ddomain[1]) & (pos_array <= ddomain[2]))[0]
				#merge ids of domain A and B
				ddd_ids = np.r_[d_ids, dd_ids]
				#index AlphaFold data for domain A+B
				coord_array_ddd = coord_array[ddd_ids, :]
				radii_array_ddd = radii_array[ddd_ids]
				pos_array_ddd = pos_array[ddd_ids]
				#calculate sasa of domain A+B
				protein_sasa = freesasa.calcCoord(np.ndarray.flatten(coord_array_ddd), radii_array_ddd)
				out_ppp = []
				for ddd_surface_id, ddd_real_id in enumerate(pos_array_ddd):
					out_ppp.append([ddd_real_id, protein_sasa.atomArea(ddd_surface_id)])
				for i in out_ppp:
					res_out_ppp.setdefault(i[0], 0)
					res_out_ppp[i[0]] += i[1]
				#calculate interface by taking the additional sasa of domain A compared to the sasa of domain A+B 
				interface = 0
				for i, j in res_out_p.items():
					interface += abs(j-res_out_ppp[i])

				domain_interface_dict.setdefault(str(p)+'_'+str(pp), 0)
				domain_interface_dict[str(p)+'_'+str(pp)] = interface

with open(path + "sasa_dict_full.pick", "wb") as handle:
	pickle.dump(sasa_dict_full, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open(path + "domain_interface_dict.pick", "wb") as handle:
	pickle.dump(domain_interface_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
