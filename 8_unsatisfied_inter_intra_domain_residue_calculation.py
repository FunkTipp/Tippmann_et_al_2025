import numpy as np
import pickle

def res_saturation(contact_map, last_emerged):
	
	'''
	calculates which residues are unsatisfied and which contacts can be formed at a certain NC length
	and distinguished between inter and intra domain contacts
 
	assumes that inter-domain contacts remains unsatisfied until translation termination
	'''
	
	stable_contact_ids = np.where((contact_map[:,0] <= last_emerged) & \
						   (contact_map[:,1] <= last_emerged))[0]
	stable_contact_ids_remains = np.where((contact_map[:,0] <= last_emerged) & \
						   (contact_map[:,1] <= last_emerged) & \
						   (contact_map[:,3] != 1))[0]
	unsat_interD_contact_ids = np.where((contact_map[:,0] <= last_emerged) & \
						   (contact_map[:,1] > last_emerged) & \
						   (contact_map[:,3] == 1))[0]
	unsat_interD_contact_ids_remains = np.where((contact_map[:,0] <= last_emerged) & \
						   (contact_map[:,1] <= last_emerged) & \
						   (contact_map[:,3] == 1))[0]
	unsat_intraD_contact_ids = np.where((contact_map[:,0] <= last_emerged) & \
						   (contact_map[:,1] > last_emerged) & \
						   (contact_map[:,3] == 2))[0]
	#sum up interaction strength for every residues in each category
	#1: all stable contacts
	#2: all stable contacts for unsatisfied inter_domain interactions
	#3: all interdomain contacts
	#4: additional interdomain conatcs for inter_domain interactions
	#5: all intradomain contacts
	res_mat = np.zeros((last_emerged, 5))
	for res in np.arange(last_emerged):
		stable_res_ids = np.where(contact_map[stable_contact_ids,:2] == res)[0]
		res_mat[res, 0] = np.sum(contact_map[stable_res_ids, 2])
		stable_res_inter_remains = np.where(contact_map[stable_contact_ids_remains,:2] == res)[0]
		res_mat[res, 1] = np.sum(contact_map[stable_res_inter_remains, 2])
		unsat_inter_res_ids = np.where(contact_map[unsat_interD_contact_ids,:2] == res)[0]
		res_mat[res, 2] = np.sum(contact_map[unsat_inter_res_ids, 2])
		unsat_inter_res_add_ids = np.where(contact_map[unsat_interD_contact_ids_remains,:2] == res)[0]
		res_mat[res, 3] = np.sum(contact_map[unsat_inter_res_add_ids, 2])
		unsat_intra_res_ids = np.where(contact_map[unsat_intraD_contact_ids,:2] == res)[0]
		res_mat[res, 4] = np.sum(contact_map[unsat_intra_res_ids, 2])
	return res_mat

#load amino acid sequence of MG1655 genes to compare gene_names and gene_length
path = "/Users/.../"
with open(path + ".../aa_seq_MG1655.pick", "rb") as handle:
	aa_seq = pickle.load(handle)
#loads the AlphaFold data base with domain annotation generated in the 5_structure_assisted_domain_annotation.py script
with open(path + ".../AF_domains.pick", "rb") as handle:
	afD_domains = pickle.load(handle)	

#iterate over contact maps
unsatisfied_domain_residues_out = {}
for g_count, (gene, contacts) in enumerate(afD_domains.items()):
	if gene not in genes_to_take:
		continue
	#get contacts and gene length
	contact_map = contacts[3]
	gene_length = len(aa_seq[gene])
	#create outfile
	#row: translation state
	#column 1: unsatisfied intra-domain contacts
	#column 2: unsatisfied intra-domain residues
	#column 3: unsatisfied inter-domain contacts
	#column 4: unsatisfied inter-domain residues

	outfile = np.zeros((gene_length, 4))
	#iterate over translation states
	for ts in np.arange(gene_length):
		#correct for the tunnel
		last_emerged = ts-30
		if last_emerged < 0:
			continue
		#get emerged residue saturation
		residue_mat = res_saturation(contact_map, last_emerged)
		if len(residue_mat[:, 0]) == 0:
			continue
		#get stabalized residues to determine compacted regions
		stable_res_A = residue_mat[:, 0]
		stable_res_B = np.amin([np.cumsum(stable_res_A), np.cumsum(stable_res_A[::-1])[::-1]], axis=0)
		stable_ids = np.where(stable_res_B > 10)[0]
		try:
			start = np.min(stable_ids)
			end = np.max(stable_ids)
		except:
			start = last_emerged
			end = last_emerged
		####intra####
		intra = residue_mat[start:end, 4]
		outfile[ts, 0] = np.sum(intra)
		outfile[ts, 1] = np.count_nonzero(intra)

		####inter####
		inter_stable_res_A = residue_mat[:, 1]
		inter_stable_res_B = np.amin([np.cumsum(inter_stable_res_A), np.cumsum(inter_stable_res_A[::-1])[::-1]], axis=0)
		inter_stable_ids = np.where(inter_stable_res_B > 10)[0]
		try:
			inter_start = np.min(inter_stable_ids)
			inter_end = np.max(inter_stable_ids)
		except:
			inter_start = last_emerged
			inter_end = last_emerged
		inter1 = residue_mat[inter_start:inter_end, 2]
		inter2 = residue_mat[inter_start:inter_end, 3]
		outfile[ts, 2] = np.sum(inter1)+np.sum(inter2)
		outfile[ts, 3] = np.count_nonzero(inter1)+np.count_nonzero(inter2)

	unsatisfied_domain_residues_out[gene] = outfile

with open(path + "unsatisfied_domain_residues_out.pick", "wb") as handle:
	pickle.dump(unsatisfied_domain_residues_out, handle, protocol=pickle.HIGHEST_PROTOCOL)

