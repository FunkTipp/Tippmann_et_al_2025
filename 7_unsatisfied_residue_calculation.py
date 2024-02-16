import numpy as np
import pickle

def res_saturation(contact_map, last_emerged):
	
	'''
	calculates which residues are unsatisfied and which contacts can be formed at a certain NC length
	'''
	
	stable_contact_ids = np.where((contact_map[:,0] <= last_emerged) & \
						   (contact_map[:,1] <= last_emerged))[0]
	unsat_contacts = np.where((contact_map[:,0] <= last_emerged) & \
						   (contact_map[:,1] > last_emerged))[0]
	#sum up interaction strength for every residues in each category
	#1: all stable contacts
	#2: all unsatisfied contacts
	res_mat = np.zeros((last_emerged, 2))
	for res in np.arange(last_emerged):
		stable_res_ids = np.where(contact_map[stable_contact_ids,:2] == res)[0]
		res_mat[res, 0] = np.sum(contact_map[stable_res_ids, 2])
		unsat_ids = np.where(contact_map[unsat_contacts,:2] == res)[0]
		res_mat[res, 1] = np.sum(contact_map[unsat_ids, 2])
	return res_mat

#load amino acid sequence of MG1655 genes to compare gene_names and gene_length
path = "/Users/.../"
with open(path + ".../aa_seq_MG1655.pick", "rb") as handle:
	aa_seq = pickle.load(handle)
#loads the AlphaFold data base generated in the 4_get_AF_data.py script
with open(path + ".../1_AF_data.pick", "rb") as handle:
	afD = pickle.load(handle)	

#iterate over contact maps
unsatisfied_residues_out = {}
for g_count, (gene, contacts) in enumerate(afD.items()):
	#get contacts and gene length
	contact_map = contacts[3]
	gene_length = len(aa_seq[gene])
	#create outfile
	#row: translation state
	#column 1: unsatisfied contacts
	#column 2: unsatisfied residues
	outfile = np.zeros((gene_length, 2))
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
		all_unsatisfied = residue_mat[start:end, 1]
		outfile[ts, 0] = np.sum(all_unsatisfied)
		outfile[ts, 1] = np.count_nonzero(all_unsatisfied)

	unsatisfied_residues_out[gene] = outfile

with open(path + "unsatisfied_residues_out.pick", "wb") as handle:
	pickle.dump(unsatisfied_residues_out, handle, protocol=pickle.HIGHEST_PROTOCOL)

