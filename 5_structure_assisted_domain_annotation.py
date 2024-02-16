import numpy as np
import pickle
from scipy import signal
import json
import os
#you need to install curl first

path = "/Users/.../"
#loads the AlphaFold data base generated in the 4_get_AF_data.py script
with open(path + ".../1_AF_data.pick", "rb") as handle:
	afD = pickle.load(handle)	
 #load amino acid sequence of MG1655 genes to compare gene_names and gene_length
with open(path + ".../aa_seq_MG1655.pick", "rb") as handle:
	aa_seq = pickle.load(handle)



####1. get AF domain data#####

# triangle function to normalize residue saturation and contact ratio data
conv = signal.triang(11)

signal_D_dict = {}
for gene, data in afD.items():
	if gene not in aa_seq:
		continue
	#load contact table
	contacts = data[3]
	#get gene length
	gene_length = len(aa_seq[gene])

	#collect data from structure
	measures = np.zeros((gene_length, 3))
	#iterate over structure from N to C terminus
	for ts in np.arange(gene_length):

		#calculate residue saturation
		stable_contact_ids = np.where((contacts[:,0] <= ts) & (contacts[:,1] <= ts))[0]
		measures[ts, 0] = np.sum(contacts[stable_contact_ids, 2])

		#calculate contact ratio
		res_A_in_ids = np.isin(contacts[:,0], np.arange(0, ts)).astype(np.int_)
		res_B_in_ids = np.isin(contacts[:,1], np.arange(0, ts)).astype(np.int_)
		conts_in = np.sum(contacts[(res_A_in_ids*res_B_in_ids).astype(np.bool_), 2])+1
		conts_out = np.sum(contacts[np.abs(res_A_in_ids-res_B_in_ids).astype(np.bool_), 2])+1
		measures[ts, 1] = np.log2(conts_in/conts_out)

		#calculate local contact density
		local_contact_ids = np.where((contacts[:,0] > ts-30) & (contacts[:,1] < ts+30))[0]
		measures[ts, 2] = np.sum(contacts[local_contact_ids, 2])
	#calculate first discrete differences and convolutions, apply normalization
	A = np.convolve(np.diff(measures[:, 0]), conv)[5:-5]
	A[A<5] = 5
	B = np.convolve(-np.diff(measures[:, 1]), conv)[5:-5]
	B[B<0] = 0
	B[B>1] = 1
	B[:20] = 0
	C = measures[1:, 2]
	D = B/(C*A)
	signal_D_dict[gene] = D

#with open(path + "signal_D_dict.pick", "wb") as handle:
#	pickle.dump(signal_D_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



####2. get CATH data#####

CATH_dict = {}
#iterate over genes
for p, (gene, seq) in enumerate(aa_seq.items()):
	#1write the sequence into a fasta file
	A = open(path + '5_temp_1.fasta', 'w')
	A.write('fasta=>'+gene+'\n'+seq)
	A.close()
	#submit the job to the CATH server
	os.system('curl -o "'+path+'5_temp_2.json" -w "\n" -s -X POST -H "Accept: application/json" --data-binary "@'+path+'5_temp_1.fasta" http://www.cathdb.info/search/by_funfhmmer')
	#extract the jobnumber form the recieved job file
	with open(path + '5_temp_2.json') as json_file:
		B = json.load(json_file)
	id_ = B['task_id']
	#write the code for checking your submission status
	os_in_B = 'curl -o "'+path+'5_temp_3.json" -w "\n" -s -X GET -H "Accept: application/json" http://www.cathdb.info/search/by_funfhmmer/check/'+id_
	print(p, gene)
	#check every second if the sumbmission is complete
	#usually a gene takes 10-20 seconds
	c=0
	while True:
		time.sleep(10 - ((time.time() - starttime) % 10))
		c+=10
		if c%10 == 0:
			print('- waiting for ' + str(c) + ' seconds -')
		os.system(os_in_B)
		#some genes don't give a job status file (for what ever reason),
		#so this try/except command is a safety that the program runs trough all genes
		try:
			with open(path + '5_temp_3.json') as json_file:
				C = json.load(json_file)
		except:
			break
		if C['success'] == 1:
			break
	#get the results
	os_in_C = 'curl -o "'+path+'5_temp_4.json" -w "\n" -s -X GET -H "Accept: application/json" http://www.cathdb.info/search/by_funfhmmer/results/'+id_
	os.system(os_in_C)
	#save the results (only the Top 50 are given)
	#note that genes with no results don't give an outputfile
	try:
		with open(path + '5_temp_4.json') as json_file:
			D = json.load(json_file)
	except:
		continue
	CATH_dict[gene] = D['funfam_scan']['results'][0]['hits']

#with open(path + "CATH_D_temp.pick", "wb") as handle:
#	pickle.dump(CATH_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



####3. merge AF and CATH data#####
domain_annotation = {}
for count, (gene, data) in enumerate(CATH_dict.items()):
	if gene not in signal_D_dict:
		continue

	# find peaks in the AlphaFold based domain signal_D
	signal_A = B[gene]
	signal_A_peaks = signal.find_peaks(signal_A, height=0.00002)[0]
	# add gene start and stop positions as peaks
	signal_A_peaks = np.append(signal_A_peaks, np.array([0, len(aa_seq[gene])]))
	if len(data) == 0 or len(signal_A) == 0:
		continue
	temp = []
	#iterate over all CATH_dict entries found in the CATH database search
	for p, i in enumerate(data[::-1]):
		for pp, j in enumerate(i['hsps']):
			start = j['query_start']
			stop = j['query_end']
			domain_name = i['match_cath_id']['id']
			#exclude domains shorter than 50 amino acids
			if stop-start < 50:
				continue
			#calculate metric that measures the distance between the CATH domain  \
   			#and the signal_D peak positions
			start_penatlty = np.amin(np.abs(signal_A_peaks-start))
			stop_penatlty = np.amin(np.abs(signal_A_peaks-stop))
			penalty = start_penatlty+stop_penatlty
			query_range = np.arange(start, stop)
			#iterate over already accepted domains in the temp list
			#take new domains that don't overlapp with already axisting domains
			#replace already accepted domains with better fitting domains
			take = True
			for tp, t in enumerate(temp):
				match_range = np.arange(t[1], t[2])
				compare = np.count_nonzero(np.isin(query_range, match_range))
				if compare > 10:
					take = False

					if penalty < t[3] and compare/len(match_range) > 0.8:
						temp[tp] = [i['match_cath_id']['id'], start, stop, penalty]
			if take == True:	
				temp.append([i['match_cath_id']['id'], start, stop, penalty])
	if len(temp) == 0:
		continue	
	# take the largest domain for every region, exclude domoains spanning the complete gene
	if len(temp) == 1:
		domain_annotation[gene] = temp
	else:
		temp = sorted(temp, key=lambda x: x[1])
		repeat_ = True		
		while repeat_ == True:
			repeat_ = False	
			for p, d1 in enumerate(temp):
				query_range = np.arange(d1[1], d1[2])
				for pp, d2 in enumerate(temp):
					if p > pp:
						if d1[1] == d2[1] and d1[2] == d2[2]:
							try:
								del temp[p]
							except:
								print(gene, '--- !!! ---')
							repeat_ = True
						match_range = np.arange(d2[1], d2[2])
						compare = np.count_nonzero(np.isin(query_range, match_range))
						if len(query_range) == len(match_range) or compare <= 25:
							continue
						else:
							if len(query_range) > len(ref[gene])*0.9:
								del temp[p]	
								repeat_ = True
							elif len(match_range) > len(ref[gene])*0.9:
								del temp[pp]	
								repeat_ = True											
							else:
								if d1[2]-d1[1] > d2[2]-d2[1]:
									del temp[pp]
								else:
									del temp[p]            
								repeat_ = True						


	# sort the annotation from N to C terminus
	# add new domains to N and C terminus if more than 75 reidues of the N/C terminus are not annotated with a domain
	temp_sort = sorted(temp, key=lambda x: x[1])
	if temp_sort[0][1] > 75:
		temp_sort.insert(0, ['FG', 0, temp_sort[0][1]-1, -1])
	if temp_sort[-1][2] < len(aa_seq[gene])-75:
		temp_sort.append(['FG', temp_sort[-1][2]+1, len(aa_seq[gene]), -1])

	domain_annotation[gene] = temp_sort

with open(path + "domain_annotation.pick", "wb") as handle:
	pickle.dump(domain_annotation, handle, protocol=pickle.HIGHEST_PROTOCOL)



####4. update AF contact data with new domain annotation#####

AF_domains = {}
for gene_count, (gene, domains) in enumerate(domain_annotation.items()):
	if gene not in afD:
		continue
	contacts = afD[gene][3]
	contacts_add = np.zeros((len(contacts), 5))
	#check if contacts are spread among two or are within one domain
	for p, d in enumerate(domains):
		for pp, dd in enumerate(domains):
			A = np.isin(contacts[:, 0], np.arange(d[1], d[2]))
			B = np.isin(contacts[:, 1], np.arange(dd[1], dd[2]))
			C = A*B
			# add domain border information to contact table
			contacts_add[C, 1] = d[1]
			contacts_add[C, 2] = d[2]
			contacts_add[C, 3] = dd[1]
			contacts_add[C, 4] = dd[2]
			# add 1 if contacts are between domains
			if p != pp:
				contacts_add[C, 0] += 1
			# add 2 if contacts are within a single domain
			if p == pp:
				contacts_add[C, 0] += 2
			cc = np.where(contacts_add[:, 3]==3)[0]
			contacts_add[cc, 0] = 2
	contacts_update = np.c_[contacts, contacts_add]
	AF_domains[gene] = [afD[gene][0], afD[gene][1], afD[gene][2], contacts_update]

with open(path + "AF_domains.pick", "wb") as handle:
	pickle.dump(AF_domains, handle, protocol=pickle.HIGHEST_PROTOCOL)
