import copy
import Candidate
import pickle

def prob_cover_genes_lst(candidate, genes_lst):
	cover_all = 1
	for gene in genes_lst:
		cover_all *= candidate.genes_score_dict[gene]
	return cover_all

def find_set_cover(best_permutations_DS, sg_genes_dict, thr, genes_sg_dict = None):
	'''for now, might won't work in a case when there is a gene that isn't covered by any of the permutations in the best_permutations_DS. not finished. can make it more readble'''
	temp_best_perm_DS = copy.copy(best_permutations_DS)
	res = list()#[temp_best_perm_DS[0]]
	if genes_sg_dict:
		for gene, targets in genes_sg_dict.items():
			if len(targets) == 0:
				print("no targets for gene " + gene)
				genes_name_lst.remove(gene)
				continue
			c = Candidate.Candidate(targets[0])
			c.fill_default_fildes(sg_genes_dict[targets[0]])
			temp_best_perm_DS.append(c)

	uncovered_genes = set()
	for sg, genesLst in sg_genes_dict.items():
		for gene in genesLst:
			uncovered_genes.add(gene)
	while(len(uncovered_genes)) > 0 and len(temp_best_perm_DS) > 0:
		#print(uncovered_genes)
		#for gene in uncovered_genes:
			#print(gene)
	##going over all the permutations, and return the permutation that cover the maximal amount of genes haven't been covered yet, in the highest probability among the maximal covered permutations
		#print('len uncovered genes', len(uncovered_genes))
		best_current_perm, best_num_of_coverd, best_prob_of_covered = None, 0,0  #best_current_perm is the hole tuple
		i = 0
		while i < (len(temp_best_perm_DS)):
			new_genes_coverd = list()#0
			for gene, score in temp_best_perm_DS[i].genes_score_dict.items():
				if gene in uncovered_genes and score >= thr:
					new_genes_coverd.append(gene)
					#uncovered_genes.remove(gene)
			if len(new_genes_coverd) == 0:
				i+=1
				continue
				#del temp_best_perm_DS[i]
			elif len(new_genes_coverd) >= best_num_of_coverd:## and temp_best_perm_DS[i][2] > best_prob_of_covered:  ##need to check if 2 is the right index, and not 1.
				#print(new_genes_coverd)
				if len(new_genes_coverd) > best_num_of_coverd or prob_cover > best_prob_of_covered: # cover more gene or cover the same amount with greater prob.
					prob_cover = prob_cover_genes_lst(temp_best_perm_DS[i], new_genes_coverd)
				#if prob_cover > best_prob_of_covered:
					best_num_of_coverd, best_prob_of_covered = len(new_genes_coverd), prob_cover
					best_current_perm = temp_best_perm_DS[i]
			i+=1
		if(best_current_perm):
			res.append(best_current_perm)
			for gene, score in best_current_perm.genes_score_dict.items():
				if gene in uncovered_genes and score >= thr: #there is a probability that this gene had already been covered bya prevuis sgRNA
					uncovered_genes.remove(gene)
	return res

def find_w_set_cover(best_permutations_DS, sg_genes_dict, thr, distance_matrix, genes_sg_dict = None):
	return find_w_set_cover_heated(best_permutations_DS, sg_genes_dict, thr, distance_matrix, alfa = 1, genes_sg_dict = None)

def find_w_set_cover_heated(best_permutations_DS, sg_genes_dict, thr, distance_matrix, alfa, genes_sg_dict = None):
	'''for now, might won't work in a case when there is a gene that isn't covered by any of the permutations in the best_permutations_DS. not finished. can make it more readble'''
	temp_best_perm_DS = copy.copy(best_permutations_DS)
	res = list()#[temp_best_perm_DS[0]]
	gene_names_lst = distance_matrix.names
	if genes_sg_dict:
		for gene, targets in genes_sg_dict.items():
			if len(targets) == 0:
				print("no targets for gene " + gene)
				genes_name_lst.remove(gene)
				continue
			c = Candidate.Candidate(targets[0])
			c.fill_default_fildes(sg_genes_dict[targets[0]])
			temp_best_perm_DS.append(c)

	uncovered_genes = set()
	for sg, genesLst in sg_genes_dict.items():
		for gene in genesLst:
			uncovered_genes.add(gene)
	while(len(uncovered_genes)) > 0 and len(temp_best_perm_DS) > 0:
		best_current_perm, best_num_of_coverd, best_prob_of_covered = None, 0,0  #best_current_perm is the hole tuple
		best_w = len(uncovered_genes)
		i = 0
		while i < (len(temp_best_perm_DS)):
			#find utility of sgRNA
			new_genes_coverd = list()
			for gene, score in temp_best_perm_DS[i].genes_score_dict.items():
				if gene in uncovered_genes and score >= thr:
					new_genes_coverd.append(gene) #new_genes_coverd.append((gene, score)) #
			if len(new_genes_coverd) == 0:
				i+=1
				continue
			#compute the weight
			price = 0 # the lower the price, the lighter the set
			deniminator = 0 # deniminator = len(new_genes_coveres)*(len(new_genes_coveres) - 1)
			for j in range(len(new_genes_coverd)):
				score_j = temp_best_perm_DS[i].genes_score_dict[new_genes_coverd[j]]
				for k in range(i, len(new_genes_coverd)):
					score_k = temp_best_perm_DS[i].genes_score_dict[new_genes_coverd[k]]
					curr_avg = (1 - score_j + 1 - score_k)/2
					deniminator += curr_avg
					index_j, index_k = gene_names_lst.index(new_genes_coverd[j]), gene_names_lst.index(new_genes_coverd[k])
					dist = distance_matrix[index_j, index_k]
					price += (dist**alfa) * curr_avg
			if deniminator == 0:
				price = 1# sure?
			else:
				price = price/deniminator
			w = price / len(new_genes_coverd)

			#del temp_best_perm_DS[i]
			#keep the best sgRNA in current iteration
			if w <= best_w:## and temp_best_perm_DS[i][2] > best_prob_of_covered:  ##need to check if 2 is the right index, and not 1.
				#print(new_genes_coverd)
				#if len(new_genes_coverd) > best_num_of_coverd or prob_cover > best_prob_of_covered: # cover more gene or cover the same amount with greater prob.
				prob_cover = prob_cover_genes_lst(temp_best_perm_DS[i], new_genes_coverd)
				if prob_cover > best_prob_of_covered:
					best_w, best_prob_of_covered = w, prob_cover
					best_current_perm = temp_best_perm_DS[i]
			i+=1
		#add the best sgRNA to res
		if(best_current_perm):
			res.append(best_current_perm)
			for gene, score in best_current_perm.genes_score_dict.items():
				if gene in uncovered_genes and score >= thr: #there is a probability that this gene had already been covered bya prevuis sgRNA
					uncovered_genes.remove(gene)
	return res


###old version###
def find_w_set_cover_old(candidates_lst ,distance_matrix, Omega):
	annealing_coef = 0
	return find_w_set_cover_heated(candidates_lst ,distance_matrix, annealing_coef, Omega)

def find_w_set_cover_sevral(candidates_lst ,distance_matrix):
	res = [] #list of lists f lists. each level 2 list represent the size of the set cover, and contains the best candidate with the lowest considering in homology and with highest consideting in homology, correspondingly
	annealing_coef = 0
	min_set_cover = find_w_set_cover_heated(candidates_lst ,distance_matrix, annealing_coef)
	res.append([0,min_set_cover])
	current_set_cover = min_set_cover
	for i in range(1,10):
		prev_set_cover = current_set_cover
		current_set_cover = find_w_set_cover_heated(candidates_lst ,distance_matrix, i/10)
		if len(current_set_cover)>len(prev_set_cover):
			res[len(res)-1].append(prev_set_cover)
			res.append([i/10,current_set_cover])
	return res[0]

	annealing_coef = 0.5
	return find_w_set_cover_heated(candidates_lst ,distance_matrix, annealing_coef)

def find_w_set_cover_heated_old(candidates_lst ,distance_matrix, annealing_coef, Omega):
	'''	 the standard greedy aproximation algorithm
	:param candidates_lst ; a toy example:  ([['ACGCACCC', 0.6, 3.7961971025899606e-06, [('gene1', 0.01608054522924407), ('gene2', 0.01424411400247827), ('gene4', 0.01657343550446999)], [['gene1', [('ACGTACGT', {3: {'C', 'T'}, 6: {'C', 'G'}, 7: {'C', 'T'}})]], ['gene2', [('ACGTAGCT', {3: {'C', 'T'}, 5: {'C', 'G'}, 7: {'C', 'T'}})]], ['gene3', [('ACGTATTG', {3: {'C', 'T'}, 5: {'C', 'T'}, 6: {'C', 'T'}, 7: {'C', 'G'}})]], ['gene4', [('ATGCACGT', {1: {'C', 'T'}, 6: {'C', 'G'}, 7: {'C', 'T'}})]], ['gene5', [('ATGCATGC', {1: {'C', 'T'}, 5: {'C', 'T'}, 6: {'C', 'G'}})]]]], ['ACGCACCT', 0.6, 4.393668496780193e-05, [('gene1', 0.036452247191011256), ('gene2', 0.03157967032967035), ('gene4', 0.03816764705882347)], [['gene1', [('ACGTACGT', {3: {'C', 'T'}, 6: {'C', 'G'}})]], ['gene2', [('ACGTAGCT', {3: {'C', 'T'}, 5: {'C', 'G'}})]], ['gene3', [('ACGTATTG', {3: {'C', 'T'}, 5: {'C', 'T'}, 6: {'C', 'T'}, 7: {'T', 'G'}})]], ['gene4', [('ATGCACGT', {1: {'C', 'T'}, 6: {'C', 'G'}})]], ['gene5', [('ATGCATGC', {1: {'C', 'T'}, 5: {'C', 'T'}, 6: {'C', 'G'}, 7: {'C', 'T'}})]]]] ...
	:param distance_matrix:
	:return:
	'''
	res = []
	names_lst = distance_matrix.names #names of genes by their  order in the distance matrix
	covered_genes = [0 for i in range(len(names_lst))] # 0 if not covered, 1 if covered
	used_candidates = [0 for i in range(len(candidates_lst))]  # the same here
	prices_lst = [0 for i in range(len(candidates_lst))]
	while(0 in covered_genes):
		print(covered_genes)
		for c in range(len(candidates_lst)): #going over all of the candidates
			single_gene_by_the_set = -1
			weight = 0 # a list of tuples to consider to the avg that makes the weight
			num_of_pairs = 0 #to be later davide the weight
			num_of_coveres_genes = 0
			if used_candidates[c] == 1:#this candidate was already chosen
				continue
			for gene_name in range(len(list(candidates_lst[c].genes_score_dict.keys()))):
				if covered_genes[names_lst.index(gene_name)] == 1:# #for knowing which to include in take_into_agv list
					continue
				if num_of_coveres_genes == 0:
					single_gene_by_the_set = i
				num_of_coveres_genes += 1
				for j in range(i+1, len(candidates_lst[c][3])):
					if covered_genes[names_lst.index(candidates_lst[c][3][j][0])] == 1:# #for knowing which to include in take_into_agv list
						continue
					num_of_pairs += 1 # also = num_of_coveres_genes(num_of_coveres_genes-1)
					single_gene_by_the_set = -1
					#still need to normalise the 1/cut prob.
					#weight += (distance_matrix[i,j]**annealing_coef) *(1/candidates_lst[c][3][i][1] + 1/candidates_lst[c][3][j][1])/2  # (1/cleaving_probobility_i+ 1/cleaving_proboblity_j)/2 # the distance times the avereged covering probobility is the weight
					#in assuminig distance matrix is similarity matrix
					if annealing_coef == 0:
						weight += (candidates_lst[c][3][i][1] + candidates_lst[c][3][j][1])/2
					else: #dosen't have to have the 'else' here
						weight += ((1-distance_matrix[i,j])**annealing_coef) *(candidates_lst[c][3][i][1] + candidates_lst[c][3][j][1])/2 #is it correct? # (1/cleaving_probobility_i+ 1/cleaving_proboblity_j)/2 # the distance times the avereged covering probobility is the weight
						print(1-distance_matrix[i,j])
			if num_of_coveres_genes == 0:  #this candidate is irelevant for now
				prices_lst[c] = 100000 #just a big number
			else:
				if num_of_coveres_genes == 1: #assuming taking an sgRNA with perfect mach to the target site
					#weight = 1 #candidates_lst[single_gene_by_the_set] #cleaving_probobility_single_gene_by_the_set  #givn by the flage single_gene_by_the_set
					weight = 2 - candidates_lst[c][3][single_gene_by_the_set][1]  #just for testing - have to find what to do with those perfect mached sg.
				else:
					weight = 2 - (weight/num_of_pairs)


				#prices_lst[c] = weight/num_of_coveres_genes #standard weight by the set cover greedy algorithm
				prices_lst[c] = 1/num_of_coveres_genes #unwighted set cover

		#stage 2: find the minimal priced set
		min_priced = [-1, max(prices_lst)] #index and price
		for c in range(len(prices_lst)):
			#print(prices_lst[c])
			if used_candidates[c] == 0 and prices_lst[c] < min_priced[1]:
				min_priced = [c, prices_lst[c]]
				#print("temp_min_price", min_priced)
		#add the result to the cover and update the covered genes:
		print("min priced",min_priced)
		res.append(candidates_lst[min_priced[0]])
		used_candidates[min_priced[0]] = 1
			#update the covered genes:
		for i in range(len(candidates_lst[min_priced[0]][3])): #the covered genes lst
			gene_index = names_lst.index(candidates_lst[min_priced[0]][3][i][0])
			covered_genes[gene_index] = 1
	return res

###
def test1(path = '/bioseq/data/results/multicrispr/1516613143' , thr = 0.45):
	best_permutations_DS, sg_genes_dict = pickle.load(open("/".join([path, "res_in_lst.p"]),'rb')), pickle.load(open("/".join([path, "sg_genes_dict.p"]), 'rb'))
	find_set_cover(best_permutations_DS, sg_genes_dict, thr)