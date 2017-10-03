import Candidate
import pickle
import copy
import os

import cplex
from cplex.exceptions import CplexError

##########cplexCovers_part##########
import CplexCovers

def call_CplexCovers(list_of_candidates, genes_lst, thr, method):

	list_of_candidates, genes_lst = pickle.load(open(list_of_candidates, 'rb')), pickle.load(open(genes_lst, 'rb'))
	if len(list_of_candidates) > 999: #if using the free version of Cplex
		list_of_candidates = list_of_candidates[:999]
	try:
		if method == "SC":
			cplex_problem_object = CplexCovers.CplexSetCover(list_of_candidates, genes_lst, thr)
		elif method == "F_SC":
			cplex_problem_object = CplexCovers.Cplex_fuzzy_set_cover(list_of_candidates, genes_lst, thr)
	except CplexError as exc:
		print("no solution")
		return
	else:
		print("which method?")
		print(method)
	return cover_from_cplex_promblem_obj(list_of_candidates, cplex_problem_object)
	

	#####partly fuzzy set cover####
def bounded_cover(lst_of_candidates_path, genes_lst_path, cover_size, maximal_thr = 0.99):
	'''if thr is lower the the maximal thr and cover size is at least as long as the genes_lst, the result is the fuzzy set cover''' 
	# first, find the fuzzy set cover
	lst_of_candidates, genes_lst = pickle.load(open(lst_of_candidates_path, 'rb')), pickle.load(open(genes_lst_path, 'rb'))
	if len(lst_of_candidates) > 900:
		lst_of_candidates = lst_of_candidates[:900]
	#maximal_thr = 0.99
	try:
		F_sc = CplexCovers.Cplex_fuzzy_set_cover(lst_of_candidates, genes_lst, maximal_thr) #return a problem object

	except CplexError as exc:
	#if failed to find a solution, go greedy.
	#if my_prob.solution.get_status() !=  101:
		#best_score, group = gready_cover(lst_of_candidates, genes_lst, cover_size) #works well as well
		best_score, group = full_cover_V0(lst_of_candidates, genes_lst, cover_size) #exact cover. possible to compute fast since the initial list is very short
		return best_score, group
	group = cover_from_cplex_promblem_obj(lst_of_candidates, F_sc)	#group here is a list of indices

	#if all the genes are covered, we are done.
	if F_sc.solution.get_objective_value() <= cover_size:
		best_score = calculate_score(lst_of_candidates, group, genes_lst)
	else:
	#cover_size is to large, need to pick from it
		best_score, group = reverse_gready_cover(lst_of_candidates, genes_lst, cover_size, group)

	return best_score, group  #score here is how many genes are not covered
	
def from_index_to_candidate_lst(full_candidates_lst,group):
	res = list()
	for index in group:
		res.append(full_candidates_lst[index])
	return res

def cover_from_cplex_promblem_obj(lst_of_candidates, cplex_obj):
	'''return the indecies of the cover'''
	res = []
	x = cplex_obj.solution.get_values()
	
	numcols = cplex_obj.variables.get_num() #sgRNAs
	numrows = cplex_obj.linear_constraints.get_num() #genes

	slack = cplex_obj.solution.get_linear_slacks()
	x = cplex_obj.solution.get_values()


	for j in range(numcols):
		if x[j] == 1:
			res.append(j)
	return res
	

def calculate_score(lst_of_candidates, lst, genes_lst):
	'''
	:param lst_of_candidates: lst of all the candidates. "res" of algorithm A.
	:param lst: list of candidates of the current group - each candidate represented as an index
	:param genes_lst: list of the genes of the family
	:return: the objective for this group: sigma pi 1-phi(sg_j, gene)*Xj
	'''
	#cdef int score, gene_score
	score = 0

	for gene in genes_lst:
		gene_score = 1
		for candidate in lst:
			if gene in lst_of_candidates[candidate].genes_score_dict:
				gene_score *= (1- lst_of_candidates[candidate].genes_score_dict[gene])
		#gene_score = 1- gene_score
		score += gene_score
	return score
#################
## brute force ##
#################
def full_cover_V0(lst_of_candidates, genes_lst, k):
	'''
	Going over all of the k-mer in O(1) space. Finds the best scored one
	:param lst_of_candidates: the res from the algorithm before the cover
	:param k: size of wanted group
	:return:
	'''
	#cdef int num_of_candidates, best_score, score, not_done
	#cdef int lst[k]
	num_of_candidates = len(lst_of_candidates)
	lst = [i for i in range(k)] #thats the first set
	not_done = True
	#best_score = len(genes_lst) #find the lowest score
	best_score = len(genes_lst)
	best_group = []
	#best_score = 0
	while(not_done):
		score = calculate_score(lst_of_candidates, lst, genes_lst)
		#update best score and best group
		if score < best_score:
			best_score, best_group = copy.deepcopy(score), copy.deepcopy(lst)
		not_done = increment(lst, num_of_candidates)
	return best_score, best_group



def upper_bound(lst, index, num_of_candidates):
	return lst[index] == num_of_candidates - 1 - (len(lst)-1 -index)

def reset(lst, index):
	'''
	reset all of the digit from index and forword.
	:param lst:
	:param index: > 0
	:return:
	'''
	#cdef int j
	lst[index -1] += 1
	for j in range(index, len(lst)):
		lst[j] = lst[j-1] + 1
	return 1 #1 is True


def increment(lst, num_of_candidates):
	#cdef int i
	i = len(lst) -1
	while upper_bound(lst, i, num_of_candidates):
		if i == 0:
			return 0
		i -= 1
	if i == len(lst) -1:
		lst[len(lst) -1] = lst[len(lst) -1] + 1
		return 1
	return reset(lst, i+1)


def increment_old(lst, num_of_candidates):
	for i in range(len(lst), -1, -1):
		if upper_bound(lst, i, num_of_candidates):
			if i == 0:
				return False
			reset(lst, i)
		else:
			lst[len(lst)] = lst[len(lst)] + 1
	return True


##########################
## gready approximation ##
##########################

def gready_cover(lst_of_candidates, genes_lst, k):
	'''aproximate minimal sigma pi(1-cleaving_prob(candidate, gene)), meaning the expected number of genes that won't be cleaved
	lst_of_candidates: sorted
	'''
	#group = []
	if len(group) >0:
		best_score = calculate_score(lst_of_candidates, group, genes_lst)
	else:
		best_score = len(genes_lst)
	for n in range(k):
		#find the that minimize the score, no repetitions.
		for candidate in range(len(lst_of_candidates)):
			if (not candidate in group):
				#break
				#continue
				#print(group)
				temp_score = calculate_score(lst_of_candidates, group + [candidate], genes_lst)
				if temp_score < best_score: #update
					best_score = copy.copy(temp_score) #in the end, the best score will be updated here
					temp_best_candidate = copy.copy(candidate)
		if temp_best_candidate not in group:
			group.append(temp_best_candidate)
			best_score = best_score
		if best_score == 0.0:
			break
	return best_score, group
	
def reverse_gready_cover(lst_of_candidates, genes_lst, k, group):
	'''aproximate minimal sigma pi(1-cleaving_prob(candidate, gene)), meaning the expected number of genes that won't be cleaved
	start with a the resulted set from the fuzzy set cover
	lst_of_candidates: sorted
	'''
	for n in range(len(group)-k):
		best_score_per_stage = len(genes_lst)
		#find the that minimize the score, no repetitions.
		for candidate in group: #which candidate to remove?
			if (candidate in group):
				temp_group =  copy.copy(group)
				temp_group.remove(candidate)
				temp_score = calculate_score(lst_of_candidates, temp_group, genes_lst)
				if temp_score < best_score_per_stage: #update
					best_score_per_stage = copy.copy(temp_score) #in the end, the best score will be updated here
					candidate_to_remove = copy.copy(candidate)
		if candidate_to_remove in group:
			group.remove(candidate_to_remove)
			best_score = best_score_per_stage #update the score
		if best_score == 0.0:
			break
	return best_score, group




########################################### thr covers ########################

def test_thr_best_only_amount(candidates_path, thr = 0.426036):
	candidates_lst = pickle.load(open(candidates_path, "rb"))
	return thr_best_candidate(candidates_lst, thr)	
	
def mul_candidate(candidates_lst, genesList, thr):
	'''
	:param candidates_lst: fi
	:param thr:
	:return: nds the best candidate by score = mult(genes)
	'''
	best_score = 0
	best_cut_expectation = 0
	best_candidate = None
	best_amount = 0
	for candidate in candidates_lst:
		amount, mult_score = 0 , 1
		for gene in genesList:
			if gene not in candidate.genes_score_dict:
				score = 0
			else:
				score = candidate.genes_score_dict[gene]
		#for gene, score in candidate.genes_score_dict.items():
			if score > thr:
				amount += 1
			mult_score *= score
		#print(mult_score)
		if mult_score > best_score or (mult_score == best_score and amount > best_amount):
			best_amount, best_cut_them_all, best_candidate = amount, mult_score, candidate
	return best_candidate, best_amount, best_cut_them_all

def test_mult(candidates_path, genesList_path ,thr = 0.426036):
	candidates_lst = pickle.load(open(candidates_path, "rb"))
	genesList = pickle.load(open(genesList_path, "rb"))
	return mul_candidate(candidates_lst, genesList,thr)

def find_fraction(candidates_path, thr):
	amount = 0
	lst = pickle.load(open(candidates_path, 'rb'))
	#count num of above thr:
	if len(lst) == 0:
		return None, None
	for gene, score in lst[0].genes_score_dict.items():
		if score >= thr:
			amount += 1
	return lst[0], amount
############################################
######           test        ###############
############################################
def test_score():
	path = "D:\Lab\\test5"
	genes_lst = pickle.load(open(path + "\\genesNames.p", 'rb'))
	#print(genes_lst)
	lst_of_candidates = pickle.load(open(path + "\\res_in_lst.p", 'rb'))
	r = calculate_score(lst_of_candidates, [0], genes_lst)
	print(r)

def test_fc(k):
	#cdef int k
	#k = l
	path = "D:\Lab\\test5"
	genes_lst = pickle.load(open(path + "\\genesNames (2).p", 'rb'))
	lst_of_candidates = pickle.load(open(path + "\\res_in_lst (2).p", 'rb'))
	r = full_cover_V0(lst_of_candidates, genes_lst, k)
	print(r)


def test_gready_cover(k):
	path = "D:\Lab\\test5"
	genes_lst = pickle.load(open(path + "\\genesNames (2).p", 'rb'))
	print(genes_lst)
	lst_of_candidates = pickle.load(open(path + "\\res_in_lst (2).p", 'rb'))
	r = gready_cover(lst_of_candidates, genes_lst, k)
	print(r)

def test_CC():
	#print(call_CplexCovers("/groups/itay_mayrose/galhyams/test6/res_in_lstHOM03D000566.p", "/groups/itay_mayrose/galhyams/test6/genesNamesHOM03D000566.p", 0.66, "SC"))
	print(bounded_cover("/groups/itay_mayrose/galhyams/test6/res_in_lstHOM03D000566.p", "/groups/itay_mayrose/galhyams/test6/genesNamesHOM03D000566.p", 5))


