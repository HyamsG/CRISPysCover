from __future__ import print_function

import sys
import numpy as np
import Candidate
import cplex
from cplex.exceptions import CplexError
import math


def populatebyrow(prob, my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows):
	prob.objective.set_sense(prob.objective.sense.minimize)

	prob.variables.add(obj=my_obj, lb=my_lb, ub=my_ub, types=my_ctype,
					   names=my_colnames)
	prob.linear_constraints.add(lin_expr=rows, senses=my_sense, rhs=my_rhs, names=my_rownames)

def mipex1(pop_method, my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows):

	try:
		my_prob = cplex.Cplex()
		if pop_method == "r":
			handle = populatebyrow(my_prob, my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows)
		my_prob.solve()
	except CplexError as exc:
		print(exc)
		return
	#numcols = my_prob.variables.get_num()
	#numrows = my_prob.linear_constraints.get_num()
	#slack = my_prob.solution.get_linear_slacks()
	#x = my_prob.solution.get_values()


	return my_prob
		

def CplexSetCover(candidatesLst, genesLst, thr):

	my_obj = [1.0 for i in range(len(candidatesLst))]
	my_ub = [1.0 for i in range(len(candidatesLst))]
	my_lb = [0.0 for i in range(len(candidatesLst))]
	my_ctype = "I" * len(candidatesLst)
	my_colnames = ["x" + str(i+1) for i in range(len(candidatesLst))]

	my_rhs = [1.0 for i in range(len(genesLst))]
	my_rownames = ["r" + str(i+1) for i in range(len(genesLst))]
	my_sense = "G" * len(genesLst) #upper bound

	rows = make_rows_set_cover(candidatesLst, genesLst, thr)

	return mipex1("r", my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows)#(pop_method, my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows)

def Cplex_fuzzy_set_cover(candidatesLst, genesLst, thr, FULL_COVER_THR = 0.99):
	#candidatesLst = candidatesLst[:500]#test

	thr = -math.log2(1-thr)#can be math.log as well
	my_obj = [1.0 for i in range(len(candidatesLst))]
	my_ub = [1.0 for i in range(len(candidatesLst))]
	my_lb = [0.0 for i in range(len(candidatesLst))]
	my_ctype = "I" * len(candidatesLst)
	my_colnames = ["x" + str(i+1) for i in range(len(candidatesLst))]
	#my_colnames = ["x1", "x2", "x3", "x4"]
	my_rhs = [thr for i in range(len(genesLst))]
	my_rownames = ["r" + str(i+1) for i in range(len(genesLst))]
	my_sense = "G" * len(genesLst) #change for upper bound somehow
	#input matrix for the problem
	rows = make_rows_fuzzy_set_cover(candidatesLst, genesLst, thr, FULL_COVER_THR)

	prob_obj = mipex1("r", my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows)#(pop_method, my_obj, my_lb, my_ub, my_colnames, my_ctype, my_rhs, my_sense, my_rownames, rows)
	return prob_obj

	
def make_rows_set_cover(candidatesLst, genesLst, thr):
	'''make a table of a_ij: does the i gene is covered by the j sgRNA'''
	rows = []
	for i in range(len(genesLst)): # a row
		row = [[],[]]#row = list(list())
		for j in range(len(candidatesLst)):
			if genesLst[i] in candidatesLst[j].genes_score_dict:
				if candidatesLst[j].genes_score_dict[genesLst[i]]> thr:
					row[0].append('x' + str(j+1))
					row[1].append(1.0)
					print("gene " + str(i)+ "is covered by candidate" + str(j))
		rows.append(row)
	return rows
	
def make_rows_fuzzy_set_cover(candidatesLst, genesLst, thr, FULL_COVER_THR):
	'''make a table of a_ij: does the i gene is covered by the j sgRNA'''
	rows = []
	for i in range(len(genesLst)): # a row
		row = [[],[]]#row = list(list())
		for j in range(len(candidatesLst)):
			if genesLst[i] in candidatesLst[j].genes_score_dict:
				row[0].append('x' + str(j+1))

				if candidatesLst[j].genes_score_dict[genesLst[i]] == 1:
					row[1].append(-math.log2(1-(FULL_COVER_THR)))
				else:
					row[1].append(-math.log2(1-(candidatesLst[j].genes_score_dict[genesLst[i]])))
				#print("gene " + str(i)+ "is covered by candidate" + str(j))
					#A[i][j] = True
		rows.append(row)
	return rows

