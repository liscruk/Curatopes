#!/usr/bin/python3

## This data is loaded here as a test, later this function should be called from the curaLib liberary

import pandas as pd
import numpy as np

bindData = pd.read_table("/mnt/e/Work/Projects/pancancercuratopes/Data/Samples/OutPut",header = None)

#Random TPM data to test f3 and f4
bindData['testTPM'] = np.random.randint(1,101,bindData.shape[0])

def calcGPIE(data):

	bindAffTotal = data[4].tolist()
	immunoGen = data[5].tolist()


	def exprNorm(Expr):
		pass
		

	#IC50
	def f1(pepIC50):

		# Calculate the score of an individual peptide by taking all predicted ones into account
		# c1 calcualtes the IC50 coefficient for the gpie score and is returned.
		c1 = (pepIC50 - max(bindAffTotal))/(min(bindAffTotal) - max(bindAffTotal))

		return c1

	#Provide ImmunoInput as list
	def f2(pepImmuno):

		# Calculate the normalized immunogenicty score of a peptide c2.
		c2 = (pepImmuno - min(immunoGen))/(max(immunoGen) - min(immunoGen))

		return c2

	# Median transcript expression of the source gen of a peptide
	def f3(pepGenExpr):

		if pepGenExpr<=100:
			pepGenExpr = pepGenExpr/100
		else:
			pepGenExpr = 1

		return pepGenExpr

	data['f1'] = data[4].apply(f1)
	data['f2'] = data[5].apply(f2)
	data['f3'] = data['testTPM'].apply(f3)
	# data['f4'] = data[].apply(f4)

	
	data['gPie'] = 100 * data['f1'] * data['f2'] * data['f3']

	data = data.drop(columns=['f1','f2','f3'])

	# Placeholder functions kommen hier für expression und expressions index
	
	return data


result = calcGPIE(bindData)

print(result)