#July 23 2017
import time
import sys
import copy
import itertools
import cPickle as pickle

def read_data(path):
	with open(path + ".pickle", "rU") as fp:
		obj = pickle.load(fp)
	print len(obj), path + " elements load over.", time.ctime()
	return obj

def which_beta(line, GSM_info):
	gsm = line[0]
	age = GSM_info[gsm][1]
	#print age
	if age < 20:
		label = "A"
	elif age < 40:
		label = "B"
	elif age < 60:
		label = "C"
	else:
		label = "D"
	return label

def which_beta_B(line, GSM_info):
	gsm = line[0]
	age = GSM_info[gsm][1]
	#print age
	if age < 40:
		label = "AB"
	else:
		label = "CD"
	return label

#convert number to position by order of pos' num
def num2pos(line):
	s = line
	gsm = s[0]
	s[0] = float("inf")
	ns = zip(s, range(len(s)))
	ns.sort(key = lambda x: x[0], reverse = True)
	#print ns[:5], ns[-5:]
	ns[0] = gsm
	for i in range(1, len(ns)):
		ns[i] = ns[i][1]
	#print ns[:5], ns[-5:]
	return ns

def convert(matrix):
	m = []
	for line in matrix:
		nl = num2pos(copy.deepcopy(line))
		m.append(nl)
	print "Convert over.", time.ctime()
	return m

##############################
# longest common subsequence
def lcs_len(s1, s2):
	# notice that the first element won't be compare!
	# s1, s2 = "02579312", "035328"
	m, n = len(s1), len(s2)

	# DP table
	dp = [[0] * n for i in range(m)]
	# pre table
	# if top is equal to left, default by top
	pre = [[0] * n for i in range(m)]

	for i in range(1, m):
		#if  i%1000 == 0:
			#print i
		for j in range(1, n):
			if s1[i] == s2[j]:
				dp[i][j] = dp[i - 1][j - 1] + 1
			else:
				dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])
				if dp[i - 1][j] >= dp[i][j - 1]:
					pre[i][j] = "T"
				else:
					pre[i][j] = "L"
	return dp[-1][-1], pre

def get_lcs_inner(pre, s1, i, j, lcs):
	#print i, j
	if i == 0 or j == 0:
		return
	if pre[i][j] == 0:
		#print 0
		get_lcs_inner(pre, s1, i - 1, j - 1, lcs)
		lcs.append(s1[i])
	elif pre[i][j] == "T":
		#print "T"
		get_lcs_inner(pre, s1, i - 1, j, lcs)
	else:
		#print "L"
		get_lcs_inner(pre, s1, i, j - 1, lcs)

def get_lcs(pre, s1, i, j):
	lcs = []
	get_lcs_inner(pre, s1, i - 1, j - 1, lcs)
	return lcs

# non-recursion version
def get_lcs_nr(pre, s1, i, j):
	lcs = []
	i -= 1
	j -= 1
	while(i > 0 and j > 0):
		if pre[i][j] == 0:
			lcs.append(s1[i])
			i -= 1
			j -= 1
		elif pre[i][j] == "T":
			i -= 1
		else:
			j -= 1
	lcs.reverse()
	return lcs

def lcs_test():
	s1, s2 = "02579312", "035328"
	s1, s2 = "02579312", "03"
	m, n = len(s1), len(s2)
	l, pre = lcs_len(s1, s2)
	
	#for e in pre:
		#print e
	lcs1 = get_lcs(pre, s1, m, n)
	lcs2 = get_lcs_nr(pre, s1, m, n)
	print l, lcs1, lcs2
# longest common subsequence
##############################

# Is s1 a subsequence of s2?
def isSubsequence(s1, s2):
	len1 = len(s1)
	len2 = len(s2)
	i, j = 0, 0
	while(i < len1 and j < len2):
		if s1[i] == s2[j]:
			i += 1
		j += 1
	return i == len1

# ubr: [beta, alpha, sup, conf]
def UBR_Patch():
	ubr = []
	#iteration = ["A", "B", "C", "D"]
	iteration = ["D"]
	for beta in iteration:
		# make data
		m1 = []
		for line in matrix:
			beta_1 = which_beta(line, GSM_info)
			if beta_1 == beta:
				m1.append(line)
		N = len(m1)
		print beta, N

		#start search
		for i in range(N):
			print i
			s1 = m1[i]
			j = i

			while(True):
				j += 1
				j = j % len(m1)
				s2 = m1[j]

				l, pre = lcs_len(s1, s2)
				lcs = get_lcs_nr(pre, s1, len(s1), len(s2))
				s1 = list(lcs)
				s1.insert(0, "x")

				if l <= 10:
					if l < 2:
						print ""
						break
					# get sup & conf of rule
					count = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
					for line in matrix:
						if isSubsequence(lcs, line):
							beta_1 =  which_beta(line, GSM_info)
							count[beta_1] += 1
					sup = count[beta]
					conf = 1.0  * sup / sum(count.values())
					ceiling = 1.0  * N / (N + sum(count.values()) - sup)
					print l, sup, conf, ceiling, count, time.ctime()
					temp = [beta, lcs, sup, conf]
					
					#ubr.append(temp)
			#break
	
	with open("ubr.AD.pickle", "w") as fp:
		pickle.dump(ubr, fp)

def DFS(m1, s1, start, N, beta):
	mini_conf = 0.05
	for j in range(start, N):
		s2 = m1[j]

		l, pre = lcs_len(s1, s2)
		lcs = get_lcs_nr(pre, s1, len(s1), len(s2))
		new =  list(lcs)
		new.insert(0, "x")


		if l > 12:
			temp = DFS(m1, new, start + 1, N, beta)
			if temp != None:
				return temp
		elif l > 1:
			# get sup & conf of rule
			count = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
			for line in matrix:
				if isSubsequence(lcs, line):
					beta_1 =  which_beta(line, GSM_info)
					count[beta_1] += 1
			sup = count[beta]
			conf = 1.0  * sup / sum(count.values())
			ceiling = 1.0  * N / (N + sum(count.values()) - sup)
			temp = [beta, lcs, sup, conf]

			if mini_conf <= conf:
				print l, sup, conf, ceiling, count, time.ctime()
				return temp
			elif mini_conf <= ceiling:
				temp = DFS(m1, new, start + 1, N, beta)
				if temp != None:
					return temp
			else:
				pass
		else:
			pass
	return None

# ubr: [beta, alpha, sup, conf]
def UBR_DFS():
	ubr = []
	#iteration = ["A", "B", "C", "D"]
	iteration = ["A"]
	for beta in iteration:
		# make data
		remain = []
		for line in matrix:
			beta_1 = which_beta(line, GSM_info)
			if beta_1 == beta:
				remain.append(line)
		N = len(remain)
		print beta, N

		m1 = copy.deepcopy(remain)
		while(True):
			N = len(m1)
			if N < 10:
				break
			print N, time.ctime()
			
			temp= DFS(m1, m1[0], 1, N, beta)

			if temp == None:
				print m1[0]
				del m1[0]
			else:
				ubr.append(temp)
				lcs = temp[1]
				remain = [line for line in remain if not isSubsequence(lcs, line)]
				m1 = copy.deepcopy(remain)
			

	with open("ubr.A.pickle", "w") as fp:
		pickle.dump(ubr, fp)

def M2B():
	matrix = read_data("age_train_topk2")
	GSM_info = read_data("GSM_info")
	matrix = convert(matrix)
	ubr = read_data("ubr.25")
	
	ubrB = []
	for rule in ubr:
		beta, lcs = rule[0], rule[1]
		if beta == "A" or beta == "B":
			beta = "AB"
		else:
			beta = "CD"

		# get sup & conf
		count = {'AB': 0, 'CD': 0}
		for line in matrix:
			if isSubsequence(lcs, line):
				real =  which_beta_B(line, GSM_info)
				count[real] += 1
		print count
		sup = sum(count.values())
		conf = 1.0  * count[beta] / sup
	
		print beta, sup, conf, count
		temp = [beta, lcs, sup, conf]
		ubrB.append(temp)
		#break
	
	with open("ubrB.pickle", "w") as fp:
		pickle.dump(ubrB, fp)

def UBR_Select():
	ubr = read_data("ubr.M.25")

	ubrs = []
	'''
	#iteration = ["A", "B", "C", "D"]
	iteration = ["AB", "CD"]
	for beta in iteration:
		m = []
		for rule in ubr:
			pass
			if rule[0] == beta:
				m.append(rule)
		m.sort(key = lambda e: (e[-1], e[-2]), reverse = True)
		
		part = len(m) / 2
		temp = m[:part]
		print temp[-1]
		for e in temp:
			ubrs.append(e)
	'''
	for e in ubr:
		if e[-1] >= 0.7:
			ubrs.append(e)
	print "Select over.", time.ctime()

	with open("ubrs.M.25.pickle", "w") as fp:
		pickle.dump(ubrs, fp)

# find LBR from top to bottom, but it will be factorial explosion
def LBR_inner(matrix, alpha, sup):
	for e in alpha:
		na = copy.copy(alpha)
		na.remove(e)
		count = 0
		for line in matrix:
			if isSubsequence(na, line):
				count += 1
		if count == sup:
			inner(matrix, na, sup)
		else:
			print alpha, e
			return
	pass

# lbr: [beta, alpha, sup, conf]
def LBR():
	matrix = read_data("age_train_topk")
	GSM_info = read_data("GSM_info")
	matrix = convert(matrix)
	
	ubrc = read_data("ubrc")
	lbr = []
	for k in range(0, 20):
		rule = ubrc[k]
		print rule
		beta, alpha, sup, conf = rule[0], rule[1], rule[2], rule[3]
		nl = 0

		#LBR_inner(matrix, alpha, sup)
		for i in range(6, len(alpha) + 1):
			if nl == 2:
				break
			print i, time.ctime()

			l = list(itertools.combinations(alpha, i))
			for e in l:
				count = 0
				for line in matrix:
					if isSubsequence(e, line):
						count += 1
				if count == sup:
					temp = [beta, e, sup, conf]
					print temp
					lbr.append(temp)

					nl += 1
					if nl == 5:
						break
	with open("lbr1.pickle", "w") as fp:
		pickle.dump(lbr, fp)	

def classifier():
	matrix = read_data("age_train_topk2")
	GSM_info = read_data("GSM_info")
	matrix = convert(matrix)
	lbr = read_data("ubrs.M.25")

	#age distribution
	ad = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
	for line in matrix:
		beta = which_beta(line, GSM_info)
		ad[beta] += 1
	print ad, sum(ad.values())

	s = 0
	no = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
	for line in matrix:
		real =  which_beta(line, GSM_info)
		score = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
		hit = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
	
		
		for rule in lbr:
			#print rule
			beta, alpha, sup, conf = rule[0], rule[1], rule[2], rule[3]
			if isSubsequence(alpha, line):
				score[beta] += conf #* sup / ad[beta]
				hit[beta] += 1
		if max(hit.values()) == 0:
			no[real] += 1
			predict = "Q"
		else:
			predict = max(score, key = score.get)
		#print real, predict, hit, score
		if real  == predict:
			s += 1
	print 1.0 * s / len(matrix), no
	s += max(no.values())
	print 1.0 * s / len(matrix)

def age_distribution():
	matrix = read_data("age_test_topk2")
	GSM_info = read_data("GSM_info")
	
	count = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
	for line in matrix:
		beta = which_beta(line, GSM_info)
		count[beta] += 1
	print count, sum(count.values())

def merge():
	keys1 = read_data("ubr.25.AB")
	keys2 = read_data("ubr.25.CD")
	
	for key in keys2:
		keys1.append(key)
	print len(keys1), keys1[-1]

	with open("ubr.25.pickle", "w") as fp:
		pickle.dump(keys1, fp)

def test():
	pass
	lbr = read_data("ubr")
	count = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
	c = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
	for rule in lbr:
		print rule
		beta, alpha, sup, conf = rule[0], rule[1], rule[2], rule[3]
		count[beta] += conf * sup / 20
		c[beta] += 1

	print count
	print c

if __name__ == "__main__":
	print "Start.", time.ctime()
	
	matrix = read_data("age_train_topk2")
	GSM_info = read_data("GSM_info")
	matrix = convert(matrix)
	
	#lcs_test()
	#merge()
	#age_distribution()
	#UBR_Patch()
	UBR_DFS()
	#LBR()
	#classifier()
	print "End.", time.ctime()