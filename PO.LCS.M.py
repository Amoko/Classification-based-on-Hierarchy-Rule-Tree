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

def DFS(m1, s1, start, N, beta, mini_conf):
	global deep
	deep += 1
	print deep, N, len(matrix)
	for j in range(start, N):
		s2 = m1[j]

		l, pre = lcs_len(s1, s2)
		lcs = get_lcs_nr(pre, s1, len(s1), len(s2))
		new =  list(lcs)
		new.insert(0, "x")

		if l > 12:
			temp = DFS(m1, new, start + 1, N, beta, mini_conf)
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
			ceiling = 1.0  * subsum / (subsum + sum(count.values()) - sup)
			temp = [beta, lcs, sup, conf]

			if mini_conf <= conf:
				print l, sup, conf, ceiling, count, time.ctime()
				return temp
			elif mini_conf <= ceiling:
				temp = DFS(m1, new, start + 1, N, beta, mini_conf)
				if temp != None:
					return temp
			else:
				pass
		else:
			pass
	return None

# ubr: [beta, alpha, sup, conf]
def UBR_DFS():
	remain = read_data("mixed.sample.A2")
	current = copy.deepcopy(remain)
	ubr, beta, mini_conf = [], "A", 0.5

	while(True):
		N = len(current)
		if N < 2:
			break
		print "current", N, time.ctime()
		
		global deep
		deep = 0
		temp= DFS(current, current[0], 1, N, beta, mini_conf)

		if temp == None:
			print "pop", time.ctime()
			del current[0]
		else:
			ubr.append(temp)
			lcs = temp[1]
			remain = [line for line in remain if not isSubsequence(lcs, line)]
			current = copy.deepcopy(remain)
		
	print len(ubr)
	with open("ubr.A2.pickle", "w") as fp:
		pickle.dump(ubr, fp)

def  mix():
	'''
	beta = "D"
	remain = []
	for line in matrix:
		beta_1 = which_beta(line, GSM_info)
		if beta_1 == beta:
			remain.append(line)
	'''
	remain = read_data("mixed.sample.CD2")
	N = len(remain)
	print N

	ms = []
	for i in range(0, N, 2):
		j = (i + 1) % N
		s1, s2= remain[i], remain[j]

		l, pre = lcs_len(s1, s2)
		lcs = get_lcs_nr(pre, s1, len(s1), len(s2))
		new =  list(lcs)
		new.insert(0, "CD4")
		ms.append(new)
		print i, j, l, time.ctime()

	print len(ms)
	with open("mixed.sample.CD4.pickle", "w") as fp:
		pickle.dump(ms, fp)

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

def age_distribution(goal):
	count = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
	for line in matrix:
		beta = which_beta(line, GSM_info)
		count[beta] += 1
	print count
	print sum(count.values()), count[goal]
	return count[goal]

def merge():
	keys1 = read_data("mixed.sample.A2")
	keys2 = read_data("mixed.sample.B2")
	keys3 = read_data("mixed.sample.C2")
	
	for key in keys2:
		keys1.append(key)
	print len(keys1), keys1[-1]

	for key in keys3:
		keys1.append(key)
	print len(keys1), keys1[-1]

	with open("mixed.sample.ABC2.pickle", "w") as fp:
		pickle.dump(keys1, fp)

def test():
	lbr = read_data("ubr.A2")

	s = 0
	for line in matrix:
		real =  which_beta(line, GSM_info)
		if real  == "A":
			for rule in lbr:
				#print rule
				beta, alpha, sup, conf = rule[0], rule[1], rule[2], rule[3]
				if isSubsequence(alpha, line):
					s += 1
					break
	print s

if __name__ == "__main__":
	print "Start.", time.ctime()
	mix()
	#merge()
	'''
	matrix = read_data("age_train_802")
	GSM_info = read_data("GSM_info")
	matrix = convert(matrix)
	subsum = age_distribution("A")
	
	deep = 0
	UBR_DFS()
	test()
	'''
	#LBR()
	#classifier()
	print "End.", time.ctime()