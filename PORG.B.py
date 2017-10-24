#July 23 2017
import time
import sys
import copy
import itertools
import cPickle as pickle
import matplotlib.pyplot as plt

def read_data(path):
	with open(path + ".pickle", "r") as fp:
		obj = pickle.load(fp)
	print len(obj), path + " elements load over.", time.ctime()
	return obj

def read_data_b(path):
	with open(path + ".pickle", "rb") as fp:
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
	if age < 60:
		label = "ABC"
	else:
		label = "D"
	return label

def which_beta_B2(line, GSM_info):
	gsm = line[0]
	age = GSM_info[gsm][1]
	#print age
	if age < 40:
		label = "AB"
	else:
		label = "CD"
	return label

def merge():
	keys1 = read_data("mixed.sample.C2")
	keys2 = read_data("mixed.sample.D2")
	#keys3 = read_data("mixed.sample.C2")
	
	for key in keys2:
		keys1.append(key)
	print len(keys1), keys1[-1]

	with open("mixed.sample.CD2.pickle", "w") as fp:
		pickle.dump(keys1, fp)

def mix():
	beta = "AB"
	remain = []
	for line in matrix:
		beta_1 = which_beta(line, GSM_info)
		if beta_1 == beta:
			remain.append(line)
	N = len(remain)
	print beta, N

	ms = []
	for i in range(0, N, 2):
		j = (i + 1) % N
		s1, s2= remain[i], remain[j]

		l, pre = lcs_len(s1, s2)
		lcs = get_lcs_nr(pre, s1, len(s1), len(s2))
		new =  list(lcs)
		new.insert(0, "D2")
		ms.append(new)
		print i, j, l, time.ctime()

	print len(ms)
	with open("mixed.sample.D2.pickle", "w") as fp:
		pickle.dump(ms, fp)

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

def DFS(m1, s1, path, N, beta, mini_conf):
	global deep
	deep += 1
	print deep, path, len(s1), time.ctime()

	for j in range(path[-1] + 1, N):
		s2 = m1[j]

		l, pre = lcs_len(s1, s2)
		lcs = get_lcs_nr(pre, s1, len(s1), len(s2))
		new =  list(lcs)
		new.insert(0, "x")
		if new == s1:
			temp = DFS(m1, new, path + [j], N, beta, mini_conf)
			if temp != None:
				return temp
		else:
			if l > 12:
				temp = DFS(m1, new, path + [j], N, beta, mini_conf)
				if temp != None:
					return temp
			elif l > 1:
				# get sup & conf of rule
				#count = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
				#count = {'ABC': 0, 'D': 0}
				count = {'AB': 0, 'CD': 0}
				for line in matrix:
					if isSubsequence(lcs, line):
						#beta_1 =  which_beta(line, GSM_info)
						#beta_1 =  which_beta_B(line, GSM_info)
						beta_1 =  which_beta_B2(line, GSM_info)
						count[beta_1] += 1
				sup = count[beta]
				conf = 1.0  * sup / sum(count.values())
				ceiling = 1.0  * subsum[beta] / (subsum[beta] + sum(count.values()) - sup)
				temp = [beta, lcs, sup, conf]
				#print temp

				if mini_conf <= conf:
					print l, sup, conf, ceiling, count, time.ctime()
					return temp
				elif mini_conf <= ceiling:
					#print conf, count
					temp = DFS(m1, new, path + [j], N, beta, mini_conf)
					if temp != None:
						return temp
				else:
					#print j, "ceiling fail", time.ctime()
					pass
			else:
				#print j, "empty", time.ctime()
				pass
	return None

# ubr: [beta, alpha, sup, conf]
def UBR_DFS():
	remain = read_data("mixed.sample.AB2")
	current = copy.deepcopy(remain)
	ubr, beta, mini_conf = [], "AB", 0.85

	'''
	for e in remain:
		print len(e), e
	return
	'''
	while(True):
		N = len(current)
		if N < 2:
			break
		print "current", N, time.ctime()
		
		global deep
		deep = 0
		temp= DFS(current, current[0], [0], N, beta, mini_conf)

		if temp == None:
			print "pop", time.ctime()
			del current[0]
		else:
			ubr.append(temp)
			lcs = temp[1]
			current = [line for line in current if not isSubsequence(lcs, line)]
			#remain = [line for line in remain if not isSubsequence(lcs, line)]
			#current = copy.deepcopy(remain)
		
	print len(ubr)
	with open("ubr.AB2.pickle", "w") as fp:
		pickle.dump(ubr, fp)

# lbr: [beta, alpha, sup, conf]
def LBR():
	ubr = read_data("ubr.D2")
	lbr = []

	for rule in ubr:
		beta, alpha, sup, conf = rule[0], rule[1], rule[2], rule[3]
		print rule
		flag = False

		for i in range(2, len(alpha) + 1):
			if flag == True:
				break
			print i, time.ctime()
			l = list(itertools.combinations(alpha, i))
			for e in l:
				# get sup & conf of rule
				count = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
				for line in matrix:
					if isSubsequence(e, line):
						beta_1 = which_beta(line, GSM_info)
						count[beta_1] += 1
				lbr_sup = count[beta]
				lbr_conf = 1.0  * lbr_sup / sum(count.values())
				if lbr_conf >= conf:
					temp = [beta, e, lbr_sup, lbr_conf]
					lbr.append(temp)
					
					print temp
					flag = True
					break
			
	with open("lbr.D2.pickle", "w") as fp:
		pickle.dump(lbr, fp)	

def classifier_builder():
	rules = read_data("ubr.ABC4.0.85")
	rules_2 = read_data_b("lbr.D2b")
	for e in rules_2:
		rules.append(e)

	#sort
	rules.sort(key = lambda x: (x[-1], x[-2]), reverse = True)
	print rules[0], "\n", rules[-1]
	print time.ctime()

	#build
	classifier = []
	hit = []
	remain = copy.deepcopy(matrix)
	for rule in rules:
		beta, alpha = rule[0], rule[1]
		#hit = [line for line in remain if isSubsequence(alpha, line) and which_beta(line, GSM_info) == beta]
		#hit = [line for line in remain if isSubsequence(alpha, line) and which_beta_B(line, GSM_info) == beta]
		hit = [line for line in remain if isSubsequence(alpha, line) and which_beta_B2(line, GSM_info) == beta]
		print len(hit), 0.05 * subsum[beta]
		if len(hit) >= 0.05 * subsum[beta]:
			for e in hit:
				remain.remove(e)
			print rule
			classifier.append(rule)
	#default
	#count = {'ABC': 0, 'D': 0}
	count = {'AB': 0, 'CD': 0}
	for line in remain:
		#beta = which_beta_B(line, GSM_info)
		beta = which_beta_B2(line, GSM_info)
		count[beta] += 1
	print len(remain), count

	with open("classifier.ABCD.pickle", "w") as fp:
		pickle.dump(classifier, fp)

def accuracy(classifier):
	tp, tn, fp, fn = 0.0, 0.0, 0.0, 0.0
	default = "D"
	positive = "ABC"
	for line in matrix:
		#real =  which_beta(line, GSM_info)
		#real =  which_beta_B(line, GSM_info)
		real =  which_beta_B2(line, GSM_info)
		hit = False
		for rule in classifier:
			beta, alpha = rule[0], rule[1]
			if isSubsequence(alpha, line):
				if real == beta:
					if real == positive:
						tp += 1
					else:
						tn += 1
				else:
					if real == positive:
						fn += 1
					else:
						fp += 1
				hit = True
				break
		if not hit:
			if real == default:
				if real == positive:
					tp += 1
				else:
					tn += 1
			else:
				if real == positive:
					fn += 1
				else:
					fp += 1
	precision, recall, f1, acc, precision2, recall2, f12 = 0,0,0,0,0,0,0
	if (tp + fp) != 0:
		precision = tp / (tp + fp)
	recall = tp / (tp + fn)
	if (precision + recall) != 0:
		f1 = 2 * precision * recall / (precision + recall)
	acc = (tp + tn) / (tp + tn + fp + fn)
	temp = (precision, recall, f1, acc)
	
	precision2 = tn / (tn + fn)
	recall2 = tn / (tn + fp)
	f12 = 2 * precision2 * recall2 / (precision2 + recall2)
	temp = (precision, recall, f1, acc, precision2, recall2, f12)
	
	return temp

def show_results():
	classifier = read_data("classifier.ABCD")
	#classifier = read_data_b("classifier.D2b")
	index = []
	for i in range(1, len(classifier) + 1):
		temp = accuracy(classifier[0:i])
		index.append(temp)
		print i, classifier[i - 1]
	print temp
	plt.plot(range(1, len(index)+1), [e[0] for e in index], "ko")
	plt.plot(range(1, len(index)+1), [e[1] for e in index], "ko")
	plt.plot(range(1, len(index)+1), [e[2] for e in index], "ko")
	plt.plot(range(1, len(index)+1), [e[3] for e in index], "ko")
	plt.plot(range(1, len(index)+1), [e[0] for e in index], "-", label="precision")
	plt.plot(range(1, len(index)+1), [e[1] for e in index], "-", label="recall")
	plt.plot(range(1, len(index)+1), [e[2] for e in index], "-", label="f1")
	plt.plot(range(1, len(index)+1), [e[3] for e in index], "-", label="acc")
	
	plt.plot(range(1, len(index)+1), [e[4] for e in index], "-", label="precision2")
	plt.plot(range(1, len(index)+1), [e[5] for e in index], "-", label="recall2")
	plt.plot(range(1, len(index)+1), [e[6] for e in index], "-", label="f12")
	
	plt.title("Age.bi.train " + time.ctime())
	plt.grid(True)
	plt.legend()
	plt.show()

def rcbt():
	rules1 = read_data("ubr.ABC4.0.85")
	rules2 = read_data_b("lbr.D2")
	
	tp, tn, fp, fn = 0.0, 0.0, 0.0, 0.0
	positive = "ABC"

	for e in rules1:
		print e
	print ""
	for e in rules2:
		print e

	for line in matrix:
		real =  which_beta_B(line, GSM_info)
		s1, count = 0, 0
		for rule in rules1:
			beta, alpha, sup, conf = rule[0], rule[1], rule[2], rule[3]
			if isSubsequence(alpha, line):
				#print conf , sup , subsum[beta]
				#s1 += conf #* sup / subsum[beta]
				#count += 1
				s1 = max(s1, conf)
		#s1 = s1 / count
		s2, count = 0, 0
		for rule in rules2:
			beta, alpha, sup, conf = rule[0], rule[1], rule[2], rule[3]
			if isSubsequence(alpha, line):
				#print conf , sup , subsum[beta]
				#s2 += conf #* sup / subsum[beta]
				#count += 1
				s2 = max(s2, conf)
		#if count != 0:
			#s2 = s2 / count
		print s1, s2
		if s1 > s2:
			pre = "ABC"
		else:
			pre = "D"
		
		if real == positive and pre == positive:
			tp += 1
		elif real == positive and pre != positive:
			fn += 1
		elif real != positive and pre == positive:
			fp += 1
		elif real != positive and pre != positive:
			tn += 1
			if real == positive:
				tp += 1
			else:
				tn += 1
		else:
			if real == positive:
				fn += 1
			else:
				fp += 1
	
	precision, recall, f1, acc, precision2, recall2, f12 = 0,0,0,0,0,0,0
	if (tp + fp) != 0:
		precision = tp / (tp + fp)
	recall = tp / (tp + fn)
	if (precision + recall) != 0:
		f1 = 2 * precision * recall / (precision + recall)
	acc = (tp + tn) / (tp + tn + fp + fn)
	temp = (precision, recall, f1, acc)
	
	precision2 = tn / (tn + fn)
	recall2 = tn / (tn + fp)
	f12 = 2 * precision2 * recall2 / (precision2 + recall2)
	temp = (precision, recall, f1, acc, precision2, recall2, f12)
	print temp

def age_distribution():
	count = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
	for line in matrix:
		beta = which_beta(line, GSM_info)
		count[beta] += 1
	count["AB"] = count["A"] + count["B"]
	count["CD"] = count["C"] + count["D"]
	count["ABC"] = count["AB"] + count["C"]
	return count

if __name__ == "__main__":
	print "Start.", time.ctime()
	#mix()
	#merge()
	
	matrix = read_data("age_train_802")
	#matrix = read_data("age_test_198")
	GSM_info = read_data("GSM_info")
	matrix = convert(matrix)
	subsum = age_distribution()
	print subsum
	
	deep = 0
	UBR_DFS()
	#LBR()
	#classifier_builder()
	#show_results()
	#rcbt()

	print "End.", time.ctime()