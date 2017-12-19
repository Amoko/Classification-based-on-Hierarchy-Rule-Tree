#July 23 2017
import time
import sys
import copy
import random
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

def sex_beta(line):
	e = line[0][-1]
	return e

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
	beta = "M"
	remain = []
	for line in matrix:
		beta_1 = sex_beta(line)
		if beta_1 == beta:
			remain.append(line)
	N = len(remain)
	print beta, N
	
	#remain = read_data("mixed.sample.F2")
	#N = len(remain)
	ms = []
	for i in range(0, N, 3):
		j = (i + 1) % N
		k = (i + 2) % N
		s1, s2 = remain[i], remain[j]
		s3 = remain[k]

		l, pre = lcs_len(s1, s2)
		lcs = get_lcs_nr(pre, s1, len(s1), len(s2))
		new =  list(lcs)
		info = [s1[0], s2[0]]
		#info = s1[0] + s2[0]
		new.insert(0, info)

		l, pre = lcs_len(new, s3)
		lcs = get_lcs_nr(pre, new, len(new), len(s3))
		new =  list(lcs)
		info = info +[s3[0]]
		new.insert(0, info)

		ms.append(new)
		print i, j, l, info, time.ctime()

	print len(ms)
	with open("mixed.sample.M3.pickle", "w") as fp:
		pickle.dump(ms, fp)



#convert number to position by order of pos' num
def num2pos(line):
	s = line
	gsm = s[0]
	s[0] = -float("inf")
	ns = zip(s, range(len(s)))
	ns.sort(key = lambda x: x[0])
	#print ns[:5], ns[-5:]
	ns[0] = gsm
	for i in range(1, len(ns)):
		ns[i] = ns[i][1]
	#print ns[:5], ns[-5:]
	return ns

def convert(matrix):
	for i in range(len(matrix)):
		matrix[i] = num2pos(matrix[i])
	print "Convert over.", time.ctime()

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
	'''
	len1 = len(s1)
	len2 = len(s2)
	i, j = 0, 0
	while(i < len1 and j < len2):
		if s1[i] == s2[j]:
			i += 1
		j += 1
	return i == len1
	'''
	a, b = s1[0], s1[1]
	if s2.index(a) < s2.index(b):
		return True
	else:
		return False

def accuracy2(matrix, classifiers):
	tp, tn, fp, fn = 0.0, 0.0, 0.0, 0.0
	positive = "M"
	for line in matrix:
		#real =  which_beta(line, GSM_info)
		#real =  which_beta_B(line, GSM_info)
		real =  sex_beta(line)
		score = {'M': 0, 'F': 0}
		for classifier in classifiers:
			default = classifier[-1]
			classifier = classifier[:-1]
			hit = False
			for rule in classifier:
				beta, alpha = rule[0], rule[1]
				if isSubsequence(alpha, line):
					pre = beta
					hit = True
					break
			if not hit:
				pre = default
			score[pre] += 1
		pre = max(score, key=score.get)
		
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

	#index
	precision, recall, f1, acc, precision2, recall2, f12 = 0,0,0,0,0,0,0
	if (tp + fp) != 0:
		precision = tp / (tp + fp)
	recall = tp / (tp + fn)
	if (precision + recall) != 0:
		f1 = 2 * precision * recall / (precision + recall)
	acc = (tp + tn) / (tp + tn + fp + fn)
	temp = (precision, recall, f1, acc)
	
	if (tn + fn) != 0:
		precision2 = tn / (tn + fn)
	recall2 = tn / (tn + fp)
	if (precision2 + recall2) != 0:
		f12 = 2 * precision2 * recall2 / (precision2 + recall2)
	temp = (precision, recall, f1, acc, precision2, recall2, f12)
	
	return temp

def origin():
	# sd [[3, 42], [1, 54], [2, 60]]
	# the first is origin index, thesecond is na count
	

	classifier = read_data("rule.nov.21th")
	for rule in classifier:
		print rule
		beta, alpha, sup, conf = rule[0], rule[1], rule[2], rule[3]
		new = []
		for e in alpha:
			new.append(sd[e][0])
		print new

def sex_distribution(matrix):
	count = {'M': 0, 'F': 0}
	for line in matrix:
		beta = sex_beta(line)
		count[beta] += 1
	return count

def get_w(matrix):
	pass
	L = len(matrix[0])
	wp = [0] * L
	wn = [0] * L
	for line in matrix:
		beta_1 = sex_beta(line)
		if beta_1 == "M":
			for i in range(1, L):
				wp[line[i]] += i
				pass
		else:
			for i in range(1, L):
				wn[line[i]] += i
				pass
	#[key, value]
	#[20, 50, 30, 40, 60]
	#[[1, 20], [2, 50], [3, 30], [4, 40], [5, 60]]
	#[[1, 20], [3, 30], [4, 40], [2, 50], [5, 60]]
	#[1, 3, 4, 2, 5]
	#[[1, 0], [3, 1], [4, 2], [2, 3], [5, 4]]
	#[[1, 0], [2, 3], [3, 1], [4, 2], [5, 4]]
	wp = zip(range(len(wp)), wp)
	del wp[0]
	wp.sort(key=lambda x: x[1])
	wp = [e[0] for e in wp]
	wp = zip(wp, range(len(wp)))
	wp.sort(key=lambda x: x[0])

	wn = zip(range(len(wn)), wn)
	del wn[0]
	wn.sort(key=lambda x: x[1])
	wn = [e[0] for e in wn]
	wn = zip(wn, range(len(wn)))
	wn.sort(key=lambda x: x[0])

	w = []
	for i in range(len(wp)):
		if wp[i][0] == wn[i][0]:
			w.append([wp[i][0], wp[i][1] - wn[i][1]])
	w.sort(key=lambda x: x[1])

	return w

def pop(matrix, k):
	matrix = list(matrix)
	matrix = [e for e in matrix if matrix.index(e)%3!=k]
	#N = int(len(matrix) * 0.618)
	#matrix = random.sample(matrix, N)
	print len(matrix)
	classifier = []

	last = [0, 0, 0, 0]
	while len(matrix) > 20:
		w = get_w(matrix)
		print w[0], w[-1], len(matrix)

		pairs = []
		for i in range(50):
			#print i
			j = - (i + 1)
			e = [w[i][0], w[j][0]]

			#hit, nohit = {'AB': 0, 'CD': 0}, {'AB': 0, 'CD': 0}
			#hit, nohit = {'ABC': 0, 'D': 0}, {'ABC': 0, 'D': 0}
			hit, nohit = {'M': 0, 'F': 0}, {'M': 0, 'F': 0}
			for line in matrix:
				beta = sex_beta(line)
				if isSubsequence(e, line):
					hit[beta] += 1
				else:
					nohit[beta] += 1
			# hit
			beta_h = max(hit, key=hit.get)
			sup_h = hit[beta_h]
			if sum(hit.values()) == 0:
				conf_h = 0
			else:
				conf_h = 1.0  * sup_h / sum(hit.values())
			z = copy.deepcopy(e)
			rule_h = [beta_h, z, sup_h, conf_h]
			# nohit
			beta_no = max(nohit, key=nohit.get)
			sup_no = nohit[beta_no]
			if sum(nohit.values()) == 0:
				conf_no = 0
			else:
				conf_no = 1.0  * sup_no / sum(nohit.values())
			e.reverse()
			rule_no = [beta_no, e, sup_no, conf_no]
			#judge
			if conf_h > conf_no:
				win, loser = rule_h, rule_no
			else:
				win, loser = rule_no, rule_h
			
			if win[-2] > 20:
				pairs.append([win, loser])

		if pairs != []:
			pairs.sort(key=lambda x:x[0][-1], reverse=True)
			if pairs[0][0][-1] > last[-1]:
				top = pairs[0][0]
				last = pairs[0][1]
				classifier.append(top)
				matrix = [line for line in matrix if not isSubsequence(top[1], line)]
				print len(classifier), top, time.ctime()
			else:
				classifier.append(last)
				break
		else:
			break

	remain = {'M': 0, 'F': 0}
	for line in matrix:
		beta = sex_beta(line)
		remain[beta] += 1
	default = max(remain, key=remain.get)
	classifier.append(default)
	print remain, default
	
	return classifier	
	#with open("rule.nov.F.20.pickle", "w") as fp:
		#pickle.dump(classifier, fp)	

def acc():
	css = read_data("pop.css")
	j = 0
	t = []
	sample = sampling("sex_test", 2000)
	convert(sample)
	for i in [1000, 2000, 4000, 6000, 8000, 10000]:
		#i = int(i * 0.2)
		for c in css[j]:
			print c
		print ""
		acc = accuracy2(sample, css[j])[3]
		t.append(acc)
		j += 1
	
	with open("acc.pop.pickle", "w") as fp:
		pickle.dump(t, fp)

	#plt.ylim((0.5, 1))
	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t, "-")
	plt.ylabel("Accuracy")
	plt.xlabel("# of samples")
	plt.grid(True)
	plt.show()

def sampling(data, n):
	matrix = read_data(data)
	subsum = sex_distribution(matrix)
	
	expect = {}
	expect["M"] = n * 1.0 * subsum["M"] / len(matrix)
	expect["F"] = n * 1.0 * subsum["F"] / len(matrix)

	count = {}
	m = []
	for line in matrix:
		beta = sex_beta(line)
		if count.get(beta):
			if count[beta] <= expect[beta]:
				m.append(line)
				count[beta] += 1
		else:
			m.append(line)
			count[beta] = 1
		if len(m) >= n:
			break
	print count, "sampling over.", time.ctime()
	return m

def timing():
	css = []
	t = []
	for  i in [1000, 2000, 4000, 6000, 8000, 10000]:
		i = int(i * 0.8)
		sample = sampling("sex_train", i)
		convert(sample)
	
		start  = time.clock()
		#rule = nov(sample)
		classifiers = bagging(sample)
		elapsed = time.clock() - start
	
		t.append(elapsed)
		css.append(classifiers)
	
	with open("pop.css.pickle", "w") as fp:
		pickle.dump(css, fp)
	with open("pop.timing.pickle", "w") as fp:
		pickle.dump(t, fp)	

	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t, "-")
	plt.ylabel("CPU time(second)")
	plt.xlabel("# of samples")
	plt.grid(True)
	plt.show()

def bagging(matrix):
	classifiers = []
	for i in range(3):
		print i, len(matrix)
		
		c = pop(matrix, i)
		classifiers.append(c)

	return classifiers
	'''
	matrix = read_data("sex_test_177")
	convert()
	acc = accuracy2(classifiers)
	print acc
	'''

def test():
	css = read_data("bagging.css.timing")
	for cs in css:
		print len(cs)
		for c in cs:
			print len(c)


if __name__ == "__main__":
	print "Start.", time.ctime()

	#timing()
	acc()
	
	#nov()
	#bagging()
	#test()

	print "End.", time.ctime()