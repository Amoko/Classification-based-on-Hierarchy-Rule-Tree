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

def which_beta_D(line, GSM_info):
	gsm = line[0]
	age = GSM_info[gsm][1]
	#print age
	if age < 60:
		label = "ABC"
	else:
		label = "D"
	return label

def which_beta_CD(line, GSM_info):
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
	beta = "D"
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

def head():
	count = {}
	for line in matrix:
		gsm = line[0]
		if gsm in count:
			count[gsm] += 1
		else:
			count[gsm] = 1
	print len(count)

	mixed = read_data("mixed.sample.D2")
	for e in mixed:
		e[0] = []
		for line in matrix:
			if isSubsequence(e[1:], line):
				e[0].append(line[0])
		print e[0], time.ctime()

	with open("mixed.sample.D2.pickle", "w") as fp:
		pickle.dump(mixed, fp)

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
				count = {'ABC': 0, 'D': 0}
				#count = {'AB': 0, 'CD': 0}
				for line in matrix:
					if isSubsequence(lcs, line):
						#beta_1 =  which_beta(line, GSM_info)
						beta_1 =  which_beta_B(line, GSM_info)
						#beta_1 =  which_beta_B2(line, GSM_info)
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
	mixed = read_data("mixed.sample.D2")
	ubr, beta, mini_conf = [], "D", 0.8

	hit ={}
	for line in matrix:
		gsm = line[0]
		hit[gsm] = 0

	while(True):
		N = len(mixed)
		if N < 2:
			break
		print "current", N, time.ctime()
		
		global deep
		deep = 0
		temp= DFS(mixed, mixed[0], [0], N, beta, mini_conf)

		if temp == None:
			print "pop", time.ctime()
			del mixed[0]
		else:
			ubr.append(temp)
			lcs = temp[1]
			for line in matrix:
				gsm = line[0]
				if isSubsequence(lcs, line):
					hit[gsm] += 1
			mixed = [line for line in mixed if (hit[line[0][0]] == 0 or hit[line[0][1]] == 0)]
		
	print len(ubr)
	with open("ubr.D2.0.8.pickle", "w") as fp:
		pickle.dump(ubr, fp)

# lbr: [beta, alpha, sup, conf]
def LBR():
	ubr = read_data("ubr.D2.0.8")
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
				#count = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
				#count = {'AB': 0, 'CD': 0}
				count = {'ABC': 0, 'D': 0}
				for line in matrix:
					if isSubsequence(e, line):
						beta_1 = which_beta_B(line, GSM_info)
						count[beta_1] += 1
				lbr_sup = count[beta]
				lbr_conf = 1.0  * lbr_sup / sum(count.values())
				if lbr_conf >= conf:
					temp = [beta, e, lbr_sup, lbr_conf]
					lbr.append(temp)
					
					print temp
					flag = True
					break
			
	with open("lbr.D2.0.8.pickle", "w") as fp:
		pickle.dump(lbr, fp)	

def classifier_builder():
	rules = read_data("lbr.D2.0.8")
	'''
	rules_2 = read_data_b("lbr.D2b")
	for e in rules_2:
		rules.append(e)
	'''
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
		hit = [line for line in remain if isSubsequence(alpha, line) and which_beta_B(line, GSM_info) == beta]
		#hit = [line for line in remain if isSubsequence(alpha, line) and which_beta_B2(line, GSM_info) == beta]
		print len(hit), 0.05 * subsum[beta]
		#if len(hit) >= 0:
		if len(hit) >= 0.05 * subsum[beta]:
			for e in hit:
				remain.remove(e)
			print rule
			classifier.append(rule)
	#default
	count = {'ABC': 0, 'D': 0}
	#count = {'AB': 0, 'CD': 0}
	for line in remain:
		beta = which_beta_B(line, GSM_info)
		#beta = which_beta_B2(line, GSM_info)
		count[beta] += 1
	print len(remain), count

	with open("classifier.D2.0.8.pickle", "w") as fp:
		pickle.dump(classifier, fp)

def accuracy(classifier):
	tp, tn, fp, fn = 0.0, 0.0, 0.0, 0.0
	default = "ABC"
	positive = "ABC"
	for line in matrix:
		#real =  which_beta(line, GSM_info)
		real =  which_beta_D(line, GSM_info)
		#real =  which_beta_B2(line, GSM_info)
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
	
	if (tn + fn) != 0:
		precision2 = tn / (tn + fn)
	recall2 = tn / (tn + fp)
	if (precision2 + recall2) != 0:
		f12 = 2 * precision2 * recall2 / (precision2 + recall2)
	temp = (precision, recall, f1, acc, precision2, recall2, f12)
	
	return temp

def show_results():
	#classifier = read_data("classifier.D2.0.8")
	classifier = read_data_b("rule.nov.D.20th")

	index = []
	for i in range(1, len(classifier) + 1):
		temp = accuracy(classifier[0:i])
		index.append(temp)
		print i, classifier[i - 1]
	print temp
	plt.subplot(121)
	plt.plot(range(1, len(index)+1), [e[0] for e in index], "ko")
	plt.plot(range(1, len(index)+1), [e[1] for e in index], "ko")
	plt.plot(range(1, len(index)+1), [e[2] for e in index], "ko")
	plt.plot(range(1, len(index)+1), [e[3] for e in index], "ko")
	plt.plot(range(1, len(index)+1), [e[0] for e in index], "-", label="precision")
	plt.plot(range(1, len(index)+1), [e[1] for e in index], "-", label="recall")
	plt.plot(range(1, len(index)+1), [e[2] for e in index], "-", label="f1")
	plt.plot(range(1, len(index)+1), [e[3] for e in index], "-", label="acc")
	plt.grid(True)
	plt.legend()

	plt.subplot(122)
	plt.plot(range(1, len(index)+1), [e[4] for e in index], "-", label="precision2")
	plt.plot(range(1, len(index)+1), [e[5] for e in index], "-", label="recall2")
	plt.plot(range(1, len(index)+1), [e[6] for e in index], "-", label="f12")
	plt.plot(range(1, len(index)+1), [e[3] for e in index], "-", label="acc")
	plt.grid(True)
	plt.legend()

	plt.suptitle("Age.bi.train " + time.ctime())
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

def origin():
	# sd [[3, 42], [1, 54], [2, 60]]
	# the first is origin index, thesecond is na count
	sd = read_data("sd")

	classifier = read_data("rule.nov.D.20th")
	for rule in classifier:
		print rule
		beta, alpha, sup, conf = rule[0], rule[1], rule[2], rule[3]
		new = []
		for e in alpha:
			new.append(sd[e][0])
		print new

def get_w():
	pass
	L = len(matrix[0])
	wp = [0] * L
	wn = [0] * L
	for line in matrix:
		beta_1 = which_beta_D(line, GSM_info)
		if beta_1 == "D":
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

def nov():
	global matrix 
	rules = []
	while len(matrix) > 10:
		w = get_w()
		print w[0], w[5000], w[-1], len(matrix)
		top = [0, 0, 0, 0]
		for i in range(200):
			j = - (i + 1)
			e = [w[i][0], w[j][0]]

			#hit, nohit = {'AB': 0, 'CD': 0}, {'AB': 0, 'CD': 0}
			hit, nohit = {'ABC': 0, 'D': 0}, {'ABC': 0, 'D': 0}
			for line in matrix:
				beta = which_beta_D(line, GSM_info)
				if isSubsequence(e, line):
					hit[beta] += 1
				else:
					nohit[beta] += 1

			beta_h = max(hit, key=hit.get)
			sup_h = hit[beta_h]
			if sum(hit.values()) == 0:
				conf_h = 0
			else:
				conf_h = 1.0  * sup_h / sum(hit.values())

			beta_no = max(nohit, key=nohit.get)
			sup_no = nohit[beta_no]
			if sum(nohit.values()) == 0:
				conf_no = 0
			else:
				conf_no = 1.0  * sup_no / sum(nohit.values())

			if conf_h > conf_no:
				temp = [beta_h, e, sup_h, conf_h]
			else:
				e.reverse()
				temp = [beta_no, e, sup_no, conf_no]
			#print temp[-1], temp[-2]
			#if temp[-1] > 0.9 and temp[-2] > 10:
			if temp[-1] > top[-1] and temp[-2] > 10:
				top = temp
		if top[0] != 0:
			rules.append(top)
			matrix = [line for line in matrix if not isSubsequence(top[1], line)]
			print len(rules), top, time.ctime()
		else:
			break
				
	with open("rule.nov.D.20th.pickle", "w") as fp:
		pickle.dump(rules, fp)

if __name__ == "__main__":
	print "Start.", time.ctime()

	#matrix = read_data("age_train_802")
	matrix = read_data("age_test_198")
	GSM_info = read_data("GSM_info")
	matrix = convert(matrix)
	subsum = age_distribution()
	print subsum
	
	#head()
	#deep = 0
	#UBR_DFS()
	#LBR()
	#classifier_builder()
	#show_results()
	#rcbt()
	origin()
	#nov()

	print "End.", time.ctime()