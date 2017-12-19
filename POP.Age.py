#July 23 2017
import time
import random
import copy
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

def age_beta_D(line, GSM_info):
	gsm = line[0]
	age = GSM_info[gsm][1]
	#print age
	if age < 60:
		label = "ABC"
	else:
		label = "D"
	return label

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
	positive = "D"
	for line in matrix:
		real =  age_beta_D(line, GSM_info)
		#real =  sex_beta(line)
		#score = {'M': 0, 'F': 0}
		score = {'ABC': 0, 'D': 0}
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

def age_distribution(matrix):
	count = {'ABC': 0, 'D': 0}
	for line in matrix:
		beta = age_beta_D(line, GSM_info)
		count[beta] += 1
	return count

def get_w(matrix):
	pass
	L = len(matrix[0])
	wp = [0] * L
	wn = [0] * L
	for line in matrix:
		#beta_1 = sex_beta(line)
		beta_1 = age_beta_D(line, GSM_info)
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

		#rules = []
		pairs = []
		for i in range(50):
			#print i
			j = - (i + 1)
			e = [w[i][0], w[j][0]]

			#hit, nohit = {'AB': 0, 'CD': 0}, {'AB': 0, 'CD': 0}
			hit, nohit = {'ABC': 0, 'D': 0}, {'ABC': 0, 'D': 0}
			#hit, nohit = {'M': 0, 'F': 0}, {'M': 0, 'F': 0}
			for line in matrix:
				#beta = sex_beta(line)
				beta = age_beta_D(line, GSM_info)
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
				#rules.append(win)
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

	#remain = {'M': 0, 'F': 0}
	remain = {'ABC': 0, 'D': 0}
	for line in matrix:
		#beta = sex_beta(line)
		beta = age_beta_D(line, GSM_info)
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
	sample = sampling("age_test", 2000)
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
	subsum = age_distribution(matrix)
	
	expect = {}
	expect["ABC"] = n * 1.0 * subsum["ABC"] / len(matrix)
	expect["D"] = n * 1.0 * subsum["D"] / len(matrix)

	count = {}
	m = []
	for line in matrix:
		#beta = sex_beta(line)
		beta = age_beta_D(line, GSM_info)
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
		sample = sampling("age_train", i)
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

def test():
	css = read_data("bagging.css.timing")


	for cs in css:
		print len(cs)
		for c in cs:
			print c
	
if __name__ == "__main__":
	print "Start.", time.ctime()
	GSM_info = read_data("GSM_info")

	#timing()
	acc()
	
	#nov()
	#bagging()
	#test()

	print "End.", time.ctime()