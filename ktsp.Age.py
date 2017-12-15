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


def age_beta_D(line, GSM_info):
	gsm = line[0]
	age = GSM_info[gsm][1]
	#print age
	if age < 60:
		label = "ABC"
	else:
		label = "D"
	return label

def age_distribution(matrix):
	count = {'ABC': 0, 'D': 0}
	for line in matrix:
		beta = age_beta_D(line, GSM_info)
		count[beta] += 1
	return count

def ktsp(matrix):
	matrix = list(matrix)
	#matrix = [e for e in matrix if matrix.index(e)%3!=k]
	print len(matrix)

	k = 5
	c0, c1 = "ABC", "D"
	tops = [[0] * 3 for i in range(k)]
	
	L = len(matrix[0])
	for i in range(1, L, 400):
		#print i, time.ctime()
		for j in range(i+1, L):
			e = [i, j]
			hit, nohit = {'ABC': 0, 'D': 0}, {'ABC': 0, 'D': 0}
			#hit, nohit = {'M': 0, 'F': 0}, {'M': 0, 'F': 0}
			for line in matrix:
				beta = age_beta_D(line, GSM_info)
				if line[i] < line[j]:
					hit[beta] += 1
				else:
					nohit[beta] += 1

			pijc0 = 1.0 * hit[c0] / (hit[c0] + nohit[c0])
			pjic0 = 1.0 * nohit[c0] / (hit[c0] + nohit[c0])
			pijc1 =  1.0 * hit[c1] / (hit[c1] + nohit[c1])
			pjic1 = 1.0 * nohit[c1] / (hit[c1] + nohit[c1])
			score = abs(pijc0 - pijc1)

			if pijc0 > pijc1:
				temp = [c0, e, score]
			else:
				temp = [c1, e, score]

			if temp[-1] > tops[-1][-1]:
				tops[-1] = temp
				tops.sort(key=lambda x:x[-1], reverse=True)
				#print top
				#print pijc0, pjic0, pijc1, pjic1
	print tops, time.ctime()
	return tops
	#with open("tsp.pickle", "w") as fp:
		#pickle.dump(top, fp)	

def sampling(data, n):
	matrix = read_data(data)
	subsum = age_distribution(matrix)

	expect = {}
	expect["ABC"] = n * 1.0 * subsum["ABC"] / len(matrix)
	expect["D"] = n * 1.0 * subsum["D"] / len(matrix)

	count = {}
	m = []
	for line in matrix:
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
	print count
	return m

def timing():
	css = []
	t = []
	for  i in [1000, 2000, 4000, 6000, 8000, 10000]:
		i = int(i * 0.8)
		sample = sampling("age_train", i)
	
		start  = time.clock()
		cs = ktsp(sample)
		elapsed = time.clock() - start
		
		t.append(elapsed * 40)
		css.append(cs)
	
	with open("ktsp.pairs.pickle", "w") as fp:
		pickle.dump(css, fp)	
	
	with open("ktsp.timing.pickle", "w") as fp:
		pickle.dump(t, fp)	

	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t, "-")
	plt.ylabel("CPU time(second)")
	plt.xlabel("# of samples")
	plt.grid(True)
	plt.show()


def accuracy(matrix, tsp):
	tc, fc = 0.0, 0.0
	for line in matrix:
		real =  age_beta_D(line, GSM_info)
		if line[tsp[1][0]] < line[tsp[1][1]]:
			predict  = tsp[0]
		else:
			if tsp[0] == "ABC":
				predict = "D"
			else:
				predict = "ABC"
	
		if real == predict:
			tc += 1
		else:
			fc += 1
	acc = tc / (tc + fc)
	return acc

def accuracy2(matrix, ktsp):
	tc, fc = 0.0, 0.0
	for line in matrix:
		real =  age_beta_D(line, GSM_info)
		#score = {'M': 0, 'F': 0}
		score = {'ABC': 0, 'D': 0}
		for tsp in ktsp:
			if line[tsp[1][0]] < line[tsp[1][1]]:
				predict  = tsp[0]
			else:
				if tsp[0] == "ABC":
					predict = "D"
				else:
					predict = "ABC"
			score[predict] += 1
		predict = max(score, key=score.get)
		if real == predict:
			tc += 1
		else:
			fc += 1
	acc = tc / (tc + fc)
	return acc

def acc():
	cs = read_data("ktsp.pairs")
	j = 0
	t1, t2 = [], []
	sample = sampling("age_test", 2000)
	for i in [1000, 2000, 4000, 6000, 8000, 10000]:
		#i = int(i * 0.2)

		print cs[j]
	
		acc2 = accuracy2(sample, cs[j])
		t2.append(acc2)

		j += 1
	
	with open("acc.ktsp.pickle", "w") as fp:
		pickle.dump(t2, fp)

	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t2, "-")
	plt.ylabel("Accuracy")
	plt.xlabel("# of samples")
	plt.grid(True)
	plt.show()


def compare_t():
	t1 = read_data("nov.timing")
	t2 = read_data("ktsp.timing")
		
	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t1, "-", label="our")
	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t2, "-", label="k-tsp")
	plt.ylabel("CPU time(second)")
	plt.xlabel("# of samples")
	plt.grid(True)
	plt.legend()
	plt.show()

def compare_a():
	t1 = read_data("acc.bagging")
	t2 = read_data("acc.ktsp")
	
	plt.ylim((0.9,1))
	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t1, "-", label="our")
	#plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t2, "-", label="k-tsp")
	plt.ylabel("Accuracy")
	plt.xlabel("# of samples")
	plt.grid(True)
	plt.legend()
	plt.show()


def test():
	cs = read_data("ktsp.pairs")
	for c in cs:
		print c

if __name__ == "__main__":
	print "Start.", time.ctime()
	GSM_info = read_data("GSM_info")
	#tsp()
	timing()
	#acc()

	#test()
	#compare_a()

	print "End.", time.ctime()