#coding: utf-8
#July 23 2017
import time
import sys
import copy
import itertools
import cPickle as pickle
import matplotlib.pyplot as plt

from matplotlib.font_manager import FontManager
from pylab import mpl
import subprocess

def get_matplot_zh_font():
    fm = FontManager()
    mat_fonts = set(f.name for f in fm.ttflist)

    output = subprocess.check_output('fc-list :lang=zh -f "%{family}\n"', shell=True)
    zh_fonts = set(f.split(',', 1)[0] for f in output.split('\n'))
    available = list(mat_fonts & zh_fonts)

    print '*' * 10, '可用的字体', '*' * 10
    for f in available:
        print f
    return available

def set_matplot_zh_font():
    available = get_matplot_zh_font()
    if len(available) > 0:
        mpl.rcParams['font.sans-serif'] = [available[0]]    # 指定默认字体
        mpl.rcParams['axes.unicode_minus'] = False          # 解决保存图像是负号'-'显示为方块的问题

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

	k = 3
	c0, c1 = "ABC", "D"
	tops = [[0] * 3 for i in range(100)]
	
	L = len(matrix[0])
	for i in range(1, L, 47):
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
	final = []
	final.append(tops[0])
	s = set(tops[0][1])

	
	tops = [e for e in tops if e[1][0] not in s and e[1][1] not in s]
	tops.sort(key=lambda x:x[-1], reverse=True)
	final.append(tops[0])
	s = s | set(tops[0][1])

	tops = [e for e in tops if e[1][0] not in s and e[1][1] not in s]
	tops.sort(key=lambda x:x[-1], reverse=True)
	final.append(tops[0])


	print final, time.ctime()
	return final

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
		#i = int(i * 0.8)
		sample = sampling("age_train", i)
	
		start  = time.clock()
		cs = ktsp(sample)
		elapsed = time.clock() - start
		
		t.append(elapsed * 47)
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
	#sample = sampling("age_test", 2000)
	sample  = read_data("age_test")
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
	t1 = read_data("pop.timing")
	t2 = read_data("ktsp.timing")
	t3 = read_data("svm_rfe.timing")
	
	'''
	fig, ax1 = plt.subplots()

	ax1.plot([1000, 2000, 4000, 6000, 8000, 10000], t1, 'b-')
	ax1.set_xlabel('# of samples')
	# Make the y-axis label, ticks and tick labels match the line color.
	ax1.set_ylabel('CPU time(second)', color='b')
	ax1.tick_params('y', colors='b')

	ax2 = ax1.twinx()
	ax2.plot([1000, 2000, 4000, 6000, 8000, 10000], [e/3600 for e in t2], 'r-')
	ax2.set_ylabel('CPU time(h)', color='r')
	ax2.tick_params('y', colors='r')

	fig.tight_layout()
	plt.show()
	'''
	fig, ax = plt.subplots()
	ax.plot([1000, 2000, 4000, 6000, 8000, 10000], t1, "-", label="k-ppt")
	ax.plot([1000, 2000, 4000, 6000, 8000, 10000], [e/3 for e in t2], "-", label="k-tsp")
	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t3, "-", label="svm-rfe")
	ax.set_ylabel(u"CPU时间/秒")
	ax.set_yscale("log")
	ax.set_xlabel(u"样本数量")
	
	plt.grid(True)
	plt.legend()
	plt.show()

def compare_a():
	t1 = read_data("acc.pop")
	t2 = read_data("acc.ktsp")
	t3 = read_data("acc.svm_rfe")

	print t1[-1], t2[-1], t3[-1]
	
	plt.ylim((0,1))
	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t1, "-", label="k-pop")
	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t2, "-", label="k-tsp")
	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t3, "-", label="svm-rfe")
	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t1, "ko")
	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t2, "ko")
	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t3, "ko")
	
	plt.ylabel(u"分类准确度")
	plt.xlabel(u"样本数量")
	plt.grid(True)
	plt.legend()
	plt.show()


def test():
	cs = read_data("ktsp.pairs")
	for c in cs:
		print c

if __name__ == "__main__":
	print "Start.", time.ctime()
	set_matplot_zh_font()
	
	GSM_info = read_data("GSM_info")
	#tsp()
	#timing()
	#acc()
	#test()
	compare_a()

	print "End.", time.ctime()