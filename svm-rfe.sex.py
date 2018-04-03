import time
import cPickle as pickle
import matplotlib.pyplot as plt
from sklearn.feature_selection import RFE
from sklearn.svm import SVC

from sklearn.datasets import load_digits

def read_data(path):
	with open(path + ".pickle", "r") as fp:
		obj = pickle.load(fp)
	print len(obj), path + " elements load over.", time.ctime()
	return obj

def sex_beta(line):
	e = line[0][-1]
	return e

def sex_distribution(matrix):
	count = {'M': 0, 'F': 0}
	for line in matrix:
		beta = sex_beta(line)
		count[beta] += 1
	return count

def svm_rfe(X_train):
	#X_train = read_data("age_test")
	#subsum = age_distribution(X_train)
	#print subsum
	Y_train = [sex_beta(line) for line in X_train]
	X_train = [line[1:] for line in X_train]
	
	
	svc = SVC(kernel="linear", C=1)
	rfe = RFE(estimator=svc, n_features_to_select=None, step=0.5)
	rfe.fit(X_train, Y_train)
	#print rfe.ranking_
	#print rfe.score(X_train, Y_train)
	return rfe

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
	print count
	return m

def timing():
	css = []
	t = []
	for  i in [1000, 2000, 4000, 6000, 8000, 10000]:
		sample = sampling("sex_train", i)
	
		start  = time.clock()
		cs = svm_rfe(sample)
		elapsed = time.clock() - start
		
		t.append(elapsed)
		css.append(cs)
	
	with open("svm_rfe.models.pickle", "w") as fp:
		pickle.dump(css, fp)	
	
	with open("svm_rfe.timing.pickle", "w") as fp:
		pickle.dump(t, fp)	

	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t, "-")
	plt.ylabel("CPU time(second)")
	plt.xlabel("# of samples")
	plt.grid(True)
	plt.show()

def acc():
	X_test = read_data("sex_test")
	subsum = sex_distribution(X_test)
	print subsum
	Y_test = [sex_beta(line) for line in X_test]
	X_test = [line[1:] for line in X_test]

	rfes = read_data("svm_rfe.models")
	t, j = [], 0
	
	for j in range(len(rfes)):
		ac = rfes[j].score(X_test, Y_test)
		t.append(ac)

	
	with open("acc.svm_rfe.pickle", "w") as fp:
		pickle.dump(t, fp)

	plt.plot([1000, 2000, 4000, 6000, 8000, 10000], t, "-")
	plt.ylabel("Accuracy")
	plt.xlabel("# of samples")
	plt.grid(True)
	plt.show()


def demo():
	digits = load_digits()
	X = digits.images.reshape((len(digits.images), -1))
	y = digits.target

	# Create the RFE object and rank each pixel
	svc = SVC(kernel="linear", C=1)
	rfe = RFE(estimator=svc, n_features_to_select=None, step=1)
	rfe.fit(X, y)
	ranking = rfe.ranking_.reshape(digits.images[0].shape)

	print rfe.score(X, y)
	# Plot pixel ranking
	plt.matshow(ranking, cmap=plt.cm.Blues)
	plt.colorbar()
	plt.title("Ranking of pixels with RFE")
	plt.show()


if __name__ == "__main__":
	print "Start.", time.ctime()
	
	GSM_info = read_data("GSM_info")
	#demo()
	#svm_rfe()
	#timing()
	acc()

	print "End.", time.ctime()