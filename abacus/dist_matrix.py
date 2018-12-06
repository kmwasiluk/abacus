import itertools
from random import randint
from numpy.polynomial import polynomial as P
import numpy as np
#For comparison/timing
def substrings(s):
	return [s[i:i+j] for i in range(0,len(s)) for j in range(0,len(s))]
def gen_strs(ba, m):
	return [''.join(x) for x in itertools.product(ba, repeat=m)]
#This is the function to find the suggested proxy for relative edit distances
def compare_strs(s1,s2, max_str_len=1, cpc_function="Square", kernel_smoothing=False, norm="Jaccard"):
	v1 = compacted_dvectors(s1,max_str_len,cpc_function,kernel_smoothing)
	v2 = compacted_dvectors(s2,max_str_len,cpc_function,kernel_smoothing)
	L = []
	for i in range(0,len(v1)):
		L.append(vnorm(v1[i][0],v2[i][0],norm))
	return L
#Performs same as above, but returns sum.
def ld_proxy(s1,s2,max_str_len=1,cpc_function="Square", kernel_smoothing=False, norm="Jaccard"):
	L = compare_strs(s1,s2,max_str_len,cpc_function,kernel_smoothing,norm)
	return sum(L)

# Distance matrix for constructing the phylogenetic tree, input arr is a list of strings from which 
# we are creating the tree.
def pd_matrix(arr,max_str_len=1,cpc_function="Square", kernel_smoothing=False, norm="Jaccard"):
	M = np.zeros(shape=(len(arr), len(arr)))
	carr = [compacted_dvectors(w,max_str_len,cpc_function,kernel_smoothing) for w in arr]
	for i in range(0,len(carr)):
		for j in range(i+1,len(carr)):
			d = 0
			#possibly a poor design decision at the beginning... 
			if(norm == "Jaccard"):
				M[i][j] = M[j][i] = jaccard(carr[i],carr[j])
				continue
			for k in range(0,len(carr[i])):
				d += vnorm(carr[i][k][0],carr[j][k][0],norm)
			M[i][j] = M[j][i] = d
	return M

def vnorm(v1,v2,norm="L1"):
    if norm == "JSD":
    	return jensen_shannon_divergence(v1,v2)
    elif norm == "L1":
    	return l1_norm(v1,v2)
    elif norm == "L2":
    	return l2_norm(v1,v2)
    elif norm == "KL":
    	return KL_divergence(v1,v2) # Does not work well at all
    else: #Can add different norms later.
        return 0

def jaccard(M1,M2):
	numerator = 0
	denominator = 0
	for k in range(0,max(len(M1),len(M2))):
		n = len(M1[k][0])
		m = len(M2[k][0])
		for p in range(0, max(len(M1[k][0]), len(M2[k][0]))):
			if(p >= n):
				denominator+= M2[k][0][p]
			elif(p >= m):
				denominator += M1[k][0][p]
			else:
				denominator += max(M1[k][0][p], M2[k][0][p])
				numerator += min(M1[k][0][p], M2[k][0][p])
	return 1 - float(numerator/denominator)
			

def l1_norm(v1,v2):
    ret = 0
    if(len(v1) < len(v2)):
        smaller = v1
        larger = v2
    else:
        smaller = v2
        larger = v1
    i = 0
    while(i < len(smaller)):
        ret += abs(smaller[i] - larger[i])
        i+=1
    while(i < len(larger)):
        ret+= abs(larger[i])
        i+=1
    return ret

def l2_norm(v1,v2):
	ret = 0
	if(len(v1) < len(v2)):
		smaller = v1
		larger = v2
	else:
		smaller = v2
		larger = v1
	i = 0
	while(i < len(smaller)):
		ret += (smaller[i] - larger[i])**2
		i+=1
	while(i < len(larger)):
		ret+= larger[i]**2
		i+=1
	return float(np.sqrt(ret))

#Not as efficient as it could be
def jensen_shannon_divergence(v1, v2):
	if(len(v1) < len(v2)):
		while(len(v1) < len(v2)):
			v1.append(0)
	if(len(v2) < len(v1)):
		while(len(v2) < len(v1)):
			v2.append(0)
	P = np.array(v1, dtype='float32')
	Q = np.array(v2, dtype = 'float32')
	M = (P+Q)/2
	sp = sum(P)
	sq = sum(Q)
	sm = sum(M)
	#Normalize distributions
	for i in range(0,np.size(P)):
		P[i] = P[i]/sp
		Q[i] = Q[i]/sq
		M[i] = M[i]/sm
	return 0.5*KL_divergence(P,M) + 0.5*KL_divergence(Q,M)

def KL_divergence(A,B):
	ret = 0
	for i in range(0,min(np.size(B),np.size(A))): #IF coming from JSD these will be equal..
		if(B[i] == 0 or A[i] == 0): #Must be how to handle this case?
			continue
		ret+= A[i] * (np.log(A[i]) - np.log(B[i]))#(np.log(A[i]))/(np.log(B[i])) #default is log_e
	return ret


# This function finds distribution of distances between ocurrences of each kmer, does not combine them.
def dist_vector(s, max_str_len=1):
	ba = 'ACGT'
	arr = []
	for i in range(1, max_str_len+1):
		arr = arr + gen_strs(ba,i)
	ret = []
	for w in arr:
		ret.append([dstr_polynomial(s,[w]),w])
	return ret


#arr being the characters/strs that we are measuring dist between.
def dstr_polynomial(s,arr):
	g1 = cpolynomial_num(s,arr)
	arr_rev = [w[::-1] for w in arr]
	g2 = cpolynomial_num(s,arr_rev,True)
	return P.polymul(g1,g2)


#Will take each dist_vector, smooth and compact it so that there are fewer comparisons. default will be O(sqrt(N)).
def compacted_dvectors(s,max_str_len=1, shape="Square", kernel_smoothing=False):
	dv = dist_vector(s,max_str_len)
	ret = dv[:]
	for i in range(0,len(dv)):
		ret[i][0] = compaction(ret[i][0], shape, kernel_smoothing)
	return ret

def compaction(p, shape="Square", kernel_smoothing=False):
	ret = p[int((np.size(p))/2):] #Coefficients will be symmetric, so just take the second half.
	ret = list(ret)
	if(kernel_smoothing):
		ret = ksmoothing(ret,shape)
	#Add different compaction functions later.
	if(shape == "Square"):
		ret = compact_square(ret)
	if(shape == "Square25"):
		ret = compact_square_delay_k(ret,25)
	if(shape == "Square50"):
		ret = compact_square_delay_k(ret,50)
	if(shape == "Square147"):
		ret = compact_square_delay_k(ret, 147)
	return ret

#Not a great method since we are converting to python lists, which have significant overhead.
def compact_square(p):
	i = 0
	j = 1
	counter = 1
	ret = []
	b = True
	while(b):
		ret.append(sum(p[i:j]))
		temp = j
		counter+=1
		j += counter
		i = temp
		if(j > np.size(p)):
			ret.append(sum(p[i:(np.size(p))]))
			b = False # Not really necessary, but w/e
			break
	return ret

def compact_square_delay_k(p,k):
	ret = p[0:k]
	i = k
	j = k+1
	counter = 1
	b = True
	while(b):
		ret.append(sum(p[i:j]))
		temp = j
		counter += 1
		j += counter
		i = temp
		if( j > np.size(p)):
			ret.append(sum(p[i:(np.size(p))]))
			b = False
			break
	return ret

#Implement this later if needed. It should smooth results along the edges of our compactions.
def ksmoothing(p,shape):
	return p

def cpolynomial_num(s,arr,rev=False):
	ret = np.zeros(len(s))
	if(not rev):
		for i in range(0,len(s)):
			for w in arr:
				k = len(w)
				if(s[i:(i+k)] == w):
					ret[i] += 1
	else:
		v = s[::-1]
		for i in range(0,len(v)):
			for w in arr:
				k = len(w)
				if(v[i:(i+k)] == w):
					ret[i] += 1
	return bring_pol_to_zero(ret)

def bring_list_to_zero(L):
	for j in range(0,len(L)):
		if(L[j] >= 1): #Leading coefficient (and every coef.) should be 1
			break
    	return L[j:]

def bring_pol_to_zero(g):
	for j in range(0,np.size(g)):
		if(g[j] >= 1):
			break
        #ret = g*x^(-j)
    	return g[j:]

def rand_actg_str(numchars):
	s = ""
	A = 'ACGT'
	while(len(s) < numchars):
		r = randint(0,3)
		s+= A[r]
	return s