#input: n_parties, scale of train of each party
import numpy as np
# import matplotlib.pyplot as plt
from sklearn import linear_model
from scipy import linalg
import sys
np.random.seed(41)
n_train = 463715
n_test = 51630
# -- global trainning data and parameters
# data = np.genfromtxt("small_data.txt",delimiter=',')
data = np.genfromtxt("../data/YearPredictionMSD.txt", delimiter=',')
data_train = data[:n_train]
# data_train = data[:463715]
data_test = data[n_train:n_train + n_test]
# data_test = data[463715:]
nParties = int(sys.argv[1])
n_train = 1000 * nParties * int(sys.argv[2])
np.random.shuffle(data_train)
data_train = data_train[:n_train]
nFeatures = 90  # TODO: should be 90
nIter = 10                      # TODO: experiment and change this
param_rho = 10                  # TODO: experiment and change this
param_lamb = 0.1                # TODO: experiment and change this
max_num = 0
# -- MPC preparation parameters
MPC_sizes = [1000, 2000, 4000, 6000, 8000,
             10000, 20000, 40000, 60000, 80000, 100000]


def divide_data(data_in, nParties, nFeatures):
    dataSize = data_in.shape[0]
    perPartySize = dataSize // nParties 
    cutOff = perPartySize * nParties
    X = data_in[:cutOff, 1:nFeatures+1]
    Y = data_in[:cutOff, 0]
    Xs = X.reshape(nParties, perPartySize, nFeatures)
    Ys = Y.reshape(nParties, perPartySize)
    return Xs, Ys


def compress_data(Xi, Yi, rho):
    # nCols = Xi.shape[1]
    # Ai = (Xi^T * Xi + rho * I)^(-1)
    _, s, Vh = linalg.svd(Xi, False)
    Ai = Vh.T @ np.diag(1 / (s * s + rho)) @ Vh
    # Ai = np.linalg.inv(Xi.T @ Xi + rho*np.identity(nCols))
    # bi = Xi^T * Yi
    bi = Xi.T @ Yi
    return Ai, bi


def soft_thresh(thresh, z):
    # s(k, x) = +(|x| - k) if  x > k
    # s(k, x) = 0          if -k < x < k
    # s(k, x) = -(|x| - k) if -x > k
    # zp = np.sign(z) * np.maximum(0, (np.abs(z) - thresh))
    zp = np.maximum(0, z - thresh) - np.maximum(0, -z-thresh)
    return zp


def compute_error(prediction, truth):
    # return np.linalg.norm(np.abs(prediction - truth) / truth)**2
    nPredictions = prediction.shape[0]
    return np.linalg.norm(prediction - truth)**2 / nPredictions

def compute_mae(prediction, truth):
    n = prediction.shape[0]
    return sum(np.abs(truth - prediction)) / n

def preprocess_data(data_train, data_test):
    X_average = np.average(data_train[:, 1:])
    Y_average = np.average(data_train[:, 0])
    data_train[:, 1:] -= X_average
    data_train[:, 0] -= Y_average
    data_test[:, 1:] -= X_average
    data_test[:, 0] -= Y_average
    return data_train, data_test


sft = 2**52
#sft = 1
# -- initialization: step 1, 2, 3
data_train, data_test = preprocess_data(data_train, data_test)
# now only test first 10 feature ?
tmp_train = data_train[:, :1+nFeatures]

Xs, Ys = divide_data(tmp_train, nParties, nFeatures)
#Xs, Ys = divide_data(data_train, nParties, 90)
As = [np.empty((nFeatures, nFeatures))]*nParties
bs = [np.empty(nFeatures)]*nParties
for i in range(nParties):
    As[i], bs[i] = compress_data(Xs[i], Ys[i], param_rho)
    fo = open("input"+str(i)+".txt","w")
    #fo = open("input"+str(i)+"_90_"+str(n_train)+".txt", "w")
    #print(As[i].shape[0],bs[i].shape[0])
    fo.write(str(As[i].shape[0])+'\n')
    for j in range(As[i].shape[0]):
        for k in range(As[i].shape[1]):
            fo.write(str(int(As[i][j][k]*sft))+' ')
        fo.write('\n')
    for j in range(bs[i].shape[0]):
        if(bs[i][j] > max_num):#if(bs[i][j] > (2**(62-sft))):
            max_num = bs[i][j] #print("Too large B!")
        #    exit()
        fo.write(str(int(bs[i][j]*sft))+' ')
    fo.write('\n')
    fo.close()
us = [np.zeros(nFeatures)] * nParties
ws = [np.zeros(nFeatures)] * nParties
z = np.zeros(nFeatures)

# -- iteration: step 4
for _ in range(nIter):
    # a: wi = Ai * (bi + rho * (z - ui))
    ws = [As[i] @ (bs[i] + param_rho * (z - us[i])) for i in range(nParties)]
    for b1 in (ws[i] for i in range(nParties)):
        for a1 in b1:
            if(a1>max_num):
                max_num = a1
    # b: z = soft_thresh(average(wi + ui))
    z = soft_thresh(param_lamb/nParties/param_rho, (sum(ws)+sum(us))/nParties)
    for z1 in z:
        if(z1>max_num):
            max_num = z1
    # c: ui = ui + wi - z
    us = [us[i] + ws[i] - z for i in range(nParties)]
    #print("value of w",ws)
    #print("value of z",z)
    #print("value of us",us)

print("Max number:",max_num)
for _ in range(z.shape[0]):
    print(z[_])
# -- test accuracy
prediction = (data_test[:, 1:1+nFeatures] @ z).round()
truth = (data_test[:, 0]).round()
l2Error = compute_error(prediction, truth)
mae = compute_mae(prediction,truth)
print("Helen prediction error: ", l2Error,mae)

fo2 = open("output.txt","r")
_z = []
for _ in range(nFeatures):
    tmp = fo2.readline()
    v1,v2 = map(int,tmp.split()) #should be a length-2 vector
    curv = v1 * (2**64) + v2
    if curv > (2**90): #negative number
        curv = (2**127) - 1 - curv
        curv = -curv
    curv = curv / sft
    _z.append(curv)
    print(curv)
#calculated accuracy
res_file = open("result_song.txt","a")
res_file.write(str(n_train / nParties)+":")
prediction = (data_test[:, 1:1+nFeatures] @ _z).round()
truth = (data_test[:, 0]).round()
l2Error = compute_error(prediction, truth)
mae = compute_mae(prediction, truth)
print("My Helen prediction error: ", l2Error, mae)
res_file.write(str(l2Error)+","+str(mae)+"/")
# -- baseline accuracy
clf = linear_model.Lasso(alpha=param_lamb)
# clf = linear_model.LinearRegression()
clf.fit(data_train[:, 1:1+nFeatures], data_train[:, 0])
sk_prediction = clf.predict(data_test[:, 1:1+nFeatures]).round()
skl2Error = compute_error(sk_prediction, truth)
mae = compute_mae(sk_prediction,truth)
print("sklearn prediction error: ", skl2Error,mae)
res_file.write(str(skl2Error)+","+str(mae)+"\n")
res_file.close()
# -- prepare input data for MPC
for curSize in MPC_sizes:
    
    pass