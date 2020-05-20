# Hanjun Li <lihanjun1212@gmail.com>

import numpy as np
# import matplotlib.pyplot as plt
from sklearn import linear_model
from scipy import linalg
from data_util import *

# -- global trainning data and parameters
data = np.genfromtxt("ethylene_CO.txt", skip_header=1)
rng = np.random.default_rng()
rng.shuffle(data)

data1 = np.append(np.expand_dims(data[:, 1], axis=1), data[:, 3:], axis=1)
data2 = data[:, 2:]
nSamples = data.shape[0]
nTests = nSamples // 5          # 80 - 20 divide
nTrain = nSamples - nTests
nParties = 5
nFeatures = 16
nIter = 10                      # TODO: experiment and change this
param_rho = 100                 # TODO: experiment and change this
param_lamb = 0.1                # TODO: experiment and change this

# -- training
data_train = data1[:nTrain, :]  # Helen only uses data1
data_test = data1[nTrain:, :]
# data_train = data2[:nTrain, :]
# data_test = data2[nTrain:, :]

data_train, data_test = process_data(data_train, data_test)
Xs, Ys = divide_data(data_train, nParties, nFeatures)
As = [np.empty((nFeatures, nFeatures))]*nParties
bs = [np.empty(nFeatures)] * nParties
us = [np.zeros(nFeatures)] * nParties
ws = [np.zeros(nFeatures)] * nParties
z = np.zeros(nFeatures)

for i in range(nParties):
    As[i], bs[i] = compress_data(Xs[i], Ys[i], param_rho)

# -- iteration: step 4
for _ in range(nIter):
    # a: wi = Ai * (bi + rho * (z - ui))
    ws = [As[i] @ (bs[i] + param_rho * (z - us[i])) for i in range(nParties)]
    # b: z = soft_thresh(average(wi + ui))
    z = soft_thresh(param_lamb/nParties/param_rho, (sum(ws)+sum(us))/nParties)
    print(max(map(abs, z)))
    # c: ui = ui + wi - z
    us = [us[i] + ws[i] - z for i in range(nParties)]

truth = (data_test[:, 0])
prediction = (data_test[:, 1:1+nFeatures] @ z)
l2Error = compute_error(prediction, truth)
print("Plaintext prediction error: ", l2Error)

clf = linear_model.Lasso(alpha=param_lamb)
# clf = linear_model.LinearRegression
clf.fit(data_train[:, 1:1+nFeatures], data_train[:, 0])
sk_prediction = clf.predict(data_test[:, 1:1+nFeatures])
skl2Error = compute_error(sk_prediction, truth)
print("sklearn prediction error: ", skl2Error)
