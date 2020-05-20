import numpy as np
from scipy import linalg

def divide_data(data_in, nParties, nFeatures):
    # rng = np.random.default_rng()
    # rng.shuffle(arr)
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


def process_data(data_train, data_test):
    X_average = np.average(data_train[:, 1:])
    Y_average = np.average(data_train[:, 0])
    data_train[:, 1:] -= X_average
    data_train[:, 0] -= Y_average
    data_test[:, 1:] -= X_average
    data_test[:, 0] -= Y_average
    return data_train, data_test

def write_fixed_pt_output(name, As, bs, sft):
    max_num = 0
    for i in range(len(As)):
        dim = As[i].shape[0]
        Ai_fpt = As[i] * sft
        bi_fpt = bs[i] * sft
        max_num = max(bi_fpt.max(), max_num)
        with open(name+str(i)+".txt", "w") as fo:
            fo.write(str(dim)+'\n')
            for row in Ai_fpt:
                fo.write(' '.join(map(str, map(int, np.rint(row)))) + '\n')
            fo.write(' '.join(map(str, map(int, np.rint(bi_fpt)))) + '\n')
