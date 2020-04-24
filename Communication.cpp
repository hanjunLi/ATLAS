#include "Communication.h"

void Communication::
reset(int N, int myId, int numThreads, string partiesFile, Dispute* disp_pt) {
  _N = N;
  _myId = myId;
  _nCommThreads = numThreads-1;
  if (_N-1 <= numThreads) {
    _nCommThreads = _N-1-1;
    _commLoad = 1;
  } else {
    _commLoad = (_N-1 + numThreads - 1) / numThreads;
  }

  MPCCommunication comm;
  _partyChannels =
    comm.setCommunication(_io_service, _myId, _N, partiesFile);

  _disp_pt = disp_pt;
};

void Communication::
allToT(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs) {
  // no relay
  for (int i = 0; i < _N; i++) {
    assert(sendBufs[i].size() == sendBufs[0].size());
    assert(sendBufs[i].size() == recBufs[i].size());
  }
  recBufs[_myId] = sendBufs[_myId];
  vector<bool> readMask, writeMask;
  _disp_pt->tMaskP(_myId, writeMask); // write to my TSet
  _disp_pt->tMaskPRev(_myId, readMask); // read from all whose TSet include me
  vector<thread> threads(_nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t] =
      thread(&Communication::rwWorkerM, this,
             ref(sendBufs), ref(recBufs), ref(readMask), ref(writeMask), t);
  }
  rwWorkerM(sendBufs, recBufs, readMask, writeMask, _nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
}

void Communication::
allToAll(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, bool relay) {
  for (int i = 0; i < _N; i++) {
    assert(sendBufs[i].size() == sendBufs[0].size());
    assert(sendBufs[i].size() == recBufs[i].size());
  }
  
  // direct send to myself
  recBufs[_myId] = sendBufs[_myId];
  vector<thread> threads(_nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
      threads[t] =
        thread(&Communication::rwWorker, this,
               ref(sendBufs), ref(recBufs), t);
  }
  rwWorker(sendBufs, recBufs, _nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
}

void Communication::
allToOne(vector<byte> &myShare, vector<vector<byte>> &recBufs, int king, bool relay) {
  int sendSize = myShare.size();
  if (_myId == king) {
    // recBufs.clear();
    recBufs.resize(_N, vector<byte>(sendSize, 0));
  }

  if (_myId == king) {          // direct receive
    recBufs[_myId] = myShare;
    vector<thread> threads(_nCommThreads);
    for (int t = 0; t < _nCommThreads; t++) {
      threads[t] =
        thread(&Communication::rWorker, this, ref(recBufs), t);
    }
    rWorker(recBufs, _nCommThreads);
    for (int t = 0; t < _nCommThreads; t++) {
      threads[t].join();
    }
  } else { // direct send
    int kIdx = getChannelIdx(king);
    _partyChannels[kIdx]->getChannel()->write(myShare.data(), sendSize);
  }
}

void Communication::
oneToAll(vector<byte> &myShare, vector<vector<byte>> &sendBufs, int king, bool relay) {
  int recSize = myShare.size();
  if (_myId == king) {          // direct send
    myShare = sendBufs[_myId];
    vector<thread> threads(_nCommThreads);
    for (int t = 0; t < _nCommThreads; t++) {
      threads[t] =
        thread(&Communication::wWorker, this, ref(sendBufs), t);
    }
    wWorker(sendBufs, _nCommThreads);
    for (int t = 0; t < _nCommThreads; t++) {
      threads[t].join();
    }
  } else { // direct receive
    int kIdx = getChannelIdx(king);
    _partyChannels[kIdx]->getChannel()->read(myShare.data(), recSize);
  }
}

void Communication::
rwWorkerM(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs,
         vector<bool> &readMask, vector<bool> &writeMask, int id) {
  for (int i = id; i < _N-1; i+= (_nCommThreads+1)) {
    int pid = getChannelId(i);
    int sendSize = sendBufs[pid].size();
    int recSize = recBufs[pid].size();
    if (_myId < pid) {
      if (writeMask[pid]) {
        _partyChannels[i]->getChannel()->write(sendBufs[pid].data(), sendSize);
      }
      if (readMask[pid]) {
        _partyChannels[i]->getChannel()->read(recBufs[pid].data(), recSize);
      }
    } else {
      if (readMask[pid]) {
        _partyChannels[i]->getChannel()->read(recBufs[pid].data(), recSize);
      }
      if (writeMask[pid]) {
        _partyChannels[i]->getChannel()->write(sendBufs[pid].data(), sendSize);
      }
    }
  }
}

void Communication::
rWorker(vector<vector<byte>> &recBufs, int id ) {
  for (int i = id; i < _N-1; i+= (_nCommThreads+1)) {
    int pid = getChannelId(i);
    _partyChannels[i]->getChannel()->read(recBufs[pid].data(),
                                          recBufs[pid].size());
  }
}

void Communication::
wWorker(vector<vector<byte>> &sendBuf, int id) {
  for (int i = id; i < _N-1; i+= (_nCommThreads+1)) {
    int pid = getChannelId(i);
    _partyChannels[i]->getChannel()->write(sendBuf[pid].data(),
                                           sendBuf[pid].size());
  }
}

void Communication::
rwWorker(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int id) {
  for (int i = id; i < _N-1; i+= (_nCommThreads+1)) {
    int pid = getChannelId(i);
    int sendSize = sendBufs[pid].size();
    int recSize = recBufs[pid].size();
    if (_myId < pid) {
        _partyChannels[i]->getChannel()->write(sendBufs[pid].data(), sendSize);
        _partyChannels[i]->getChannel()->read(recBufs[pid].data(), recSize);
    } else {
        _partyChannels[i]->getChannel()->read(recBufs[pid].data(), recSize);
        _partyChannels[i]->getChannel()->write(sendBufs[pid].data(), sendSize);
    }
  }
}

