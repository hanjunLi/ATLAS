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
singleBroadcast(int king, vector<byte> &sendBuf, vector<byte> &recBuf) {
  int recSize = recBuf.size();
  if (_myId != king) {          // non-king receives
    int kIdx = getChannelIdx(king);
    _partyChannels[kIdx]->getChannel()->read(recBuf.data(), recSize);
    return;
  } // king sends

  recBuf = sendBuf;
  vector<vector<byte>> sendBufs(_N, sendBuf);
  vector<thread> threads(_nCommThreads);
  vector<bool> activeVec;
  _disp_pt->nonCorrMask(activeVec);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t] =
      thread(&Communication::wWorker, this, ref(sendBufs), ref(activeVec), t);
  }
  wWorker(sendBufs, activeVec, _nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
}

void Communication::
allBroadcast(vector<byte> &sendBuf, vector<vector<byte> > &recBufs) {
  for (int i = 0; i < _N; i++) {
    assert(sendBuf.size() == recBufs[i].size());
  }
  recBufs[_myId] = sendBuf;
  vector<vector<byte>> sendBufs(_N, sendBuf);
  vector<thread> threads(_nCommThreads);
  vector<bool> activeVec;
  _disp_pt->nonCorrMask(activeVec);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t] =
      thread(&Communication::rwWorker, this,
             ref(sendBufs), ref(recBufs), ref(activeVec), ref(activeVec), t);
  }
  rwWorker(sendBufs, recBufs, activeVec, activeVec, _nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
}

void Communication::
kingToT(int king, vector<byte> &myShare, vector<vector<byte>> &sendBufs) {
  // NOTE: no relay needed
  vector<bool> TVec;
  _disp_pt->tMaskP(king, TVec);
  int recSize = myShare.size();
  if (!TVec[_myId]) {           // non-T party gets zeros
    myShare.clear();
    myShare.resize(recSize, 0);
    return;
  } else if (_myId != king) {   // T party receives
    int kIdx = getChannelIdx(king);
    _partyChannels[kIdx]->getChannel()->read(myShare.data(), recSize);
    return;
  } // else: king sends

  myShare = sendBufs[_myId];
  vector<thread> threads(_nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t] =
      thread(&Communication::wWorker, this, ref(sendBufs), ref(TVec), t);
  }
  wWorker(sendBufs, TVec, _nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
}

void Communication::
TToKing(int king, vector<byte> &myShare, vector<vector<byte>> &recBufs) {
  // NOTE: no relay needed
  vector<bool> TMask;
  _disp_pt->tMask(TMask);
  int sendSize = myShare.size();
  if (!TMask[_myId]) {          // non-T party sends zeros
    return;
  } else if (_myId != king) {   // T party sends
    int kIdx = getChannelIdx(king);
    _partyChannels[kIdx]->getChannel()->write(myShare.data(), sendSize);
    return;
  } // else: king receives

  recBufs.clear();
  recBufs.resize(_N, vector<byte>(sendSize, 0));
  recBufs[_myId] = myShare;
  vector<thread> threads(_nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t] =
      thread(&Communication::rWorker, this, ref(recBufs), ref(TMask), t);
  }
  rWorker(recBufs, TMask, _nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
}

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
      thread(&Communication::rwWorker, this,
             ref(sendBufs), ref(recBufs), ref(readMask), ref(writeMask), t);
  }
  rwWorker(sendBufs, recBufs, readMask, writeMask, _nCommThreads);
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
  if (relay && _disp_pt->hasDisp()) {        // do relay
    allToAllRelay(sendBufs, recBufs);
  }
  
  // direct send to myself
  recBufs[_myId] = sendBufs[_myId];
  vector<bool> myMask;
  _disp_pt->nonDispMask(_myId, myMask);
  vector<thread> threads(_nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
      threads[t] =
        thread(&Communication::rwWorker, this,
               ref(sendBufs), ref(recBufs), ref(myMask), ref(myMask), t);
  }
  rwWorker(sendBufs, recBufs, myMask, myMask, _nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
}

void Communication::
allToOne(vector<byte> &myShare, vector<vector<byte>> &recBufs, int king, bool relay) {
  int sendSize = myShare.size();
  if (_myId == king) {
    recBufs.clear();
    recBufs.resize(_N, vector<byte>(sendSize, 0));
  }
  
  if (relay && _disp_pt->hasDisp(king)) {    // do relay
    vector<int> relay;
    vector<bool> relayerMask;
    vector<vector<int>> relayerLoad;
    _disp_pt->rlyeeVecs(king, relay, relayerMask, relayerLoad);
    if (_myId == king) {        // receive from relayers
      recFromRelay(recBufs, relayerMask, relayerLoad, sendSize);
    } else if (relay[_myId] > -1) { // send to relayers
      int relayerIdx = getChannelIdx(relay[_myId]);
      _partyChannels[relayerIdx]->getChannel()->write(myShare.data(), sendSize);
    } else if (relayerMask[_myId]) { // receive from relayees, send to king
      allToOneRelay(relayerLoad[_myId], king, sendSize);
    } // else: send directly to king later
  }

  if (_myId == king) {          // direct receive
    recBufs[_myId] = myShare;
    vector<bool> myMask;
    _disp_pt->nonDispMask(_myId, myMask);
    vector<thread> threads(_nCommThreads);
    for (int t = 0; t < _nCommThreads; t++) {
      threads[t] =
        thread(&Communication::rWorker, this, ref(recBufs), ref(myMask), t);
    }
    rWorker(recBufs, myMask, _nCommThreads);
    for (int t = 0; t < _nCommThreads; t++) {
      threads[t].join();
    }
  } else if (!_disp_pt->isDisp(_myId, king)) { // direct send
    int kIdx = getChannelIdx(king);
    _partyChannels[kIdx]->getChannel()->write(myShare.data(), sendSize);
  }
}

void Communication::
oneToAll(vector<byte> &myShare, vector<vector<byte>> &sendBufs, int king, bool relay) {
  int recSize = myShare.size();
  if (relay && _disp_pt->hasDisp(king)) {    // do relay
    vector<int> relay;
    vector<bool> relayerMask;
    vector<vector<int>> relayerLoad;
    _disp_pt->rlyeeVecs(king, relay, relayerMask, relayerLoad);
    if (_myId == king) {        // send to relayers
      sendToRelay(sendBufs, relayerMask, relayerLoad, recSize);
    } else if (relay[_myId] > -1) { // receive from relayers
      int relayerIdx = getChannelIdx(relay[_myId]);
      _partyChannels[relayerIdx]->getChannel()->read(myShare.data(), recSize);
    } else if (relayerMask[_myId]) { // receive from king, send to relayees
      oneToAllRelay(relayerLoad[_myId], king, recSize);
    } // else: receive directly from king later
  }

  if (_myId == king) {          // direct send
    myShare = sendBufs[_myId];
    vector<bool> myMask;
    _disp_pt->nonDispMask(_myId, myMask);
    vector<thread> threads(_nCommThreads);
    for (int t = 0; t < _nCommThreads; t++) {
      threads[t] =
        thread(&Communication::wWorker, this, ref(sendBufs), ref(myMask), t);
    }
    wWorker(sendBufs, myMask, _nCommThreads);
    for (int t = 0; t < _nCommThreads; t++) {
      threads[t].join();
    }
  } else if (!_disp_pt->isDisp(_myId, king)) { // direct receive
    int kIdx = getChannelIdx(king);
    _partyChannels[kIdx]->getChannel()->read(myShare.data(), recSize);
  }
}

void Communication::
recFromRelay(vector<vector<byte>>& recBufs, vector<bool>& relayerMask,
             vector<vector<int>>& relayerLoad, int sendSize) {
  vector<vector<byte>> tmpBufs(_N);
  for (int i=0; i<_N; i++) {    // combined size for each relayer
    tmpBufs[i].resize(relayerLoad[i].size() * sendSize);
  }
  vector<thread> threads(_nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t] =
      thread(&Communication::rWorker, this, ref(tmpBufs), ref(relayerMask), t);
  }
  rWorker(tmpBufs, relayerMask, _nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
  for (int i=0; i<_N; i++) {    // scatter relayed messages
    auto segStart = tmpBufs[i].begin();
    for (int recIdx : relayerLoad[i]) {
      copy(segStart, segStart + sendSize, recBufs[recIdx].begin());
      segStart += sendSize;
    }
  }
}

void Communication::
sendToRelay(vector<vector<byte>>& sendBufs, vector<bool>& relayerMask,
            vector<vector<int>>& relayerLoad, int recSize) {
  vector<vector<byte>> tmpBufs(_N);
  for (int i=0; i<_N; i++) { // combined size for each relayer
    tmpBufs[i].resize(relayerLoad[i].size() * recSize);
  }
  for (int i=0; i<_N; i++) { // combine relay messages
    auto segStart = tmpBufs[i].begin();
    for (int sendIdx : relayerLoad[i]) {
      copy(sendBufs[sendIdx].begin(), sendBufs[sendIdx].end(), segStart);
      segStart += recSize;
    }
  }
  vector<thread> threads(_nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t] =
      thread(&Communication::wWorker, this, ref(tmpBufs), ref(relayerMask), t);
  }
  wWorker(tmpBufs, relayerMask, _nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
}

void Communication::
allToOneRelay(vector<int>& relayLoad, int king, int sendSize) {
  // -- gather relayees
  vector<bool> relayee(_N, false);
  for (int rlyeeIdx : relayLoad) {
    relayee[rlyeeIdx] = true;
  }
  // -- receive from relayees
  vector<vector<byte>> recBufs(_N, vector<byte>(sendSize));
  vector<thread> threads(_nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t] =
      thread(&Communication::rWorker, this, ref(recBufs), ref(relayee), t);
  }
  rWorker(recBufs, relayee, _nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
  // -- compbine relayee messages
  vector<byte> sendBuf(relayLoad.size() * sendSize);
  auto segStart = sendBuf.begin();
  for (int rlyeeIdx : relayLoad) {
    copy(recBufs[rlyeeIdx].begin(), recBufs[rlyeeIdx].end(), segStart);
    segStart += sendSize;
  }
  // -- send to king
  int kIdx = getChannelIdx(king);
  _partyChannels[kIdx]->getChannel()->write(sendBuf.data(),
                                            sendSize * relayLoad.size());
}


void Communication::
oneToAllRelay(vector<int>& relayLoad, int king, int recSize) {
  // -- receive from king
  vector<byte> recBuf(relayLoad.size() * recSize);
  int kIdx = getChannelIdx(king);
  _partyChannels[kIdx]->getChannel()->read(recBuf.data(), recBuf.size());
  // -- scatter relayee messages
  vector<vector<byte>> sendBufs(_N, vector<byte>(recSize));
  auto segStart = recBuf.begin();
  for (int rlyeeIdx : relayLoad) {
    copy(segStart, segStart + recSize, sendBufs[rlyeeIdx].begin());
    segStart += recSize;
  }
  // -- gather relayees
  vector<bool> relayee(_N, false);
  for (int rlyeeIdx : relayLoad) {
    relayee[rlyeeIdx] = true;
  }
  // -- send to relayees
  vector<thread> threads(_nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t] =
      thread(&Communication::wWorker, this, ref(sendBufs), ref(relayee), t);
  }
  wWorker(sendBufs, relayee, _nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }  
}

void Communication::
allToAllRelay(vector<vector<byte> > &sendBufs, vector<vector<byte> > &recBufs) {  
  int msgSize = sendBufs[0].size();
  vector<int> relay;
  vector<bool> rlyerMask, rlyeeMask;
  vector<vector<int>> rlyerLoad, rlyeeLoad;
  bool isRlyee = _disp_pt->rlyeeVecs(_myId, relay, rlyerMask, rlyerLoad);
  bool isRlyer = _disp_pt->rlyerVecs(_myId, rlyeeMask, rlyeeLoad);

  // Note: this only works with uniform msgSize.
  vector<vector<byte>> rlyeeSendBufs(_N), rlyeeRecBufs(_N), 
    rlyerRecBufs(_N), rlyerSendBufs(_N);
  if (isRlyee) {            // gather relay messages to send
    for (int i=0; i<_N; i++) {
      rlyeeSendBufs[i].resize(rlyerLoad[i].size() * msgSize);
      auto segStart = rlyeeSendBufs[i].begin();
      for (int sIdx : rlyerLoad[i]) {
        copy(sendBufs[sIdx].begin(), sendBufs[sIdx].end(), segStart);
        segStart += msgSize;
      } // allocate buffer to receive relay messages
      rlyeeRecBufs[i].resize(rlyerLoad[i].size() * msgSize);
    }
  }
  if (isRlyer) {            // allocate relay buffer to receive
    for (int i=0; i<_N; i++) {
      rlyerRecBufs[i].resize(rlyeeLoad[i].size() * msgSize);
    }
  }
  allToRelayer(rlyeeSendBufs, rlyerRecBufs,
               rlyerMask, rlyeeMask, isRlyer, isRlyee);
  
  if (isRlyer) {           // scatter and recombine relay messages
    for (int i=0; i<_N; i++) {
      auto segStart = rlyerRecBufs[i].begin();
      for (int rIdx : rlyeeLoad[i]) {
        rlyerSendBufs[rIdx].insert(rlyerSendBufs[rIdx].end(),
                              segStart, segStart + msgSize);
        segStart += msgSize;
      }
    }
  }
  relayerToAll(rlyerSendBufs, rlyeeRecBufs,
               rlyerMask, rlyeeMask, isRlyer, isRlyee);
  if (isRlyee) {            // scatter relayed messages
    for (int i=0; i<_N; i++) {
      auto segStart = rlyeeRecBufs[i].begin();
      for (int rIdx : rlyerLoad[i]) {
        copy(segStart, segStart + msgSize, recBufs[rIdx].begin());
        segStart += msgSize;
      }
    }
  }
}

void Communication::
allToRelayer(vector<vector<byte>>& rSendBufs, vector<vector<byte>>& relayBufs,
             vector<bool>& rlyerMask, vector<bool>& rlyeeMask,
             bool isRelayer, bool isRelayee) {
  // -- send / receive round
  vector<thread> threads(_nCommThreads);
  if (isRelayee) {
    if (isRelayer) {           // write to relayer and read from relayee
      for (int t = 0; t < _nCommThreads; t++) {
        threads[t] =
          thread(&Communication::rwWorker, this, ref(rSendBufs), ref(relayBufs),
                 ref(rlyeeMask), ref(rlyerMask), t);
      }
      rwWorker(rSendBufs, relayBufs, rlyeeMask, rlyerMask,
                          _nCommThreads);
    } else {                  // write to relayer
      for (int t = 0; t < _nCommThreads; t++) {
        threads[t] = thread(&Communication::wWorker, this,
                            ref(rSendBufs), ref(rlyerMask), t);
      }
      wWorker(rSendBufs, rlyerMask, _nCommThreads);
    }
  } else if (isRelayer) {     // read from relayee
    for (int t = 0; t < _nCommThreads; t++) {
      threads[t] = thread(&Communication::rWorker, this,
                          ref(relayBufs), ref(rlyeeMask), t);
    }
    rWorker(relayBufs, rlyeeMask, _nCommThreads);
  } else {                    // nothing to do
    return;
  }
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
}

void Communication::
relayerToAll(vector<vector<byte>>& rSendBufs, vector<vector<byte>>& rRecBufs,
             vector<bool>& rlyerMask, vector<bool>& rlyeeMask,
             bool isRelayer, bool isRelayee) {

  vector<thread> threads(_nCommThreads);
  if (isRelayee) {
    if (isRelayer) {          // read from relayer, write to relayee
      for (int t = 0; t < _nCommThreads; t++) {
        threads[t] =
          thread(&Communication::rwWorker, this, ref(rSendBufs), ref(rRecBufs),
                 ref(rlyerMask), ref(rlyeeMask), t);
      }
      rwWorker(rSendBufs, rRecBufs, rlyerMask, rlyeeMask, _nCommThreads);
    } else {                  // read from relayer
      for (int t = 0; t < _nCommThreads; t++) {
        threads[t] = thread(&Communication::rWorker, this,
                            ref(rRecBufs), ref(rlyerMask), t);
      }
      rWorker(rRecBufs, rlyerMask, _nCommThreads);
    }
  } else if (isRelayer) {     // write to relayee
    for (int t = 0; t < _nCommThreads; t++) {
      threads[t] = thread(&Communication::wWorker, this,
                          ref(rSendBufs), ref(rlyeeMask), t);
    }
    wWorker(rSendBufs, rlyeeMask, _nCommThreads);
  } else {
    return;
  }
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
}

void Communication::
rWorker(vector<vector<byte>> &recBufs, vector<bool>& mask, int id ) {
  for (int i = id; i < _N-1; i+= (_nCommThreads+1)) {
    int pid = getChannelId(i);
    if (mask[pid]) {
      _partyChannels[i]->getChannel()->read(recBufs[pid].data(),
                                            recBufs[pid].size());
    } // else: reads nothing from non-T party
  }
}

void Communication::
wWorker(vector<vector<byte>> &sendBuf, vector<bool>& mask, int id) {
  for (int i = id; i < _N-1; i+= (_nCommThreads+1)) {
    int pid = getChannelId(i);
    if (mask[pid]) {
      _partyChannels[i]->getChannel()->write(sendBuf[pid].data(),
                                             sendBuf[pid].size());
    } // else: writes nothing to non-T party
  }
}

void Communication::
rwWorker(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs,
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
allToOneStore(vector<byte> &myShare, vector<vector<byte>> &myRelay,
              vector<vector<byte>> &recBufs, int king) {
  int sendSize = myShare.size();
  if (_myId == king) {
    recBufs.clear();
    recBufs.resize(_N, vector<byte>(sendSize, 0));
  }

  myRelay.assign(_N, vector<byte>());
  if (_disp_pt->hasDisp(king)) {
    vector<int> relay;
    vector<bool> relayerMask;
    vector<vector<int>> relayerLoad;
    _disp_pt->rlyeeVecs(king, relay, relayerMask, relayerLoad);
    if (_myId == king) {        // receive from relayers
      recFromRelay(recBufs, relayerMask, relayerLoad, sendSize);
    } else if (relay[_myId] > -1) { // send to relayers
      int relayerIdx = getChannelIdx(relay[_myId]);
      _partyChannels[relayerIdx]->getChannel()->write(myShare.data(), sendSize);
    } else if (relayerMask[_myId]) { // receive from relayees, send to king
      allToOneRelayStore(relayerLoad[_myId], myRelay, king, sendSize);
    } // else: send directly to king later
  }

  if (_myId == king) {          // direct receive
    recBufs[_myId] = myShare;
    vector<bool> myMask;
    _disp_pt->nonDispMask(_myId, myMask);
    vector<thread> threads(_nCommThreads);
    for (int t = 0; t < _nCommThreads; t++) {
      threads[t] =
        thread(&Communication::rWorker, this, ref(recBufs), ref(myMask), t);
    }
    rWorker(recBufs, myMask, _nCommThreads);
    for (int t = 0; t < _nCommThreads; t++) {
      threads[t].join();
    }
  } else if (!_disp_pt->isDisp(_myId, king)) { // direct send
    int kIdx = getChannelIdx(king);
    _partyChannels[kIdx]->getChannel()->write(myShare.data(), sendSize);
  }
}

void Communication::
allToOneRelayStore(vector<int>& relayLoad, vector<vector<byte>>& relayBufs,
                   int king, int sendSize) {
  // -- gather relayees
  vector<bool> relayee(_N, false);
  for (int rlyeeIdx : relayLoad) {
    relayee[rlyeeIdx] = true;
    relayBufs[rlyeeIdx].resize(sendSize);
  }
  // -- receive from relayees
  vector<thread> threads(_nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t] =
      thread(&Communication::rWorker, this, ref(relayBufs), ref(relayee), t);
  }
  rWorker(relayBufs, relayee, _nCommThreads);
  for (int t = 0; t < _nCommThreads; t++) {
    threads[t].join();
  }
  // -- compbine relayee messages
  vector<byte> sendBuf(relayLoad.size() * sendSize);
  auto segStart = sendBuf.begin();
  for (int rlyeeIdx : relayLoad) {
    copy(relayBufs[rlyeeIdx].begin(), relayBufs[rlyeeIdx].end(), segStart);
    segStart += sendSize;
  }
  // -- send to king
  int kIdx = getChannelIdx(king);
  _partyChannels[kIdx]->getChannel()->write(sendBuf.data(),
                                            sendSize * relayLoad.size());
}

void Communication::
oneToOneStore(int fromId, int toId, vector<byte> &sendBuf,
              vector<byte> &relayBuf, vector<byte> &recBuf) {
  relayBuf.clear();
  int recSize = recBuf.size();
  if (_disp_pt->isDisp(fromId, toId)) {
    // use relay
    int relayer = _disp_pt->relayer(toId, fromId);
    if (_myId == fromId) {
      int idx = getChannelIdx(relayer);
      _partyChannels[idx]->getChannel()->write(sendBuf.data(), recSize);
    }
    if (_myId == relayer) {
      relayBuf.resize(recSize);
      int fromIdx = getChannelIdx(fromId);
      _partyChannels[fromIdx]->getChannel()->read(relayBuf.data(), recSize);
      int toIdx = getChannelIdx(toId);
      _partyChannels[toIdx]->getChannel()->write(relayBuf.data(), recSize);
    }
    if (_myId == toId) {
      int idx = getChannelIdx(relayer);
      _partyChannels[idx]->getChannel()->read(recBuf.data(), recSize);
    }
  } else {
    if (_myId == fromId) {
      int idx = getChannelIdx(toId);
      _partyChannels[idx]->getChannel()->write(sendBuf.data(), recSize);
    }
    if (_myId == toId) {
      int idx = getChannelIdx(fromId);
      _partyChannels[idx]->getChannel()->read(recBuf.data(), recSize);
    }
  }
}

void Communication::
write(vector<byte> &sendBuf, int id) {
  int idx = getChannelIdx(id);
  _partyChannels[idx]->getChannel()->write(sendBuf.data(), sendBuf.size());
}

void Communication::
read(vector<byte> &recBuf, int id) {
  int idx = getChannelIdx(id);
  _partyChannels[idx]->getChannel()->read(recBuf.data(), recBuf.size());
}
