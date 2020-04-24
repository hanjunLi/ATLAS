// Hanjun Li <lihanjun1212@gmail.com>

#ifndef COMMUNICATION_H__
#define COMMUNICATION_H__

#include <libscapi/include/comm/MPCCommunication.hpp>
#include <vector>
#include <thread>
#include <fstream>
#include "Dispute.h"

class Communication {
private:
  // N: number of parties;
  int _N, _myId, _nCommThreads, _commLoad; 
  boost::asio::io_service _io_service;
  vector<shared_ptr<ProtocolPartyData>> _partyChannels;
  Dispute* _disp_pt;            // updated in outer protocol

  // local helper
  __inline__
  int getChannelIdx(int id) {
    assert(id != _myId);
    return (id < _myId) ? id : id-1;
  };

  __inline__
  int getChannelId(int idx) {
    return (idx < _myId) ? idx : idx+1;
  }
  
public:
  Communication(){};
  ~Communication(){};           // _disp_pt deleted in outer protocol

  void reset(int N, int myId, int numThreads,
             string partiesFile, Dispute* disp_pt);
  // w/o relay round functions
  void allToT(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs);
  // w/ relay round functions
  void allToAll(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, bool relay = true);
  void allToOne(vector<byte> &myShare, vector<vector<byte>> &recBufs, int king, bool relay = true);
  void oneToAll(vector<byte> &myShare, vector<vector<byte>> &sendBufs, int king, bool relay = true);

  // thread functions
  void rwWorkerM(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs,
                vector<bool>& readMask, vector<bool>& writeMask, int id);
  void rWorker(vector<vector<byte>> &recBufs, int id);
  void wWorker(vector<vector<byte>> &sendBufs, int id);
  void rwWorker(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs,
                int id);
  
};



#endif /* COMMUNICATION_H__ */
