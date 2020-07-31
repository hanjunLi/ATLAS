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
  // improve: consistent way for resizing recBuf
  void allBroadcast(vector<byte>& sendBuf, vector<vector<byte>>& recBufs);
  void singleBroadcast(int king, vector<byte>& sendBuf, vector<byte>& recBuf);
  // -- w/o relay round functions
  void kingToT(int king, vector<byte> &myShare, vector<vector<byte>> &sendBufs);
  void TToKing(int king, vector<byte> &myShare, vector<vector<byte>> &recBufs);
  void allToT(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs);
  // -- w/ relay round functions
  void allToAll(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, bool relay = true);
  // NOTE: receive size calculated by the size of sendBuf
  void allToOne(vector<byte> &myShare, vector<vector<byte>> &recBufs, int king, bool relay = true);
  void oneToAll(vector<byte> &myShare, vector<vector<byte>> &sendBufs, int king, bool relay = true);

  // -- relayee helpers
  void recFromRelay(vector<vector<byte>>& recBufs, vector<bool>& relayerMask,
                    vector<vector<int>>& relayerLoad, int sendSize);
  void sendToRelay(vector<vector<byte>>& sendBufs, vector<bool>& relayerMask,
                   vector<vector<int>>& relayerLoad, int recSize);
  // -- relayer helpers
  void allToOneRelay(vector<int>& relayLoad, int king, int sendSize);
  void oneToAllRelay(vector<int>& relayLoad, int king, int recSize);
  // -- all-to-all relay helper (for both relayer and relayee)
  void allToAllRelay(vector<vector<byte>>& sendBufs, vector<vector<byte>>& recBufs);
  void allToRelayer(vector<vector<byte>>& rSendBufs, vector<vector<byte>>& relayBufs,
                    vector<bool>& rlyerMask, vector<bool>& rlyeeMask,
                    bool isRelayer, bool isRelayee);
  void relayerToAll(vector<vector<byte>>& rSendBufs, vector<vector<byte>>& rRecBufs,
                    vector<bool>& rlyerMask, vector<bool>& rlyeeMask,
                    bool isRelayer, bool isRelayee);
  // -- thread functions
  void rWorker(vector<vector<byte>> &recBufs, vector<bool>& mask, int id);
  void wWorker(vector<vector<byte>> &sendBufs, vector<bool>& mask, int id);
  void rwWorker(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs,
                vector<bool>& readMask, vector<bool>& writeMask, int id);

  // -- special version for dispute
  // improve: avoid writing these. parametrize above ones
  void allToOneStore(vector<byte> &myShare, vector<vector<byte>> &myRelay,
                     vector<vector<byte>> &recBufs, int king);
  void oneToOneStore(int fromId, int toId, vector<byte>& sendBuf, vector<byte>& relayBuf, vector<byte>& recBuf);
  
  
  void allToOneRelayStore(vector<int>& relayLoad, vector<vector<byte>>& relayBufs, int king, int sendSize);
  void write(vector<byte>& sendBuf, int id);
  void read(vector<byte>& recBuf, int id);
};



#endif /* COMMUNICATION_H__ */
