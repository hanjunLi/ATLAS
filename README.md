### ATLAS: Efficient and Scalable MPCin the Honest Majority Setting

The repository implements the ATLAS MPC protocol in Section 4.

##### Installation

The protocol written in c++ and uses c++11 standard. It uses [libscapi](https://github.com/cryptobiu/libscapi).  
For `libscapi` installation instructions, visit [here](https://github.com/cryptobiu/libscapi/blob/master/build_scripts/INSTALL.md).  
After you installed `libscapi`, run `cmake . && make`

##### Usage

The protocol is designed for odd number (>= 3) of parties.

To run the the protocol open a terminal and run: 
`run_protocol.sh <min_party_id> <max_party_id> <number_of_parties> <input_file> <circuit_file> <filed type> <parties_file> <number_of_iterations>` 

* field_type can be one of this values:
    * ZpMersenne31
    * ZpMersenne61
    * ZpKaratsuba
    * GF2m
    * Zp

* parties_file - a file that contains the ip addresses and the port of all the parties. An example file can be found [here](../master/Parties.txt).
    
