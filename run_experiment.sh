#! /bin/bash
n=5
inp=input.txt
typ=ZpMersenne127
cir=circuit.txt
pty=ore_party.txt
itr=1
i=0
./MPCHonestMajorityNoTriples -partyID $i -partiesNumber $n -inputFile $inp -outputFile output.txt -circuitFile $cir -fieldType $typ -partiesFile $pty -internalIterationsNumber $itr &
