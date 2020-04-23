#! /bin/bash
n=3
inp=inputs.txt
typ=Zp
cir=circuit.txt #comp_circuit2.txt
pty=Parties.txt
itr=1
for i in `seq 0 1 2`;
do	
	./MPCHonestMajorityNoTriples -partyID $i -partiesNumber $n -inputFile $inp -outputFile output.txt -circuitFile $cir -fieldType $typ -partiesFile $pty -internalIterationsNumber $itr &
	echo "Running $i..."
done
