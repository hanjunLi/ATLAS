#! /bin/bash
n=3
#inp=(input0.txt input1.txt input2.txt input3.txt)
inp=(input0.txt input1.txt input2.txt)
#inp=(input.txt input.txt input.txt)
typ=ZpMersenne127
cir=circuit.txt
pty=Parties.txt
itr=1
for i in `seq 0 1 2`;
do	
	./MPCHonestMajorityNoTriples -partyID $i -partiesNumber $n -inputFile ${inp[$i]} -outputFile output.txt -circuitFile $cir -fieldType $typ -partiesFile $pty -internalIterationsNumber $itr &
	echo "Running $i..."
done
