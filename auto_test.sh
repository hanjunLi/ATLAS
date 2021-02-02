#!/bin/bash
n=5
ips=(172.31.5.61 172.31.14.20 172.31.1.153 172.31.8.119 172.31.2.241 172.31.5.80 172.31.10.211 172.31.4.242 172.31.0.201 172.31.4.232 172.31.0.173 172.31.0.203 172.31.7.69 172.31.9.197 172.31.9.220 172.31.4.118 172.31.5.101)
for j in 1 2 4 6 8 10 20 40 60 80
do
echo $j
python3 preproc.py $n $j
for ((i=1; i<$n; i++))
do
	echo $i
	#ssh -i GBcalifornia.pem ubuntu@${ips[0]} << flag0
	scp -i GBcalifornia.pem ~/ATLAS/input$i.txt ubuntu@${ips[$i]}:~/ATLAS/input$i.txt
	ssh -i GBcalifornia.pem ubuntu@${ips[$i]} << flag
	cd ATLAS
	git pull
	make
 	nohup ./MPCHonestMajorityNoTriples -partyID $i -partiesNumber $n -inputFile input$i.txt -outputFile output.txt -circuitFile circuit.txt -fieldType ZpMersenne127 -partiesFile Parties.txt -internalIterationsNumber 1 > log.txt 2>&1 &
	#nohup ./run_experiment.sh > log.txt 2>&1 &
	exit
flag
#flag0
done
#ssh -i GBcalifornia.pem ubuntu@${ips[0]} << flag2
#cd ATLAS
git pull
make
./run_experiment.sh #start running local
#flag2
python3 preproc.py $n $j
echo "1 running done"
done
echo "Program Done"
