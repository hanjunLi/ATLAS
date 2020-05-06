dim = 20
sft = 36
for i in range(10):
	print("output #",i)
	print("value of w[0]:")
	fo = open("output"+str(i*2)+".txt","r");
	out = []
	for j in range(dim):
		b1,b2 = map(int,fo.readline().split())
		b1 = b1 * (2**64) + b2
		if(b1 > 2**90):
			b1 = -(2**127 - 1 - b1)
		out.append(b1 / (2**sft))
	print(out)
	fo.close()
	print("value of z:")
	fo = open("output"+str(i*2+1)+".txt","r");
	out = []
	for j in range(dim):
		b1,b2 = map(int,fo.readline().split())
		b1 = b1 * (2**64) + b2
		if(b1 > 2**90):
			b1 = -(2**127 - 1 - b1)
		out.append(b1 / (2**sft))
	print(out)
