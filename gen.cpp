#include<cstdio>
#include<cstdlib>
int main()
{
	int n=3;
	int n_input = 8;
	int n_layer = 100;
	int gate_type = 2;
	int g_per_layer = n * n * n_input * n_input;
	int g = n_layer * g_per_layer;
	freopen("input30.txt","w",stdout);
	for(int i=1; i<=n_input; i++)
		printf("%d\n",rand());
	fclose(stdout);
	freopen("comp_circuit2.txt","w",stdout);
	printf("%d %d\n",g,n);
	int tot=0;
	for(int i=1; i<=n; i++)
	{
		printf("%d %d\n",i,n_input);
		for(int j=1; j<=n_input; j++)
			printf("%d\n",tot++);
	}
	int tmp = (n_layer-1) * g_per_layer + n * n_input;
	for(int i=1; i<=n; i++)
	{
		printf("%d %d\n",i,g_per_layer / n);
		for(int j=1; j<=g_per_layer / n; j++)
			printf("%d\n",tmp++);
	}
	//deal with first layer
	tmp = n * n_input;
	tot = tmp;
	for(int i=0; i<tmp; i++)
		for(int j=0; j<tmp; j++)
			printf("2 1 %d %d %d %d\n",i,j,tot++,gate_type);
	//other layers
	for(int i=2; i<=n_layer; i++)
	{
		int bas = tot - g_per_layer;
		for(int k=0; k<g_per_layer; k++)
		{
			int j1 = rand()%g_per_layer, j2 = rand()%g_per_layer;
			printf("2 1 %d %d %d %d\n",j1+bas,j2+bas,tot++,gate_type);
		}
	}
	fclose(stdout);
	return 0;
}
