#include<cstdio>
#include<iostream>
using namespace std;
char tmp[100010];
int main()
{
	freopen("YearPredictionMSD.txt","r",stdin);
	freopen("small_data.txt","w",stdout);
	int train_cnt = 100;
	int test_cnt = 10;
	int tot_train = 461525;
	int cnt=0;
	for(; cnt<train_cnt; cnt++)
	{
		scanf("%s",tmp);
		printf("%s\n",tmp);
	}
	for(; cnt<tot_train; cnt++)
		scanf("%s",tmp);
	for(; cnt<tot_train + test_cnt; cnt++)
	{
		scanf("%s",tmp);
		printf("%s\n",tmp);
	}
	fclose(stdin);
	fclose(stdout);
	return 0;
}
