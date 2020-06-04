#include "uncertain_clique.h"

#include <time.h>

bool fc(int a, int b)
{
	return a < b;
}

int main(int argc, char** argv)
{
	//clock_t cs, ct;
	struct timeval start_tm, end_tm;

	//get_graph("as20000102.txt", "datasets/as20000102.txt");
	//return 0 ;
	if (argc != 6){
		printf("Please input parameters: (file, k-constraint, p-probabilistic, algorithm, prune)\n");
		exit(1);
	}
	int k = atoi(argv[2]);
	double p = atof(argv[3]);
	int alg = atoi(argv[4]);
	int prune = atoi(argv[5]);
	Uncertain_Clique uc;
	//cs = clock();
	uc.read_graph(argv[1]);
	gettimeofday(&start_tm, NULL);
	uc.compute(k, p, alg, prune);
	//ct = clock();
	gettimeofday(&end_tm, NULL);
	//printf("time %d\n", int(ct - cs));
	printf("All time %.3f ms\n", double(end_tm.tv_sec - start_tm.tv_sec) * 1000.0f + 
		double(end_tm.tv_usec - start_tm.tv_usec) / 1000.0f);
	return 0;
}
