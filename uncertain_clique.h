#pragma once
#ifndef UNCERTAIN_CLIQUE
#define UNCERTAIN_CLIQUE

#include <iostream>
#include <algorithm>
#include <queue>
#include <stack>
#include <map>

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>

using namespace std;
#define pairs pair<int, double>

#define doube_deviation 10000000
#define precision 0.00000000001f

#define Stack stack<int *>
#define Stack_N stack<pairs *>

#define GENERAL_CLIQUE 1
#define CORE_ALG 2
#define EFF_CLIQUE 3
#define MAXIMUM_CLIQUE 4
#define ORIGINAL_MAXIMUM_CLQ 5

class Uncertain_Clique
{
private:
	int V, E;
	int *ver; // offset
	pairs *edges, **adj;
	int algorithm;

	int *degrees;

	int *colors; // MAXIMUM_CLIQUE
	//pairs *edges_in_degrees;
	int *pos;
	int *sub_set;
	int MAX_SIZE;
	int coloring_nums;

	int max_nei_nums;
	long maximal_cliques;

	//memory
	long graph_memory;
	long max_memory;
	long temp_memory;
	long max_temp_memory;

	int k_core;
	double eta;

	int *omega;

	double *max_clr_p;

	double **prob; // CORE_ALG
	double **eff_pro; // EFF_CLIQUE
	pairs **vert_nei_pro; //prune2 

	FILE *outc;

	//old GENERAL_CLIQUE
	void general_enumerate(int *R,double q, pairs *I, int I_size, pairs *C, int C_size);
	void general_max_enumerate(int *R,double q, pairs *I, int I_size, pairs *C, int C_size);
	int generateI(pairs *I, pairs *Is, int I_size, double q_n, int u, pairs * In);
	int generateX(pairs *C, int C_size, double q_n, int u, pairs * Cn);
	void general_algorithm_of_cliques();

	//core CORE_ALG
	void get_cores();
	int get_eta_degrees(int v, double *d1, double *d2);

	//EFF_CLIQUE
	//prune1
	int core_prune1(int *C);
	int is_k_core_or_not(int v, double *d1, double *d2);
	// enhanced prune
	int core_prune2(int *C);

	void free_Dynamic_pointers();

	//MAXIMUM_CLIQUE
	void maximun_clique(int *R, double q, pairs *I, int I_size, pairs *C, int C_size);
	int color_vertices();
	int color_vertices(pairs *Can ,int c_size, int set_NO);
	int count_colors(pairs *Can, int c_size, int *color_class);
	int provide_degrees(pairs * Can, int c_size, double q);
	int require_degrees(int *R, double pr, int color_nums, int set_NO, int max_v);
	int probability_of_vR(double pr,pairs * Can, int size, int color_nums);

	//an exact maximum probabilistic algorithm (2014)
	void Max_Pclq();
	int Get_Can(int i, int v, pairs *Can, double p);
	void search(pairs *Can, int c_size, int *R, int r_size,double p, bool &found);
	int build_Can(int v, pairs *Can, int st, int lt, double p, pairs *res);

	//colorful top-k probability core
	void colorful_topk_core(const int &k, int *left_vertices, int &left_size);
	void vertices_sort(int *vertices, const int vt_size);
	int heuristic_clique(int *vertices, const int vt_size);

public:
	Uncertain_Clique();
	~Uncertain_Clique();

	void read_graph(const char *str);
	void init_parameters(int prunes);

	void compute(int k, double n, int alg,int prunes);
};

extern int myfunc(const void *a,const void *b);
extern bool myfunc1(double &a, double &b);
extern bool fun_v(pairs &a, pairs &b);

#endif // !UNCERTAIN_CLIQUE