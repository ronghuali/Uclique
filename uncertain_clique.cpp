#include "uncertain_clique.h"

int myfunc(const void *a, const void *b)
{
	return int((*(double *)b - *(double *)a) * doube_deviation);
}

bool myfunc1(double &a, double &b)
{
	return a > b;
}

bool fun_v(pairs &a, pairs &b)
{
	if (a.second != b.second)
		return a.second > b.second;
	else
		return a.first > b.first;
}

Uncertain_Clique::Uncertain_Clique()
{
	V = 0; E = 0;
	algorithm = 0;

	k_core = 0;
	eta = 0;

	graph_memory = 0;
	max_memory = 0;
	temp_memory = 0;
	max_temp_memory = 0;

	max_nei_nums = 0;
	maximal_cliques = 0;
	MAX_SIZE = 0;
	coloring_nums = 0;

	edges = NULL;
	ver = NULL;
	degrees = NULL;
	adj = NULL;

	prob = NULL;
	eff_pro = NULL;
	vert_nei_pro = NULL;

	colors = NULL;
	pos = NULL;
	sub_set = NULL;

	max_clr_p = NULL;
	omega = NULL;

	outc = NULL;
}

Uncertain_Clique::~Uncertain_Clique()
{
	printf("destruction\n");
	//if (edges != NULL) delete[] edges;
	//if (adj != NULL) delete[] adj;
	//if (ver != NULL) delete[] ver;
	//if (degrees != NULL) delete[] degrees;
	// if (prob != NULL){
	// 	for (int i = 0; i <= V; ++i)
	// 		if (prob[i] != NULL)
	// 			delete[] prob[i];
	// 	delete[] prob;
	// }
	// if (eff_pro != NULL){
	// 	for (int i = 0; i <= V; ++i)
	// 		if (eff_pro[i] != NULL)
	// 			delete[] eff_pro[i];
	// 	delete[] eff_pro;
	// }

	// if (vert_nei_pro != NULL){
	// 	for (int i = 0; i <= V; ++i)
	// 		if (vert_nei_pro[i] != NULL)
	// 			delete[] vert_nei_pro[i];
	// 	delete[] vert_nei_pro;
	// }

	// if (colors != NULL) delete[] colors;
	// if (pos != NULL) delete[] pos;
	// if (sub_set != NULL) delete[] sub_set;
	// if(max_clr_p != NULL) delete[] max_clr_p;
	// if(omega != NULL) delete[] omega;
}

void Uncertain_Clique::read_graph(const char * str)
{
	FILE *in = fopen(str, "r");

	if (in == NULL)
	{
		printf("No such the input file\n");
		exit(1);
	}
	printf("InFile: %s\n", str);
	if (fscanf(in, "%d %d", &V, &E) != 2)
		exit(1);
	assert(V > 0);
	assert(E > 0);

	int x, y, counts = 0;
	double p;

	ver = new int[V + 2];
	memset(ver, 0, sizeof(int) * (V + 2));

	//图文件为双向边
	/*graph_memory += sizeof(pairs) * E;
	graph_memory += sizeof(int) * (V + 1);
	memset(ver, 0, sizeof(int) * (V + 2));
	if (edges == NULL)
		edges = new pair<int, double>[E + 1];
	memset(edges, 0, sizeof(pairs) * (E + 1));
	int i = 0;
	while (fscanf(in, "%d %d %lf", &x, &y, &p) == 3)
	{
		assert(x <= V && x >= 0);
		assert(y <= V && y >= 0);
		assert(p > 0 && p <= 1);
		if (x == y)
			continue;
		ver[x]++;
		edges[i].first = y;
		edges[i].second = p;
		++i;
		if (feof(in)) break;
	}
	fclose(in);
	printf("|V| = %d\t|E| = %d\n", V, E / 2);
	for (int i = 0; i <= V; ++i){
		adj[i] = edges + counts;
		counts += ver[i];
		max_nei_nums = max(max_nei_nums, ver[i])
	}
	adj[V+1] = adj[V] + ver[V];
	printf("max nbr \t%d\n", max_nei_nums);
	exit(1);*/
	//for (int i = 0; i < V + 2; ++i)
	//	cout << i <<" "<< ver[i] << endl;

	//图文件为单向边
	E = 0;
	while (fscanf(in, "%d %d %lf", &x, &y, &p) == 3)
	{
		assert(x <= V && x >= 0);
		assert(y <= V && y >= 0);
		assert(p > 0 && p <= 1);
		if (x >= y) continue;
		++ver[x], ++ver[y];
		++E;
		if (feof(in)) break;
	}
	fclose(in);
	graph_memory += (sizeof(pairs)) * E * 2;
	graph_memory += sizeof(int) * (V + 1);
	//return ;

	counts = 0;
	if (edges == NULL) edges = new pair<int, double>[E * 2 + 1];
	if (adj == NULL) adj = new pairs*[V+2];
	memset(edges, 0, sizeof(pairs) * (E * 2 + 1));
	for (int i = 0; i <= V; ++i)
	{
		int d = ver[i];
		adj[i] = edges + counts;
		counts += d;
		ver[i] = 0;
		max_nei_nums = max(max_nei_nums, d);
		
	}
	adj[V+1] = adj[V] + ver[V];
	in = fopen(str, "r");
	if (fscanf(in, "%d %d", &V, &counts) != 2) exit(1);
	while (fscanf(in, "%d %d %lf", &x, &y, &p) == 3)
	{
		if (x >= y) continue;
		int id = ver[x]++;
		adj[x][id].first = y;
		adj[x][id].second = p;

		id = ver[y]++;
		adj[y][id].first = x;
		adj[y][id].second = p;
		if (feof(in)) break;
	}
	fclose(in);
	for (int i = 1; i <= V; ++i){
		assert(ver[i]>=0 && ver[i] <= 2*E);
	}
	printf("V=%d, E=%d, max_deg=%d\n", V, E, max_nei_nums);
	// for (int i = 0 ; i <= V; ++i)
	// {
	// 	int d = ver[i];
	// 	for (int j = 0; j < d; ++j)
	// 		printf("%d %d %.2f\n", i, adj[i][j].first, adj[i][j].second);
	// 	if (d > 0) break;
	// }
	// exit(1);
}

void Uncertain_Clique::init_parameters(int prunes)
{
	int *neis = NULL;
	pairs *s, *t, *e , *ed ;
	switch (algorithm)
	{
	case CORE_ALG:
		if (degrees == NULL){
			degrees = new int[V + 2];
			max_memory += sizeof(int) * (V+2);
		}
		memset(degrees, 0, sizeof(int) *(V + 2));
		if (prob == NULL){
			prob = new double*[V + 1];
			max_memory += sizeof(double) * (V + 1);
			int n_size;
			for (int i = 0; i <= V; ++i)
			{
				n_size = ver[i] + 1;
				prob[i] = new double[n_size];
				max_memory += sizeof(double) * (n_size);
				memset(prob[i], 0, sizeof(double)*(n_size));
			}
		}
		if (pos == NULL)
		{
			pos = new int[V + 1];
			memset(pos, 0, sizeof(int) *(V + 1));
			max_memory += sizeof(int) * (V+1);
		}
		break;
	case EFF_CLIQUE:
		break;
	case MAXIMUM_CLIQUE:

		if (sub_set == NULL)
		{
			sub_set = new int[V + 1];
			memset(sub_set, -1, sizeof(int) *(V + 1));
			max_memory += sizeof(int) * (V+1);
		}
		if (colors == NULL)
		{
			colors = new int[V + 1];
			memset(colors, -1, sizeof(int) *(V + 1));
			max_memory += sizeof(int) * (V+1);
		}
		if (pos == NULL)
		{
			pos = new int[V + 1];
			memset(pos, 0, sizeof(int) *(V + 1));
			max_memory += sizeof(int) * (V+1);
		}
		if (degrees == NULL){
			degrees = new int[V + 2];
			memset(degrees, 0, sizeof(int) * (V + 2));
		}
		break;
	default:
		break;
	}
	
	if (algorithm == MAXIMUM_CLIQUE || algorithm == EFF_CLIQUE)
	{
		if (eff_pro == NULL && prunes == 1){
			eff_pro = new double*[V + 1];
			max_memory += sizeof(double) * (V + 1);
			for (int i = 0; i <= V; ++i)
			{
				eff_pro[i] = new double[k_core + 1];
				max_memory += sizeof(double) * (k_core + 1);
				memset(eff_pro[i], 0, sizeof(double)*(k_core + 1));
			}
		}
		
	}
}

void Uncertain_Clique::general_enumerate(int *R, double q, pairs *I, int I_size, pairs *C, int C_size)
{
	if (I_size == 0 && C_size == 0 && *R > k_core)
	{
		//for (int i = 1; i <= *R; ++i)
		//	fprintf(outc, "%d\t%ld\n",R[i], maximal_cliques );
		++maximal_cliques;
		if (*R > MAX_SIZE)
			MAX_SIZE = *R;
		return;
	}
	if (I_size == 0)
		return;
	pairs *I_n, *C_n, *Is, *It, *Cs;
	int *R_n, *R_ns;
	int sizes = I_size > max_nei_nums ? max_nei_nums : I_size;
	R_n = new int[*R + 2];
	I_n = new pairs[sizes];
	C_n = new pairs[C_size + sizes];
	
	R_ns = R_n;
	for (int i = 0; i <= *R; ++i)
		*R_ns++ = R[i];
	++(*R_n);

	int heap_size = C_size + sizes;
	temp_memory += sizeof(pairs) * (heap_size + sizes);
	temp_memory += sizeof(int) * (*R + 2);
	if (temp_memory > max_temp_memory)
		max_temp_memory = temp_memory;

	int u, In_size = 0, Cn_size = 0;
	double r, q_n = 0;
	Is = I; It = I + I_size;
	Cs = C + C_size;
	while (Is != It)
	{
		u = (*Is).first;
		r = (*Is).second; ++Is;
		*R_ns = u;
		q_n = q * r;

		In_size = generateI(I, Is, I_size, q_n, u, I_n);
		if (*R_n + In_size <= k_core)
			continue;
		Cn_size = generateX(C, C_size, q_n, u, C_n);

		general_enumerate(R_n, q_n, I_n, In_size, C_n, Cn_size);

		(*Cs).first = u;
		(*Cs).second = r;
		++Cs; ++C_size;
	}
	temp_memory -= sizeof(pairs) * (heap_size + sizes);
	temp_memory -= sizeof(int) * (*R + 2);
	delete[] I_n;
	delete[] C_n;
	delete[] R_n;
}

void Uncertain_Clique::general_max_enumerate(int *R, double q, pairs *I, int I_size, pairs *C, int C_size)
{
	if (I_size == 0 && C_size == 0 && *R > MAX_SIZE)
	{
		//++maximal_cliques;
		MAX_SIZE = *R;
		return;
	}
	if (I_size == 0) return;
	pairs *I_n, *C_n, *Is, *It, *Cs;
	int *R_n, *R_ns;
	int sizes = I_size > max_nei_nums ? max_nei_nums : I_size;
	R_n = new int[*R + 2];
	I_n = new pairs[sizes];
	C_n = new pairs[C_size + sizes];

	R_ns = R_n;
	for (int i = 0; i <= *R; ++i)
		*R_ns++ = R[i];
	++(*R_n);

	int heap_size = C_size + sizes;
	temp_memory += sizeof(pairs) * (heap_size + sizes);
	temp_memory += sizeof(int) * (*R + 2);
	if (temp_memory > max_temp_memory)
		max_temp_memory = temp_memory;

	int u, In_size = 0, Cn_size = 0;
	double r, q_n = 0;
	Is = I; It = I + I_size;
	Cs = C + C_size;
	while (Is != It)
	{
		u = (*Is).first;
		r = (*Is).second; ++Is;
		*R_ns = u;
		q_n = q * r;
		if (*R_n + It - Is <= MAX_SIZE)
			break;
		In_size = generateI(I, Is, I_size, q_n, u, I_n);
		if (*R_n + In_size <= MAX_SIZE)
			continue;
		Cn_size = generateX(C, C_size, q_n, u, C_n);

		general_enumerate(R_n, q_n, I_n, In_size, C_n, Cn_size);

		(*Cs).first = u;
		(*Cs).second = r;
		++Cs; ++C_size;
	}
	temp_memory -= sizeof(pairs) * (heap_size + sizes);
	temp_memory -= sizeof(int) * (*R + 2);
	delete[] I_n;
	delete[] C_n;
	delete[] R_n;
}

int Uncertain_Clique::generateI(pairs * I, pairs * Is, int I_size, double q_n, int u, pairs * In)
{
	pairs *Isg, *It, *Ins;
	pairs *us, *ut;
	int a, b, c_nei = 0;
	double r, p;
	Isg = Is; It = I + I_size;
	us = adj[u];
	ut = adj[u + 1];
	Ins = In;
	while (Isg < It && us < ut)
	{
		a = (*Isg).first;
		b = (*us).first;
		if (a < b)
			++Isg;
		else if (a > b)
			++us;
		else
		{
			r = (*Isg++).second;
			p = (*us++).second;
			if (q_n * r * p >= eta)
			{
				(*Ins).first = a;
				(*Ins++).second = r * p;
				++c_nei;
				//++Ins;
			}
			//++Isg; ++us;
		}
	}
	return c_nei;
}

int Uncertain_Clique::generateX(pairs * C, int C_size, double q_n, int u, pairs * Cn)
{
	pairs *Cs, *Ct, *Cns;
	pairs *us, *ut;
	int a, b, c_nei = 0;
	double r, p;
	Cs = C; Ct = C + C_size;
	us = adj[u];
	ut = adj[u + 1];
	Cns = Cn;
	while (Cs < Ct && us < ut)
	{
		a = (*Cs).first;
		b = (*us).first;
		if (a < b)
			++Cs;
		else if (a > b)
			++us;
		else
		{
			r = (*Cs++).second;
			p = (*us++).second;
			if (q_n * r * p >= eta)
			{
				(*Cns).first = a;
				(*Cns++).second = r * p;
				++c_nei;
				//++Cns;
			}
			//++Cs; ++us;
		}
	}
	return c_nei;
}

void Uncertain_Clique::general_algorithm_of_cliques()
{
	pairs * I, *C, *Is;
	I = new pairs[V + 2];
	C = new pairs[V + 2];
	//memset(I, 0, sizeof(pairs) * (V + 1));
	memset(C, 0, sizeof(pairs) * (V + 2));
	max_memory += sizeof(pairs) * (V + 2) * 2;
	temp_memory = 0;

	Is = I;
	for (int i = 0; i <= V; ++i)
	{
		(*Is).first = i;
		(*Is).second = 1.0f;
		++Is;
	}
	int R[1] = { 0 };
	if (algorithm == 1)
		general_enumerate(R, 1.0f, I, V + 1, C, 0);
	else if (algorithm == 6)
		general_max_enumerate(R, 1.0f, I, V + 1, C, 0);
	delete[] I;
	delete[] C;
	//max_memory += max_temp_memory;
	//max_temp_memory = 0;
	temp_memory = 0;
	if (algorithm == 1){
		cout << "Maximal cliques\t" << maximal_cliques << endl;
		cout << "Maximum clique\t" << MAX_SIZE << endl;
		//cout << "max_temp_memory "<< max_temp_memory /1024 <<endl;
	}
	max_temp_memory = 0;
}

void Uncertain_Clique::get_cores()
{
	int *ds = degrees, max_d = 0;
	int *pos = new int[V + 2], i = 0;
	double *d1, *d2, p, pr_n = 0.0f;
	memset(pos, 0, sizeof(int) * (V + 2));

	d1 = new double[max_nei_nums + 2];
	d2 = new double[max_nei_nums + 2];

	for (i = 0; i <= V; ++i)
	{
		//cout << "ver " << i << endl;
		*ds = get_eta_degrees(i, d1, d2);
		if (*ds > max_d)
			max_d = *ds;
		++ds;
	}
	delete[] d1;
	delete[] d2;
	/*ds = degrees;
	for (int i = 0; i <= V; ++i)
	cout << i << " degrees " << *ds++ << endl;*/

	int *D = new int[max_d + 1];
	int *D_t = new int[max_d + 1];
	int *sort_d = new int[V + 1];
	memset(D, 0, sizeof(int)*(max_d + 1));
	ds = degrees;
	for (i = 0; i <= V; ++i)
		++D[*ds++];
	i = 0;
	for (int counts = 0, c; i <= max_d; ++i)
	{
		c = D[i];
		D[i] = counts;
		counts += c;
		D_t[i] = D[i];
		//	cout << D_t[i] << endl;
	}

	ds = degrees;
	for (i = 0; i <= V; ++i)
	{
		*pos++ = D_t[*ds];
		sort_d[D_t[*ds++]++] = i;
	}
	pos -= (V + 1);

	int u, v, v_d, v_dn;
	pairs *neis, *neit;
	ds = sort_d; i = 0;
	for (int vd_pos, vdn_pos, x; i <= V; ++i)
	{
		u = ds[i];
		if (degrees[u] >= k_core) break;

		neis = adj[u];
		neit = adj[u + 1];
		while (neis != neit)
		{
			v = (*neis).first;
			v_d = degrees[v];
			if (v_d >= k_core)
			{
				p = (*neis).second;
				prob[v][0] /= (1.0f - p);
				for (int j = 1; j <= v_d; ++j)
				{
					prob[v][j] = (prob[v][j] - p * prob[v][j - 1]) / (1.0f - p);
				}
				pr_n = 0.0f;
				v_dn = v_d;
				for (int j = 0; j <= v_d; ++j)
				{
					pr_n += prob[v][j];
					if (pr_n > 1 - eta)
					{
						v_dn = j;
						break;
					}
				}
				degrees[v] = v_dn;
				while (v_dn < v_d && v_d >= k_core)
				{
					vd_pos = pos[v];
					vdn_pos = D[v_d];
					if (vd_pos != vdn_pos)
					{
						x = sort_d[vdn_pos];
						sort_d[vd_pos] = sort_d[vdn_pos];
						sort_d[vdn_pos] = v;
						pos[v] = vdn_pos;
						pos[x] = vd_pos;
					}
					++D[v_d--];
				}
			}
			++neis;
		}
	}
	cout << "After k-core, remaining nodes are " << V + 1 - i << endl;
	/*--ds;
	for (; i <= V; ++i)
	cout << *ds++ << " ";
	cout << endl;*/
	ds = degrees;
	/*for (i = 0; i <= V; ++i)
	cout << i << " degrees " << *ds++ << endl;*/
	delete[] D;
	delete[] D_t;
	delete[] sort_d;
	delete[] pos;
}

int Uncertain_Clique::get_eta_degrees(int v, double *d1, double *d2)
{
	int v_size = ver[v], c_k = 0;
	if (v_size <= 0)
		return 0;
	double *ds1, *ds2, *ds_temp, p = 0.0f;
	ds1 = d1;
	ds2 = d2;
	pairs *es = adj[v];
	*ds2 = 1.0f;

	for (int i = 1; i < v_size; ++i)
	{
		p = (*es++).second;
		for (int j = 0; j <= i; ++j)
		{
			if (i == j)
				ds1[j] = p * ds2[j - 1];
			else if (j == 0)
				ds1[j] = (1 - p) * ds2[j];
			else
				ds1[j] = p * ds2[j - 1] + (1 - p) * ds2[j];
			//cout << i << " " << j << " " << ds[i][j] << endl;
		}
		/*for (int j = 0; j <= i; ++j)
			ds2[j] = ds1[j];*/
		ds_temp = ds1;
		ds1 = ds2;
		ds2 = ds_temp;
	}
	p = (*es++).second;
	for (int j = 0; j <= v_size; ++j)
	{
		if (v_size == j)
			prob[v][j] = p * ds2[j - 1];
		else if (j == 0)
			prob[v][j] = (1 - p) * ds2[j];
		else
			prob[v][j] = p * ds2[j - 1] + (1 - p) * ds2[j];
		//cout << v_size << " " << j << " " << prob[v][j] << endl;
	}

	p = 0.0f;
	for (int i = v_size; i >= 0; --i)
	{
		p += prob[v][i];
		if (p >= eta)
		{
			c_k = i;
			break;
		}
	}

	return c_k;
}

int Uncertain_Clique::core_prune1(int * PV)
{
	temp_memory = 0;
	queue<int> Q;
	double *d1, *d2 , q, p = 0.0f;
	d1 = new double[k_core + 1];
	d2 = new double[k_core + 1];
	int double_size = sizeof(double);
	temp_memory += double_size * (k_core + 1) * 2;

	for (int i = 0, k; i <= V; ++i)
	{
		//cout << "ver " << i << endl;
		k = is_k_core_or_not(i, d1, d2);
		if (k == 0)
			Q.emplace(i);
		PV[i] = k;
		/*if (i > 10)
		return 0;*/
	}
	delete[] d1;
	delete[] d2;

	int u = 0, v = 0, counts = 0;
	pairs *s, *t;
	while (!Q.empty())
	{
		u = Q.front(); Q.pop();
		s = adj[u];
		t = adj[u + 1];
		while (s < t)
		{
			v = (*s).first;
			if (PV[v] != 0)
			{
				p = (*s).second;
				eff_pro[v][0] /= (1.0f - p);
				q = eff_pro[v][0];
				for (int j = 1; j < k_core; ++j)
				{
					eff_pro[v][j] = (eff_pro[v][j] - p * eff_pro[v][j - 1]) / (1.0f - p);
					q += eff_pro[v][j];
				}
				if (q > 1 - eta)
				{
					Q.emplace(v);
					PV[v] = 0;
				}
			}
			++s;
		}
		++counts;
	}

	int *x = PV, *y = PV;
	for (int i = 0; i <= V; ++i)
	{
		if (*y++ != 0)
			*x++ = i;
	}
	cout << "Prune1: the left nodes " << V + 1 - counts << endl;
	
	if (max_temp_memory < temp_memory)
		max_temp_memory = temp_memory;
	return V + 1 - counts;
}

int Uncertain_Clique::is_k_core_or_not(int v, double *d1, double *d2)
{
	int v_size = ver[v], c_k = 0;
	if (v_size <= 0)
		return 0;
	if (v_size < k_core)
		return 0;

	pairs *es = adj[v];
	double *ds1, *ds2, *ds_temp, p = 0.0f;
	ds1 = d1;
	ds2 = d2;
	*ds2 = 1.0f;

	for (int i = 1; i < v_size; ++i)
	{
		p = (*es++).second;
		for (int j = 0; j <= i && j <= k_core; ++j)
		{
			if (i == j)
				ds1[j] = p * ds2[j - 1];
			else if (j == 0)
				ds1[j] = (1 - p) * ds2[j];
			else
				ds1[j] = p * ds2[j - 1] + (1 - p) * ds2[j];
			//cout << i << " " << j << " " << d[i][j] << endl;
		}
		/*for (int j = 0; j <= i && j <= k_core; ++j)
			ds2[j] = ds1[j];*/
		ds_temp = ds1;
		ds1 = ds2;
		ds2 = ds_temp;
	}
	p = (*es++).second;
	double q = 0.0f;
	for (int j = 0; j <= v_size && j < k_core; ++j)
	{
		if (v_size == j)
			eff_pro[v][j] = p * ds2[j - 1];
		else if (j == 0)
			eff_pro[v][j] = (1 - p) * ds2[j];
		else
			eff_pro[v][j] = p * ds2[j - 1] + (1 - p) * ds2[j];
		q += eff_pro[v][j];
		//cout << v_size << " " << j << " " << eff_pro[v][j] << endl;
	}
	if (q > 1 - eta)
		return 0;
	else
		return k_core;
}

int Uncertain_Clique::core_prune2(int * PV)
{
	temp_memory = 0;
	pairs *prn = new pairs[max_nei_nums], *p;
	double *pn = new double[V + 1];
	int pairs_size = sizeof(pairs);
	if (vert_nei_pro == NULL){
		vert_nei_pro = new pairs*[V + 1];
		max_memory += pairs_size * (V+1);
	}
	temp_memory += (pairs_size * max_nei_nums + sizeof(double) * (V+1));
	memset(pn, 0, sizeof(double) *(V + 1));
	pairs *s, *t;
	double q = 1.0, qn;
	queue<int> Q;
	for (int i = 0, n_size; i <= V; ++i)
	{
		vert_nei_pro[i] = NULL;
		n_size = ver[i];
		if (n_size < k_core)
		{
			PV[i] = 0; Q.emplace(i);
		}
		else
		{
			s = adj[i];
			t = adj[i + 1];
			p = prn;
			while (s < t)
			{
				p->first = s->first;
				p->second = s->second;
				++p;
				//*p++ = (*s).second;
				++s;
			}
			sort(prn, prn + n_size, fun_v);
			//qsort(prn, n_size, sizeof(double), myfunc);
			q = 1.0;
			p = prn;
			for (int x = 0; x < k_core; ++x)
			{
				//cout << "p " << p->second <<" " << p->first << endl;
				q *= p->second; ++p;
				//q *= *p++;
			}
			//return 0;
			if (q < eta)
			{
				PV[i] = 0; Q.emplace(i);
			}
			else
			{
				vert_nei_pro[i] = new pairs[n_size];
				max_memory += pairs_size * n_size;
				p = prn;
				for (int x = 0; x < n_size; ++x){
					vert_nei_pro[i][x].first = p->first;
					vert_nei_pro[i][x].second = p->second;
					++p;
				}
				PV[i] = k_core;
				pn[i] = q;
			}
		}
	}

	int u, v, v_size, v_pos, counts = 0;
	p = prn;
	while (!Q.empty())
	{
		u = Q.front(); Q.pop();
		s = adj[u];
		t = adj[u + 1];
		PV[u] = -1;
		++counts;
		while (s < t)
		{
			v = (*s).first;
			q = (*s).second;
			++s;
			v_size = ver[v];
			v_pos = PV[v];
			if (v_pos <= 0)
				continue;

			if (v_pos < v_size)
			{
				qn = vert_nei_pro[v][v_pos - 1].second;
				if (qn >= q ){
					if (qn > q)
						continue;
					//当qn = q考虑u被计算了还是没有被计算，如果没被计算跳过，计算了则更新
					int n = vert_nei_pro[v][v_pos - 1].first;
					if (n > u)
						continue;
				}
			
				int n = vert_nei_pro[v][v_pos].first;
				while (PV[n] == -1)
				{
					++v_pos;
					if (v_pos == v_size)
						break;
					n = vert_nei_pro[v][v_pos].first;
				}
				if (v_pos == v_size)
				{
					PV[v] = 0;
					Q.emplace(v);
					continue;
				}
				pn[v] *= vert_nei_pro[v][v_pos].second / q;
				PV[v] = ++v_pos;
				if (pn[v] < eta)
				{
					PV[v] = 0;
					Q.emplace(v);
				}
			}
			else
			{
				PV[v] = 0;
				Q.emplace(v);
			}
		}
	}


	int *x = PV, *y = PV;
	for (int i = 0; i <= V; ++i)
	{
		if (*y++ > 0)
			*x++ = i;
	}
	delete[] prn;
	delete[] pn;
	if (max_temp_memory < temp_memory)
		max_temp_memory = temp_memory;
	cout << "Prune2: the left nodes " << V + 1 - counts << endl;
	return V + 1 - counts;
}

void Uncertain_Clique::free_Dynamic_pointers()
{
	int double_size = sizeof(double);
	int pairs_size = sizeof(pairs);
	if (prob != NULL)
	{
		for (int i = 0; i <= V; ++i)
		{
			if (prob[i] != NULL){
				delete[] prob[i];
				prob[i] = NULL;
				max_memory -= double_size * (ver[i] + 1);
			}
		}
		delete[] prob;
		prob = NULL;
		max_memory -= double_size *(V+1);
	}
	if (eff_pro != NULL)
	{
		for (int i = 0; i <= V; ++i)
		{
			if (eff_pro[i] != NULL)
				delete[] eff_pro[i];
		}
		delete[] eff_pro;
		eff_pro = NULL;
	}

	if (vert_nei_pro != NULL)
	{
		for (int i = 0; i <= V; ++i)
		{
			if (vert_nei_pro[i] != NULL)
				delete[] vert_nei_pro[i];
		}
		delete[] vert_nei_pro;
		vert_nei_pro = NULL;
	}
}

void Uncertain_Clique::maximun_clique(int * R, double q, pairs * I, int I_size, pairs * C, int C_size)
{
	if (I_size == 0 && C_size == 0 && *R > MAX_SIZE)
	{
		MAX_SIZE = *R;
		return;
	}
	if (I_size == 0)
		return;
	int nums = 0, set_NO = *R + 1;
	/*int col = color_vertices(I, I_size, set_NO);
	if (*R + col <= MAX_SIZE)
	return;*/

	pairs *I_n, *C_n, *Is, *It, *Cs;
	int *R_n, *R_ns;
	int sizes = I_size > max_nei_nums ? max_nei_nums : I_size + 1;
	R_n = new int[*R + 2];
	I_n = new pairs[sizes];
	C_n = new pairs[C_size + sizes];
	int heap_size = C_size + sizes;

	R_ns = R_n;
	for (int i = 0; i <= *R; ++i)
		*R_ns++ = R[i];
	++(*R_n);

	int u, In_size = 0, Cn_size = 0, max_v;
	double r, q_n = 0;
	Is = I; It = I + I_size;
	Cs = C + C_size;
	max_v = (*(It - 1)).first;

	bool flag = true;
	/*int col = count_colors(Is, I_size - nums);
	if (*R + col <= MAX_SIZE)
	flag = false;
	if (flag)
	{
	int theta = probability_of_vR(q, Is, I_size - nums, coloring_nums);
	if (*R + theta <= MAX_SIZE)
	flag = false;
	}*/
	while (Is != It)
	{
		sub_set[(*Is).first] = set_NO;
		++Is;
	}
	Is = I;
	int col = 0;
	int *test = new int[coloring_nums];
	memset(test, 0, sizeof(int) * coloring_nums);
	col = count_colors(Is, I_size - nums, test);

	while (Is != It && flag)
	{
		//col = count_colors(Is, I_size - nums);
		//int col = color_vertices(Is, I_size - nums, set_NO);
		if (*R + col <= MAX_SIZE)
			break;
		if (*R != 0)
		{
			int theta = probability_of_vR(q, Is, I_size - nums, coloring_nums);
			if (*R + theta <= MAX_SIZE)
				break;
			/*int min = require_degrees(R, q, coloring_nums, set_NO, max_v);
			if (*R + min <= MAX_SIZE)
				break;*/
		}
		u = (*Is).first;
		r = (*Is).second;
		++Is; ++nums;
		*R_ns = u;
		q_n = q * r;
		sub_set[u] = set_NO - 1;
		if (--test[colors[u]] <= 0)
			--col;
		In_size = generateI(I, Is, I_size, q_n, u, I_n);

		int dertas = provide_degrees(I_n, In_size,  q_n);
		if (*R_n + dertas > MAX_SIZE)
		{
			//In_size = generateI(I, Is, I_size, q_n, u, I_n);
			Cn_size = generateX(C, C_size, q_n, u, C_n);

			maximun_clique(R_n, q_n, I_n, In_size, C_n, Cn_size);
		}

		/*int dertas = provide_degrees(Is, I_size - ++nums, u, col, q_n, set_NO);
		if (*R_n + dertas <= MAX_SIZE)
		return;*/

		/*In_size = generateI(I, Is, I_size, q_n, u, I_n);
		Cn_size = generateX(C, C_size, q_n, u, C_n);

		maximun_clique(R_n, q_n, I_n, In_size, C_n, Cn_size);*/

		(*Cs).first = u;
		(*Cs).second = r;
		++Cs; ++C_size;
	}
	while (Is != It)
	{
		sub_set[(*Is++).first] = set_NO - 1;
	}

	delete[] test;
	delete[] C_n;
	delete[] I_n;
	delete[] R_n;
}

int Uncertain_Clique::color_vertices()
{
	temp_memory = 0;
	int *D = new int[max_nei_nums + 1];
	int *sorted_ver = new int[V + 1], *q;
	if (colors == NULL) colors = new int[V+1];
	memset(colors, -1, sizeof(int) * (V + 1));
	memset(D, 0, sizeof(int) * (max_nei_nums + 1));
	temp_memory = sizeof(int) * ((max_nei_nums + 1)  + V + 1);
	int max_dg = 0;
	for (int i = 0, d; i <= V; ++i)
	{
		d = ver[i];
		++D[d];
		max_dg = max(max_dg, d);
	}
	//printf("max_dg=%d\n", max_dg);

	for (int i = max_nei_nums, d = 0, counts = 0; i >= 0; --i)
	{
		d = D[i];
		D[i] = counts;
		counts += d;
	}
	//printf("max_dg=%d\n", max_dg);

	for (int i = 0, d; i <= V; ++i)
	{
		d = ver[i];
		sorted_ver[D[d]++] = i;
	}
	printf("max_dg=%d\n", max_dg);

	int size_int = sizeof(int);
	memset(D, 0, size_int * (max_nei_nums + 1));
	pairs *s, *t;
	int color = 0, color_nums = 1;
	for (int i = 0, v, u ; i <= V; ++i)
	{
		v = sorted_ver[i];
		s = adj[v];
		t = adj[v + 1];
		if (s == t) {
			colors[v] = 0;
			continue;
		}
		
		int max_c = 0, c, j;
		while (s < t)
		{
			u = (*s++).first;
			c = colors[u];
			if (c == -1) continue;
			max_c = max(max_c, c);
			++D[c];
		}
		for (j = 0; j <= max_c; ++j) {
			if (D[j] == 0)
				break;
		}
		colors[v] = j;
		if (j >= color_nums)
			++color_nums;
		memset(D, 0, size_int * (max_c + 1));
	}
	if (temp_memory > max_temp_memory)
		max_temp_memory = temp_memory;

	cout << "Color nums " << color_nums << endl;
	delete[] D;
	delete[] sorted_ver;
	return color_nums;
}

int Uncertain_Clique::color_vertices(pairs * Can, int c_size, int set_NO)
{
	int *D = new int[max_nei_nums + 1];
	int *sorted_ver = new int[c_size], *q;
	memset(D, 0, sizeof(int) * (max_nei_nums + 1));

	for (int i = 0, u; i < c_size; ++i)
	{
		u = Can[i].first;
		sub_set[u] = set_NO;
		pos[u] = i;
		degrees[u] = 0;
		colors[u] = -1;
	}

	pairs *s, *t;
	int max_d = 0;
	for (int i = 0, u, v, d; i < c_size; ++i)
	{
		u = Can[i].first;

		s = adj[u];
		t = adj[u + 1];
		while (s < t)
		{
			v = (*s++).first;
			if (sub_set[v] == set_NO)
				++degrees[u];
		}
		d = degrees[u];
		if (max_d < d)
			max_d = d;
	}
	for (int i = 0, d; i < c_size; ++i)
	{
		d = degrees[Can[i].first];
		++D[d];
	}

	for (int i = max_d, d = 0, counts = 0; i >= 0; --i)
	{
		d = D[i];
		//counts += D[i];
		D[i] = counts;
		counts += d;
	}

	for (int i = 0, u, d; i < c_size; ++i)
	{
		u = Can[i].first;
		d = degrees[u];
		sorted_ver[D[d]++] = u;
	}

	int size_int = sizeof(int);
	memset(D, 0, size_int * (max_nei_nums + 1));
	q = sorted_ver;
	int color = 0, color_nums = 1;
	for (int i = 0, v, u; i < c_size; ++i)
	{
		v = *q++;
		s = adj[v];
		t = adj[v + 1];
		if (s == t){
			colors[v] = 0;
			continue;
		}
	
		int max_c = 0, c, j;
		while (s <t)
		{
			u = (*s).first;
			if (sub_set[u] == set_NO)
			{
				c = colors[u];
				if (c != -1)
				{
					if (c > max_c)
						max_c = c;
					++D[c];
				}
			}
			++s;
		}
		for (j = 0; j <= max_c; ++j)
			if (D[j] == 0)
				break;
		colors[v] = j;
		if (j >= color_nums)
			++color_nums;
		memset(D, 0, size_int * (max_c + 1));

	}
	
	//cout << "color nums " << color_nums << endl;
	delete[] D;
	delete[] sorted_ver;
	return color_nums;
}

int Uncertain_Clique::count_colors(pairs *Can, int c_size, int *color_class)
{
	if (c_size <= 0)
		return 0;
	int i, c, u, nums;
	i = 0; nums = 0;
	for (; i < c_size; ++i)
	{
		u = Can[i].first;
		c = colors[u];
		if (color_class[c] == 0)
			++nums;
		++color_class[c];
	}
	return nums;
}

int Uncertain_Clique::provide_degrees(pairs * Can, int c_size, double q)
{
	if (c_size <= 0)
		return 0;
	int derta = 0;
	memset(max_clr_p, 0, sizeof(double) * (coloring_nums));

	pairs *s, *t;
	s = Can;
	t = Can + c_size;
	int v, c, max_c;
	double p;
	max_c = 0;
	while (s < t)
	{
		v = (*s).first;
		c = colors[v];
		p = (*s).second;
		if (max_clr_p[c] < p)
			max_clr_p[c] = p;
		if (max_c < c)
			max_c = c;
		++s;
	}
	sort(max_clr_p, max_clr_p + max_c + 1, myfunc1);
	p = q;
	for (int i = 0; i <= max_c; ++i)
	{
		p *= max_clr_p[i];
		if (p < eta)
			break;
		++derta;
	}
	return derta;
}

int Uncertain_Clique::require_degrees(int *R, double pr,int color_nums,int set_NO, int max_v)
{
	int r_size = *R;
	if (r_size <= 0)
		return color_nums;
	double p = 1.0f;
	memset(max_clr_p, 0, sizeof(double) * (color_nums));
	pairs *s, *t;
	int min = color_nums, max_c = -1;
	for (int i = 1, u, v, c, add_in; i <= r_size; ++i)
	{
		u = R[i];
		max_c = -1;
		s = adj[u];
		t = adj[u + 1];
		while (s < t)
		{
			v = (*s).first;
			if (v > max_v)
				break;
			if (sub_set[v] == set_NO)
			{
				p = (*s).second;
				c = colors[v];
				if (max_clr_p[c] < p)
					max_clr_p[c] = p;
				if (c > max_c)
					max_c = c;
			}
			++s;
		}
		++max_c;
		sort(max_clr_p, max_clr_p + (max_c), myfunc1);
		p = pr;
		add_in = max_c;
		for (int j = 0; j < max_c; ++j)
		{
			p *= max_clr_p[j];
			if (p < eta)
			{
				add_in = j;
				break;
			}
		}
		if (min > add_in)
			min = add_in;
		memset(max_clr_p, 0, sizeof(double) * (max_c));
	}
	return min;
}

int Uncertain_Clique::probability_of_vR(double pr, pairs * Can, int size, int color_nums)
{
	if (size <= 0)
		return 0;
	int u, c, max_c = -1;
	double p;
	memset(max_clr_p, 0, sizeof(double) * (color_nums));

	pairs *s, *t;
	s = Can;
	t = Can + size;
	while ( s < t)
	{
		u = (*s).first;
		p = (*s).second;
		c = colors[u];
		if (max_c < c)
			max_c = c;
		if (c >= color_nums || c < 0)
			cout << c << endl;
		++s;

		if (max_clr_p[c] < p)
			max_clr_p[c] = p;
	}

	sort(max_clr_p, max_clr_p + (++max_c), myfunc1);
	int theta = max_c + 1;
	p = pr;
	for (int i = 0; i < max_c; ++i)
	{
		pr *= max_clr_p[i];
		if (pr < eta)
		{
			theta = i + 1;
			break;
		}
	}
	return theta;
}

void Uncertain_Clique::Max_Pclq()
{
	//Order();
	int c_size;
	bool found = false;
	int *R = new int[max_nei_nums];
	memset(R, 0, sizeof(int) * max_nei_nums);
	if (omega == NULL)
		omega = new int[V + 1];
	memset(omega, 0, sizeof(int) * (V + 1));
	double p = 1.0f;
	pairs * Can = new pairs[max_nei_nums];
	max_memory = 0;
	max_temp_memory = 0;
	max_memory += (sizeof(pairs) * (max_nei_nums) * 2 + sizeof(int) * (V+1));
	for (int i = V; i >= 0; --i)
	{
		found = false;
		c_size = Get_Can(i, i, Can, 1.0);
		*R = i;
		temp_memory = 0;
		//cout << "Max_Pclq " << i << endl;
		search(Can, c_size, R, 1, p ,found);
		omega[i] = MAX_SIZE;
	}
	delete[] Can;
	delete[] R;
}

int Uncertain_Clique::Get_Can(int i, int v, pairs * Can, double p)
{
	int vs, vt, c_size;
	pairs *ns, *nt;
	vs = i; vt = V;
	ns = adj[v];
	nt = adj[v+1] - 1;
	if (nt->first < i) return 0;
	c_size = 0;
	while (vs <= vt && ns <= nt)
	{
		if (vs < ns->first)
			++vs;
		else if (vs > ns->first)
			++ns;
		else
		{
			if (p * ns->second >= eta) {
				Can[c_size].first = vs;
				Can[c_size].second = p * ns->second;
				++c_size;
			}
			++vs; ++ns;
		}
	}
	return c_size;
}

void Uncertain_Clique::search(pairs * Can, int c_size, int * R, int r_size, double p, bool & found)
{
	if (c_size == 0)
	{
		if (r_size > MAX_SIZE)
		{
			MAX_SIZE = r_size;
			found = true;
		}
		return;
	}
	int i, u, rn_size, cn_size, n_size;
	double pn;
	n_size = c_size > max_nei_nums ? max_nei_nums : c_size;
	pairs *res = new pairs[n_size];
	temp_memory += sizeof(pairs) * n_size;
	if (temp_memory > max_temp_memory)
		max_temp_memory = temp_memory;
	i = 0;
	rn_size = r_size + 1;
	while (i < c_size)
	{
		u = Can[i].first;
		if (c_size + r_size < MAX_SIZE + i)
			break;
		if (omega[u] + r_size < MAX_SIZE)
			break;
		pn = p * Can[i].second;
		R[r_size] = u; 
		++i;

		cn_size = build_Can(u, Can, i, c_size, pn, res);

		search(res, cn_size, R, rn_size, pn, found);
		if (found) break;
	}
	delete[] res;
	temp_memory -= sizeof(pairs) * n_size;
}

int Uncertain_Clique::build_Can(int v, pairs * Can, int st,int lt, double p, pairs * res)
{
	pairs *vs, *vt, *cs, *ct;
	vs = adj[v];
	vt = adj[v + 1];
	cs = Can + st;
	ct = Can + lt;
	int res_size = 0;
	while (vs < vt && cs < ct)
	{
		if (vs->first > cs->first)
			++cs;
		else if (vs->first < cs->first)
			++vs;
		else 
		{
			if (vs->second * p * cs->second >= eta) {
				res[res_size].first = vs->first;
				res[res_size].second = cs->second * vs->second;
				++res_size;
			} ++vs; ++cs;
		}
	}
	return res_size;
}

void Uncertain_Clique::compute(int k, double p, int alg, int prunes)
{
	assert(prunes >= 1 && prunes <= 2);
	assert(alg >= 1 && alg <= 6);
	assert(p >= 0 && p <= 1.0);
	assert(k > 0);
	k_core = k;
	eta = p;
	algorithm = alg;
	init_parameters(prunes);

	printf("k=%d, eta=%.4f, alg=%d, prune=%d\n", k, eta, alg, prunes);

	struct timeval start_tm, end_tm, prune_tm;
	double prune_time;
	gettimeofday(&start_tm, NULL);
	int *PV = NULL, pv_size = 0;
	int R[1] = {0};
	int *cores = new int[V]();
	switch (alg) {
	case 1:
		general_algorithm_of_cliques();
		break;
	case 2:
		//get_cores();
		colorful_topk_core(k_core, cores, pv_size);
		break;
	case 3:
		//输出团
		//outc = fopen("out_cliques.txt", "w");

		PV = new int[V + 2];
		memset(PV, 0, sizeof(int) * (V + 2));
		if (prunes == 1)
			pv_size = core_prune1(PV);
		else if (prunes == 2)
			pv_size = core_prune2(PV);
		
		gettimeofday(&prune_tm, NULL);
		prune_time = double(prune_tm.tv_sec - start_tm.tv_sec) * 1000.0f + 
			double(prune_tm.tv_usec - start_tm.tv_usec) / 1000.0f;
		printf("Prune%d, time\t%.3f ms\n", prunes, prune_time);

		free_Dynamic_pointers();
		gettimeofday(&prune_tm, NULL);
		max_memory += sizeof(int) * (V+2);
		if (pv_size != 0)
		{
			pairs * I, *C, *Is;
			I = new pairs[pv_size + 2];
			C = new pairs[pv_size + 2];
			//memset(I, 0, sizeof(I) * (c_size + 1));
			memset(C, 0, sizeof(pairs) * (pv_size + 2));

			Is = I;
			for (int i = 0; i < pv_size; ++i)
			{
				(*Is).first = PV[i];
				(*Is).second = 1.0;
				++Is;
			}
			temp_memory = 0;
			general_enumerate(R, 1.0, I, pv_size, C, 0);
			delete[] I;
			delete[] C;
			cout << "Maximal clique nums\t" << maximal_cliques << endl;
			cout << "The maximum clique size\t" << MAX_SIZE << endl;
		}
		if (PV != NULL)
		{
			delete[] PV;
			PV = NULL;
		}
		gettimeofday(&end_tm, NULL);
		printf("Enumerating time\t%.3f ms\n", double(end_tm.tv_sec - prune_tm.tv_sec) * 1000.0f + 
			double(end_tm.tv_usec - prune_tm.tv_usec) / 1000.0f);
		
		//fclose(outc);
		break;
	case 4:
		coloring_nums = color_vertices();
		if (max_clr_p == NULL)
			max_clr_p = new double[coloring_nums];
		PV = new int[V + 2];
		memset(PV, 0, sizeof(int) * (V + 2));
		if (prunes == 1)
			pv_size = core_prune1(PV);
		else if (prunes == 2)
			pv_size = core_prune2(PV);
		
		gettimeofday(&prune_tm, NULL);
		prune_time = double(prune_tm.tv_sec - start_tm.tv_sec) * 1000.0f + 
			double(prune_tm.tv_usec - start_tm.tv_usec) / 1000.0f;
		printf("Prune%d, time\t%.3f ms\n", prunes, prune_time);

		free_Dynamic_pointers();
	
		gettimeofday(&prune_tm, NULL);
		max_memory += sizeof(int) * (V+2);

		if (pv_size != 0)
		{
			pairs * I, *C, *Is;
			I = new pairs[pv_size + 2];
			C = new pairs[pv_size + 2];
			//memset(I, 0, sizeof(pairs) * (V + 1));
			memset(C, 0, sizeof(pairs) * (pv_size + 2));

			Is = I;
			for (int i = 0; i < pv_size; ++i)
			{
				//assert(PV[i] >= 0 && PV[i] <=V);
				(*Is).first = PV[i];
				(*Is).second = 1.0f;
				++Is;
			}
			delete[] PV; PV = NULL;
			temp_memory = 0;
			MAX_SIZE = k_core;
			maximun_clique(R, 1.0f, I, pv_size , C, 0);
			delete[] I;
			delete[] C;
			delete[] PV; PV = NULL;
		}
		gettimeofday(&end_tm, NULL);
		printf("Enumerating time\t%.3f ms\n", double(end_tm.tv_sec - prune_tm.tv_sec) * 1000.0f + 
			double(end_tm.tv_usec - prune_tm.tv_usec) / 1000.0f);
		cout << "The maximum clique size\t" << MAX_SIZE << endl;
		break;
	case 5:
		MAX_SIZE = 0;
		cout << "Max_Pclq" << endl;
		Max_Pclq();
		cout << "The maximum clique size\t" << MAX_SIZE << endl;
		break;
	case 6:
		MAX_SIZE = k_core;
		general_algorithm_of_cliques();
		cout << "The maximum clique size\t" << MAX_SIZE << endl;
		break;
	default:
		printf("The parameter 'alg' errors! \n");
		break;
	}

	gettimeofday(&end_tm, NULL);
	printf("Alg %d, time\t%.3f ms\n", alg, double(end_tm.tv_sec - start_tm.tv_sec) * 1000.0f + 
		double(end_tm.tv_usec - start_tm.tv_usec) / 1000.0f);

	graph_memory /= 1024;
	max_memory /= 1024;
	max_temp_memory /= 1024;
	printf("Graph memory\t%ld kb\n",graph_memory);
	printf("Max memory\t%ld kb\n", graph_memory + max_memory + max_temp_memory);
}

void Uncertain_Clique::colorful_topk_core(const int &k, int *left_vertices, int &left_size)
{
	printf("Colorful k-core decomposition\n");
	int color_nums = color_vertices();
	vector< vector<int> > nbr_index(V+1);
	vector< int > nbr_index_cnts(V+1);
	pair<double, int> * edges_new = new pair<double, int>[E*2];
	pair<double, int> ** adj_new = new pair<double, int>*[V+1];
	int nbr_max_clr = 0, cnt = 0;
	left_size = 0;

	int *vertices = new int[V+1]();
	for (int i=0; i <= V; ++i) vertices[i] = i;
	vertices_sort(vertices, V+1);
	printf("heuristic clique: %d\n", heuristic_clique(vertices, V+1));
	return ;

	for (int i = 0; i <= V; ++i){
		int u, c, d ;
		d = ver[i];
		adj_new[i] = edges_new + cnt;
		cnt += d;
		if (d <= 0) continue;
		nbr_max_clr = 0;
		for (int j = 0; j < d; ++j)
		{
			u = adj[i][j].first;
			c = colors[u];
			nbr_max_clr = max(nbr_max_clr, c);
			double p = adj[i][j].second;
			adj_new[i][j].first = p;
			adj_new[i][j].second = u;
		}
		sort(adj_new[i], adj_new[i] + d, greater< pair<double, int> >());
		nbr_index[i].resize(nbr_max_clr+1, -1);
		nbr_index_cnts[i] = nbr_max_clr;
	}

	/* for (int i = 0; i <= V; ++i){
		int u, c, d ;
		d = ver[i];
		adj_new[i] = edges_new + cnt;
		cnt += d;
		if (d <= 0) continue;
		nbr_max_clr = 0;
		for (int j = 0; j < d; ++j)
		{
			u = adj[i][j].first;
			c = colors[u];
			nbr_max_clr = max(nbr_max_clr, c);
		}
		nbr_index[i].resize(nbr_max_clr + 2, 0);
		nbr_index_cnts[i] = nbr_max_clr + 1;
		for (int j = 0; j < d; ++j)
		{
			u = adj[i][j].first;
			c = colors[u];
			nbr_index[i][c]++;
		}
		nbr_index[i][nbr_max_clr+1] = d;
		for (int j = 0, cnt = 0; j <= nbr_max_clr; ++j){
			c = nbr_index[i][j];
			nbr_index[i][j] = cnt;
			cnt += c;
		}
		for (int j = 0; j < d; ++j)
		{
			u = adj[i][j].first;
			double p = adj[i][j].second;;
			int pos ;
			c = colors[u];
			pos = nbr_index[i][c]++;
			adj_new[i][pos].first = p;
			adj_new[i][pos].second = u;
		}
		for (int j = nbr_max_clr; j > 0; --j)
			nbr_index[i][j] = nbr_index[i][j-1];
		nbr_index[i][0] = 0;
		for (int j = 0; j < nbr_max_clr; ++j)
			sort(adj_new[i] + nbr_index[i][j], adj_new[i] + nbr_index[i][j+1], greater< pair<double, int> >());
	} */

	// for (int i = 0; i <= V; ++i){
	// 	int d = ver[i];
	// 	for (int j = 0; j < d; ++j){
	// 		int u = adj_new[i][j].second;
	// 		double p = adj_new[i][j].first;
	// 		printf("%d\t%.4lf\t%d\n",colors[u], p, u);
	// 	}

	// 	if (d > 0) break;
	// }
	
	vector<long double> prob(V+1);
	vector<int> D; D.reserve(V/2);
	int max_cr = 0;
	for (int i = 0; i <= V; ++i){
		int u, c, d, cr = 0;
		long double q = 1.0, p = 1.0;
		d = ver[i];
		if (d < k) { 
			prob[i] = 0;
			nbr_index_cnts[i] = 0;
			D.push_back(i);
			continue;
		}
		for (int j = 0; j < d; ++j){
			u = adj_new[i][j].second;
			c = colors[u];
			if (nbr_index[i][c] == -1) {
				q = p * adj_new[i][j].first;
				if (q + 1e-16 < eta) break;
				nbr_index[i][c] = j;
				p = q;
				cr++;
			}
		}
		prob[i] = p;
		nbr_index_cnts[i] = cr;
		if (cr < k_core) D.push_back(i);
		max_cr = max(max_cr, cr);
	}
	//printf("v=%d, p=%.5Lf, k=%d, D_size=%ld\n", 1, prob[1], nbr_index_cnts[1], D.size());
	printf("max_cr=%d, k_core=%d\n", max_cr, k_core);
	for (size_t i = 0; i < D.size(); ++i){
		int v, u, c, d;
		v = D[i];
		d = ver[v];
		c = colors[v];
		for (int j = 0; j < d; ++j){
			u = adj_new[v][j].second;
			double p = adj_new[v][j].first;
			int cr = nbr_index_cnts[u];
			//printf("v=%d, u=%d, cru=%d\n", v, u, cr);
			if ( cr >= k_core) {
				//int c = colors[v];
				int du = ver[u];
				int index_c = nbr_index[u][c];
				// for (int x = 0; x < du; ++x){
				// 	int w = adj_new[u][x].second;
				// 	int cw = colors[w];
				// 	//if (nbr_index[u][cw] > -1)
				// 	printf("  u=%d, w=%d, cl=%d, index=%d, prob=%.5lf\n", 
				// 		u, w, cw, nbr_index[u][cw], adj_new[u][x].first);
				// }
				if (adj_new[u][index_c].second != v)
					continue;
				index_c++;
				nbr_index[u][c] = -1;
				int w = 0, cw = 0;
				while (index_c < du) {
					w = adj_new[u][index_c].second;
					cw = colors[w];
					if (nbr_index[u][cw] == -1){
						nbr_index[u][cw] = index_c;
						break;
					}
					++index_c;
				}
				//printf("index=%d\n", nbr_index[u][cw]);
				if (index_c < du){
					long double pro = prob[u];
					pro /= p;
					pro *= adj_new[u][index_c].first;
					if (pro + 1e-16 < eta) {
						nbr_index[u][cw] = -1;
						cr--;
					}
					else prob[u] = pro;
				}
				else cr--;
				
				//nbr_index_cnts[u] = max(cr, k_core);
				nbr_index_cnts[u] = cr;
				//printf("prob=%Lf, cr=%d, index=%d\n",prob[u], cr, nbr_index[u][cw]);

				if (cr < k_core) {
					//cout << "Left vertices: " << V+1 - D.size() << endl;
					D.push_back(u);
				}
			}
		}
		//break;
	}
	//cout << "D vertices: " << D.size() << endl; 
	cout << "Left vertices: " << V + 1 - (int) D.size() << endl; 
}

void Uncertain_Clique::vertices_sort(int *vertices, const int vt_size)
{
	printf("Starting heuristic\n");
	if (vt_size <= 0) return;
	int *bin = new int[max_nei_nums+1]();
	for (int i = 0; i < vt_size; ++i) {
		int v = vertices[i];
		bin[ver[v]]++;
	}
	for (int i = 0, cnt = 0; i <= max_nei_nums; ++i)
	{
		int d = bin[i];
		bin[i] = cnt;
		cnt += d;
	}
	// for (int i = 0; i <= max_nei_nums; ++i) {
	// 	//int v = vertices[i];
	// 	//printf("v=%d, d=%d\n", v, ver[v]);
	// 	printf("bin[%d]=%d\n", i, bin[i]);
	// }
	for (int i = 0; i <= V; ++i)
		vertices[bin[ver[i]]++] = i;

	// for (int i = 0; i < vt_size; ++i) {
	// 	int v = vertices[i];
	// 	printf("v=%d, d=%d\n", v, ver[v]);
	// }
	delete[] bin;
}

int Uncertain_Clique::heuristic_clique(int *vertices, const int vt_size)
{
	int max_clq = 0;
	int *visited = new int[V+1]();
	int can_size = 0;
	for (int i = vt_size - 1; i >= 0; --i){
		int v = vertices[i];
		int d = ver[v];
		int u, du, id, w, max_du = 0;
		int len = 1;
		can_size = d;
		if (d <= max_clq) continue;
		for (int j = 0; j < d; ++j){
			u = adj[v][j].first;
			du = ver[u];
			visited[u] = len;
			if (du > max_du){
				max_du = du;
				id = u;
			}
		}
		//printf("can size: %d\n", can_size);
		while (can_size > 0) {
			if (can_size + len <= max_clq) break;
			u = id;
			du = ver[u];
			can_size = 0;
			max_du = 0;
			//printf("u=%d, d=%d, len=%d\n",u, du, len);
			for (int j = 0; j < du; ++j){
				w = adj[u][j].first;
				if (visited[w] == len){
					visited[w] = len+1;
					int dw = ver[w];
					if (dw > max_du){
						max_du = dw;
						id = w;
					}
					can_size++;
				}
			}
			//cout << "Iterator " << len << endl;
			//printf("can size: %d\n", can_size);
			len++;
			//if (len > 10)
			//	break;
		}
		for (int j = 0; j < d; ++j){
			w = adj[v][j].first;
			du = ver[w];
			visited[w] = 0;
		}
		max_clq = max(max_clq, len);
		//break;
	}
	return max_clq;
}