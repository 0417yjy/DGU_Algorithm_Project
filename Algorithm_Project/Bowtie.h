#pragma once

typedef struct str_with_index {
	char * str;
	int idx;
} str_with_index;

typedef struct nucleo_cnt {
	int a_cnt;
	int c_cnt;
	int g_cnt;
	int t_cnt;
} nucleo_cnt;

typedef struct reconstruct_struct {
	char more_ch;
	int earlier_cnt;
	int later_cnt;
} reconstruct_struct;

void merge(str_with_index *list, int left, int mid, int right);
void mergeSort(str_with_index *list, int left, int right);
char* bwt(char *t);
int lf(idx, qc); //last to first function
void backtracking_search(reconstruct_struct* nc_arr, char *codes, char *ref_str, char *pattern, int mismatch_cnt, int top_range, int bot_range, int steps);
void compute_sp(char* p, int* sp, int pl);
void kmp_matcher(reconstruct_struct* nc_arr, char *codes, char* t, char* p, int tl, int pl);