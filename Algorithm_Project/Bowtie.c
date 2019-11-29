#include "Bowtie.h"
#include "FileHandler.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef D
#define N 50000
#define L 50
#define M 1000000
#define D 2
#endif // !D

//number of characters lexically smaller than qc in BWT(T)
int occ[4] = { 1 }; // ACGT는 각각 인덱스 0,1,2,3을 의미 ($다음에 A가 나오므로 A의 값은 1)
char nucleobase[4] = { 'A', 'C', 'G', 'T' };
str_with_index *arr;
nucleo_cnt *count_table; //number of qc characters before position idx in BWT(T)
char *res_str;

//두 배열 합병
void merge(str_with_index *list, int left, int mid, int right) {
	int i, j;
	str_with_index *buffer; //합병하기 위한 임시 버퍼
	int buffer_length, left_ptr, right_ptr, buf_ptr;
	/*
	buffer_length : 임시 버퍼의 길이
	left_ptr : 왼쪽 배열의 현재 인덱스
	right_ptr : 오른쪽 배열의 현재 인덱스
	buf_ptr : 임시 버퍼의 현재 인덱스
	*/
	buffer_length = right - left + 1; //인덱스 값들 초기화
	left_ptr = left;
	right_ptr = mid + 1;
	buf_ptr = 0;

	//버퍼 동적할당
	buffer = (str_with_index*)malloc(sizeof(str_with_index)*buffer_length);
	if (buffer == NULL) {
		printf("Memory allocation failed!\n");
		exit(1);
	}

	//정렬된 리스트 합침
	while (left_ptr <= mid && right_ptr <= right) { //한쪽 배열을 모두 배치할 때까지 정렬작업 반복
		if(strcmp(list[left_ptr].str, list[right_ptr].str) <= 0) {
			buffer[buf_ptr++] = list[left_ptr++];
		}
		else {
			buffer[buf_ptr++] = list[right_ptr++];
		}
	}

	//한쪽 배열의 원소를 모두 사용하면 남은 배열을 모두 버퍼에 복사 - 이미 정렬되어 있으므로
	if (left_ptr > mid) {
		for (i = right_ptr; i <= right; i++)
			buffer[buf_ptr++] = list[i];
	}
	else {
		for (i = left_ptr; i <= mid; i++)
			buffer[buf_ptr++] = list[i];
	}

	//원본 리스트에 정렬된 버퍼 복사
	for (i = left, j = 0; i <= right; i++, j++)
		list[i] = buffer[j];

	//임시 버퍼 메모리 삭제
	free(buffer);
}

void mergeSort(str_with_index *list, int left, int right) {
	int mid;
	if (left < right) {
		mid = (left + right) / 2;
		//왼쪽 절반과 오른쪽 절반에 대하여 merge sort
		mergeSort(list, left, mid);
		mergeSort(list, mid + 1, right);

		//정렬된 두 부분을 합병
		merge(list, left, mid, right);
	}
}

char* bwt(char *t) {
	int length_src = strlen(t);
	int length_res = length_src + 1;
	int i, j; //used in loop
	
	res_str = (char*)malloc(sizeof(char)*(length_res + 1)); //'$' and '\0' should be added, allocate size of length+2
	strcpy(res_str, t);
	res_str[length_src] = '$';
	res_str[length_res] = '\0'; //make res_str	

	//initialize count table
	count_table = (nucleo_cnt*)malloc(sizeof(nucleo_cnt)*(length_res));
	for (i = 0; i < length_src + 1; i++) {
		count_table[i].a_cnt = 0;
		count_table[i].c_cnt = 0;
		count_table[i].g_cnt = 0;
		count_table[i].t_cnt = 0;
	}

	//dynamic allocation of structure array
	arr = (str_with_index*)malloc(sizeof(str_with_index)*(length_res));
	//fill first index
	arr[0].idx = 0;
	arr[0].str = (char*)malloc(sizeof(char)*(length_src + 2));
	strcpy(arr[0].str, res_str);

	//fill other indicies
	for (i = 1; i < length_res; i++) {
		arr[i].idx = i;
		arr[i].str = (char*)malloc(sizeof(char)*(length_src + 2 - i));
		//save suffix
		strcpy(arr[i].str, res_str + i);
	}
	//filling structure array completed

	//merge sort with the suffix string
	mergeSort(arr, 0, length_src);

	//make res_str, occ table, and count table
	for (i = 0; i < length_res; i++) {
		//making occ table
		if (!occ[1] && arr[i].str[0] == 'C')
			occ[1] = i;
		else if (!occ[2] && arr[i].str[0] == 'G')
			occ[2] = i;
		else if (!occ[3] && arr[i].str[0] == 'T')
			occ[3] = i;
		//free(arr[i].str);

		//making count table
		int j = arr[i].idx - 1;
		count_table[i + 1] = count_table[i];

		if (j < 0)
			res_str[i] = '$';
		else {
			switch (t[j]) {
			case 'A':
				count_table[i + 1].a_cnt++;
				break;
			case 'T':
				count_table[i + 1].t_cnt++;
				break;
			case 'G':
				count_table[i + 1].g_cnt++;
				break;
			case 'C':
				count_table[i + 1].c_cnt++;
				break;
			}
			res_str[i] = t[j];
		}
	}
	return res_str;
}

//int search_occ_idx(char qc) {
//	int qc_idx;
//	switch (qc)
//	{
//	case 'A':
//		qc_idx = 0;
//		break;
//	case 'C':
//		qc_idx = 1;
//		break;
//	case 'G':
//		qc_idx = 2;
//		break;
//	case 'T':
//		qc_idx = 3;
//		break;
//	}
//	return qc_idx;
//}

//int count(int idx, char qc) { //인덱싱해서 빨리 할 필요성
//	int i, cnt = 0;
//	for (i = 0; i < idx; i++) {
//		if (res_str[i] == qc)
//			cnt++;
//	}
//	return cnt;
//}

int lf(int idx, char qc) {
	//lf(idx, qc) = occ(qc) + count(idx, qc)
	switch (qc) {
	case 'A':
		return occ[0] + count_table[idx].a_cnt;
	case 'C':
		return occ[1] + count_table[idx].c_cnt;
	case 'G':
		return occ[2] + count_table[idx].g_cnt;
	case 'T':
		return occ[3] + count_table[idx].t_cnt;
	}
	//return occ[search_occ_idx(qc)] + count(idx, qc);
}

void backtracking_search(reconstruct_struct* nc_arr, char *codes, char *ref_str, char *pattern, int mismatch_cnt, int top_range, int bot_range, int steps) {
	int length_pattern = strlen(pattern);
	//int top_range = 0, bot_range = length_res;
	int top_temp , bot_temp; //flag if top or bot pointer is available
	int original_steps = steps;
	//int mismatch_cnt = 0;
	int i, j;
	char target_ch;
	char *sub_pattern;

	for (i = length_pattern - 1; i >= 0; i--) {
		target_ch = pattern[i];

		top_temp = lf(top_range, target_ch);
		bot_temp = lf(bot_range, target_ch);

		if (top_temp >= bot_temp) {//perfect match does not occur
			//printf("mismatch occured in idx %d\n", arr[top_range].idx);
			//save current range and nucleobase
 			if (mismatch_cnt > 0) {
				//if the mismatched position makes another mismatches, don't count it and just terminate
				if (steps == original_steps)
					return;
			}
			mismatch_cnt++;
			if (mismatch_cnt <= D) { //if mismatch count is in acceptable range
				/* try again with a different character */
				//make substring of pattern string
				sub_pattern = malloc(sizeof(char)*(i + 2));
				strncpy(sub_pattern, pattern, i + 1);
				sub_pattern[i + 1] = '\0';
				//0:AC, 1:AG, 2:AT, 3:CG, 4:CT, 5:GT
				switch (target_ch) {
				case 'A':
					switch (codes[arr[top_temp].idx]) { //read code of the index and try the other code
					case '0':
						sub_pattern[i] = 'C';
						backtracking_search(nc_arr, codes, ref_str, sub_pattern, mismatch_cnt, top_range, bot_range, steps);
						break;
					case '1':
						sub_pattern[i] = 'G';
						backtracking_search(nc_arr, codes, ref_str, sub_pattern, mismatch_cnt, top_range, bot_range, steps);
						break;
					case '2':
						sub_pattern[i] = 'T';
						backtracking_search(nc_arr, codes, ref_str, sub_pattern, mismatch_cnt, top_range, bot_range, steps);
						break;
					}
					break;
				case 'C':
					switch (codes[arr[top_temp].idx]) {
					case '0':
						sub_pattern[i] = 'A';
						backtracking_search(nc_arr, codes, ref_str, sub_pattern, mismatch_cnt, top_range, bot_range, steps);
						break;
					case '3':
						sub_pattern[i] = 'G';
						backtracking_search(nc_arr, codes, ref_str, sub_pattern, mismatch_cnt, top_range, bot_range, steps);
						break;
					case '4':
						sub_pattern[i] = 'T';
						backtracking_search(nc_arr, codes, ref_str, sub_pattern, mismatch_cnt, top_range, bot_range, steps);
						break;
					}
					break;
				case 'G':
					switch (codes[arr[top_temp].idx]) {
					case '1':
						sub_pattern[i] = 'A';
						backtracking_search(nc_arr, codes, ref_str, sub_pattern, mismatch_cnt, top_range, bot_range, steps);
						break;
					case '3':
						sub_pattern[i] = 'C';
						backtracking_search(nc_arr, codes, ref_str, sub_pattern, mismatch_cnt, top_range, bot_range, steps);
						break;
					case '5':
						sub_pattern[i] = 'T';
						backtracking_search(nc_arr, codes, ref_str, sub_pattern, mismatch_cnt, top_range, bot_range, steps);
						break;
					}
					break;
				case 'T':
					switch (codes[arr[top_temp].idx]) {
					case '4':
						sub_pattern[i] = 'C';
						backtracking_search(nc_arr, codes, ref_str, sub_pattern, mismatch_cnt, top_range, bot_range, steps);
						break;
					case '5':
						sub_pattern[i] = 'G';
						backtracking_search(nc_arr, codes, ref_str, sub_pattern, mismatch_cnt, top_range, bot_range, steps);
						break;
					case '2':
						sub_pattern[i] = 'A';
						backtracking_search(nc_arr, codes, ref_str, sub_pattern, mismatch_cnt, top_range, bot_range, steps);
						break;
					}
					break;
				}
				free(sub_pattern);
				for (j = 0; j <= length_pattern - i; j++) { //save matched indecies
					//insert_node(&(ll_arr[arr[top_range].idx + i + j]), pattern[i + j]);
					switch (pattern[i + j])
					{
						//0:AC, 1:AG, 2:AT, 3:CG, 4:CT, 5:GT
					case 'A':
						nc_arr[arr[top_temp].idx + i + j].earlier_cnt++;
						if (nc_arr[arr[top_temp].idx + i + j].earlier_cnt > nc_arr[arr[top_temp].idx + i + j].later_cnt && nc_arr[arr[top_temp].idx + i + j].more_ch != 'A')
							nc_arr[arr[top_temp].idx + i + j].more_ch = 'A';
						break;
					case 'C':
						if (codes[arr[top_temp].idx + i + j] == '3' || codes[arr[top_temp].idx + i + j] == '4') {
							nc_arr[arr[top_temp].idx + i + j].earlier_cnt++;
							if (nc_arr[arr[top_temp].idx + i + j].earlier_cnt > nc_arr[arr[top_temp].idx + i + j].later_cnt && nc_arr[arr[top_temp].idx + i + j].more_ch != 'C')
								nc_arr[arr[top_temp].idx + i + j].more_ch = 'C';
						}
						else {
							nc_arr[arr[top_temp].idx + i + j].later_cnt++;
							if (nc_arr[arr[top_temp].idx + i + j].later_cnt > nc_arr[arr[top_temp].idx + i + j].earlier_cnt && nc_arr[arr[top_temp].idx + i + j].more_ch != 'C')
								nc_arr[arr[top_temp].idx + i + j].more_ch = 'C';
						}
						break;
					case 'G':
						if (codes[arr[top_temp].idx + i + j] == '5') {
							nc_arr[arr[top_temp].idx + i + j].earlier_cnt++;
							if (nc_arr[arr[top_temp].idx + i + j].earlier_cnt > nc_arr[arr[top_temp].idx + i + j].later_cnt && nc_arr[arr[top_temp].idx + i + j].more_ch != 'G')
								nc_arr[arr[top_temp].idx + i + j].more_ch = 'G';
						}
						else {
							nc_arr[arr[top_temp].idx + i + j].later_cnt++;
							if (nc_arr[arr[top_temp].idx + i + j].later_cnt > nc_arr[arr[top_temp].idx + i + j].earlier_cnt && nc_arr[arr[top_temp].idx + i + j].more_ch != 'G')
								nc_arr[arr[top_temp].idx + i + j].more_ch = 'G';
						}
						break;
					case 'T':
						nc_arr[arr[top_temp].idx + i + j].later_cnt++;
						if (nc_arr[arr[top_temp].idx + i + j].later_cnt > nc_arr[arr[top_temp].idx + i + j].earlier_cnt && nc_arr[arr[top_temp].idx + i + j].more_ch != 'T')
							nc_arr[arr[top_temp].idx + i + j].more_ch = 'T';
						break;
					}
				}
				return;
			}
		}
		top_range = top_temp;
		bot_range = bot_temp;
		steps++;
	}
	//if pattern matches, run codes below
	//weigh the nucleobase of searched index
	for (i = 0; i < length_pattern; i++) {
		for (j = 0; j < bot_range - top_range; j++) {
			switch (pattern[i])
			{
			case 'A':
				nc_arr[arr[top_range].idx + i].earlier_cnt++;
				if (nc_arr[arr[top_range].idx + i].earlier_cnt > nc_arr[arr[top_range].idx + i].later_cnt && nc_arr[arr[top_range].idx + i].more_ch != 'A')
					nc_arr[arr[top_range].idx + i].more_ch = 'A';
				break;
			case 'C':
				if (codes[arr[top_range].idx + i] == '3' || codes[arr[top_range].idx + i] == '4') {
					nc_arr[arr[top_range].idx + i].earlier_cnt++;
					if (nc_arr[arr[top_range].idx + i].earlier_cnt > nc_arr[arr[top_range].idx + i].later_cnt && nc_arr[arr[top_range].idx + i].more_ch != 'C')
						nc_arr[arr[top_range].idx + i].more_ch = 'C';
				}
				else {
					nc_arr[arr[top_range].idx + i].later_cnt++;
					if (nc_arr[arr[top_range].idx + i].later_cnt > nc_arr[arr[top_range].idx + i].earlier_cnt && nc_arr[arr[top_range].idx + i].more_ch != 'C')
						nc_arr[arr[top_range].idx + i].more_ch = 'C';
				}
				break;
			case 'G':
				if (codes[arr[top_range].idx + i] == '5') {
					nc_arr[arr[top_range].idx + i].earlier_cnt++;
					if (nc_arr[arr[top_range].idx + i].earlier_cnt > nc_arr[arr[top_range].idx + i].later_cnt && nc_arr[arr[top_range].idx + i].more_ch != 'G')
						nc_arr[arr[top_range].idx + i].more_ch = 'G';
				}
				else {
					nc_arr[arr[top_range].idx + i].later_cnt++;
					if (nc_arr[arr[top_range].idx + i].later_cnt > nc_arr[arr[top_range].idx + i].earlier_cnt && nc_arr[arr[top_range].idx + i].more_ch != 'G')
						nc_arr[arr[top_range].idx + i].more_ch = 'G';
				}
				break;
			case 'T':
				nc_arr[arr[top_range].idx + i].later_cnt++;
				if (nc_arr[arr[top_range].idx + i].later_cnt > nc_arr[arr[top_range].idx + i].earlier_cnt && nc_arr[arr[top_range].idx + i].more_ch != 'T')
					nc_arr[arr[top_range].idx + i].more_ch = 'T';
				break;
			}
		}
	}
}

//kmp functions below
void compute_sp(char* p, int* sp, int pl) {
	int i;
	int k = -1;
	sp[0] = -1;
	for (i = 1; i < pl; i++) {
		while (k > 0 && p[k + 1] != p[i]) //character doesn't match
			k = sp[k];
		if (p[k + 1] == p[i]) //if suffix character matches prefix character, count
			k++;
		sp[i] = k; //save counting value
	}
}

void kmp_matcher(reconstruct_struct* nc_arr, char *codes, char* t, char* p, int tl, int pl) {
	int *sp = (int*)malloc(sizeof(int)*pl); //suffix-prefix table array
	int q = 0; //numbers of characters matched
	int i, j; //used in loop
	compute_sp(p, sp, pl); //preprocess maximum prefix table

	for (i = 1; i <= tl; i++) {
		while (q > 0 && p[q + 1] != t[i]) {//next character does not match
			q = sp[q];
		}
		if (p[q + 1] == t[i]) //next character matches
			q++;
		if (q == pl - 1) {
			//if pattern matches, run codes below
			//weigh the nucleobase of searched index
			for (j = 0; j < pl; j++) {
				switch (p[j])
				{
				case 'A':
					nc_arr[i + j].earlier_cnt++;
					if (nc_arr[i + j].earlier_cnt > nc_arr[i + j].later_cnt && nc_arr[i + j].more_ch != 'A')
						nc_arr[i + j].more_ch = 'A';
					break;
				case 'C':
					if (codes[i + j] == '3' || codes[i + j] == '4') {
						nc_arr[i + j].earlier_cnt++;
						if (nc_arr[i + j].earlier_cnt > nc_arr[i + j].later_cnt && nc_arr[i + j].more_ch != 'C')
							nc_arr[i + j].more_ch = 'C';
					}
					else {
						nc_arr[i + j].later_cnt++;
						if (nc_arr[i + j].later_cnt > nc_arr[i + j].earlier_cnt && nc_arr[i + j].more_ch != 'C')
							nc_arr[i + j].more_ch = 'C';
					}
					break;
				case 'G':
					if (codes[i + j] == '5') {
						nc_arr[i + j].earlier_cnt++;
						if (nc_arr[i + j].earlier_cnt > nc_arr[i + j].later_cnt && nc_arr[i + j].more_ch != 'G')
							nc_arr[i + j].more_ch = 'G';
					}
					else {
						nc_arr[i + j].later_cnt++;
						if (nc_arr[i + j].later_cnt > nc_arr[i + j].earlier_cnt && nc_arr[i + j].more_ch != 'G')
							nc_arr[i + j].more_ch = 'G';
					}
					break;
				case 'T':
					nc_arr[i + j].later_cnt++;
					if (nc_arr[i + j].later_cnt > nc_arr[i + j].earlier_cnt && nc_arr[i + j].more_ch != 'T')
						nc_arr[i + j].more_ch = 'T';
					break;
				}
			}
			q = sp[q];
		}
	}
	free(sp);
}