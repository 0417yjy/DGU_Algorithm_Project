#include "FileHandler.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef N
//modify line 116
#define N 50000
#define L 50
#define M 1000000
#endif

int makeFileofString(char *filename, char *str) {
	FILE *f;
	f = fopen(filename, "w");
	if (f == NULL) {
		printf("fopen failed!\n");
		return 1;
	}

	fputs(str, f);
	fclose(f);
	return 0;
}
int readFiletoString(char *filename, char *str) {
	FILE *f;
	f = fopen(filename, "r");
	if (f == NULL) {
		printf("fopen failed!\n");
		return 1;
	}

	fscanf(f, "%s", str);
	fclose(f);
	return 0;
}

void generateSequence(char *ref, char *str_name) {
	int i, tmp, tmp2;
	char *bases_code = (char *)malloc(sizeof(int)*(N + 1));
	printf("\nFile \"%s\" doesn't exist. Make a new file.\n", str_name);
	printf("Generating random sequence...");
	for (i = 0; i < N; i++) {
		tmp = rand() % 4;
		tmp2 = rand() % 3;
		//0:AC, 1:AG, 2:AT, 3:CG, 4:CT, 5:GT
		switch (tmp) {
		case 0:
			ref[i] = 'A';
			switch (tmp2) {
			case 0:
				bases_code[i] = '0';
				break;
			case 1:
				bases_code[i] = '1';
				break;
			case 2:
				bases_code[i] = '2';
				break;
			}
			break;
		case 1:
			ref[i] = 'T';
			switch (tmp2) {
			case 0:
				bases_code[i] = '2';
				break;
			case 1:
				bases_code[i] = '4';
				break;
			case 2:
				bases_code[i] = '5';
				break;
			}
			break;
		case 2:
			ref[i] = 'G';
			switch (tmp2) {
			case 0:
				bases_code[i] = '1';
				break;
			case 1:
				bases_code[i] = '3';
				break;
			case 2:
				bases_code[i] = '5';
				break;
			}
			break;
		case 3:
			ref[i] = 'C';
			switch (tmp2) {
			case 0:
				bases_code[i] = '0';
				break;
			case 1:
				bases_code[i] = '3';
				break;
			case 2:
				bases_code[i] = '4';
				break;
			}
			break;
		}
	}
	ref[N] = '\0';
	bases_code[N] = '\0';
	printf("Completed\n");
	makeFileofString(str_name, ref);
	makeFileofString("bases_code.txt", bases_code);
	free(bases_code);
}

void readSequence(char *ref, char *myseq) {
	char str_name[] = "input_50000.txt";
	int i, tmp;

	//read input str
	printf("Reading %s...", str_name);
	if (readFiletoString(str_name, ref)) { //if there isn't  "input.txt", make new one
		generateSequence(ref, str_name);
		strcpy(myseq, ref);
	}
	else { //if opening input.txt file succeeded, save the string in ReferenceSeq
		strcpy(myseq, ref);
		printf("Completed\n");
	}
}

void makeSnp(char *str) {
	int i;
	char *codes = (char*)malloc(sizeof(char)*(N + 1));
	int *snps = (int*)malloc(sizeof(int)*(N/100 + 1));
	FILE *f;
	printf("Reading base codes...");
	readFiletoString("bases_code.txt", codes);
	printf("Completed\n");

	//make snps (1% of input string)
	printf("Making SNPs...");
	for (i = 0; i < N / 100; i++) {
		int rnd_index = getRandomNumber(0, N - L);
		//0:AC, 1:AG, 2:AT, 3:CG, 4:CT, 5:GT
		switch (codes[rnd_index]) {
		case '0':
			switch (str[rnd_index]) {
			case 'A':
				str[rnd_index] = 'C';
				break;
			case 'C':
				str[rnd_index] = 'A';
				break;
			}
			break;
		case '1':
			switch (str[rnd_index]) {
			case 'A':
				str[rnd_index] = 'G';
				break;
			case 'G':
				str[rnd_index] = 'A';
				break;
			}
			break;
		case '2':
			switch (str[rnd_index]) {
			case 'A':
				str[rnd_index] = 'T';
				break;
			case 'T':
				str[rnd_index] = 'A';
				break;
			}
			break;
		case '3':
			switch (str[rnd_index]) {
			case 'G':
				str[rnd_index] = 'C';
				break;
			case 'C':
				str[rnd_index] = 'G';
				break;
			}
			break;
		case '4':
			switch (str[rnd_index]) {
			case 'T':
				str[rnd_index] = 'C';
				break;
			case 'C':
				str[rnd_index] = 'T';
				break;
			}
			break;
		case '5':
			switch (str[rnd_index]) {
			case 'G':
				str[rnd_index] = 'T';
				break;
			case 'T':
				str[rnd_index] = 'G';
				break;
			}
			break;
		}
		snps[i] = rnd_index;
	}
	snps[N / 100] = '\0';
	printf("Completed\n");
	printf("Making myGenome.txt...");
	if (makeFileofString("myGenome.txt", str))
		exit(0);

	f = fopen("SNPs_indicies.txt", "w");
	for (i = 0; i < N / 100; i++) {
		fprintf(f, "%d %c\n", snps[i], codes[snps[i]]);
	}
	fclose(f);
	free(codes);
	free(snps);
	printf("Completed\n");
}

void makeShortReads(char *str) {
	char short_name[] = "reads.txt";
	char tmp_str[L + 1];
	int i, j, cnt = 0;
	int *coverages = (int*)malloc(sizeof(int)*N);
	if (coverages == NULL) {
		printf("Memory allocation failed!\n");
		exit(0);
	}
	FILE *fr, *fc;

	for (i = 0; i < N; i++)
		*(coverages+i) = 0;

	//make a new file named "read.txt"
	fr = fopen(short_name, "w");
	if (fr == NULL) {
		printf("fopen failed!\n");
		exit(1);
	}

	//make a new file named "coverages.txt"
	fc = fopen("coverages.txt", "w");
	if (fc == NULL) {
		printf("fopen failed!\n");
		exit(1);
	}

	//make short read
	printf("Making short strings...");
	for (i = 0; i < M; i++) {
		int start_idx = getRandomNumber(0, N - L); //get random int number for index
		strncpy(tmp_str, &str[start_idx], L);
		tmp_str[L] = '\0'; //make tmp_str by adding '\0' in last index
		for (j = 0; j < L; j++)
			coverages[start_idx + j]++;
		fprintf(fr, "%s\n", tmp_str);
	}

	for (i = 0; i < N; i++)
		fprintf(fc, "%d\n", coverages[i]);

	printf("Completed\n");
	fclose(fc);
	fclose(fr);
}