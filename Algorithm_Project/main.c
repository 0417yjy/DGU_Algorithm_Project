#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "FileHandler.h"
#include "Bowtie.h"
#define N 50000
#define L 50
#define M 1000000

int getRandomNumber(int min, int max) {
	//return (int)(min + (double)rand()*rand() / (RAND_MAX*RAND_MAX) * (max - min + 1));
	return min + (rand()*rand()) % (max - min + 1);
}

int main(void) {
	char *ReferenceSeq = (char*)malloc(sizeof(char)*(N + 1));
	char *MyGenomeSeq = (char*)malloc(sizeof(char)*(N + 1));
	char *bwt_reference;
	char short_read_buf[L + 1];
	int mismatch_cnt = 0;
	clock_t start_time, end_time;
	double time_taken;
	FILE *f;
	if (!MyGenomeSeq || !ReferenceSeq) {
		printf("Memory allocation failed. close the program..\n");
		exit(0);
	}

	srand(time(NULL));
	
	//read sequence and make files
	readSequence(ReferenceSeq, MyGenomeSeq);
	makeSnp(MyGenomeSeq);
	makeShortReads(MyGenomeSeq);

	//start counting time (bwt)
	printf("\n*****Counting time starts here (bwt)*****\n\n");
	start_time = clock();

	//read short reads and reconstruct myGenome
	printf("Transforming Reference Sequence...");
	//bwt_reference = bwt(ReferenceSeq);
	bwt_reference = bwt(MyGenomeSeq);
	printf("Completed\n");

	//printf("Making file of transformed string...");
	//makeFileofString("bwt_ref.txt", bwt_reference);
	//printf("Completed\n");

	//free(ReferenceSeq);
	free(MyGenomeSeq);
	printf("Initializing arrays...");
	reconstruct_struct *reconstruct_arr = (reconstruct_struct*)malloc(sizeof(reconstruct_struct)*N);
	for (int i = 0; i < N; i++) {
		reconstruct_arr[i].more_ch = ' ';
		reconstruct_arr[i].earlier_cnt = 0;
		reconstruct_arr[i].later_cnt = 0;
	}
	char *codes = (char*)malloc(sizeof(char)*(N + 1));
	readFiletoString("bases_code.txt", codes);
	printf("Completed\n");

	printf("Reading short reads and matching string...\n");
	f = fopen("reads.txt", "r");
	if (f == NULL) {
		printf("fopen failed!\n");
		exit(-1);
	}
	for (int i = 0; i < M; i++) {
		fscanf(f, "%s\n", short_read_buf);
		backtracking_search(reconstruct_arr, codes, bwt_reference, short_read_buf, 0, 0, N + 1, 0);
		//printf("A pattern matched completed %d / %d\n", i + 1, M);
		if (i == M / 10)
			printf("10%% achieved\n");
		else if (i == M*2 / 10)
			printf("20%% achieved\n");
		else if (i == M*3 / 10)
			printf("30%% achieved\n");
		else if (i == M*4 / 10)
			printf("40%% achieved\n");
		else if (i == M*5 / 10)
			printf("50%% achieved\n");
		else if (i == M*6 / 10)
			printf("60%% achieved\n");
		else if (i == M*7 / 10)
			printf("70%% achieved\n");
		else if (i == M*8 / 10)
			printf("80%% achieved\n");
		else if (i == M*9 / 10)
			printf("90%% achieved\n");
		else if (i == M - 1)
			printf("100%% achieved\n");
	}
	fclose(f);
	printf("Completed\n");

	printf("Reconstructing myGenome string...");
	MyGenomeSeq = (char*)malloc(sizeof(char)*(N + 1));
	for (int i = 0; i < N; i++)
		MyGenomeSeq[i] = reconstruct_arr[i].more_ch;
	MyGenomeSeq[N] = '\0';
	printf("Completed\n");

	//end counting time
	printf("\n*****Counting time ends here***** (bwt)\n");
	end_time = clock();
	time_taken = ((double)end_time - start_time) / CLOCKS_PER_SEC; //get execution time
	printf("Execution time = %.3lf\n\n", time_taken);

	free(reconstruct_arr);
	free(codes);
	free(bwt_reference);

	printf("Make a file of myGenome string...");
	makeFileofString("reconstructed_myGenome.txt", MyGenomeSeq);
	printf("Completed\n");

	printf("Caculating accuracy of reconstructed string...\n");
	readFiletoString("input_50000.txt", ReferenceSeq);
	for (int i = 0; i < N; i++)
		if (MyGenomeSeq[i] != ReferenceSeq[i])
			mismatch_cnt++;
	printf("Accuracy with reference: %lf%%\n", (1 - ((double)mismatch_cnt / N)) * 100);
	mismatch_cnt = 0;
	readFiletoString("myGenome.txt", ReferenceSeq);
	for (int i = 0; i < N; i++)
		if (MyGenomeSeq[i] != ReferenceSeq[i])
			mismatch_cnt++;
	printf("Accuracy with myGenome: %lf%%\n", (1 - ((double)mismatch_cnt / N)) * 100);

	//start counting time (KMP match)
	printf("\n*****Counting time starts here (kmp)*****\n\n");
	start_time = clock();

	printf("Initializing arrays...");
	reconstruct_arr = (reconstruct_struct*)malloc(sizeof(reconstruct_struct)*N);
	for (int i = 0; i < N; i++) {
		reconstruct_arr[i].more_ch = ' ';
		reconstruct_arr[i].earlier_cnt = 0;
		reconstruct_arr[i].later_cnt = 0;
	}
	*codes = (char*)malloc(sizeof(char)*(N + 1));
	readFiletoString("bases_code.txt", codes);
	printf("Completed\n");

	printf("Reading short reads and matching string...\n");
	f = fopen("reads.txt", "r");
	if (f == NULL) {
		printf("fopen failed!\n");
		exit(-1);
	}
	for (int i = 0; i < M; i++) {
		fscanf(f, "%s\n", short_read_buf);
		//backtracking_search(reconstruct_arr, codes, bwt_reference, short_read_buf, 0, 0, N + 1, 0);
		kmp_matcher(reconstruct_arr, codes, ReferenceSeq, short_read_buf, N, L);
		//printf("A pattern matched completed %d / %d\n", i + 1, M);
		if (i == M / 10)
			printf("10%% achieved\n");
		else if (i == M * 2 / 10)
			printf("20%% achieved\n");
		else if (i == M * 3 / 10)
			printf("30%% achieved\n");
		else if (i == M * 4 / 10)
			printf("40%% achieved\n");
		else if (i == M * 5 / 10)
			printf("50%% achieved\n");
		else if (i == M * 6 / 10)
			printf("60%% achieved\n");
		else if (i == M * 7 / 10)
			printf("70%% achieved\n");
		else if (i == M * 8 / 10)
			printf("80%% achieved\n");
		else if (i == M * 9 / 10)
			printf("90%% achieved\n");
		else if (i == M - 1)
			printf("100%% achieved\n");
	}
	fclose(f);
	printf("Completed\n");

	printf("Reconstructing myGenome string...");
	MyGenomeSeq = (char*)malloc(sizeof(char)*(N + 1));
	for (int i = 0; i < N; i++)
		MyGenomeSeq[i] = reconstruct_arr[i].more_ch;
	MyGenomeSeq[N] = '\0';
	printf("Completed\n");

	//end counting time
	printf("\n*****Counting time ends here***** (kmp)\n");
	end_time = clock();
	time_taken = ((double)end_time - start_time) / CLOCKS_PER_SEC; //get execution time
	printf("Execution time = %.3lf\n\n", time_taken);

	/*free(reconstruct_arr);
	free(codes);*/

	printf("Make a file of myGenome string...");
	makeFileofString("reconstructed_myGenome_kmp.txt", MyGenomeSeq);
	printf("Completed\n");

	printf("Caculating accuracy of reconstructed string...\n");
	readFiletoString("input_50000.txt", ReferenceSeq);
	for (int i = 0; i < N; i++)
		if (MyGenomeSeq[i] != ReferenceSeq[i])
			mismatch_cnt++;
	printf("Accuracy with reference: %lf%%\n", (1 - ((double)mismatch_cnt / N)) * 100);
	mismatch_cnt = 0;
	readFiletoString("myGenome.txt", ReferenceSeq);
	for (int i = 0; i < N; i++)
		if (MyGenomeSeq[i] != ReferenceSeq[i])
			mismatch_cnt++;
	printf("Accuracy with myGenome: %lf%%\n", (1 - ((double)mismatch_cnt / N)) * 100);

	return 0;
}