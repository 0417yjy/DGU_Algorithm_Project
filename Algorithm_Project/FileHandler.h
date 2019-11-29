#pragma once
#include <stdio.h>

int makeFileofString(char *filename, char *str);
int readFiletoString(char *filename, char *str);

void generateSequence(char *ref, char *str_name);
void readSequence(char *ref, char *myseq);
void makeSnp(char *str);
void makeShortReads(char*);