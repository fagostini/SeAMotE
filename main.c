#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <pthread.h>

#include "main.h"
#include "my_library.o"
#include "RNA_lib.o"
#include "DNA_lib.o"

/* This functions need to be here because of the threading */
/* Probably it is possible to move then to another file but I do not know how to properly do that */
struct thread_data{
	int  thread_id;
	char **arrSeq;
	char ***arrRnd;
	char **arrMot;
	int depth;
	int arrNum;
	int motNum;
	int motLen;
	float *cov;
	int *count;
	int **motDistr;
	int ***pVal;
};
struct thread_data thread_data_array[NUM_THREADS*3];
void *set_coverage(void *threadarg){

	int taskid, seqNum, motNum, motLen;
	float *coverage;
	char **seqSet, **motSet;
	
	sleep(1);
	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;
	taskid = my_data->thread_id;
	seqNum = my_data->arrNum;
	motNum = my_data->motNum;
   	motLen = my_data->motLen;
   	coverage = my_data->cov;
   	seqSet = my_data->arrSeq;
   	motSet = my_data->arrMot;
   	
	int m, s, pos, i, br, count;
	for( m=0; m<motNum; m++ ){
		count = 0;
		for( s=0; s<seqNum; s++ ){
			pos = br = 0;
			do{
				i = 0;
				do{
					if( motSet[m][i] != '-' & motSet[m][i] != seqSet[s][pos+i]){
						break;
					}
					i++;
					if( i==motLen ){
						count++;
						br = 1;
						break;
					}
				}while( i<motLen );
				if( br == 1 ){
					break;
				}
				pos++;
			}while( pos<(strlen(seqSet[s])-motLen+1));
		}
		coverage[m] += count > 0 ? (float)count/seqNum : 0;
	}
	pthread_exit(NULL);
}
void *set_coverage_avg(void *threadarg){

	int taskid, seqNum, motNum, motLen, arrDepth;
	float *coverage;
	char ***seqSet, **motSet;
	
	sleep(1);
	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;
	taskid = my_data->thread_id;
	arrDepth = my_data->depth;
	seqNum = my_data->arrNum;
	motNum = my_data->motNum;
   	motLen = my_data->motLen;
   	coverage = my_data->cov;
   	seqSet = my_data->arrRnd;
   	motSet = my_data->arrMot;
   	
	int r, m, s, pos, i, br, count;
	for( m=0; m<motNum; m++ ){
		for( r=0; r<arrDepth; r++ ){
			count = 0;
			for( s=0; s<seqNum; s++ ){
				pos = br = 0;
				do{
					i = 0;
					do{
						if( motSet[m][i] != '-' & motSet[m][i] != seqSet[r][s][pos+i]){
							break;
						}
						i++;
						if( i==motLen ){
							count++;
							br = 1;
							break;
						}
					}while( i<motLen );
					if( br == 1 ){
						break;
					}
					pos++;
				}while( pos<(strlen(seqSet[r][s])-motLen+1));
			}
			coverage[m] += count > 0 ? (float)count/seqNum : 0;
		}
		coverage[m] /= arrDepth;
	}
	pthread_exit(NULL);
}
void *set_count(void *threadarg){

	int taskid, seqNum, motNum, motLen, *count;
	float *coverage;
	char **seqSet, **motSet;
	int **mDist;
	
	sleep(1);
	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;
	taskid = my_data->thread_id;
	seqNum = my_data->arrNum;
	motNum = my_data->motNum;
   	motLen = my_data->motLen;
   	count = my_data->count;
   	seqSet = my_data->arrSeq;
   	motSet = my_data->arrMot;
   	mDist = my_data->motDistr;
  	
	int m, s, pos, i, c, sc;
	for( m=0; m<motNum; m++ ){
		c = 0;
		for( s=0; s<seqNum; s++ ){
			pos = sc = 0;
			do{
				i = 0;
				do{
					if( motSet[m][i] != '-' & motSet[m][i] != seqSet[s][pos+i]){
						break;	/* Escape when a mismatch is found */
					}
					i++;
					if( i==motLen ){
						c++;
						sc++;	/* Put the condition sc += sc > 0 ? 0 : 1; if you want to consider only one motif per sequence   */
						break;	/* Escape when the motif is completely covered */
					}
				}while( i<motLen );
				pos++;
			}while( pos<(strlen(seqSet[s])-motLen+1));
			mDist[m][s] = sc > 0 ? sc : 0; /* Multiple motifs in the same sequence are considered */
		}
		count[m] = c > 0 ? c : 0;
	}
	pthread_exit(NULL);
}
void *set_count_rand(void *threadarg){

	int  *count;
	int taskid, seqNum, motNum, motLen, arrDepth, ***distSet;
	float *coverage;
	char ***seqSet, **motSet;

	sleep(1);
	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;
	taskid = my_data->thread_id;
	arrDepth = my_data->depth;
	seqNum = my_data->arrNum;
	motNum = my_data->motNum;
   	motLen = my_data->motLen;
	count = my_data->count;
   	seqSet = my_data->arrRnd;
   	motSet = my_data->arrMot;
   	distSet = my_data->pVal;

	int r, m, s, pos, i, c, ctmp;
	for( m=0; m<motNum; m++ ){
		ctmp = 0;
		for( r=0; r<arrDepth; r++ ){
			for( s=0; s<seqNum; s++ ){
				c = pos = 0;
				do{
					i = 0;
					do{
						if( motSet[m][i] != '-' & motSet[m][i] != seqSet[r][s][pos+i]){
							break;
						}
						i++;
						if( i==motLen ){
							c++;
							break;
						}
					}while( i<motLen );
					pos++;
				}while( pos<(strlen(seqSet[r][s])-motLen+1));
				distSet[m][r][s] = c > 0 ? c : 0;
				ctmp += c;
			}
		/*	count[m] += c > 0 ? c : 0; */
		}
		count[m] = (int)(ctmp+(arrDepth/2))/arrDepth;
	}
	pthread_exit(NULL);

}

int main(int argc, char* argv[]){
	
/* 	srand(time(NULL)); */
	srand(1);
	
/*	char *my_seq = NULL;
	random_sequence(100000, &my_seq);
	puts(my_seq);
	exit(0); */
	
	int i, numPos, numNeg, dna, rna, prot;
	numPos = numNeg = 0;
	int ms = MIN_MOT;
	int mn = pow(4,ms);
	float th = 0.95;
	char **motifs, **posID, **posSeq, **negID, **negSeq, **shuSeq;
	motifs = posID = posSeq = negID = negSeq = shuSeq = NULL;
	size_t size_mot;
	FILE *log, *Fpositive, *Fnegative, *Fshuffle, *Fposrand, *Fnegrand, *Fposshuf, *Fnegshuf;
	char *fileM = (char *)calloc(100, sizeof(char));
	char *fileP = (char *)calloc(100, sizeof(char));
	char *fileN = (char *)calloc(100, sizeof(char));
	char *type = (char *)calloc(10, sizeof(char));
/* 	log = open_file(log, "log.txt", "w"); */
	printf("Reading args... ");  fflush(stdout);
/* 	fprintf(log, "Reading args... "); */
	read_args(argc, argv, fileM, fileP, fileN, &th, type);
	printf("done\n"); fflush(stdout);
/* 	fprintf(log, "done\n"); */
	if( type[0] != '\0' ){
		if((strcmp(type, "rna") == 0) | (strcmp(type, "RNA") == 0)){
			rna = 1;
		}	
		else if((strcmp(type, "dna") == 0) | (strcmp(type, "DNA") == 0)){
			dna = 1;
		}
		else{
			printf("The molecule type specified cannot be recognized!\n"); fflush(stdout);
			exit(1);			
		}
	}
	if( fileM[0] != '\0' ){
		printf("Motifs file: %s\n", fileM); fflush(stdout);
/* 		fprintf(log, "Motifs file: %s\n", fileM); */
	}
	else{
		printf("Motifs file: Not specified\nGenerating the seed (%d nt) motifs... ", ms); fflush(stdout);
/* 		fprintf(log, "Motifs file: Not specified\nGenerating the seed (%d nt) motifs... ", ms); */
		motifs = calloc2Dchar(ms, mn);
		size_mot = (ms)*mn*sizeof(char)*4;
		create_motifs_nt(motifs, ms);
		printf("done\n"); fflush(stdout);
/* 		fprintf(log, "done\n"); */
	}
	if( fileP[0] != '\0' ){
		printf("Positive file: %s\nReading the Positive file... ", fileP); fflush(stdout);
/* 		fprintf(log, "Positive file: %s\nReading the Positive file... ", fileP); */
		Fpositive = open_file(Fpositive, fileP, "r");
		numPos = read_lines(Fpositive);
		posID = dna == 0 ? calloc2Dchar(MAX_ID, numPos) : calloc2Dchar(MAX_ID, numPos*2);
		posSeq = calloc2Dchar(MAX_SEQ, numPos);
		read_oneline(Fpositive, posID, posSeq, MAX_ID, MAX_SEQ, "%s %[^\n]%*c");
		fclose(Fpositive);
		if( dna == 1 ){
			compl_inverse(numPos, &posSeq);
		}
		printf("done\n"); fflush(stdout);
/* 		fprintf(log, "done\n"); */
	}
	else{
		printf("Positive file: Not specified\nThe positive dataset must be specified!\n"); fflush(stdout);
		exit(1);
	}
	if( fileN[0] != '\0' ){
		printf("Negative file: %s\nReading the Negative file... ", fileN); fflush(stdout);	
/* 		fprintf(log, "Negative file: %s\nReading the Negative file... ", fileN); */
		Fnegative = open_file(Fnegative, fileN, "r");
		numNeg = read_lines(Fnegative);
		negID = dna == 0 ? calloc2Dchar(MAX_ID, numNeg) : calloc2Dchar(MAX_ID, numNeg*2);
		negSeq = calloc2Dchar(MAX_SEQ, numNeg);
		read_oneline(Fnegative, negID, negSeq, MAX_ID, MAX_SEQ, "%s %[^\n]%*c");
		fclose(Fnegative);
		if( dna == 1 ){
			compl_inverse(numNeg, &negSeq);
		}
		printf("done\n"); fflush(stdout);
/* 		fprintf(log, "done\n"); */
	}
	else{
		printf("Negative file: Not specified\nGenerating a shuffled set from the Positive file... "); fflush(stdout);
/* 		fprintf(log, "Negative file: Not specified\nGenerating a shuffled set from the Positive file... "); */
		Fshuffle = open_file(Fshuffle, "Ref_shuffle.txt", "w");
		shuSeq = calloc2Dchar(MAX_SEQ, numPos*NUM_SHUFFLE);
		for( i=0; i<numPos; i++ ){
			char *sequence = (char *)calloc(MAX_SEQ, sizeof(char));
			strcpy(sequence, posSeq[i]);
			int c;
			for( c=0; c<NUM_SHUFFLE; c++ ){
				str_shuffle(sequence, strlen(sequence));
				strcpy(shuSeq[i*NUM_SHUFFLE+c], sequence);
				fprintf(Fshuffle, "%s\n", sequence);
			}
			free(sequence);
		}
		fclose(Fshuffle);
		numNeg = numPos*NUM_SHUFFLE;
		negSeq = shuSeq;
		printf("done\n"); fflush(stdout);
/* 		fprintf(log, "done\n"); */
	}
	printf("Coverage threshold: %.2f\n", th); fflush(stdout);
/* 	fprintf(log, "Coverage threshold: %.2f\n", th); */

/* 	GENERATION OF THE RANDOMIZED REFERENCE DATASETS */	
	int r;
	Fposrand = open_file(Fposrand, "Ref_PosRand.txt", "w");
	char ***PrndSeq = (char ***)calloc(NUM_RAND, sizeof(**PrndSeq));
	for( r=0; r<NUM_RAND; r++){
		PrndSeq[r] = (char **)calloc(numPos, sizeof(*PrndSeq[r]));
		for( i=0; i<numPos; i++){
			int len = strlen(posSeq[i]);
			random_sequence(len, &PrndSeq[r][i]);
			fprintf(Fposrand, "%s\n", PrndSeq[r][i]);
		}
	}
	fclose(Fposrand);
	Fnegrand = open_file(Fnegrand, "Ref_NegRand.txt", "w");
	char ***NrndSeq = (char ***)calloc(NUM_RAND, sizeof(**NrndSeq));
	for( r=0; r<NUM_RAND; r++){
		NrndSeq[r] = (char **)calloc(numNeg, sizeof(*NrndSeq[r]));
		for( i=0; i<numNeg; i++){
			int len = strlen(negSeq[i]);
			random_sequence(len, &NrndSeq[r][i]);
			fprintf(Fnegrand, "%s\n", NrndSeq[r][i]);
		}
	}
	fclose(Fnegrand);

/* 	GENERATION OF THE SHUFFLED REFERENCE DATASETS */	
	int s;
	Fposshuf = open_file(Fposshuf, "Ref_PosShuf.txt", "w");
	char ***PshuSeq = (char ***)calloc(NUM_SHUFFLE, sizeof(**PshuSeq));
	for( s=0; s<NUM_SHUFFLE; s++){
		PshuSeq[s] = (char **)calloc(numPos, sizeof(*PshuSeq[s]));
		for( i=0; i<numPos; i++){
			int len = strlen(posSeq[i]);
			PshuSeq[s][i] = (char *)calloc(len, sizeof(char));
			strncpy(PshuSeq[s][i], posSeq[i], len);
			str_shuffle_arr(len, &PshuSeq[s][i]);
			fprintf(Fposshuf, "%s\n", PshuSeq[s][i]);
		}
	}
	fclose(Fposshuf);
	Fnegshuf = open_file(Fnegshuf, "Ref_NegShuf.txt", "w");
	char ***NshuSeq = (char ***)calloc(NUM_SHUFFLE, sizeof(**NshuSeq));
	for( s=0; s<NUM_SHUFFLE; s++){
		NshuSeq[s] = (char **)calloc(numNeg, sizeof(*NshuSeq[s]));
		for( i=0; i<numNeg; i++){
			int len = strlen(negSeq[i]);
			NshuSeq[s][i] = (char *)calloc(len, sizeof(char));
			strncpy(NshuSeq[s][i], negSeq[i], len);
			str_shuffle_arr(len, &NshuSeq[s][i]);
			fprintf(Fnegshuf, "%s\n", NshuSeq[s][i]);
		}
	}
	fclose(Fnegshuf);

/* 	GENERATION OF THE BOOTSTRAPPED REFERENCE DATASETS */	
	int b;
	int numTot = numPos+numNeg;
	int **newOrder = (int **)calloc(NUM_BOOT, sizeof(*newOrder));
	for( b=0; b<NUM_BOOT; b++ ){
		newOrder[b] = (int *)calloc(numTot, sizeof(int));
		for( i=0; i<numTot; i++ ){
			newOrder[b][i] = i;
		}
	}
	bootstrap_sampling(numTot, newOrder);

	printf("Calculating the seed coverages... "); fflush(stdout);
/* 	fprintf(log, "Calculating the seed coverages... "); */
	size_t size_cov = mn*sizeof(float);
	float *posCov = (float *)calloc(mn, sizeof(posCov));
	float *negCov = (float *)calloc(mn, sizeof(negCov));
	float *PrndCov = (float *)calloc(mn, sizeof(PrndCov));
	float *NrndCov = (float *)calloc(mn, sizeof(NrndCov));
	float *PshuCov = (float *)calloc(mn, sizeof(PshuCov));
	float *NshuCov = (float *)calloc(mn, sizeof(NshuCov));

	/* ----- MULTI-THREADING COVERAGE CALCULATION ----- BEGINS ----- */	
	int rc;
	void *status;
	pthread_t threads[NUM_THREADS*3];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	thread_data_array[0].thread_id = 0;
	thread_data_array[0].arrSeq = posSeq;
	thread_data_array[0].arrMot = motifs;
	thread_data_array[0].depth = 1;
	thread_data_array[0].arrNum = numPos;
	thread_data_array[0].motNum = mn;
	thread_data_array[0].motLen = ms;
	thread_data_array[0].cov = posCov;
	thread_data_array[1].thread_id = 1;
	thread_data_array[1].depth = 1;
	thread_data_array[1].arrSeq = negSeq;
	thread_data_array[1].arrMot = motifs;
	thread_data_array[1].arrNum = numNeg;
	thread_data_array[1].motNum = mn;
	thread_data_array[1].motLen = ms;
	thread_data_array[1].cov = negCov;
	thread_data_array[2].thread_id = 2;
	thread_data_array[2].depth = NUM_RAND;
	thread_data_array[2].arrRnd = PrndSeq;
	thread_data_array[2].arrMot = motifs;
	thread_data_array[2].arrNum = numPos;
	thread_data_array[2].motNum = mn;
	thread_data_array[2].motLen = ms;
	thread_data_array[2].cov = PrndCov;
	thread_data_array[3].thread_id = 3;
	thread_data_array[3].depth = NUM_RAND;
	thread_data_array[3].arrRnd = NrndSeq;
	thread_data_array[3].arrMot = motifs;
	thread_data_array[3].arrNum = numNeg;
	thread_data_array[3].motNum = mn;
	thread_data_array[3].motLen = ms;
	thread_data_array[3].cov = NrndCov;
	thread_data_array[4].thread_id = 4;
	thread_data_array[4].depth = NUM_SHUFFLE;
	thread_data_array[4].arrRnd = PshuSeq;
	thread_data_array[4].arrMot = motifs;
	thread_data_array[4].arrNum = numPos;
	thread_data_array[4].motNum = mn;
	thread_data_array[4].motLen = ms;
	thread_data_array[4].cov = PshuCov;
	thread_data_array[5].thread_id = 5;
	thread_data_array[5].depth = NUM_SHUFFLE;
	thread_data_array[5].arrRnd = NshuSeq;
	thread_data_array[5].arrMot = motifs;
	thread_data_array[5].arrNum = numNeg;
	thread_data_array[5].motNum = mn;
	thread_data_array[5].motLen = ms;
	thread_data_array[5].cov = NshuCov;

	for( i=0; i<NUM_THREADS; i++ ){
		rc = pthread_create(&threads[i], &attr, set_coverage, (void *) &thread_data_array[i]);
		if (rc) {
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
		rc = pthread_create(&threads[i+2], &attr, set_coverage_avg, (void *) &thread_data_array[i+2]);
		if (rc) {
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}

	}
	pthread_attr_destroy(&attr);
	for( i=0; i<NUM_THREADS*2; i++ ){
		rc = pthread_join(threads[i], &status);
		if (rc) {
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
	}
	/* ----- MULTI-THREADING COVERAGE CALCULATION ----- ENDS ----- */

	char **backup_mot = (char **)calloc(mn, sizeof(*backup_mot));
	for( i=0; i<mn; i++ ){
		backup_mot[i] = (char *)calloc(ms, sizeof(char));
	}
	float *backup_pcov = (float *)calloc(mn, sizeof(float));
	float *backup_ncov = (float *)calloc(mn, sizeof(float));
	float *backup_prcov = (float *)calloc(mn, sizeof(float));
	float *backup_nrcov = (float *)calloc(mn, sizeof(float));
	
	memmove(backup_mot, motifs, size_mot);
	memmove(backup_pcov, thread_data_array[0].cov, size_cov);
	memmove(backup_ncov, thread_data_array[1].cov, size_cov);
	memmove(backup_prcov, thread_data_array[2].cov, size_cov);
	memmove(backup_nrcov, thread_data_array[3].cov, size_cov);
	
	printf("done\nBeginnig the motif coverage optimization:\n"); fflush(stdout);
/* 	fprintf(log, "done\nBeginnig the motif coverage optimization:\n"); */
	int old_mn, tmp_mn;
	char **old_motifs = motifs;
	char **new_motifs = NULL;
	int new_mn = tmp_mn = mn;
	int loop = 0;
 	do{
		old_mn = new_mn;
		new_mn = above_threshold(thread_data_array[0].cov, thread_data_array[1].cov, old_mn, th)*10;
		if( new_mn != 0 ){
			/* MOTIFS */
			size_mot = old_mn*(ms)*sizeof(char)*2;
			backup_mot = (char **)realloc(backup_mot, size_mot);
			memcpy(backup_mot, old_motifs, size_mot);
			ms++;
			new_mn = filter_and_expand_nt(old_motifs, thread_data_array[0].cov, thread_data_array[1].cov, old_mn, ms, th, &new_motifs);
		 	printf("   %d: Testing %d motifs (%d nt)... ", loop+1, new_mn, ms); fflush(stdout);
/* 			fprintf(log, "   %d: Testing %d motifs (%d nt)... ", loop+1, new_mn, ms); */
 			/* POSITIVE COVERAGE */
 			size_cov = old_mn*sizeof(float);
 			backup_pcov = (float *)realloc(backup_pcov, size_cov);
			memmove(backup_pcov, thread_data_array[0].cov, size_cov);
			posCov = realloc(posCov, new_mn*sizeof(float));
			bzero(posCov, new_mn*sizeof(float));
 			/* NEGATIVE COVERAGE */
 			backup_ncov = (float *)realloc(backup_ncov, size_cov);
			memmove(backup_ncov, thread_data_array[1].cov, size_cov);
			negCov = realloc(negCov, new_mn*sizeof(float));
			bzero(negCov, new_mn*sizeof(float));
 			/* POSITIVE RANDOM COVERAGE */
 			backup_prcov = (float *)realloc(backup_prcov, size_cov);
			memmove(backup_prcov, thread_data_array[2].cov, size_cov);
			PrndCov = realloc(PrndCov, new_mn*sizeof(float));
			bzero(PrndCov, new_mn*sizeof(float));
 			/* NEGATIVE RANDOM COVERAGE */
 			backup_nrcov = (float *)realloc(backup_nrcov, size_cov);
			memmove(backup_nrcov, thread_data_array[3].cov, size_cov);
			NrndCov = realloc(NrndCov, new_mn*sizeof(float));
			bzero(NrndCov, new_mn*sizeof(float));

	/* ----- MULTI-THREADING COVERAGE CALCULATION ----- BEGINS ----- */
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

			thread_data_array[0].thread_id = 0;
			thread_data_array[0].arrSeq = posSeq;
			thread_data_array[0].arrMot = new_motifs;
			thread_data_array[0].arrNum = numPos;
			thread_data_array[0].motNum = new_mn;
			thread_data_array[0].motLen = ms;
			thread_data_array[0].cov = posCov;
			thread_data_array[1].thread_id = 1;
			thread_data_array[1].arrSeq = negSeq;
			thread_data_array[1].arrMot = new_motifs;
			thread_data_array[1].arrNum = numNeg;
			thread_data_array[1].motNum = new_mn;
			thread_data_array[1].motLen = ms;
			thread_data_array[1].cov = negCov;
			thread_data_array[2].thread_id = 2;
			thread_data_array[2].arrRnd = PrndSeq;
			thread_data_array[2].arrMot = new_motifs;
			thread_data_array[2].arrNum = numPos;
			thread_data_array[2].motNum = new_mn;
			thread_data_array[2].motLen = ms;
			thread_data_array[2].cov = PrndCov;
			thread_data_array[3].thread_id = 3;
			thread_data_array[3].arrRnd = NrndSeq;
			thread_data_array[3].arrMot = new_motifs;
			thread_data_array[3].arrNum = numNeg;
			thread_data_array[3].motNum = new_mn;
			thread_data_array[3].motLen = ms;
			thread_data_array[3].cov = NrndCov;

			for( i=0; i<NUM_THREADS; i++ ){
				rc = pthread_create(&threads[i], &attr, set_coverage, (void *) &thread_data_array[i]);
				if (rc) {
					printf("ERROR; return code from pthread_create() is %d\n", rc);
					exit(-1);
				}
/*				rc = pthread_create(&threads[i+2], &attr, set_coverage_avg, (void *) &thread_data_array[i+2]);
				if (rc) {
					printf("ERROR; return code from pthread_create() is %d\n", rc);
					exit(-1);
				} */
			}
			pthread_attr_destroy(&attr);
			for( i=0; i<NUM_THREADS; i++ ){
				rc = pthread_join(threads[i], &status);
				if (rc) {
					printf("ERROR; return code from pthread_join() is %d\n", rc);
					exit(-1);
				}
			}
	/* ----- MULTI-THREADING COVERAGE CALCULATION ----- ENDS ----- */
			tmp_mn = old_mn;
	  		old_motifs = new_motifs;
	  		printf("done\n"); fflush(stdout);
/* 			fprintf(log, "done\n"); */
		}
		else{
			free2Dchar(new_motifs, new_mn);		  	
			break;
		}
		loop++;
 	}while( loop < MAX_MOT-MIN_MOT ); /* CHECK THE OLD_MN, TMP_MN AND NEW_MN */

	printf("Calculating the motif distribution in the sets... "); fflush(stdout);
/* 	fprintf(log, "Calculating the motif distribution in the sets.. "); */

	int *posCount, *negCount, *prandCount, *nrandCount, ***prand_matches, ***nrand_matches, *pshufCount, *nshufCount, ***pshuf_matches, ***nshuf_matches, **pmotDist, **nmotDist;
	float *posPval, *negPval;
	posCount = (int *)calloc(tmp_mn, sizeof(int));
	negCount = (int *)calloc(tmp_mn, sizeof(int));
	prandCount = (int *)calloc(tmp_mn, sizeof(int));
	nrandCount = (int *)calloc(tmp_mn, sizeof(int));
	pshufCount = (int *)calloc(tmp_mn, sizeof(int));
	nshufCount = (int *)calloc(tmp_mn, sizeof(int));
	posPval = (float *)calloc(tmp_mn, sizeof(float));
	negPval = (float *)calloc(tmp_mn, sizeof(float));
	
	int m; 
	prand_matches = (int ***)calloc(tmp_mn, sizeof(**prand_matches));
	for( m=0; m<tmp_mn; m++ ){
		prand_matches[m] = (int **)calloc(NUM_RAND, sizeof(*prand_matches[m]));
		for( i=0; i<NUM_RAND; i++ ){
			prand_matches[m][i] = (int *)calloc(numPos, sizeof(int));
		}
	}
	nrand_matches = (int ***)calloc(tmp_mn, sizeof(**nrand_matches));
	for( m=0; m<tmp_mn; m++ ){
		nrand_matches[m] = (int **)calloc(NUM_RAND, sizeof(*nrand_matches[m]));
		for( i=0; i<NUM_RAND; i++ ){
			nrand_matches[m][i] = (int *)calloc(numNeg, sizeof(int));
		}
	}
	pshuf_matches = (int ***)calloc(tmp_mn, sizeof(**pshuf_matches));
	for( m=0; m<tmp_mn; m++ ){
		pshuf_matches[m] = (int **)calloc(NUM_SHUFFLE, sizeof(*pshuf_matches[m]));
		for( i=0; i<NUM_SHUFFLE; i++ ){
			pshuf_matches[m][i] = (int *)calloc(numPos, sizeof(int));
		}
	}
	nshuf_matches = (int ***)calloc(tmp_mn, sizeof(**nshuf_matches));
	for( m=0; m<tmp_mn; m++ ){
		nshuf_matches[m] = (int **)calloc(NUM_SHUFFLE, sizeof(*nshuf_matches[m]));
		for( i=0; i<NUM_SHUFFLE; i++ ){
			nshuf_matches[m][i] = (int *)calloc(numNeg, sizeof(int));
		}
	}
	pmotDist = (int **)calloc(tmp_mn, sizeof(*pmotDist));
	for( b=0; b<tmp_mn; b++ ){
		pmotDist[b] = (int *)calloc(numPos, sizeof(int));
	}
	nmotDist = (int **)calloc(tmp_mn, sizeof(*nmotDist));
	for( b=0; b<tmp_mn; b++ ){
		nmotDist[b] = (int *)calloc(numNeg, sizeof(int));
	}

/* ----- MULTI-THREADING COUNT ----- BEGINS ----- */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	thread_data_array[0].thread_id = 0;
	thread_data_array[0].motNum = tmp_mn;
	thread_data_array[0].arrMot = backup_mot;
	thread_data_array[0].motLen = ms-1;
	thread_data_array[0].count = posCount;
	thread_data_array[0].motDistr = pmotDist;
	
	thread_data_array[1].thread_id = 1;
	thread_data_array[1].motNum = tmp_mn;
	thread_data_array[1].arrMot = backup_mot;
	thread_data_array[1].motLen = ms-1;
	thread_data_array[1].count = negCount;
	thread_data_array[1].motDistr = nmotDist;
	
	thread_data_array[2].thread_id = 2;
	thread_data_array[2].motNum = tmp_mn;
	thread_data_array[2].arrMot = backup_mot;
	thread_data_array[2].motLen = ms-1;
	thread_data_array[2].count = prandCount;
	thread_data_array[2].pVal = prand_matches;

	thread_data_array[3].thread_id = 3;
	thread_data_array[3].motNum = tmp_mn;
	thread_data_array[3].arrMot = backup_mot;
	thread_data_array[3].motLen = ms-1;
	thread_data_array[3].count = nrandCount;
	thread_data_array[3].pVal = nrand_matches;
	
	thread_data_array[4].thread_id = 4;
	thread_data_array[4].motNum = tmp_mn;
	thread_data_array[4].arrMot = backup_mot;
	thread_data_array[4].motLen = ms-1;
	thread_data_array[4].count = pshufCount;
	thread_data_array[4].pVal = pshuf_matches;

	thread_data_array[5].thread_id = 5;
	thread_data_array[5].motNum = tmp_mn;
	thread_data_array[5].arrMot = backup_mot;
	thread_data_array[5].motLen = ms-1;
	thread_data_array[5].count = nshufCount;
	thread_data_array[5].pVal = nshuf_matches;

	for( i=0; i<NUM_THREADS; i++ ){
		rc = pthread_create(&threads[i], &attr, set_count, (void *) &thread_data_array[i]);
		if (rc) {
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
 		rc = pthread_create(&threads[i+2], &attr, set_count_rand, (void *) &thread_data_array[i+2]);
		if (rc) {
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
		rc = pthread_create(&threads[i+4], &attr, set_count_rand, (void *) &thread_data_array[i+4]);
		if (rc) {
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}
	pthread_attr_destroy(&attr);
	for( i=0; i<NUM_THREADS*3; i++ ){
		rc = pthread_join(threads[i], &status);
		if (rc) {
			printf("ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
	}
	
	int p, cp, cn;
	float ppval, npval;
	for(i=0;i<tmp_mn;i++){
		ppval = npval = 0;
		for(b=0;b<NUM_BOOT;b++){
			cp = cn = 0;
			for(p=0;p<numPos;p++){
				cp += newOrder[b][p]<numPos ? pmotDist[i][newOrder[b][p]] : nmotDist[i][newOrder[b][p]-numPos];
			}
			ppval += posCount[i] > cp ? 1 : 0;
			for(p=numPos;p<numTot;p++){
				cn += newOrder[b][p]<numPos ? pmotDist[i][newOrder[b][p]] : nmotDist[i][newOrder[b][p]-numPos];
			}
			npval += negCount[i] > cn ? 1 : 0;
		}
		posPval[i] = (1+ppval)/(NUM_BOOT+1);
		negPval[i] = (1+npval)/(NUM_BOOT+1);
	}
	
/* ----- MULTI-THREADING COVERAGE CALCULATION ----- ENDS ----- */
	
	printf("done\nCalculating the optimal distance between the sets... "); fflush(stdout);	  	
/* 	fprintf(log, "done\nCalculating the optimal distance between the sets... "); */
			  	
	FILE *FoutAll = open_file(FoutAll, "motifs_last_loop.txt", "w" );

/*for( i=0; i<tmp_mn; i++){
	int pi;
	for( pi=0; pi<10; pi++){
		int ci;
		for( pi=0; pi<numPos; pi++){
			printf("%d ", prand_matches[i][pi][ci]);
		}
		printf("\n");
	}
}*/

	float RposP, RnegP, SposP, SnegP;
	int max_distance = 0;
	for( i=0; i<tmp_mn; i++ ){
		RposP = RnegP = SposP = SnegP = 0;
		int distance = (int)(abs((backup_pcov[i]-backup_ncov[i])*P_TEST));
		max_distance = distance > max_distance ? distance : max_distance;
		float *Rpos_vec = subsampling(prand_matches, numPos, NUM_RAND, P_TEST, i);
		float *Spos_vec = subsampling(pshuf_matches, numPos, NUM_SHUFFLE, P_TEST, i);
		for( p=0; p<P_TEST; p++){
			RposP += Rpos_vec[p] > backup_pcov[i] ? 1 : 0;
			SposP += Spos_vec[p] > backup_pcov[i] ? 1 : 0;
		}
		RposP /= P_TEST;
		SposP /= P_TEST;
		float *Rneg_vec = subsampling(nrand_matches, numNeg, NUM_RAND, P_TEST, i);		
		float *Sneg_vec = subsampling(nshuf_matches, numNeg, NUM_SHUFFLE, P_TEST, i);		
		for( p=0; p<P_TEST; p++){
			RnegP += Rneg_vec[p] > backup_ncov[i] ? 1 : 0;
			SnegP += Sneg_vec[p] > backup_ncov[i] ? 1 : 0;
		}
		RnegP /= P_TEST;
		SnegP /= P_TEST;
/*		fprintf(FoutAll, "%s %.2f %.2f %.2f %.2f %d %d %d %d %.3f %.3f\n", backup_mot[i], backup_pcov[i], backup_ncov[i], backup_prcov[i], backup_nrcov[i], posCount[i], negCount[i], prandCount[i], nrandCount[i], posP, negP); */
/* 		fprintf(FoutAll, "%s %.2f %.2f %d %d %d %d %d %d %.3f %.3f %.3f %.3f\n", backup_mot[i], backup_pcov[i], backup_ncov[i], posCount[i], negCount[i], prandCount[i], nrandCount[i], pshufCount[i], nshufCount[i], RposP, RnegP, SposP, SnegP); */
		fprintf(FoutAll, "%s %.2f %.2f %d %d %.5f %.5f\n", backup_mot[i], backup_pcov[i], backup_ncov[i], posCount[i], negCount[i], 1-posPval[i], 1-negPval[i]);
	}
	fclose(FoutAll);
	
/* FOR THE OUTPUT OF THE BEST MOTIFS I CAN RE-READ THE OUTALL AND DECIDE WHAT VARIABLE TO FILTER - SUCH AS THE P-VALUES OR THE GAP */	
/* 	printf("done\nWriting the final report... "); fflush(stdout);
	FILE *FoutBest = open_file(FoutBest, "best_motifs.txt", "w" );
	for( i=0; i<tmp_mn; i++ ){
		int distance = (int)(abs((backup_pcov[i]-backup_ncov[i])*100));
		if( distance >= max_distance-10 ){
			fprintf(FoutBest, "%s %.2f %.2f %.2f %.2f %d %d %d %d\n", backup_mot[i], backup_pcov[i], backup_ncov[i], backup_prcov[i], backup_nrcov[i], posCount[i], negCount[i], prandCount[i], nrandCount[i] );
		}
	}
	fclose(FoutBest);
*/	
	printf("done\nFreeing the memory... "); fflush(stdout);
/*  	fprintf(log, "done\nFreeing the memory... "); */

	free2Dchar(posID, numPos);
	free2Dchar(posSeq, numPos);
	if( negID != NULL ){
		free2Dchar(negID, numNeg);
		free2Dchar(negSeq, numNeg);
	}
	else{
		free2Dchar(negSeq, numPos);
	}
	free3Dchar(PrndSeq, numPos, NUM_RAND);
	free3Dchar(NrndSeq, numNeg, NUM_RAND);
	free2Dchar(backup_mot, tmp_mn);
	free(backup_pcov);
	free(backup_ncov);
	
	free2Dint(pmotDist, tmp_mn);
	free2Dint(nmotDist, tmp_mn);

	printf("done\nThe script executed successfully!\n"); fflush(stdout);
/* 	fprintf(log, "done\nThe script executed successfully!\n"); */
	
	pthread_exit(NULL);
		
}