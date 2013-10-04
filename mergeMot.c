#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>

#include <pthread.h>
#define NUM_THREADS 2

#define MAX_ID 100
#define MAX_SEQ 20000
#define MAX_MOT 21
#define NUM_BOOT 200

#include "my_library.lib"

struct thread_data{
	int  thread_id;
	char **arrSeq;
	char **arrMot;
	int depth;
	int arrNum;
	int motNum;
	int motLen;
	int **motDistr;
};
struct thread_data thread_data_array[NUM_THREADS];

void *set_coverage(void *threadarg){

	int taskid, seqNum, motNum, motLen;
	int **distribution;
	char **seqSet, **motSet;
	
	sleep(1);
	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;
	taskid = my_data->thread_id;
	seqNum = my_data->arrNum;
	motNum = my_data->motNum;
   	motLen = my_data->motLen;
   	distribution = my_data->motDistr;
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
						distribution[m][s] = 1;
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
	}
	pthread_exit(NULL);
}

int main(int argc, char* argv[]){

	int i;
	FILE *log, *Fpositive, *Fnegative, *Fmotifs;
	
	Fpositive = open_file(Fpositive, argv[1], "r");
	int numPos = read_lines(Fpositive);
	char **posID = (char **)malloc(numPos*sizeof(*posID));
	for(i=0;i<numPos;i++){
		posID[i] = (char *)malloc(MAX_ID*sizeof(char));
	}
	char **posSeq = (char **)malloc(numPos*sizeof(*posSeq));
	for(i=0;i<numPos;i++){
		posSeq[i] = (char *)malloc(MAX_SEQ*sizeof(char));
	}
	read_oneline(Fpositive, posID, posSeq, MAX_ID, MAX_SEQ, "%s %[^\n]%*c");
	fclose(Fpositive);

	Fnegative = open_file(Fpositive, argv[2], "r");
	int numNeg = read_lines(Fnegative);
	char **negID = (char **)malloc(numNeg*sizeof(*negID));
	for(i=0;i<numNeg;i++){
		negID[i] = (char *)malloc(MAX_ID*sizeof(char));
	}
	char **negSeq = (char **)malloc(numNeg*sizeof(*negSeq));
	for(i=0;i<numNeg;i++){
		negSeq[i] = (char *)malloc(MAX_SEQ*sizeof(char));
	}
	read_oneline(Fnegative, negID, negSeq, MAX_ID, MAX_SEQ, "%s %[^\n]%*c");
	fclose(Fnegative);

	Fmotifs = open_file(Fmotifs, argv[3], "r");
	int numMot = read_lines(Fmotifs);
	char **mot = (char **)malloc(numMot*sizeof(*mot));
	for(i=0;i<numMot;i++){
		mot[i] = (char *)malloc(MAX_MOT*sizeof(char));
	}
	float pvalues[numMot][2];
	int pos = 0;
	int lenMot = 0;
	while(!feof(Fmotifs)){
		int tr1, tr2;
		float tr3, tr4;
		fscanf(Fmotifs, "%s %f %f %d %d %f %f\n", mot[pos], &tr3, &tr4, &tr1, &tr1, &pvalues[pos][0], &pvalues[pos][1] );
/* 
		printf("%s %.2f %.2f %d\n", mot[pos], pvalues[pos][0], pvalues[pos][1], (int)strlen(mot[pos]));
 */
		lenMot = strlen(mot[pos]);
		pos++;
	}
	fclose(Fmotifs);
	
	int **posDist = (int **)calloc(numMot, sizeof(*posDist));
	for(i=0;i<numMot;i++){
		posDist[i] = (int *)calloc(numPos, sizeof(int));
	}
	int **negDist = (int **)calloc(numMot, sizeof(*posDist));
	for(i=0;i<numMot;i++){
		negDist[i] = (int *)calloc(numNeg, sizeof(int));
	}
	
	/* ----- MULTI-THREADING COVERAGE CALCULATION ----- BEGINS ----- */	
	int rc;
	void *status;
	pthread_t threads[NUM_THREADS];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	thread_data_array[0].thread_id = 0;
	thread_data_array[0].arrSeq = posSeq;
	thread_data_array[0].arrMot = mot;
	thread_data_array[0].depth = 1;
	thread_data_array[0].arrNum = numPos;
	thread_data_array[0].motNum = numMot;
	thread_data_array[0].motLen = lenMot;
	thread_data_array[0].motDistr = posDist;
	thread_data_array[1].thread_id = 1;
	thread_data_array[1].arrSeq = negSeq;
	thread_data_array[1].arrMot = mot;
	thread_data_array[1].depth = 1;
	thread_data_array[1].arrNum = numNeg;
	thread_data_array[1].motNum = numMot;
	thread_data_array[1].motLen = lenMot;
	thread_data_array[1].motDistr = negDist;

	for( i=0; i<NUM_THREADS; i++ ){
		rc = pthread_create(&threads[i], &attr, set_coverage, (void *) &thread_data_array[i]);
		if (rc) {
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
		rc = pthread_create(&threads[i], &attr, set_coverage, (void *) &thread_data_array[i]);
		if (rc) {
			printf("ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
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
	
	int m1, m2;
	for(m1=0;m1<numMot;m1++){
		for(m2=0;m2<numMot;m2++){
			if( m1!=m2 & pvalues[m1][0]<0.05 & pvalues[m2][0]<0.05 ){
				int tot = 0;
				printf("%s", mot[m1]);
				for(i=0;i<numPos;i++){
					if( posDist[m1][i] == 1 | posDist[m2][i] == 1 ){
						printf(" %d", 1);
						tot++;
					}
					else{
						printf(" %d", 0);
					}
				}
				printf(" %s %.2f\n", mot[m2], (float)tot/numPos);
			}
		}
	}

	return 0;

}