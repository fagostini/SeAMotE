#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include <pthread.h>

#include "motifSearch.h"
#include "my_library.h"
#include "RNA_lib.h"
#include "DNA_lib.h"

#define SEED_NOTA 14
#define NOTATION 14
#define print_and_continue(a) {printf("%s ", a); fflush(stdout);}
#define print_and_exit(a) {printf("%d\n", a); fflush(stdout); exit(0);}

int cmp( const void* a, const void* b){
	double* p1 = *(double**) a; 
	double* p2 = *(double**) b; 
	/*printf("Compare %f %x <=> %f %x\n", p1[0], a, p2[0], b);*/
	if( p1[1] > p2[1]){
		return -1;
	}
	else if (p1[1] < p2[1]){
		return 1;
	}
	else{
		return 0;
	}
}

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
	double *cov;
	int *count;
	int **motDistr;
	int ***pVal;
};

struct thread_data thread_data_array[NUM_THREADS];

void *set_coverage(void *threadarg){

	int taskid, seqNum, motNum, motLen, **motDistr;
	double *coverage;
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
   	motDistr = my_data->motDistr;
   	
	int m, s, pos, i, br, count;
	for( m=0; m<motNum; m++ ){
		count = 0;
		for( s=0; s<seqNum; s++ ){
			pos = br = 0;
			do{
				i = 0;
				do{
/* 
					if( motSet[m][i] != '-' && motSet[m][i] != seqSet[s][pos+i]){
						break;
					}
 */
 					if( seqSet[s][pos+i] == 'A' && ( motSet[m][i] != 'A' && motSet[m][i] != 'R' && motSet[m][i] != 'W' && motSet[m][i] != 'M' && motSet[m][i] != 'D' && motSet[m][i] != 'H' && motSet[m][i] != 'V') ){
 						break;
 					}
 					else if( seqSet[s][pos+i] == 'C' && ( motSet[m][i] != 'C' && motSet[m][i] != 'Y' && motSet[m][i] != 'S' && motSet[m][i] != 'M' && motSet[m][i] != 'B' && motSet[m][i] != 'H' && motSet[m][i] != 'V') ){
 						break;
 					}
 					else if( seqSet[s][pos+i] == 'G' && ( motSet[m][i] != 'G' && motSet[m][i] != 'R' && motSet[m][i] != 'S' && motSet[m][i] != 'K' && motSet[m][i] != 'B' && motSet[m][i] != 'D' && motSet[m][i] != 'V') ){
 						break;
 					}
 					else if( seqSet[s][pos+i] == 'T' && ( motSet[m][i] != 'T' && motSet[m][i] != 'Y' && motSet[m][i] != 'W' && motSet[m][i] != 'K' && motSet[m][i] != 'B' && motSet[m][i] != 'D' && motSet[m][i] != 'H') ){
 						break;
 					}
					i++;
					if( i == motLen ){
						count++;
						br = 1;
						break;
					}
				}while( i < motLen );
				if( br == 1 ){
					motDistr[m][s] = 1;
					break;
				}
				pos++;
			}while( pos<(strlen(seqSet[s])-motLen+1));
		}
		coverage[m] += count > 0 ? (double)count/seqNum : 0;
	}
	
	pthread_exit(NULL);
}

void *set_count(void *threadarg){

	int taskid, seqNum, motNum, motLen, *count;
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
/* 
					if( motSet[m][i] != '-' && motSet[m][i] != seqSet[s][pos+i]){
						break;
					}
*/
 					if( seqSet[s][pos+i] == 'A' && ( motSet[m][i] != 'A' && motSet[m][i] != 'R' && motSet[m][i] != 'W' && motSet[m][i] != 'M' && motSet[m][i] != 'D' && motSet[m][i] != 'H' && motSet[m][i] != 'V') ){
 						break;
 					}
 					else if( seqSet[s][pos+i] == 'C' && ( motSet[m][i] != 'C' && motSet[m][i] != 'Y' && motSet[m][i] != 'S' && motSet[m][i] != 'M' && motSet[m][i] != 'B' && motSet[m][i] != 'H' && motSet[m][i] != 'V') ){
 						break;
 					}
 					else if( seqSet[s][pos+i] == 'G' && ( motSet[m][i] != 'G' && motSet[m][i] != 'R' && motSet[m][i] != 'S' && motSet[m][i] != 'K' && motSet[m][i] != 'B' && motSet[m][i] != 'D' && motSet[m][i] != 'V') ){
 						break;
 					}
 					else if( seqSet[s][pos+i] == 'T' && ( motSet[m][i] != 'T' && motSet[m][i] != 'Y' && motSet[m][i] != 'W' && motSet[m][i] != 'K' && motSet[m][i] != 'B' && motSet[m][i] != 'D' && motSet[m][i] != 'H') ){
 						break;
 					}
					i++;
					if( i == motLen-1 ){
						c++;
						sc++;	/* Put the condition sc += sc > 0 ? 0 : 1; if you want to consider only one motif per sequence   */
						break;	/* Escape when the motif is completely covered */
					}
				}while( i < motLen );
				pos++;
			}while( pos < (strlen(seqSet[s])-motLen+1) );
			mDist[m][s] = sc > 0 ? sc : 0; /* Multiple motifs in the same sequence are considered */
		}
		count[m] = c > 0 ? c : 0;
	}
	pthread_exit(NULL);
}

void select_motifs(int *numMot, int motLen, char ***motSet, double **posCov, double **negCov, double threshold){
	int perc = (int)((*numMot)*threshold);
	int m;
	double **index = calloc((*numMot), sizeof(double *));
	for( m=0; m<(*numMot); m++ ){
		index[m] = calloc(4, sizeof(double));
	}
	/* The first selection is on the coverage over the positive dataset */
	for( m=0; m<(*numMot); m++ ){
		index[m][0] = (*posCov)[m];
		index[m][1] = (*negCov)[m];
		index[m][2] = m;
	}
	qsort(index, (*numMot), sizeof(double), cmp);
	int *tmpIndex = calloc(perc*3,sizeof(int));
	for( m=0; m<perc; m++ ){
		tmpIndex[m] = (int)index[m][2];
	}
	/* The second selection is on the coverage over the negative dataset */
	for( m=0; m<(*numMot); m++ ){
		index[m][0] = (*negCov)[m];
		index[m][1] = (*posCov)[m];
		index[m][2] = m;
	}
	qsort(index, (*numMot), sizeof(double), cmp);
	for( m=perc; m<perc*2; m++ ){
		tmpIndex[m] = (int)index[m][2];
/* 		printf("%s %.3lf\n", (*motSet)[tmpIndex[m]], (*posCov)[tmpIndex[m]]); */
	}
	/* The third selection is on the maximization of the coverage between dataset */
	for( m=0; m<(*numMot); m++ ){
		index[m][0] = (*posCov)[m]-(*negCov)[m];
		index[m][1] = (*posCov)[m];
		index[m][2] = (*negCov)[m];
		index[m][3] = m;
	}
	qsort(index, (*numMot), sizeof(double), cmp);
	for( m=perc*2; m<perc*3; m++ ){
		tmpIndex[m] = (int)index[m][3];
	}
	
/* 	printf("%p %p %p\n", (*motSet), (*posCov), (*negCov)); fflush(stdout); */
	int realNum = 0;
	int *selIndex = calloc(perc*3,sizeof(int));
	for( m=0; m<perc*3; m++ ){
		int i;
		int exist = 0;
		for( i=0; i<m; i++ ){
			if( tmpIndex[m] == tmpIndex[i] ){
				exist++;
				break;
			}
		}
		if( exist == 0 ){
			selIndex[realNum] = tmpIndex[m];
			realNum++;	
		}
	}
	
	double *PselCov = calloc(realNum, sizeof(double));
	double *NselCov = calloc(realNum, sizeof(double));
	char **selMot = malloc(realNum*sizeof(char *));
	for( m=0; m<realNum; m++ ){
		selMot[m] = malloc((motLen+1)*sizeof(char));
		memset(selMot[m], '\0', (motLen+1)*sizeof(char));
		memmove(selMot[m], (*motSet)[selIndex[m]], motLen+1 );
		PselCov[m] = (*posCov)[selIndex[m]];
		NselCov[m] = (*negCov)[selIndex[m]];
/* 		printf("%s %.3lf %.3lf\n", selMot[m], PselCov[m], NselCov[m]); fflush(stdout); */
	}
		
	for( m=0; m<(*numMot); m++ ){
		free(index[m]);
		free((*motSet)[m]);
	}
	free(index);
	free((*motSet));
	free((*posCov));
	free((*negCov));
	free(tmpIndex);
	free(selIndex);

	(*motSet) = selMot;
	(*posCov) = PselCov;
	(*negCov) = NselCov;
	
	(*numMot) = realNum;

}

void call_R(int n, double th, int o){
	
	FILE *fp;
	int status;
/* 	char path[1035]; */

	/* Open the command for reading. */
	char command[150];
/* 	sprintf(command, "cat pvalues.R | R --slave --vanilla --args %s %.2f", filename, th); */
	sprintf(command, "cat significant.R | R --slave --vanilla --args %d %.2f %d", n, th, o);
	fp = popen(command, "r");
	if (fp == NULL) {
		printf("Failed to run command\n" );
		exit(1);
	}

	/* Read the output a line at a time - output it. */
/* 
	while( fgets(path, sizeof(path)-1, fp) != NULL) {
		printf("%s", path);
	}
 */
 	int c = 0;
 	char *line = malloc(200*sizeof(char));
 	while( !feof(fp) ){
 		memset(line, '\0', 200*sizeof(char));
 		fscanf(fp, "%[^\n]%*c", line);
 		if( *line != '\0' && c != 0 ){
 			char mot[21];
 			int tp, fp, fn, tn;
 			double pval;
 			sscanf(line, "%s %d %d %d %d %lf", mot, &tp, &fp, &fn, &tn, &pval);
/*  			printf("%s %lf\n", mot, pval); */
 		}
 		c++;
 	}
 	free(line);

	/* close */
	pclose(fp);

}

int main(int argc, char* argv[]){
	
/* 	srand(time(NULL)); */
	srand(1);	
	
	int i, numPos, numNeg, dna, rna, prot, ref;
	numPos = numNeg = ref = 0;
	int ms = MIN_MOT;
	int mn = pow(SEED_NOTA,ms);
	double th = 0.8;
	char **motifs, **posID, **posSeq, **negID, **negSeq, **shuSeq;
	motifs = posID = posSeq = negID = negSeq = shuSeq = NULL;
	FILE *log, *Fpositive, *Fnegative, *Fshuffle, *Fposshuf, *Fnegshuf;
	char *fileM = malloc(100*sizeof(char));
	memset(fileM, '\0', 100*sizeof(char));
	char *fileP = malloc(100*sizeof(char));
	memset(fileP, '\0', 100*sizeof(char));
	char *fileN = malloc(100*sizeof(char));
	memset(fileN, '\0', 100*sizeof(char));
	char *type = malloc(100*sizeof(char));
	memset(type, '\0', 100*sizeof(char));
/* 	log = open_file(log, "log.txt", "w"); */
/* 	printf("Reading args... ");  fflush(stdout); */
/* 	fprintf(log, "Reading args... "); */
	read_args(argc, argv, fileM, fileP, fileN, &th, type, &ref);
/* 	printf("done\n"); fflush(stdout); */
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
		motifs = malloc2Dchar(ms, mn);
		create_motifs_nt(motifs, ms);
		printf("done\n"); fflush(stdout);
/* 		fprintf(log, "done\n"); */
	}
	if( fileP[0] != '\0' ){
		printf("Positive file: %s\nReading the Positive file... ", fileP); fflush(stdout);
/* 		fprintf(log, "Positive file: %s\nReading the Positive file... ", fileP); */
		Fpositive = open_file(Fpositive, fileP, "r");
		numPos = read_lines(Fpositive);
		posID = dna == 0 ? malloc2Dchar(MAX_ID, numPos) : malloc2Dchar(MAX_ID, numPos*2);
		posSeq = malloc2Dchar(MAX_SEQ, numPos);
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
		negID = dna == 0 ? malloc2Dchar(MAX_ID, numNeg) : malloc2Dchar(MAX_ID, numNeg*2);
		negSeq = malloc2Dchar(MAX_SEQ, numNeg);
		read_oneline(Fnegative, negID, negSeq, MAX_ID, MAX_SEQ, "%s %[^\n]%*c");
		fclose(Fnegative);
		if( dna == 1 ){
			compl_inverse(numNeg, &negSeq);
		}
		printf("done\n"); fflush(stdout);
/* 		fprintf(log, "done\n"); */
	}
	else{
		if( ref == 1 ){
			printf("Negative file: Not specified\nReference: Random\nGenerating reference set... "); fflush(stdout);
		}
		else{
			printf("Negative file: Not specified\nReference: Shuffle\nGenerating reference set from the Positive file... "); fflush(stdout);			
		}
/* 		fprintf(log, "Negative file: Not specified\nGenerating a shuffled set from the Positive file... "); */
		Fshuffle = open_file(Fshuffle, "Ref_shuffle.txt", "w");
		shuSeq = malloc2Dchar(MAX_SEQ, numPos*NO_NEGA);
		char *sequence = malloc(MAX_SEQ*sizeof(char));
		for( i=0; i<numPos; i++ ){
			memset(sequence, '\0', MAX_SEQ*sizeof(char) );
			memmove(sequence, posSeq[i], strlen(posSeq[i]));
			int c;
			for( c=0; c<NO_NEGA; c++ ){
				if( ref == 1 ){
					str_random(sequence, strlen(sequence));
				}
				else{
					str_shuffle(sequence, strlen(sequence));
				}
				memmove(shuSeq[i*NO_NEGA+c], sequence, strlen(sequence));
				fprintf(Fshuffle, "%s\n", sequence);
			}
		}
		fclose(Fshuffle);
		free(sequence);
		numNeg = numPos*NO_NEGA;
		negSeq = shuSeq;
		printf("done\n"); fflush(stdout);
/* 		fprintf(log, "done\n"); */
	}
	printf("Coverage threshold: %.2f\n", th); fflush(stdout);
/* 	fprintf(log, "Coverage threshold: %.2f\n", th); */
/* 	GENERATION OF THE SHUFFLED REFERENCE DATASETS (NOT IN USE AT THE MOMENT) */	
/* 
	int s;
	Fposshuf = open_file(Fposshuf, "Ref_PosShuf.txt", "w");
	char ***PshuSeq = malloc(NUM_SHUFFLE*sizeof(**PshuSeq));
	for( s=0; s<NUM_SHUFFLE; s++){
		PshuSeq[s] = malloc(numPos*sizeof(*PshuSeq[s]));
		for( i=0; i<numPos; i++){
			int len = strlen(posSeq[i]);
			PshuSeq[s][i] = malloc((len+1)*sizeof(char));
			memset(PshuSeq[s][i], '\0', (len+1)*sizeof(char));
			memmove(PshuSeq[s][i], posSeq[i], len);
			str_shuffle_arr(len, &PshuSeq[s][i]);
			fprintf(Fposshuf, "%s\n", PshuSeq[s][i]);
		}
	}
	fclose(Fposshuf);
	Fnegshuf = open_file(Fnegshuf, "Ref_NegShuf.txt", "w");
	char ***NshuSeq = malloc(NUM_SHUFFLE*sizeof(**NshuSeq));
	for( s=0; s<NUM_SHUFFLE; s++){
		NshuSeq[s] = malloc(numNeg*sizeof(*NshuSeq[s]));
		for( i=0; i<numNeg; i++){
			int len = strlen(negSeq[i]);
			NshuSeq[s][i] = malloc((len+1)*sizeof(char));
			memset(NshuSeq[s][i], '\0', (len+1)*sizeof(char));
			memmove(NshuSeq[s][i], negSeq[i], len);
			str_shuffle_arr(len, &NshuSeq[s][i]);
			fprintf(Fnegshuf, "%s\n", NshuSeq[s][i]);
		}
	}
	fclose(Fnegshuf);
 */

/* 	GENERATION OF THE BOOTSTRAPPED REFERENCE DATASETS */	
	int b;
	int numTot = numPos+numNeg;
	int **newOrder = calloc(NUM_BOOT, sizeof(*newOrder));
	for( b=0; b<NUM_BOOT; b++ ){
		newOrder[b] = calloc(numTot, sizeof(int));
		for( i=0; i<numTot; i++ ){
			newOrder[b][i] = i;
		}
	}
	bootstrap_sampling(numTot, newOrder);

	printf("Calculating the coverage for %d seeds... ", mn); fflush(stdout);
/* 	fprintf(log, "Calculating the seed coverages... "); */
	double *posCov = calloc(mn, sizeof(posCov));
	double *negCov = calloc(mn, sizeof(negCov));
	int **pmotDist = calloc(mn, sizeof(int *));
	for( i=0; i<mn; i++ ){
		pmotDist[i] = calloc(numPos,sizeof(int));
	}
	int **nmotDist = calloc(mn, sizeof(int *));
	for( i=0; i<mn; i++ ){
		nmotDist[i] = calloc(numNeg, sizeof(int));
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
	thread_data_array[0].arrMot = motifs;
	thread_data_array[0].depth = 1;
	thread_data_array[0].arrNum = numPos;
	thread_data_array[0].motNum = mn;
	thread_data_array[0].motLen = ms;
	thread_data_array[0].cov = posCov;
	thread_data_array[0].motDistr = pmotDist;
	thread_data_array[1].thread_id = 1;
	thread_data_array[1].depth = 1;
	thread_data_array[1].arrSeq = negSeq;
	thread_data_array[1].arrMot = motifs;
	thread_data_array[1].arrNum = numNeg;
	thread_data_array[1].motNum = mn;
	thread_data_array[1].motLen = ms;
	thread_data_array[1].cov = negCov;
	thread_data_array[1].motDistr = nmotDist;

	for( i=0; i<NUM_THREADS; i++ ){
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
/* 
	select_motifs(&mn, ms, &motifs, &posCov, &negCov, 0.1);
	thread_data_array[0].arrMot = motifs;
	thread_data_array[1].arrMot = motifs;
	thread_data_array[0].cov = posCov;
	thread_data_array[1].cov = negCov;
	thread_data_array[0].motNum = mn;
	thread_data_array[1].motNum = mn;
 */

	char *fname = malloc(50*sizeof(char));
	memset(fname, '\0', 50*sizeof(char));
	sprintf(fname,"tmp/motifs_%dnt.dat", ms);
	char *fname2 = malloc(50*sizeof(char));
	memset(fname2, '\0', 50*sizeof(char));
	sprintf(fname2,"tmp/motifs_%dnt_distribution.dat", ms);
	FILE *myTMP = open_file(myTMP, fname, "w" );
	for( i=0; i<mn; i++ ){
		fprintf(myTMP, "%s %d %d %d %d\n", motifs[i], (int)(numPos*posCov[i]), numPos, (int)(numNeg*negCov[i]), numNeg);
	}
	fclose(myTMP);
	FILE *myTMP2 = open_file(myTMP2, fname2, "w" );
	int ii;
	for( i=0; i<mn; i++ ){
		fprintf(myTMP2, "%s ", motifs[i]);
		for( ii=0; ii<numPos; ii++ ){
			fprintf(myTMP2, " %d", pmotDist[i][ii]);
		}
		fprintf(myTMP2, " \n");
	}
	fclose(myTMP2);
	
	call_R(ms, th, 0);

	for( i=0; i<mn ; i++ ){
		free(motifs[i]);
		free(pmotDist[i]);
		free(nmotDist[i]);
	}
	free(motifs);
	free(posCov);
	free(negCov);
	free(pmotDist);
	free(nmotDist);


	myTMP = open_file(myTMP, fname, "r" );
	mn = read_lines(myTMP);
	motifs = malloc(mn*sizeof(char *));
	for( i=0; i<mn ; i++ ){
		motifs[i] = malloc((ms+1)*sizeof(char));
		memset(motifs[i], '\0', (ms+1)*sizeof(char));
	}
	posCov = calloc(mn, sizeof(double));
	negCov = calloc(mn, sizeof(double));
	
	i = 0;
 	char *mline = malloc(200*sizeof(char));
 	char *mot = malloc(50*sizeof(char));
 	while( !feof(myTMP) ){
 		memset(mline, '\0', 200*sizeof(char));
 		fscanf(myTMP, "%[^\n]%*c", mline);
 		if( *mline != '\0' ){
 			memset(mot, '\0', 50*sizeof(char));
 			int tp, fp, fn, tn;
 			double pval, nval;
 			sscanf(mline, "%s %d %d %d %d %lf", mot, &tp, &fp, &fn, &tn, &pval, &nval);
 			memmove(motifs[i], mot, ms+1);
 			posCov[i] = (double)tp/numPos;
 			negCov[i] = (double)fn/numNeg;
 			i++;
 		}
 	}
 	free(mot);
 	free(mline);
	fclose(myTMP);
	free(fname);
	free(fname2);

	printf("done\nBeginnig the motif coverage optimization:\n"); fflush(stdout);
/* 	fprintf(log, "done\nBeginnig the motif coverage optimization:\n"); */

/* NEW SECTION STARTS HERE... */
	int loop = 0;
	char **old_motifs = motifs;
	char **new_motifs = NULL;
	int old_mn, tmp_mn;
	int new_mn = tmp_mn = mn;
	do{
		old_mn = new_mn;
		if( new_mn != 0 ){
			ms++;
			printf("%d", new_mn);
			new_mn = filter_and_expand_nt(old_motifs, thread_data_array[0].cov, thread_data_array[1].cov, old_mn, ms, th, &new_motifs);
		 	printf("   %d: Testing %d motifs (%d nt)", loop+1, new_mn, ms); fflush(stdout);

			posCov = realloc(posCov, new_mn*sizeof(double));
			memset(posCov, 0, new_mn*sizeof(double));
			negCov = realloc(negCov, new_mn*sizeof(double));
			memset(negCov, 0, new_mn*sizeof(double));
			
			printf("."); fflush(stdout);
			pmotDist = calloc(new_mn, sizeof(int *));
			for( i=0; i<new_mn; i++ ){
				pmotDist[i] = calloc(numPos, sizeof(int));
			}
			nmotDist = calloc(new_mn, sizeof(int *));
			for( i=0; i<new_mn; i++ ){
				nmotDist[i] = calloc(numNeg, sizeof(int));
			}
			printf("."); fflush(stdout);
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
			thread_data_array[0].motDistr = pmotDist;
			thread_data_array[1].thread_id = 1;
			thread_data_array[1].arrSeq = negSeq;
			thread_data_array[1].arrMot = new_motifs;
			thread_data_array[1].arrNum = numNeg;
			thread_data_array[1].motNum = new_mn;
			thread_data_array[1].motLen = ms;
			thread_data_array[1].cov = negCov;
			thread_data_array[1].motDistr = pmotDist;
			printf("."); fflush(stdout);
			for( i=0; i<NUM_THREADS; i++ ){
				rc = pthread_create(&threads[i], &attr, set_coverage, (void *) &thread_data_array[i]);
				if (rc) {
					printf("ERROR; return code from pthread_create() is %d\n", rc);
					exit(-1);
				}
			}
			pthread_attr_destroy(&attr);
			printf(". "); fflush(stdout);
			for( i=0; i<NUM_THREADS; i++ ){
				rc = pthread_join(threads[i], &status);
				if (rc) {
					printf("ERROR; return code from pthread_join() is %d\n", rc);
					exit(-1);
				}
			}
	/* ----- MULTI-THREADING COVERAGE CALCULATION ----- ENDS ----- */
			int out;
			if( new_mn <= 10*NOTATION ){
				FILE *bestMotifs = open_file(bestMotifs, "tmp/best_motifs.dat", "aw" );
				fprintf(bestMotifs, "%.2lf\n", th);
				fclose(bestMotifs);
				printf("- "); fflush(stdout);
				th -= 0.1;
				out = 1;
				if( th < 0.5 ){
					printf("   Breaking the loop!\n");
	/* 				free2Dchar(new_motifs, new_mn);	 */
					break;
				}
			}
			else{
				out = 0;
			}
			
			char *fname = malloc(50*sizeof(char));
			memset(fname, '\0', 50*sizeof(char));
			sprintf(fname,"tmp/motifs_%dnt.dat", ms);
			char *fname2 = malloc(50*sizeof(char));
			memset(fname2, '\0', 50*sizeof(char));
			sprintf(fname2,"tmp/motifs_%dnt_distribution.dat", ms);
			FILE *myTMP = open_file(myTMP, fname, "w" );
			for( i=0; i<new_mn; i++ ){
				fprintf(myTMP, "%s %d %d %d %d\n", new_motifs[i], (int)(numPos*posCov[i]), numPos, (int)(numNeg*negCov[i]), numNeg);
			}
			fclose(myTMP);
			FILE *myTMP2 = open_file(myTMP2, fname2, "w" );
			for( i=0; i<new_mn; i++ ){
				fprintf(myTMP2, "%s ", new_motifs[i]);
				for( ii=0; ii<numPos; ii++ ){
					fprintf(myTMP2, " %d", pmotDist[i][ii]);
				}
				fprintf(myTMP2, " \n");
			}
			fclose(myTMP2);

			call_R(ms, th, out);

			for( i=0; i<new_mn ; i++ ){
				free(new_motifs[i]);
				free(pmotDist[i]);
				free(nmotDist[i]);
			}
			free(new_motifs);
			free(pmotDist);
			free(nmotDist);
			free(posCov);
			free(negCov);

			myTMP = open_file(myTMP, fname, "r" );
			new_mn = read_lines(myTMP);
			new_motifs = malloc(new_mn*sizeof(char *));
			for( i=0; i<new_mn ; i++ ){
				new_motifs[i] = malloc((ms+1)*sizeof(char));
				memset(new_motifs[i], '\0', (ms+1)*sizeof(char));
			}
			posCov = calloc(new_mn, sizeof(double));
			negCov = calloc(new_mn, sizeof(double));

			i = 0;
			char *mline = malloc(200*sizeof(char));
			char *mot = malloc(50*sizeof(char));
			while( !feof(myTMP) ){
				memset(mline, '\0', 200*sizeof(char));
				fscanf(myTMP, "%[^\n]%*c", mline);
				if( *mline != '\0' ){
					memset(mot, '\0', 50*sizeof(char));
					int tp, fp, fn, tn;
					double pval, nval;
					sscanf(mline, "%s %d %d %d %d %lf", mot, &tp, &fp, &fn, &tn, &pval, &nval);
					memmove(new_motifs[i], mot, ms+1);
					posCov[i] = (double)tp/numPos;
					negCov[i] = (double)fn/numNeg;
					i++;
				}
			}
			free(mot);
			free(mline);
			fclose(myTMP);
			free(fname);
			free(fname2);

			tmp_mn = old_mn;
	  		old_motifs = new_motifs;
	  		printf("done\n"); fflush(stdout);

		}
		else{
			printf("   Breaking the loop!\n");
/* 			free2Dchar(new_motifs, new_mn);	 */
			break;
		}
	
		loop++;
	}while( loop < MAX_MOT-MIN_MOT );

	goto here;

	printf("Calculating the motif distribution in the sets... "); fflush(stdout);
/* 	fprintf(log, "Calculating the motif distribution in the sets.. "); */
	tmp_mn = new_mn;

	int *posCount, *negCount, *pshufCount, *nshufCount, ***pshuf_matches, ***nshuf_matches;
	double *posPval, *negPval;
	posCount = calloc(tmp_mn, sizeof(int));
	negCount = calloc(tmp_mn, sizeof(int));
	posPval = calloc(tmp_mn, sizeof(double));
	negPval = calloc(tmp_mn, sizeof(double));
	print_and_continue("0");
	pmotDist = realloc(pmotDist, tmp_mn*numPos*sizeof(int *)*sizeof(int));
	nmotDist = realloc(nmotDist, tmp_mn*numNeg*sizeof(int *)*sizeof(int));
	memset(pmotDist, 0, tmp_mn*numPos*sizeof(int *)*sizeof(int));
	memset(nmotDist, 0, tmp_mn*numNeg*sizeof(int *)*sizeof(int));
	print_and_continue("1");
/* ----- MULTI-THREADING COUNT ----- BEGINS ----- */
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	thread_data_array[0].thread_id = 0;
	thread_data_array[0].motNum = new_mn;
	thread_data_array[0].arrMot = new_motifs;
	thread_data_array[0].motLen = ms;
	thread_data_array[0].count = posCount;
	thread_data_array[0].motDistr = pmotDist;
	
	thread_data_array[1].thread_id = 1;
	thread_data_array[1].motNum = new_mn;
	thread_data_array[1].arrMot = new_motifs;
	thread_data_array[1].motLen = ms;
	thread_data_array[1].count = negCount;
	thread_data_array[1].motDistr = nmotDist;
	
	for( i=0; i<NUM_THREADS; i++ ){
		rc = pthread_create(&threads[i], &attr, set_count, (void *) &thread_data_array[i]);
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
	print_and_continue("2");
	int p, cp, cn;
	double ppval, npval;
	for( i=0; i<tmp_mn; i++ ){
		ppval = npval = 0;
		for( b=0; b<NUM_BOOT; b++ ){
			cp = cn = 0;
			for( p=0; p<numPos; p++ ){
				cp += newOrder[b][p]<numPos ? pmotDist[i][newOrder[b][p]] : nmotDist[i][newOrder[b][p]-numPos];
			}
			ppval += posCount[i] >= cp ? 1 : 0;
			for( p=numPos; p<numTot; p++ ){
				cn += newOrder[b][p]<numPos ? pmotDist[i][newOrder[b][p]] : nmotDist[i][newOrder[b][p]-numPos];
			}
			npval += negCount[i] >= cn ? 1 : 0;
		}
		posPval[i] = (1+ppval)/(NUM_BOOT+1); 
		negPval[i] = (1+npval)/(NUM_BOOT+1);
	}
	
/* ----- MULTI-THREADING COVERAGE CALCULATION ----- ENDS ----- */
	
	printf("done\nCalculating the optimal distance between the sets... "); fflush(stdout);	  	
/* 	fprintf(log, "done\nCalculating the optimal distance between the sets... "); */
			  	
	FILE *FoutAll = open_file(FoutAll, "motifs_last_loop.txt", "w" );
 	for( i=0; i<new_mn; i++ ){
		if( new_motifs[i][0] != '-' && new_motifs[i][ms-2] != '-' ){
			fprintf(FoutAll, "%s %.3f %.3f %d %d %.4f %.4f\n", new_motifs[i], posCov[i], negCov[i], posCount[i], negCount[i], 1-posPval[i], 1-negPval[i]);
		}
	}

	fclose(FoutAll);
	
	here:
	
	printf("done\nFreeing the memory... "); fflush(stdout);
/*  	fprintf(log, "done\nFreeing the memory... "); */

	free(fileM);
	free(fileP);
	free(fileN);
	free(type);
print_and_continue("0");
/* 	free2Dchar(motifs, mn); */
	free2Dchar(posID, numPos);
	free2Dchar(posSeq, numPos);
	if( negID ){
		print_and_continue("1a");
		free2Dchar(negID, numNeg);
		free2Dchar(negSeq, numNeg);
	}
	else{
		print_and_continue("1b");
		free2Dchar(negSeq, numPos);
	}
/* 	free3Dchar(PshuSeq, numPos, NUM_SHUFFLE); */
/* 	free3Dchar(NshuSeq, numNeg, NUM_SHUFFLE); */
	free2Dint(newOrder, NUM_BOOT);
print_and_continue("2");
	free(posCov);
	free(negCov);
print_and_continue("3");
/* 
	free(backup_pcov);
	free(backup_ncov);
 */
print_and_continue("4");
	if( loop != 0 ){
		printf("%d", loop);
		print_and_continue("a");
/* 		free2Dchar(backup_mot, tmp_mn); */
		free2Dchar(new_motifs, new_mn);
	}
print_and_continue("6");
	free(posCount);
	free(negCount);
	free(posPval);
	free(negPval);
print_and_continue("7");
	free2Dint(pmotDist, tmp_mn);
	free2Dint(nmotDist, tmp_mn);

	printf("done\nThe script executed successfully!\n"); fflush(stdout);
/* 	fprintf(log, "done\nThe script executed successfully!\n"); */
	
	pthread_exit(NULL);
		
}