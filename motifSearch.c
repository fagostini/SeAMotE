#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include <pthread.h>

#include "motifSearch.h"
#include "my_library.lib"
#include "RNA_lib.lib"
#include "DNA_lib.lib"

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

	int taskid, seqNum, motNum, motLen;
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
   	
	int m, s, pos, i, br, count;
	for( m=0; m<motNum; m++ ){
		count = 0;
		for( s=0; s<seqNum; s++ ){
			pos = br = 0;
			do{
				i = 0;
				do{
					if( motSet[m][i] != '-' && motSet[m][i] != seqSet[s][pos+i]){
						break;
					}
					i++;
					if( i == motLen ){
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
					if( motSet[m][i] != '-' && motSet[m][i] != seqSet[s][pos+i]){
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
void select_motifs(int *numMot, int motLen, char ***motSet, double **posCov, double **negCov, double threshold){
	int perc = (int)((*numMot)*threshold);
/* 
	double **selCov = calloc(2, sizeof(double *));
	selCov[0] = calloc(perc*3,sizeof(double));
	selCov[1] = calloc(perc*3,sizeof(double));
 */
	int m;
/* 
	char **selMot = malloc(perc*3*sizeof(char *));
	for( m=0; m<perc*3; m++ ){
		selMot[m] = malloc((motLen+1)*sizeof(char));
		memset(selMot[m], '\0', (motLen+1)*sizeof(char));
	}
 */
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
/* 
		memmove(selMot[m], (*motSet)[(int)index[m][2]], motLen+1);
		selCov[0][m] = index[m][0];
		selCov[1][m] = index[m][1];
 */
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
/* 
		memmove(selMot[m], (*motSet)[(int)index[m][2]], motLen+1);
		selCov[0][m] = index[m][1];
		selCov[1][m] = index[m][0];
 */
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
/* 
		memmove(selMot[m], (*motSet)[(int)index[m][3]], motLen+1);
		selCov[0][m] = index[m][1];
		selCov[1][m] = index[m][2];
 */
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

int main(int argc, char* argv[]){
	
/* 	srand(time(NULL)); */
	srand(1);	
	
	int i, numPos, numNeg, dna, rna, prot, ref;
	numPos = numNeg = ref = 0;
	int ms = MIN_MOT;
	int mn = pow(4,ms);
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

	printf("Calculating the seed coverages... "); fflush(stdout);
/* 	fprintf(log, "Calculating the seed coverages... "); */
	double *posCov = calloc(mn, sizeof(posCov));
	double *negCov = calloc(mn, sizeof(negCov));

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
	thread_data_array[1].thread_id = 1;
	thread_data_array[1].depth = 1;
	thread_data_array[1].arrSeq = negSeq;
	thread_data_array[1].arrMot = motifs;
	thread_data_array[1].arrNum = numNeg;
	thread_data_array[1].motNum = mn;
	thread_data_array[1].motLen = ms;
	thread_data_array[1].cov = negCov;

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

	char **backup_mot = malloc(mn*sizeof(char *));
	for( i=0; i<mn; i++ ){
		backup_mot[i] = malloc((ms+1)*sizeof(char));
		memset(backup_mot[i], '\0', (ms+1)*sizeof(char));
	}
	double *backup_pcov = calloc(mn, sizeof(double));
	double *backup_ncov = calloc(mn, sizeof(double));
	size_t size_mot = mn*(ms+1)*sizeof(char)*sizeof(char *);
	size_t size_cov = mn*sizeof(double);
	memmove(backup_mot, motifs, size_mot);
	memmove(backup_pcov, posCov, size_cov);
	memmove(backup_ncov, negCov, size_cov);
	
/* 
	for( i=0; i<mn; i++ ){
		printf("%s %.2lf %.2lf\n", backup_mot[i], backup_pcov[i], backup_ncov[i] ); fflush(stdout);
	}
 */

	printf("done\nBeginnig the motif coverage optimization:\n"); fflush(stdout);
/* 	fprintf(log, "done\nBeginnig the motif coverage optimization:\n"); */
	int old_mn, tmp_mn;
	char **old_motifs = motifs;
	char **new_motifs = NULL;
	int new_mn = tmp_mn = mn;
	int loop = 0;
 	do{
		old_mn = new_mn;
		new_mn = above_threshold(thread_data_array[0].cov, thread_data_array[1].cov, old_mn, th)*16;
/* 
		if( new_mn == 0 && th > 0.5 ){
			th -= TH_STEP;
			printf("   Lowering the threshold to %.2f of coverage\n", th);
			new_mn = above_threshold(thread_data_array[0].cov, thread_data_array[1].cov, old_mn, th)*16;
		}
 */
		if( new_mn != 0 ){
			/* MOTIFS */
			size_mot = old_mn*(ms)*sizeof(char)*sizeof(char*);
			backup_mot = realloc(backup_mot, size_mot);
			memmove(backup_mot, old_motifs, size_mot);
			ms++;
			new_mn = filter_and_expand_nt(old_motifs, thread_data_array[0].cov, thread_data_array[1].cov, old_mn, ms, th, &new_motifs);
		 	printf("   %d: Testing %d motifs (%d nt) ", loop+1, new_mn, ms); fflush(stdout);
/* 			fprintf(log, "   %d: Testing %d motifs (%d nt)... ", loop+1, new_mn, ms); */
 			/* POSITIVE COVERAGE */
 			size_cov = old_mn*sizeof(double);
 			backup_pcov = realloc(backup_pcov, size_cov);
			memmove(backup_pcov, thread_data_array[0].cov, size_cov);
			free(posCov);
			posCov = calloc(new_mn, sizeof(double));
			printf("."); fflush(stdout);
 			/* NEGATIVE COVERAGE */
 			backup_ncov = realloc(backup_ncov, size_cov);
			memmove(backup_ncov, thread_data_array[1].cov, size_cov);
			free(negCov);
			negCov = calloc(new_mn, sizeof(double));
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
			printf("."); fflush(stdout);
			for( i=0; i<NUM_THREADS; i++ ){
				rc = pthread_create(&threads[i], &attr, set_coverage, (void *) &thread_data_array[i]);
				if (rc) {
					printf("ERROR; return code from pthread_create() is %d\n", rc);
					exit(-1);
				}
			}
			pthread_attr_destroy(&attr);
			printf("."); fflush(stdout);
			for( i=0; i<NUM_THREADS; i++ ){
				rc = pthread_join(threads[i], &status);
				if (rc) {
					printf("ERROR; return code from pthread_join() is %d\n", rc);
					exit(-1);
				}
			}
	/* ----- MULTI-THREADING COVERAGE CALCULATION ----- ENDS ----- */
			printf("."); fflush(stdout);
			select_motifs(&new_mn, ms, &new_motifs, &posCov, &negCov, 0.05);
			thread_data_array[0].arrMot = new_motifs;
			thread_data_array[1].arrMot = new_motifs;
			thread_data_array[0].cov = posCov;
			thread_data_array[1].cov = negCov;
			thread_data_array[0].motNum = new_mn;
			thread_data_array[1].motNum = new_mn;
			
/* 
			for(i=0;i<new_mn;i++){
				printf("%s %.2lf %.2lf\n", new_motifs[i], posCov[i], negCov[i]);
			}
 */
			printf(". "); fflush(stdout);
			tmp_mn = old_mn;
	  		old_motifs = new_motifs;
	  		printf("done\n"); fflush(stdout);
	  		char fname[50];
			sprintf(fname,"tmp/motifs_%dnt.dat", ms-1);
	  		FILE *myTMP = open_file(myTMP, fname, "w" );
	  		for( i=0; i<tmp_mn; i++ ){
	  			if( backup_mot[i][0] != '-' && backup_mot[i][ms-2] != '-' ){
	  				fprintf(myTMP, "%s %.3f %.3f\n", backup_mot[i], backup_pcov[i], backup_ncov[i]);
	  			}
	  		}
	  		fclose(myTMP);
/* 			fprintf(log, "done\n"); */
		}
		else{
			printf("   Breaking the loop!\n");
			free2Dchar(new_motifs, new_mn);	
			break;
		}
		loop++;
 	}while( loop < MAX_MOT-MIN_MOT ); /* CHECK THE OLD_MN, TMP_MN AND NEW_MN */

	printf("Calculating the motif distribution in the sets... "); fflush(stdout);
/* 	fprintf(log, "Calculating the motif distribution in the sets.. "); */

	int *posCount, *negCount, *pshufCount, *nshufCount, ***pshuf_matches, ***nshuf_matches, **pmotDist, **nmotDist;
	double *posPval, *negPval;
	posCount = calloc(tmp_mn, sizeof(int));
	negCount = calloc(tmp_mn, sizeof(int));
	posPval = calloc(tmp_mn, sizeof(double));
	negPval = calloc(tmp_mn, sizeof(double));
	
	pmotDist = calloc(tmp_mn, sizeof(*pmotDist));
	for( b=0; b<tmp_mn; b++ ){
		pmotDist[b] = calloc(numPos, sizeof(int));
	}
	nmotDist = calloc(tmp_mn, sizeof(*nmotDist));
	for( b=0; b<tmp_mn; b++ ){
		nmotDist[b] = calloc(numNeg, sizeof(int));
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
	for( i=0; i<tmp_mn; i++ ){
/* 		if( backup_pcov[i] >= th || backup_ncov[i] >= th ){ */
		if( backup_mot[i][0] != '-' && backup_mot[i][ms-2] != '-' ){
			fprintf(FoutAll, "%s %.3f %.3f %d %d %.4f %.4f\n", backup_mot[i], backup_pcov[i], backup_ncov[i], posCount[i], negCount[i], 1-posPval[i], 1-negPval[i]);
		}
/* 		} */
	}
	fclose(FoutAll);
	printf("done\nFreeing the memory... "); fflush(stdout);
/*  	fprintf(log, "done\nFreeing the memory... "); */

	free(fileM);
	free(fileP);
	free(fileN);
	free(type);
print_and_continue("0");
	free2Dchar(motifs, mn);
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
	free(backup_pcov);
	free(backup_ncov);
print_and_continue("4");
	if( loop != 0 ){
	printf("%d", loop);
	print_and_continue("a");
		free2Dchar(backup_mot, tmp_mn);
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