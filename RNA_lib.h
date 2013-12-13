#define pA 27
#define pC 24
#define pG 24
#define pT 25

int rand_int(int n) {
  int limit = RAND_MAX - RAND_MAX % n;
  int rnd;

  do {
    rnd = rand();
  } while (rnd >= limit);
  return rnd % n;
}

static int str_shuffle(char *array, int n) { 
	int i, j;
	char tmp;
	for ( i=n-1; i>=0; i--) {
 	  j = rand_int(i + 1);
	  tmp = array[j];
	  array[j] = array[i];
	  array[i] = tmp;
	}
  
  	return 0;
}

static int str_random(char *array, int n) { 
	int i, p;
	for( i=0; i<n; i++){
		p = rand() % 100;
		if( p < pA ){
			array[i] = 'A';
		}
		else if( p < pA+pC ){
			array[i] = 'C';
		}
		else if( p < pA+pC+pG ){
			array[i] = 'G';
		}
		else{
			array[i] = 'T';
		}
	}
  
  	return 0;
}

static void str_shuffle_arr(int n, char **array) { 
	int i, j;
	char tmp;
	for ( i=n-1; i>=0; i--) {
 	  j = rand_int(i + 1);
	  tmp = (*array)[j];
	  (*array)[j] = (*array)[i];
	  (*array)[i] = tmp;
	}
}

static void random_sequence(int seqLen, char **sequence){
	int i, p;
	(*sequence) = malloc((seqLen+1)*sizeof(char));
	for( i=0; i<seqLen; i++){
		p = rand() % 100;
		if( p < pA ){
			strncat((*sequence), "A", 1);
		}
		else if( p < pA+pC ){
			strncat((*sequence), "C", 1);
		}
		else if( p < pA+pC+pG ){
			strncat((*sequence), "G", 1);
		}
		else{
			strncat((*sequence), "T", 1);
		}
	}
	(*sequence)[seqLen] = '\0';
}

static void create_motifs_nt(char **array, char mlen){
	int i, ii, c, cc;
	int n = SEED_NOTA;
	char pool[] = "ACGT";
 	int n2 = NOTATION;
	char pool2[] = "ACGTRYSWKMBDHV";
	int i1, i2, i3;
	c = 0;
	for( i1=0; i1<SEED_NOTA; i1++ ){
		for( i2=0; i2<NOTATION; i2++ ){
			for( i3=0; i3<SEED_NOTA; i3++ ){
				array[c][0] = pool[i1];
				array[c][1] = pool2[i2];
				array[c][2] = pool[i3];
				array[c++][3] = '\0';
			}
		}
	}
	
/* 
	for( i=0; i<mlen; i++ ){
		int step = pow(n, mlen-i-1);
		c = ii = cc = 0;
		do{
			array[c++][i] = pool[ii];
			cc++;
			ii += cc < step ? 0 : 1;
			cc = cc < step ? cc : 0;
			ii = ii == n ? 0 : ii;
		}while( c<pow(n, mlen));
	} 	
 */
}

static int check_motif(char *mstr, char mlen, double th){
	int i, dash, no;
	dash = no = 0;
	for( i=0; i<mlen; i++ ){
		dash += mstr[i] == '-' ? 1 : 0;
	}
	if( (double)dash/mlen > th || (mstr[0] == '-' && mstr[mlen-1] == '-') ){
		return 1;	
	}
	
	return 0;
}

/* This is a better version of the function */
static int filter_and_expand_nt(char **oldSet, double *pcov, double *ncov, int numMot, int newLen, double threshold, char ***newSet){
	
	printf("   Expanding %d motifs... ", numMot); fflush(stdout); 
	int i, ii, u;
	int n = SEED_NOTA*2;
	char pool[] = "ACGT";
/* 	int n = NOTATION*2; */
	char pool2[] = "ACGTRYSWKMBDHV";
	char **tmpSet = malloc(numMot*n*14*sizeof(char*));
	for( i=0; i<=numMot*n*14; i++ ){
		tmpSet[i] = malloc((newLen+1)*sizeof(char));
		memset(tmpSet[i], '\0', (newLen+1)*sizeof(char));
	}
	
/* 	FILE *tmp_out, *tmp_in; */
/* 	if( (tmp_out = fopen("tmpMot.txt", "w")) == NULL ){ */\
	FILE *temp;
	if( (temp = tmpfile()) == NULL ) {
		printf("Fail\nCannot open the temporary file to store motifs.\n");
		exit(1);
	}
	int c = 0;
	char *motA = malloc((newLen+1)*sizeof(char));
	memset(motA, '\0', (newLen+1)*sizeof(char));

/* 
	char *motB = malloc((newLen+1)*sizeof(char));
	memset(motB, '\0', (newLen+1)*sizeof(char));
 */
/* 
	char *motC = malloc((newLen+1)*sizeof(char));
	char *motD = malloc((newLen+1)*sizeof(char));
	memset(motC, '\0', (newLen+1)*sizeof(char));
	memset(motD, '\0', (newLen+1)*sizeof(char));
 */
	for( i=0; i<numMot; i++ ){
		for( ii=0; ii<n; ii++ ){
			/* Add the nucleotide on the right side */
			memset(motA, '\0', (newLen+1)*sizeof(char));
			memmove(motA, oldSet[i], newLen-1);
			strncat(motA, &pool[ii], 1);
			if( motA[newLen-3] == 'N' ){
				motA[newLen-2] = 'N';
				for( u=0; u<14; u++ ){
					motA[newLen-3] = pool2[u];
					fprintf(temp, "%s\n", motA);
				}
			}
			else{
				fprintf(temp, "%s\n", motA);
				motA[newLen-2] = 'N';
				fprintf(temp, "%s\n", motA);
			}
/* 
			if( check_motif(motA, newLen, 0.35) == 0 ){
				fprintf(tmp_out, "%s\n", motA);
			}
*/
			/* Add the nucleotide on the left side */
/* 
			memset(motB, '\0', (newLen+1)*sizeof(char));
			memcpy(motB, &pool[ii], 1);
			strncat(motB, oldSet[i], newLen-1);
			fprintf(tmp_out, "%s\n", motB);
*/
/* 
			if( check_motif(motB, newLen, 0.35) == 0 ){
				fprintf(tmp_out, "%s\n", motB);
			}
*/
			/* Add the nucleotide on the left and right side with a gap */
/* 
			memset(motC, '\0', (newLen+1)*sizeof(char));
			memset(motD, '\0', (newLen+1)*sizeof(char));
			memcpy(motC, oldSet[i], newLen-1);
			motC[newLen-2] = '-';
			strncat(motC, &pool[ii], 1);
			memcpy(motD, &pool[ii], 1);
			strncat(motD, oldSet[i], newLen-1);
			motD[1] = '-';
			if( check_motif(motC, newLen, 0.35) == 0 ){
				fprintf(tmp_out, "%s\n", motC);
			}
			if( check_motif(motD, newLen, 0.35) == 0 ){
				fprintf(tmp_out, "%s\n", motD);
			}
*/
		}
	}
/* 	fclose(tmp_out); */
	rewind(temp);
	free(motA);
/* 	free(motB); */
/* 
	free(motC);
	free(motD);
 */
/* 
	if( (tmp_in = fopen("tmpMot.txt", "r")) == NULL ){
		printf("Fail\nCannot open the temporary file to store motifs.\n");
		exit(1);
	}
 */

	int numUniq = 0;
	char *mot = malloc((newLen+1)*sizeof(char));
	while(!feof(temp)){
		memset(mot, '\0', (newLen+1)*sizeof(char));
		fscanf(temp,"%[^\n]%*c", mot);
		if( *mot != '\0' ){
			int found = c = 0;
			do{
				if( strncmp(tmpSet[c++], mot, newLen+1) == 0 ){
					found = 1;
					break;
				}
			}while( c <= numUniq );
			if( found == 0 ){
				memmove(tmpSet[numUniq++], mot, newLen+1);
			}
		}
	}
/* 	fclose(tmp_in); */
	fclose(temp);
	free(mot);
/* 	remove("tmpMot.txt"); */

	for( i=0; i<numMot; i++ ){
		free(oldSet[i]);
	}
	free(oldSet);
	
	(*newSet) = malloc(numUniq*sizeof(char *));
	for( i=0; i<numUniq; i++ ){
		(*newSet)[i] = malloc((newLen+1)*sizeof(char));
		memset((*newSet)[i], '\0', (newLen+1)*sizeof(char));
		memmove((*newSet)[i], tmpSet[i], newLen);
	}
	
	for( i=0; i<numMot*n*14; i++ ){
		free(tmpSet[i]);
	}
	free(tmpSet);
	
	printf("done\n"); fflush(stdout); 
	
	return numUniq;
}

static int filter_and_expand_iupac(char **oldSet, int numMot, int newLen, char ***newSet){
	
/* 	printf("   Expanding %d motifs... ", numMot); fflush(stdout);  */
	int i, ii, u;
	int n = SEED_NOTA*2;
	char pool[] = "ACGT";
	int n2 = NOTATION;
	char pool2[] = "ACGTRYSWKMBDHV";
	char **tmpSet = malloc(numMot*NOTATION*NOTATION*sizeof(char*));
	for( i=0; i<=numMot*NOTATION*NOTATION; i++ ){
		tmpSet[i] = malloc((newLen)*sizeof(char));
		memset(tmpSet[i], '\0', (newLen)*sizeof(char));
	}
	
	FILE *temp;
	if( (temp = tmpfile()) == NULL ) {
		printf("Fail\nCannot open the temporary file to store motifs.\n");
		exit(1);
	}
	
	char *motA = malloc((newLen)*sizeof(char));
	memset(motA, '\0', (newLen)*sizeof(char));
	int c = 0;
	for( i=0; i<numMot; i++ ){
		memset(motA, '\0', (newLen)*sizeof(char));
		memmove(motA, oldSet[i], newLen);
		if( motA[newLen-2] == 'N' ){
			for( u=0; u<14; u++ ){
				for( c=0; c<14; c++ ){
					motA[newLen-2] = pool2[u];
					motA[newLen-1] = pool2[c];
					fprintf(temp, "%s\n", motA);
				}
			}
		}
		else{
			fprintf(temp, "%s\n", motA);
		}
	}
	rewind(temp);
	free(motA);

	int numUniq = 0;
	char *mot = malloc((newLen)*sizeof(char));
	while(!feof(temp)){
		memset(mot, '\0', (newLen)*sizeof(char));
		fscanf(temp,"%[^\n]%*c", mot);
		if( *mot != '\0' ){
			int found = c = 0;
			do{
				if( strncmp(tmpSet[c++], mot, newLen) == 0 ){
					found = 1;
					break;
				}
			}while( c <= numUniq );
			if( found == 0 ){
				memmove(tmpSet[numUniq++], mot, newLen);
			}
		}
	}
	fclose(temp);
	free(mot);
	
	(*newSet) = malloc(numUniq*sizeof(char *));
	for( i=0; i<numUniq; i++ ){
		(*newSet)[i] = malloc((newLen)*sizeof(char));
		memset((*newSet)[i], '\0', (newLen)*sizeof(char));
		memmove((*newSet)[i], tmpSet[i], newLen);
	}
	
	for( i=0; i<numMot*NOTATION*NOTATION; i++ ){
		free(tmpSet[i]);
	}
	free(tmpSet);
	
/* 	printf("done\n"); fflush(stdout);  */
	
	return numUniq;
}
