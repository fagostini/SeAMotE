struct data{
		char *cmd;
		char *des;
		char *err;
};

static void call_help(int n, struct data *options){
	printf("\n%s\nUsage: %s [-options]\n\nCommon options are:\n", options[0].des, options[0].cmd);
	printf("\t%s : %s\n\n", options[1].cmd, options[1].des);
	printf("   options controlling inputs and parameters:\n");
	printf("\t%s <file> : %s\n", options[2].cmd, options[2].des);
	printf("\t%s <file> : %s\n", options[3].cmd, options[3].des);
	printf("\t%s <file> : %s\n", options[4].cmd, options[4].des);
	printf("\t%s <double> : %s\n", options[5].cmd, options[5].des);
	printf("\t%s <string> : %s\n", options[6].cmd, options[6].des);
	printf("\t%s <string> : %s\n", options[7].cmd, options[7].des);
	printf("\n");
	exit(n);
}

static void read_args(int argc, char **argv, char *fileM, char *fileP, char *fileN, double *th, char *type, int *ref){
	FILE *Ftest;
	struct data options[8];
	options[0].cmd = argv[0];
	options[0].des = "Version 0.1 - Description of the application (under development)";
	options[0].err = "";
	options[1].cmd = "-h";
	options[1].des = "show brief help on version and usage";
	options[1].err = "";
	options[2].cmd = "-m";
	options[2].des = "the file containing a set of pre-generated motifs (under development)";
	options[2].err = "Error! The motifs file does not exists or has not been specified!";
	options[3].cmd = "-p";
	options[3].des = "the file containing a set of positive sequences";
	options[3].err = "Error! The positive set does not exists or has not been specified!";
	options[4].cmd = "-n";
	options[4].des = "the file containing a set of negative sequences (optional)";
	options[4].err = "Error! The negative set does not exists or has not been specified!";
	options[5].cmd = "-t";
	options[5].des = "set a threshold for each selection loop (default: 0.9)";
	options[5].err = "Error! The threshold has not been specified!";
	options[6].cmd = "-a";
	options[6].des = "you can choose whether to use RNA (the set as it is) or DNA (the set and its reverse complementary)";
	options[6].err = "Error! The molecule type has not been specified!";
	options[7].cmd = "-r";
	options[7].des = "the reference set can consist of a random of shuffle dataset (default: shuffle)";
	options[7].err = "Error! The reference has not been specified correctly!";

	int i = 1;
	while( i < argc ){
		/* HELP */
		if(strcmp(argv[i], options[1].cmd) == 0 ){
			call_help(0, options);
		}
		else if(argv[i][0] == '-'){
			if(strcmp(argv[i], options[2].cmd) == 0){
				if(argv[i+1] != NULL && argv[i+1][0] != '-'){
					if( (Ftest = fopen(argv[i+1], "r")) == NULL ){
						printf("Fail\nCannot open %s file.\n", argv[i+1]);
						exit(1);
					}
					else{
						
						strcpy(fileM, argv[i+1]);
						fclose(Ftest);
					}
				}
				else{
					printf("%s\n", options[2].err);
					exit(1);
				}
			}
			else if(strcmp(argv[i], options[3].cmd) == 0){
				if(argv[i+1] != NULL && argv[i+1][0] != '-'){
					if( (Ftest = fopen(argv[i+1], "r")) == NULL ){
						printf("Fail\nCannot open %s file.\n", argv[i+1]);
						exit(1);
					}
					else{
						strcpy(fileP, argv[i+1]);
						fclose(Ftest);
					}
				}
				else{
					printf("%s\n", options[3].err);
					exit(1);
				}
			}
			else if(strcmp(argv[i], options[4].cmd) == 0){
				if(argv[i+1] != NULL && argv[i+1][0] != '-'){
					if( (Ftest = fopen(argv[i+1], "r")) == NULL ){
						printf("Fail\nCannot open %s file.\n", argv[i+1]);
						exit(1);
					}
					else{
						strcpy(fileN, argv[i+1]);
						fclose(Ftest);
					}
				}
				else{
					printf("%s\n", options[4].err);
					exit(1);
				}
			}
			else if(strcmp(argv[i], options[5].cmd) == 0){
				if(argv[i+1] != NULL && argv[i+1][0] != '-'){
					*th = atof(argv[i+1]);
				}
				else{
					printf("%s\n", options[5].err);
					exit(1);
				}
			}
			else if(strcmp(argv[i], options[6].cmd) == 0){
				if(argv[i+1] != NULL && argv[i+1][0] != '-' && ((strcmp(argv[i+1], "rna") == 0) || (strcmp(argv[i+1], "RNA") == 0) || (strcmp(argv[i+1], "dna") == 0) || (strcmp(argv[i+1], "DNA") == 0)) ){
					*type = atof(argv[i+1]);
				}
				else{
					printf("%s\n", options[6].err);
					exit(1);
				}
			}
			else if(strcmp(argv[i], options[7].cmd) == 0){
				if(argv[i+1] != NULL && argv[i+1][0] != '-' && (strcmp(argv[i+1], "random") == 0) ){
					*ref = 1;
				}
				else if(argv[i+1] != NULL && argv[i+1][0] != '-' && (strcmp(argv[i+1], "shuffle") == 0) ){
					*ref = 0;
				}
				else{
					printf("%s\n", options[7].err);
					exit(1);
				}
			}
			else{
				printf("Failed to parse command line: No such option \"%s\".\n", argv[i]);
				call_help(1, options);
			}
		}
		if( ++i >= argc ){
			break;
		}
		
	}	
}

static FILE *open_file(FILE *myFile, char *filename, char *mode ){
	if( (myFile = fopen(filename, mode)) == NULL ){
		printf("Fail\nCannot open %s file.\n", filename);
		exit(1);
	}
	return myFile;
}

static char **malloc2Dchar(int x, int y){
	int i;
	char **array = malloc(y*sizeof(char *));
	for( i=0; i<y; i++){
		array[i] = malloc((x+1)*sizeof(char));
		memset(array[i], '\0', (x+1)*sizeof(char));
	}
	return array;
}

static char ***malloc3Dchar(int x, int y, int z){
	int i, ii;
	char ***array = malloc(z*sizeof(char **));
	for( i=0; i<z; i++){
		array[i] = malloc(y*sizeof(char *));
		for( ii=0; ii<y; ii++){
			array[i][ii] = malloc((x+1)*sizeof(char));
			memset(array[i][ii], '\0', (x+1)*sizeof(char));
		}
	}
	return array;
}

static void print2Dchar(char **array, int y){
	int i;
	for( i=0; i<y; i++){
		printf("%s\n", array[i]);
	}
}

static void print3Dchar(char ***array, int y, int z){
	int i, ii;
	for( i=0; i<z; i++){
		for( ii=0; ii<y; ii++){
			printf("%s ", array[i][ii]);
		}
		printf("\n");
	}
}

static void free2Dchar(char **array, int y){
	int i;
	for( i=0; i<y; i++){
		free(array[i]);
	}
	free(array);
}

static void free3Dchar(char ***array, int y, int z){
	int i, ii;
	for( i=0; i<z; i++){
		for( ii=0; ii<y; ii++){
			free(array[i][ii]);
		}
		free(array[i]);
	}
	free(array);
}

static void free2Dint(int **array, int y){
	int i;
	for( i=0; i<y; i++){
		free(array[i]);
	}
	free(array);
}

static void free3Dint(int ***array, int y, int z){
	int i, ii;
	for( i=0; i<z; i++){
		for( ii=0; ii<y; ii++){
			free(array[i][ii]);
		}
		free(array[i]);
	}
	free(array);
}

static int bootstrap_sampling(int n, int **order ) {
	int i, j, t;
	int tmp;
	for( t=0; t<NUM_BOOT; t++ ){
		for ( i=n-1; i>=0; i--) {
		  j = rand_int(i + 1);
		  tmp = order[t][j];
		  order[t][j] = order[t][i];
		  order[t][i] = tmp;
		}
	}
	return 0;
}

static int read_lines(FILE *inFile){
	int c = 0;
	while(!feof(inFile)){
		c += getc(inFile) == '\n' ? 1 : 0;
	}
	rewind(inFile);
	
	return c;
}

static int read_fasta(FILE *inFile, int max_seq){
	int num_seq;
	num_seq = 0;
	char *line = malloc(1000*sizeof(char));
	char *id = malloc(100*sizeof(char));
	char *sequence = malloc(max_seq*sizeof(char));
	while(!feof(inFile)){
		memset(line, '\0', 1000*sizeof(char)) ;
		fscanf(inFile, "%[^\n]%*c", line);
		if( line[0] != '\0' ){
			if( line[0] == '>'){
				num_seq++;
				if( (id[0]!='\0') && (sequence[0]!='\0')){
					printf("%s %s\n", id, sequence);
				}
				memset(id, '\0', 100*sizeof(char)) ;
				memmove(id, line, 100);
				memset(sequence, '\0', max_seq*sizeof(char)) ;
			}
			else{
				strncat(sequence, line, strlen(line));
			}
		}
	}
	if( (id[0]!='\0') && (sequence[0]!='\0')){
		printf("%s %s\n", id, sequence);
	}
	free(line);
	free(id);
	free(sequence);

	return num_seq;
	
} /* Can be optimized */

static int read_oneline(FILE *inFile, char **ID, char **SEQ, int max_id, int max_seq, char *format){
	int c = 0;
	char *id = malloc(max_id*sizeof(char));
	char *sequence = malloc(max_seq*sizeof(char));
	while(!feof(inFile)){
    	memset(id, '\0', max_id*sizeof(char)) ;
   		memset(sequence, '\0', max_seq*sizeof(char)) ;
    	fscanf(inFile, format, id, sequence);
    	if( *id != '\0' && *sequence != '\0'){
    		int i;
    		for( i=0; i<strlen(sequence); i++ ){
    			sequence[i] = toupper(sequence[i]);
    			sequence[i] = sequence[i] == 'U' ? 'T' : sequence[i];
    		}
			memmove(ID[c], id, max_id);
			memmove(SEQ[c], sequence, max_seq);
    		c++;
    	}
    }
	free(id);
	free(sequence);
    
    return 0;
}

static void loadBar(int LBx, int LBn, int LBr, int LBw){
	/* x is the current position, n is the final position, r is the number of steps and w is the width of the progress bar */
    if ( LBx % (LBn/LBr) != 0 ) return;
 
    double LBratio = LBx/(double)LBn;
    int   LBc     = LBratio * LBw;
 
    printf("%3d%% [", (int)(LBratio*100) );
 
    for (LBx=0; LBx<LBc; LBx++)
       printf("=");
 
    for (LBx=LBc; LBx<LBw; LBx++)
       printf(" ");
 
     printf("]\n\033[F\033[J");
}

/* This function has been moved in the main file because of the threading implementation */
/* static void set_coverage(char **seqSet, char **motSet, int seqNum, int motNum, int motLen, double *coverage){
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
		coverage[m] = count > 0 ? (double)count/seqNum : 0;
	}
} */

static int above_threshold(double *pcov, double *ncov, int numMot, double threshold){
	int i;
	int c = 0;
	for( i=0; i<numMot; i++ ){
		if( pcov[i]>=threshold || ncov[i]>=threshold ){
			c++;
		}
	}
	return c;
}

static double *subsampling(int ***pool, int setNum, int setLen, int samples, int motNum){
	int i, s, cnt;
	double *cov = (double *)calloc(samples, sizeof(double));
	for( s=0; s<samples; s++ ){
		cnt = 0;
		for( i=0; i<setNum; i++ ){
			int rnd = rand() % setLen;
			cnt += pool[motNum][rnd][i] > 0 ? 1 : 0;
		}
		cov[s] = (double)cnt/setNum;
	}

	return cov;
}