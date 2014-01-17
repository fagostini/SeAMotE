static void compl_inverse(int n, char ***array) { 
	int i, p;
/* 	char *c = malloc(2*sizeof(char)); */
	for( i=0; i<n; i++ ){
		int len = strlen((*array)[i]);
		for( p=0; p<len; p++){
			if((*array)[i][len-p-1] == 'A'){
				(*array)[i][len+p] = 'T';
			}
			else if((*array)[i][len-p-1] == 'C'){
				(*array)[i][len+p] = 'G';
			} 
			else if((*array)[i][len-p-1] == 'G'){
				(*array)[i][len+p] = 'C';
			} 
			else if((*array)[i][len-p-1] == 'T'){
				(*array)[i][len+p] = 'A';
			} 
		}
		(*array)[i][len+p] = '\0';
	}
}