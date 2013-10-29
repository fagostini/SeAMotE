static void compl_inverse(int n, char ***array) { 
	int i, p;
	char *c = malloc(2*sizeof(char));
	for( i=0; i<n; i++ ){
		int len = strlen((*array)[i]);
		for( p=0; p<len; p++){
			memset(c, '\0', 2*sizeof(char));
			if((*array)[i][len-p-1] == 'A'){
				memcpy(c, "T", 1);
			} 
			else if((*array)[i][len-p-1] == 'C'){
				memcpy(c, "G", 1);
			} 
			else if((*array)[i][len-p-1] == 'G'){
				memcpy(c, "C", 1);
			} 
			else if((*array)[i][len-p-1] == 'T'){
				memcpy(c, "A", 1);
			} 
			strncat((*array)[i+n], c, 1);
			free(c);
		}
		(*array)[i+n][len] = '\0';
	}
}