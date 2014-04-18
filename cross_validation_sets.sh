for j in `ls meme_seq/* | grep posi.seq`; do 
	na=`echo "$j" | sed 's/\.posi.seq//g'`;
	awk 'BEGIN{srand(1)} {print $1, $2, rand()}' "$na".posi.seq | sort -gk3 | awk '{a[NR]=$1; b[NR]=$2}END{si=NR/3; c=0; for(i=1;i<=NR;i++){if(c<=si){print a[i], b[i] >> "'$na'.CV1.posi.seq"; c++} if(c>si && c<=(si*2)){print a[i], b[i] >> "'$na'.CV2.posi.seq"; c++} if(c>(si*2)){print a[i], b[i] >> "'$na'.CV3.posi.seq"; c++}}}';
	awk 'BEGIN{srand(1)} {print $1, $2, rand()}' "$na".nega.seq | sort -gk3 | awk '{a[NR]=$1; b[NR]=$2}END{si=NR/3; c=0; for(i=1;i<=NR;i++){if(c<=si){print a[i], b[i] >> "'$na'.CV1.nega.seq"; c++} if(c>si && c<=(si*2)){print a[i], b[i] >> "'$na'.CV2.nega.seq"; c++} if(c>(si*2)){print a[i], b[i] >> "'$na'.CV3.nega.seq"; c++}}}';
done
