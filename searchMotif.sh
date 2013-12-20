POS="meme_seq/TAF15.posi.seq"
NEG="meme_seq/TAF15.nega.seq"

for m1 in A C U G; do for m2 in A C U G; do for m3 in A C U G; do grep -c "$m1$m2$m3" $POS | awk '{print "'$m1$m2$m3'", $1}'; done; done; done > a

for n in 4 ; do
	echo $n
	for m in `cat a | awk '($2!=0){print $1}'`; do for m1 in A C U G; do grep -c "$m$m1" $POS | awk '{print "'$m$m1'", $1}'; done; done > b; mv b a
done

for j in `cat a | awk '($2!=0){print $1}'`; do awk 'BEGIN{printf "%s", "'$j'"; c=0}($2~/'$j'/){printf " 1"; c++} ($2!~/'$j'/){printf " 0"}END{printf " %d\n", c}' $POS; done > b

for j in `cat a | awk '($2!=0){print $1}'`; do awk 'BEGIN{printf "%s", "'$j'"; c=0}($2~/'$j'/){printf " 1"; c++} ($2!~/'$j'/){printf " 0"}END{printf " %d\n", c}' $NEG; done > c

paste b c | awk '($1002>$2004)' > d

cat mergeMot.R | R --slave --vanilla