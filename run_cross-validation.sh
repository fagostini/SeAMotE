# for d in `ls meme_seq/ | grep CV3.posi.seq`; do 
# 	na=`echo "$d" | sed 's/\.CV3.posi.seq//g'`; 
# 	cat meme_seq/"$na".CV1.posi.seq meme_seq/"$na".CV2.posi.seq > posi.cv.txt; 
# 	cat meme_seq/"$na".CV1.nega.seq meme_seq/"$na".CV2.nega.seq > nega.cv.txt; 
# 	./motifSearch -p posi.cv.txt -n nega.cv.txt -t 0.75 -a rna -w 7 >> cross_validation_results/CV1CV2_log.txt;
# 	mv tmp/best_motifs.dat cross_validation_results/"$na".CV1CV2.txt
# 	rm tmp/*.dat
# 	
# 	cat meme_seq/"$na".CV1.posi.seq meme_seq/"$na".CV3.posi.seq > posi.cv.txt; 
# 	cat meme_seq/"$na".CV1.nega.seq meme_seq/"$na".CV3.nega.seq > nega.cv.txt; 
# 	./motifSearch -p posi.cv.txt -n nega.cv.txt -t 0.75 -a rna -w 7 >> cross_validation_results/CV1CV3_log.txt;
# 	mv tmp/best_motifs.dat cross_validation_results/"$na".CV1CV3.txt
# 	rm tmp/*.dat
# 	
# 	cat meme_seq/"$na".CV3.posi.seq meme_seq/"$na".CV2.posi.seq > posi.cv.txt; 
# 	cat meme_seq/"$na".CV3.nega.seq meme_seq/"$na".CV2.nega.seq > nega.cv.txt; 
# 	./motifSearch -p posi.cv.txt -n nega.cv.txt -t 0.75 -a rna -w 7 >> cross_validation_results/CV3CV2_log.txt;
# 	mv tmp/best_motifs.dat cross_validation_results/"$na".CV3CV2.txt
# 	rm tmp/*.dat
# 
# 	./motifSearch -p meme_seq/"$na".CV1.posi.seq -n meme_seq/"$na".CV1.nega.seq -t 0.75 -a rna -w 7 >> cross_validation_results/CV1_log.txt;
# 	mv tmp/best_motifs.dat cross_validation_results/"$na".CV1.txt;
# 	./motifSearch -p meme_seq/"$na".CV2.posi.seq -n meme_seq/"$na".CV2.nega.seq -t 0.75 -a rna -w 7 >> cross_validation_results/CV2_log.txt;
# 	mv tmp/best_motifs.dat cross_validation_results/"$na".CV2.txt;
# 	./motifSearch -p meme_seq/"$na".CV3.posi.seq -n meme_seq/"$na".CV3.nega.seq -t 0.75 -a rna -w 7 >> cross_validation_results/CV3_log.txt;
# 	mv tmp/best_motifs.dat cross_validation_results/"$na".CV3.txt;
# 
# done

# for d in `ls meme_seq/ | grep CV3.posi.seq`; do 
# 	na=`echo "$d" | sed 's/\.CV3.posi.seq//g'`; 
# 	echo "$na"
# #  	sort -gk7 results/best_motifs_"$na".dat | head -5 | awk '{print $1}' > mt
# #  	sort -gk7 results/best_motifs_"$na".dat | head -5 | awk '{for(i=1;i<=length($1);i++){if(i==1){printf ".%s\n", substr($1,i+1,length($1))} if(i==length($1)){printf "%s.\n", substr($1,1,length($1)-1)} if(i>1 && i<length($1)){printf "%s.%s\n", substr($1,1,i-1), substr($1,i,length($1)-i)}}}' > mt
# 	sort -gk7 cross_validation_results/"$na".CV1CV3.txt | head -5 | awk '{print $1}' > mt
# 	sort -gk7 cross_validation_results/"$na".CV1CV2.txt | head -5 | awk '{print $1}' >> mt
# 	sort -gk7 cross_validation_results/"$na".CV3CV2.txt | head -5 | awk '{print $1}' >> mt
# # 	sort -gk7 cross_validation_results/"$na".CV1.txt | head -6 | awk '{print $1}' >> mt
# # 	sort -gk7 cross_validation_results/"$na".CV2.txt | head -6 | awk '{print $1}' >> mt
# # 	sort -gk7 cross_validation_results/"$na".CV3.txt | head -6 | awk '{print $1}' >> mt
# 	sort -u mt > mt2; mv mt2 mt
# 	for j in `cat mt`; do
# 		grep -w "$j" cross_validation_results/"$na".CV*CV* | awk 'END{printf "%d ", NR}';
# # 		grep -w "$j" cross_validation_results/"$na".CV* | awk 'END{printf "%d ", NR}';
# 	done
# 	echo ""
# done
# | awk '($0 ~/ 2/ || $0 ~/ 3/)' | wc -l

# | awk '($0 ~/ 4/ || $0 ~/ 5/ || $0 ~/ 6/)' | wc -l


# for d in `ls meme_seq/ | grep CV3.posi.seq`; do 
# 	na=`echo "$d" | sed 's/\.CV3.posi.seq//g'`; 
# 	echo "$na"
# 
# 	sort -gk7 cross_validation_results/"$na".CV1CV2.txt | head -10 | awk '{print $1}' > mt
# 	for j in `cat mt | awk '(length($1)>3)'`; do
# 		grep -w "$j" cross_validation_results/"$na".CV3.txt | awk 'END{printf "%d ", NR}';
# 	done
# 	echo "A"
# 	
# 	sort -gk7 cross_validation_results/"$na".CV1CV3.txt | head -10 | awk '{print $1}' > mt
# 	for j in `cat mt | awk '(length($1)>3)'`; do
# 		grep -w "$j" cross_validation_results/"$na".CV2.txt | awk 'END{printf "%d ", NR}';
# 	done
# 	echo "B"
# 	
# 	sort -gk7 cross_validation_results/"$na".CV3CV2.txt | head -10 | awk '{print $1}' > mt
# 	for j in `cat mt | awk '(length($1)>3)'`; do
# 		grep -w "$j" cross_validation_results/"$na".CV1.txt | awk 'END{printf "%d ", NR}';
# 	done
# 	echo "C"
# 	
# done

# USING TOMTOM
for d in `ls meme_seq/ | grep CV3.posi.seq`; do 
	na=`echo "$d" | sed 's/\.CV3.posi.seq//g'`; 
	echo "$na"

#	echo "CV1"
	iupac2meme -alpha DNA -pseudo 1 `sort -gk7 cross_validation_results/"$na".CV1CV2.txt | head -3 | awk '{printf "%s ", $1}END{printf "\n"}'` > mt
	iupac2meme -alpha DNA -pseudo 1 `sort -gk7 cross_validation_results/"$na".CV3.txt    | awk '{printf "%s ", $1}END{printf "\n"}'` > mf
	tomtom -verbosity 1 -oc tomtom_out mt mf
	head -3 tomtom_out/tomtom.txt | awk '($4<=0.05)'
	wc -l cross_validation_results/"$na".CV1CV2.txt cross_validation_results/"$na".CV3.txt
# 	for j in `head -3 tomtom_out/tomtom.txt | awk '($4<=0.05){print $1}'`; do grep -w "$j" cross_validation_results/"$na".CV1CV2.txt; done
# 	for j in `head -3 tomtom_out/tomtom.txt | awk '($4<=0.05){print $2}'`; do grep -w "$j" cross_validation_results/"$na".CV3.txt; done
		
#	echo "CV2"	
	iupac2meme -alpha DNA -pseudo 1 `sort -gk7 cross_validation_results/"$na".CV1CV3.txt | head -3 | awk '{printf "%s ", $1}END{printf "\n"}'` > mt
	iupac2meme -alpha DNA -pseudo 1 `sort -gk7 cross_validation_results/"$na".CV2.txt    | awk '{printf "%s ", $1}END{printf "\n"}'` > mf
	tomtom -verbosity 1 -oc tomtom_out mt mf
	head -3 tomtom_out/tomtom.txt | awk '($4<=0.05)'
	wc -l cross_validation_results/"$na".CV1CV3.txt cross_validation_results/"$na".CV2.txt
# 	for j in `head -3 tomtom_out/tomtom.txt | awk '($4<=0.05){print $1}'`; do grep -w "$j" cross_validation_results/"$na".CV1CV3.txt; done
# 	for j in `head -3 tomtom_out/tomtom.txt | awk '($4<=0.05){print $2}'`; do grep -w "$j" cross_validation_results/"$na".CV2.txt; done

#	echo "CV3"
	iupac2meme -alpha DNA -pseudo 1 `sort -gk7 cross_validation_results/"$na".CV3CV2.txt | head -3 | awk '{printf "%s ", $1}END{printf "\n"}'` > mt
	iupac2meme -alpha DNA -pseudo 1 `sort -gk7 cross_validation_results/"$na".CV1.txt    | awk '{printf "%s ", $1}END{printf "\n"}'` > mf
	tomtom -verbosity 1 -oc tomtom_out mt mf
	head -3 tomtom_out/tomtom.txt | awk '($4<=0.05)'
	wc -l cross_validation_results/"$na".CV3CV2.txt cross_validation_results/"$na".CV1.txt
# 	for j in `head -3 tomtom_out/tomtom.txt | awk '($4<=0.05){print $1}'`; do grep -w "$j" cross_validation_results/"$na".CV3CV2.txt; done
# 	for j in `head -3 tomtom_out/tomtom.txt | awk '($4<=0.05){print $2}'`; do grep -w "$j" cross_validation_results/"$na".CV1.txt; done
	
done
