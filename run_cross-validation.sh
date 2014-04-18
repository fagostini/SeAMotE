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

for d in `ls meme_seq/ | grep CV3.posi.seq`; do 
	na=`echo "$d" | sed 's/\.CV3.posi.seq//g'`; 
	echo "$na"
	sort -gk7 cross_validation_results/"$na".CV1CV3.txt | head -6 | awk '{print $1}' > mt
	sort -gk7 cross_validation_results/"$na".CV1CV2.txt | head -6 | awk '{print $1}' >> mt
	sort -gk7 cross_validation_results/"$na".CV3CV2.txt | head -6 | awk '{print $1}' >> mt
# 	sort -gk7 cross_validation_results/"$na".CV1.txt | head -6 | awk '{print $1}' >> mt
# 	sort -gk7 cross_validation_results/"$na".CV2.txt | head -6 | awk '{print $1}' >> mt
# 	sort -gk7 cross_validation_results/"$na".CV3.txt | head -6 | awk '{print $1}' >> mt
	sort -u mt > mt2; mv mt2 mt
	for j in `cat mt`; do
		grep -w "$j" cross_validation_results/"$na".CV*CV* | awk 'END{printf "%d ", NR}';
# 		grep -w "$j" cross_validation_results/"$na".CV* | awk 'END{printf "%d ", NR}';
	done
	echo ""
done | awk '($0 ~/ 2/ || $0 ~/ 3/)' | wc -l

# | awk '($0 ~/ 4/ || $0 ~/ 5/ || $0 ~/ 6/)' | wc -l