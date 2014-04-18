for d in `ls meme_seq/ | grep CV3.posi.seq`; do 
	na=`echo "$d" | sed 's/\.CV3.posi.seq//g'`; 
	cat meme_seq/"$na".CV1.posi.seq meme_seq/"$na".CV2.posi.seq > posi.cv.txt; 
	cat meme_seq/"$na".CV1.nega.seq meme_seq/"$na".CV2.nega.seq > nega.cv.txt; 
	./motifSearch -p posi.cv.txt -n nega.cv.txt -t 0.75 -a rna -w 7 >> cross_validation_results/CV1CV2_log.txt;
	mv tmp/best_motifs.dat cross_validation_results/"$na".CV1CV2.txt
	rm tmp/*.dat
	
	cat meme_seq/"$na".CV1.posi.seq meme_seq/"$na".CV3.posi.seq > posi.cv.txt; 
	cat meme_seq/"$na".CV1.nega.seq meme_seq/"$na".CV3.nega.seq > nega.cv.txt; 
	./motifSearch -p posi.cv.txt -n nega.cv.txt -t 0.75 -a rna -w 7 >> cross_validation_results/CV1CV3_log.txt;
	mv tmp/best_motifs.dat cross_validation_results/"$na".CV1CV3.txt
	tmp/*.dat
	
	cat meme_seq/"$na".CV3.posi.seq meme_seq/"$na".CV2.posi.seq > posi.cv.txt; 
	cat meme_seq/"$na".CV3.nega.seq meme_seq/"$na".CV2.nega.seq > nega.cv.txt; 
	./motifSearch -p posi.cv.txt -n nega.cv.txt -t 0.75 -a rna -w 7 >> cross_validation_results/CV3CV2_log.txt;
	mv tmp/best_motifs.dat cross_validation_results/"$na".CV3CV2.txt
	tmp/*.dat
done

