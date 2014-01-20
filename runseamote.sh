#!/bin/bash
# bash runseamote.sh "random_number" "motiFile" "use_as_reference" "threshold" "reverse"

cd   tmp/$1

	if [[ "$3" == "negative" && -s "negative.oneline" ]]; then
			options="-n negative.oneline "
	else
			options="-r $3 "
	fi

	options+="-t $4 "  

	./motifSearch -p positive.oneline $options -a $5
	
	mv positive.oneline outputs/positive.seq
	if [[ -s "negative.oneline" ]]; then
		mv negative.oneline outputs/reference.seq
	else
		mv reference.seq outputs/reference.seq
	fi

	if [[ -s "tmp/best_motifs.dat" ]]; then
		awk '{print $1}' tmp/best_motifs.dat > iupac.txt
		sed 's/R/[AG]/g;s/Y/[CT]/g;s/M/[CA]/g;s/K/[TG]/g;s/W/[TA]/g;s/S/[CG]/g;s/B/[CTG]/g;s/D/[ATG]/g;s/H/[ATC]/g;s/V/[ACG]/g' tmp/best_motifs.dat | awk '{print $1, $2, $4, $2/($2+$3), $4/($4+$5), $6, $NF}' > best_motifs.txt
		paste iupac.txt best_motifs.txt > outputs/best_motifs.txt
		
		for((l=1;l<=`wc -l outputs/best_motifs.txt | awk '{print $1}'`;l++)); do
			IUmot=`awk '(NR=='$l'){print $1}' outputs/best_motifs.txt`;
			REmot=`awk '(NR=='$l'){print $2}' outputs/best_motifs.txt`;
			grep -o -e "$REmot" outputs/positive.seq > matches.txt
			python createLogo.py "$IUmot" 2&> python.log
		done
		rm matches.txt
	
		awk 'BEGIN{printf "<tbody>\n"}{printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td>%s</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%.3E</td>\n\t</tr>\n", NR, $1, $2, $3, $4, $5*100, $6*100, $7*100, $8 }END{printf "</tbody>\n"}' ./outputs/best_motifs.txt > ./outputs/table.html
		mv tmp/*.dat outputs/
	else
		echo "The script has failed and no file has been created!"
		exit 1
	fi
	
# 	awk '(($10<=0.05 && $12<=0.05) || ($11<=0.05 && $13<=0.05))' motifs_last_loop.txt > ./outputs/best_motifs.txt
	
# 	awk 'BEGIN{printf "<tbody>\n"}{printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%.3f</td>\n\t\t<td>%.3f</td>\n\t</tr>\n", NR, $1, $2*100, $3*100, $4*100, $5*100, $6, $7, $10, $11}END{printf "</tbody>\n"}' ./outputs/best_motifs.txt > ./outputs/table.html
# 	awk 'BEGIN{printf "<tbody>\n"}{printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%.3f</td>\n\t\t<td>%.3f</td>\n\t\t<td>%.3f</td>\n\t\t<td>%.3f</td>\n\t</tr>\n", NR, $1, $2*100, $3*100, $4, $5, $10, $11, $12, $13}END{printf "</tbody>\n"}' ./outputs/best_motifs.txt > ./outputs/table.html
# 	mv motifs_last_loop.txt ./outputs/motifs_last_loop.txt
	
cd ..
