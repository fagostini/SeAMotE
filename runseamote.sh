#!/bin/bash
# bash runseamote.sh "random_number" "motiFile" "use_as_reference" "threshold" "reverse"
echo "bash runseamote.sh" "$1" "$2" "$3" "$4" "$5"

cd   tmp/$1
	# Set the type of reference to be used
	if [[ "$3" == "negative" && -s "negative.oneline" ]]; then
			options="-n negative.oneline "
	else
			options="-r $3 "
	fi
	# Set the threshold and the molecule type
	options+="-t $4 -a $5 -w 7"

	# Run the C script with the given parameters
	./motifSearch -p positive.oneline $options

	mv positive.oneline outputs/positive.seq
	if [[ -s "negative.oneline" ]]; then
		mv negative.oneline outputs/reference.seq
	else
		mv reference.seq outputs/reference.seq
	fi

	# Generate the output files and html tables
	if [[ -s "tmp/best_motifs.dat" ]]; then
		awk '{print $1}' tmp/best_motifs.dat > iupac.txt
		sed 's/R/[AG]/g;s/Y/[CT]/g;s/M/[CA]/g;s/K/[TG]/g;s/W/[TA]/g;s/S/[CG]/g;s/B/[CTG]/g;s/D/[ATG]/g;s/H/[ATC]/g;s/V/[ACG]/g' tmp/best_motifs.dat | awk '{print $1, $2, $4, $2/($2+$3), $4/($4+$5), $6, $NF}' > best_motifs.txt
		paste iupac.txt best_motifs.txt > outputs/best_motifs.txt

		# Generate the reverse complement for the DNA
		mv outputs/positive.seq outputs/positive.tmp
		if [[ "$5" == "dna" ]]; then
			awk '{gsub(/U/, "T", $2); print }' outputs/positive.tmp | awk 'BEGIN{ j=n=split("A C G T", t); for(i=0;++i<=n;){ map[t[i]] = t[j--] }} { printf "%s ", $0; for(i=length($2);i>0;i--){ printf "%s", map[substr($2, i, 1)]}; print x }' > outputs/positive.seq
		else
			awk '{gsub(/U/, "T", $2); print }' outputs/positive.tmp > outputs/positive.seq
		fi

		Rscript motif_filter_disc.R --vanilla --slave --args outputs/best_motifs.txt

		# Generate the logo and other files for each motif and the collected zip file
		Rscript createLogos.R --vanilla --slave 2&>R_log
# 		for((l=1;l<=`wc -l outputs/best_motifs.txt | awk '{print $1}'`;l++)); do
# 			IUmot=`awk '(NR=='$l'){print $1}' outputs/best_motifs.txt`;
# 			REmot=`awk '(NR=='$l'){print $2}' outputs/best_motifs.txt`;
# 			grep -o -e "$REmot" outputs/positive.seq > matches.txt
# 			python createLogo.py "$IUmot" 2&> python.log
# 		done
		zip -r outputs/results.zip logos/*
		rm matches.txt outputs/positive.tmp

		# Write the html files
# 		awk 'BEGIN{printf "<tbody>\n"}{printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td>%s</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%.3E</td>\n\t</tr>\n", NR, $1, $2, $3, $4, $5*100, $6*100, $7*100, $8 }END{printf "</tbody>\n"}' ./outputs/best_motifs.txt > ./outputs/table.html
# 		awk 'BEGIN{printf "<tbody>\n"}{printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td>%s</td>\n\t\t<td><img src=\"logos/%s_logo.png\" width=\"%d\" height=\"50\" style=\"solid white; display:block; margin-top:-1px; margin-left:1px; margin-right:-6px; margin-bottom:-1px;\"/></td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%.3E</td>\n\t</tr>\n", NR, $1, $2, $1, length($1)*30, $5*100, $6*100, $7*100, $8 }END{printf "</tbody>\n"}' ./outputs/best_motifs.txt > ./outputs/table.html
 		awk 'BEGIN{printf "<tbody>\n"}{if($8>0){printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td class=\"regex\">%s</td>\n\t\t<td class=\"logo\"><a href=\"logos/%s_logo.png\"><img src=\"logos/%s_logo.png\" width=\"%d\" /></a></td>\n\t\t<td class=\"pwm\"><a href=\"logos/%s_pwm.txt\"><img src=\"download.gif\" /></a></td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%.3E</td>\n\t</tr>\n", NR, $1, $2, $1, $1, length($1)*25, $1, $5*100, $6*100, $7*100, $8} else{printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td class=\"regex\">%s</td>\n\t\t<td class=\"logo\"><a href=\"logos/%s_logo.png\"><img src=\"logos/%s_logo.png\" width=\"%d\" /></a></td>\n\t\t<td class=\"pwm\"><a href=\"logos/%s_pwm.txt\"><img src=\"download.gif\" /></a></td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td> \< 5.0E-324</td>\n\t</tr>\n", NR, $1, $2, $1, $1, length($1)*25, $1, $5*100, $6*100, $7*100} }END{printf "</tbody>\n"}' ./outputs/best_motifs.txt > ./outputs/table.html
 		awk 'BEGIN{printf "<tbody>\n"}{printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td class=\"regex\">%s</td>\n\t\t<td class=\"logo\"><a href=\"logos/%s_logo.png\"><img src=\"logos/%s_logo.png\" width=\"%d\" /></a></td>\n\t\t<td class=\"pwm\"><a href=\"logos/%s_pwm.txt\">Download</a></td>\n\t\t<td><a href=\"logos/%s_transfac.txt\">Download</a></td>\n\t\t<td>%d</td>\n\t\t<td>%.3E</td>\n\t</tr>\n", NR, $1, $2, $1, $1, length($1)*25, $1, $1, $7*100, $8 }END{printf "</tbody>\n"}' ./outputs/best_motifs.txt > ./outputs/logos.html
		awk 'BEGIN{printf "#  IUPAC\t\t\t\t\tRegEX\tPos Count\tRef Count\tPos Coverage\tRef Coverage\tDiscrimination\tP-value\n"} {printf "%8s\t%37s\t%d\t\t%d\t\t%.5f\t\t%.5f\t\t%.5f\t\t%e\n", $1, $2, $3, $4, $5, $6, $7, $8}' ./outputs/best_motifs.txt > ./outputs/motifs.txt
		sed 's/0.000000e+00/< 5.0000e-324/g' ./outputs/motifs.txt > ./outputs/best_motifs.txt
		mv tmp/*.dat outputs/
	else
		echo "The script has failed and no file has been created!"
		exit 1
	fi

cd ..
