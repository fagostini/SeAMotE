#!/bin/bash

cd   tmp/$1

	./main -p positive.oneline -n negative.oneline -t 0.95
	
	awk '(($10<=0.05 && $12<=0.05) || ($11<=0.05 && $13<=0.05))' motifs_last_loop.txt > ./outputs/best_motifs.txt
	
# 	awk 'BEGIN{printf "<tbody>\n"}{printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%.3f</td>\n\t\t<td>%.3f</td>\n\t</tr>\n", NR, $1, $2*100, $3*100, $4*100, $5*100, $6, $7, $10, $11}END{printf "</tbody>\n"}' ./outputs/best_motifs.txt > ./outputs/table.html
	awk 'BEGIN{printf "<tbody>\n"}{printf "\t<tr>\n\t\t<td>%d</td>\n\t\t<td>%s</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%d</td>\n\t\t<td>%.3f</td>\n\t\t<td>%.3f</td>\n\t\t<td>%.3f</td>\n\t\t<td>%.3f</td>\n\t</tr>\n", NR, $1, $2*100, $3*100, $4, $5, $10, $11, $12, $13}END{printf "</tbody>\n"}' ./outputs/best_motifs.txt > ./outputs/table.html
	mv motifs_last_loop.txt ./outputs/motifs_last_loop.txt
	
cd ..
