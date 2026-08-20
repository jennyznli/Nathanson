
awk -v OFS='\t' 'NR==FNR{if(FNR==1){next}; newid[FNR-1]=$3; next} {$2=newid[FNR]; print} ../emilio-score.txt ../discovery/pennwtcga-66.bim > ../discovery/discovery_prs.bim
# copy others
