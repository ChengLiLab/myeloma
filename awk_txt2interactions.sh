
awk 'BEGIN{FS=","} $3==7 && $4>=21928240 && $4<=21948240 {print $3,$4-1,$4"\t"$6,$7-1,$7"\t"1}' cis.filter.txt | head > /lustre/user/liclab/wupz/dosageEffect/hic_data/RPMI8226_HindIII/bed_file/U266.HindIII.rs4487645.10k.cis.filter.txt