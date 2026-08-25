
bsub -J prs_gwas2021_z \
     -o /project/knathans_tecac/jenny/prs/logs/gwas2021_z.%J.o \
     -e /project/knathans_tecac/jenny/prs/logs/gwas2021_z.%J.e \
     bash prs.sh gwas2021_z /home/jennyzli/jenny/prs/scorefiles/gwas2021_hits_z.txt

bsub -J prs_gwas2021_effect \
     -o /project/knathans_tecac/jenny/prs/logs/gwas2021_effect.%J.o \
     -e /project/knathans_tecac/jenny/prs/logs/gwas2021_effect.%J.e \
     bash prs.sh gwas2021_effect /home/jennyzli/jenny/prs/scorefiles/gwas2021_hits_effect.txt

bsub -J prs_gwas2026_z \
     -o /project/knathans_tecac/jenny/prs/logs/gwas2026_z.%J.o \
     -e /project/knathans_tecac/jenny/prs/logs/gwas2026_z.%J.e \
     bash prs.sh gwas2026_z /home/jennyzli/jenny/prs/scorefiles/gwas2026_hits_z.txt

bsub -J prs_gwas2026_effect \
     -o /project/knathans_tecac/jenny/prs/logs/gwas2026_effect.%J.o \
     -e /project/knathans_tecac/jenny/prs/logs/gwas2026_effect.%J.e \
     bash prs.sh gwas2026_effect /home/jennyzli/jenny/prs/scorefiles/gwas2026_hits_effect.txt

