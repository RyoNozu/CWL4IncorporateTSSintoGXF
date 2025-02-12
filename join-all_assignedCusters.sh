#! /bin/zsh

# This script is used to join the .assignedClusters files, key columns are each first column of the files
# each file has header 
# specified header option 
# each file has 13 columns
# each file included 13 fields, the first field is the key column (number)
# Empty fields are filled with NA


#ova=ovary.assignedClusters.txt
#tes=testis.assignedClusters.txt
#ova_tes=joined-ovary-testis.assignedClusters.txt
#Ftissue=Ftissue.assignedClusters.txt
#ova_tes_Ftissue=joined-ovary-testis-Ftissue.assignedClusters.txt
#Mtissue=Mtissue.assignedClusters.txt
#joined=all-joined.assignedClusters.tsv

#sort $ova >sort.$ova
#sort $tes >sort.$tes
#sort $Ftissue >sort.$Ftissue
#sort $Mtissue >sort.$Mtissue
#
#join -t $'\t' --header --check-order -j 1 -a 1 -a 2 -e NA sort.$ova sort.$tes > ${ova_tes}
#join -t $'\t' --header --check-order -j 1 -a 1 -a 2 -e NA -o auto ${ova_tes} sort.$Ftissue > ${ova_tes_Ftissue}
#join -t $'\t' --header --check-order -j 1 -a 1 -a 2 -e NA -o auto ${ova_tes_Ftissue} sort.$Mtissue > $joined
#
#sort -n -k1,1 $joined >sort.$joined

d1=1day.assignedClusters.txt
d2=2days.assignedClusters.txt
joined=all-joined.assignedClusters.tsv

sort $d1 >sort.$d1
sort $d2 >sort.$d2

join -t $'\t' --header --check-order -j 1 -a 1 -a 2 -e NA -o auto sort.$d1 sort.$d2 > $joined
sort -n -k1,1 $joined >sort.$joined