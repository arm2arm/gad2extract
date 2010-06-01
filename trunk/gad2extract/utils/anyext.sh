#!/bin/bash
#example:
#sh ~/bin/any.sh "snap_gal_???" 2 40 | sh
mask=$1
snap=`ls -t ${mask}`
 
for k in $snap 
do 

# cc=`printf "${mask_format}" $c`
#echo ${k%_*}_${k##*_} 
#echo ${k%_*}0${k##*_} 
#echo  mv  $k ${k%_*}_0${k##*_}
FILE="extract_"${k##*_}".ini"
FOUT="4096_R300__2_dump_${k##*_}"
if [ -f $FOUT ];
then

echo  "file exist: ${FILE}"

else

echo "
[EXTRACTREGION]
# Radius in kpc
R = 300
#file to split 
FILE =${k}
#file format to bumped 
DUMPFILE = 4096_R%d__%d_dump_%s 
#ID of particles to be traced 
# the format is following :
# nhalo 
# np IDs  ...
# np IDs  ...
IDFILE = /home/hlrb2/h009z/h009zac/tmp/id.idx
" > ${FILE}

~/bin/gad2extract.x -p /ptmp1/h009z/h009zac/BOX64/DUMPS/ --inifile=${FILE}
fi

done
