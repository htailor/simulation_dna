#!/bin/sh

FILE_TYPE_LIST="*.log *.LAMBDA *.ROE *.T *.T00 *.T00A *.T11 *.data core.* *.status *.Output *.out Parameters *.*~"


echo "\n\nClearing Raw Data...";
sleep 2
for FILE_TYPE in $FILE_TYPE_LIST
do
	echo "\nSearching for ${FILE_TYPE}...found `find . \( ! -regex '.*/\..*' \) -name ${FILE_TYPE} | wc -l` results";
   echo "Deleting...`find . \( ! -regex '.*/\..*' \) -name ${FILE_TYPE} -exec rm -rf {} \;`";
done

echo "\n"

