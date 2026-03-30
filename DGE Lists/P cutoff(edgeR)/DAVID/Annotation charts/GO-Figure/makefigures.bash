for FILE in *.tsv
do
	echo Processing $FILE
	gofigure.py -i $FILE -o ${FILE%.*}-GoPlots -j standard-plus -v 0.2 -a 25 -e 100 -m 25 -l 'description'
	echo Completed plotting $FILE
done
