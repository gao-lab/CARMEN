rsync -aP rsync://hgdownload.cse.ucsc.edu/gbdb/hg19/multiz46way/*.wib .
rsync -aP rsync://hgdownload.cse.ucsc.edu/gbdb/hg19/multiz100way/*.wib .

wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/phastCons100way.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/phyloP100wayAll.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/phyloP46wayAll.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/phyloP46wayPlacental.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/phyloP46wayPrimates.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/phastCons46way.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/phastCons46wayPlacental.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/phastCons46wayPrimates.txt.gz

gunzip *.txt.gz
ln -s *.txt *.wig
for file in `ls *.txt | cut -d '.' -f 1`
do
ln -s ${file}.txt ${file}.wig
done
