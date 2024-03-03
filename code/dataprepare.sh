tar -xvf GSE218300_RAW.tar
for file in *.gz; do
    gunzip "$file"
done
