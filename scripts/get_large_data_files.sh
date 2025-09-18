
SCRIPT_PATH=$(dirname $(realpath -s $0))
echo $SCRIPT_PATH

mkdir -p $SCRIPT_PATH/../data/uniprot
cd  $SCRIPT_PATH/../data/uniprot

#echo "Downloading uniprot and swissprot files"
lftp https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/ <<EOF
pget -c uniprot_sprot.dat.gz
#pget -c uniprot_trembl.dat.gz
bye
EOF

mkdir -p $SCRIPT_PATH/../data/go
cd  $SCRIPT_PATH/../data/go

echo "Downloading go annotations file"
lftp https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/ <<EOF
pget -c goa_uniprot_all.gpa.gz
bye
EOF

echo "Preparing ingestion file"
zcat goa_uniprot_all.gpa.gz | tail -n +2 | cut -f2,3,4,6 | tr '\t' ','| sed 's/GO://g' | sed 's/ECO://g' | grep -v PRO_ | grep -v "!" | sed -E 's/-[0-9]+,/,/g' | sort | uniq | gzip > for_ingestion.csv.gz

#echo "Downloading interpro file"
#mkdir -p $SCRIPT_PATH/../data/interpro
#cd  $SCRIPT_PATH/../data/interpro
#lftp https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/ <<EOF
#pget -c protein2ipr.dat.gz
#bye
#EOF

#echo "filtering ipr file for gene3d annotations"
#zgrep G3DSA protein2ipr.dat.gz > g3d_trembl.txt
