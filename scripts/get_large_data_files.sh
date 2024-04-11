
SCRIPT_PATH=$(dirname $(realpath -s $0))
echo $SCRIPT_PATH

mkdir -p $SCRIPT_PATH/../data/uniprot
cd  $SCRIPT_PATH/../data/uniprot
wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt
wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz

# create a file with all the swissprot acs to filter the interpro file
zgrep "^AC" uniprot_sprot.dat | awk -F "   " '{print $2}' | tr -d ";" | tr " " "\n" > acs.txt

mkdir -p $SCRIPT_PATH/../data/go
cd  $SCRIPT_PATH/../data/go
wget -c ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpa.gz

mkdir -p $SCRIPT_PATH/../data/interpro
cd  $SCRIPT_PATH/../data/interpro
wget -c https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz
# create file with swissprot and gene 3d annotations only
zgrep -f ../uniprot/acs.txt protein2ipr.dat.gz  | grep G3DSA > g3d_swissprot_only.txt
