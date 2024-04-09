
# TODO mv this to python and fix paths
basedir=$(dirname $0)
#mkdir -p $basedir/../data/uniprot
cd  $basedir/../data/uniprot

wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt

wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz

cd  $basedir/../data/go

wget -c ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpa.gz
