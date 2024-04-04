
basedir=$(dirname $0)
mkdir -p $basedir/../data/uniprot
cd  $basedir/../data/uniprot

wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt

wget -c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz

