
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

# tried with regex to select only patterns at the start of the line but 
# it used to0 much ram. this g3d_swissprot_only will have more entries than swissprot
zgrep G3DSA protein2ipr.dat.gz | grep -f ../uniprot/acs.txt > g3d_swissprot_only.txt


# downloading uniprot data for entries associated with pdb")
#mkdir -p $SCRIPT_PATH/../data/pdb
#cd  $SCRIPT_PATH/../data/pdb
#wget -c "https://rest.uniprot.org/uniprotkb/stream?query=((structure_3d:true)+AND+(reviewed:false))&format=txt&compressed=true" > test.dat.gz
