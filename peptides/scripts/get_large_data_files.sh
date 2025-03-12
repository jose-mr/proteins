
SCRIPT_PATH=$(dirname $(realpath -s $0))
echo $SCRIPT_PATH

mkdir -p $SCRIPT_PATH/../../data/peptides/uniprot
cd  $SCRIPT_PATH/../../data/peptides/uniprot

# download peptide information
echo "downloading peptide information from Uniprot"
wget -O uniprot_peptides.dat.gz "https://rest.uniprot.org/uniprotkb/stream?format=txt&compressed=true&query=%28length%3A%5B10+TO+20%5D%29"

# create a file with all the swissprot acs to filter the big dat files
echo "create files with all peptide accession codes in uniprot"
zgrep "^AC" uniprot_peptides.dat | awk -F "   " '{print $2}' | tr -d ";" | tr " " "\n" > peptides_acs.txt

# filger go annotations
mkdir -p $SCRIPT_PATH/../../data/peptides/go
cd  $SCRIPT_PATH/../../data/peptides/go
echo "filtering go annotations file"
zgrep ../../go/goa_uniprot_all.gpa.gz -f ../uniprot/peptides_acs.txt > goa_uniprot_peptides_filter.txt

# filter interpro annotations

mkdir -p $SCRIPT_PATH/../../data/peptides/interpro
cd  $SCRIPT_PATH/../../data/peptides/interpro

# tried with regex to select only patterns at the start of the line but 
# it used to0 much ram. this g3d_swissprot_only will have more entries than swissprot
echo "filtering ipr file"
zgrep G3DSA ../../interpro/protein2ipr.dat.gz | grep -f ../uniprot/peptides_acs.txt > g3d_peptides_only.txt

