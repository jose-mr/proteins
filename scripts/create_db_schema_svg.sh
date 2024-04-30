
SCRIPT_PATH=$(dirname $(realpath -s $0))
OUT_PATH=$SCRIPT_PATH/../out/

mkdir -p $OUT_PATH


python manage.py graph_models cath ec eco go uniprot -g -o $OUT_PATH/db_schema.svg
