VERSION="1.22.0"
REFERENCE_FOLDER="igblast_databases_nov24"

cd ~/my_software/

wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/${VERSION}/ncbi-igblast-${VERSION}-x64-linux.tar.gz

tar -zxf ncbi-igblast-${VERSION}-x64-linux.tar.gz

# Download reference databases and setup IGDATA directory

module load HGI/softpack/users/cs54/immcantation_igblast/1

fetch_igblastdb.sh -o ~/${REFERENCE_FOLDER}/igblast
cp -r ~/my_software/ncbi-igblast-${VERSION}/internal_data ~/${REFERENCE_FOLDER}/igblast
cp -r ~/my_software/ncbi-igblast-${VERSION}/optional_file ~/${REFERENCE_FOLDER}/igblast

# Build IgBLAST database from IMGT reference sequences

fetch_imgtdb.sh -o ~/${REFERENCE_FOLDER}/germlines/imgt
imgt2igblast.sh -i ~/${REFERENCE_FOLDER}/germlines/imgt -o ~/${REFERENCE_FOLDER}/igblast #must have igblast in your path for makeblastdb




# imgt2igblast.sh -i ~/igblast_databases_nov24/germlines/imgt -o ~/igblast_databases_nov24/igblast #must have igblast in your path for makeblastdb
