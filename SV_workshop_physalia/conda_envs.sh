# Conda installation
curl -O  https://repo.anaconda.com/miniconda/Miniconda3-py314_26.9.4-1-Linux-x86_64.sh

BASE_PATH=""  # <== this need to be added, where are you installing conda
bash Miniconda3-py314_26.9.4-1-Linux-x86_64.sh -b -p ${BASE_PATH}/miniconda3
# uncomment the line below to auto initiate conda every time
# source ${BASE_PATH}/miniconda3/etc/profile.d/conda.sh
# source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
# conda activate base

# ENV 1: manta uses python2
conda create --name manta samtools bcftools manta=1.6 libboost=1.85.0 -vv
# conda install -c conda-forge --force-reinstall gsl=2.6

# ENV 2: mostr of the tools needed for the workshop
conda create --name svw python=3.12 sniffles=2.8.1 mummer=3.23 bwa-mem2=2.3 survivor=1.0.7 samtools=1.24 bcftools=1.24 dipcall=0.3 samplot=1.3 minimap2=2.31 annotsv=3.5.10 assemblytics=2.0.1 gsl=2.6 openblas libblas -vv
# conda install -c conda-forge --force-reinstall gsl=2.6

# ENV 3
conda create --name svafotate python=3.8 -vv
conda activate svafotate
wget https://raw.githubusercontent.com/fakedrtom/SVAFotate/master/requirements.txt
conda install --file requirements.txt -vv
pip install git+https://github.com/fakedrtom/SVAFotate.git

# for AnnotSV
# conda activate svw
cd miniconda3/envs/svw/share/AnnotSV
curl -O https://cstb-icube.fr/~geoffroy/Annotations/Annotations_Human_3.5.tar.gz
tar -xzf Annotations_Human_3.4.6.tar.gz
