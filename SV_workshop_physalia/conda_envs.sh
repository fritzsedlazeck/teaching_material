# Conda installation
curl -O  https://repo.anaconda.com/miniconda/Miniconda3-py312_25.9.1-3-Linux-x86_64.sh

BASE_PATH=""  # <== this need to be added, where are you installing conda
bash Miniconda3-py312_25.9.1-3-Linux-x86_64.sh -b -p ${BASE_PATH}/miniconda3
# uncomment the line below to auto initiate conda every time
# source ${BASE_PATH}/miniconda3/etc/profile.d/conda.sh

# ENV 1
conda create --name sr samtools=1.19 bcftools=1.19 samplot=1.3 manta=1.6 mummer=3.23 bwa-mem2=2.2.1 survivor=1.0.7 -vv

# ENV 2
conda create --name lr python=3.12 sniffles=2.7.1 samtools=1.21 bcftools=1.21 dipcall=0.3 samplot=1.3 minimap2=2.30 annotsv=3.5 -vv

# ENV 3
conda create --name svafotate python=3.8 -vv
conda activate svafotate
conda install --file https://raw.githubusercontent.com/fakedrtom/SVAFotate/master/requirements.txt -vv
pip install git+https://github.com/fakedrtom/SVAFotate.git

# for AnnotSV
cd ${BASE_PATH}/miniconda3/envs/SVW_lr/share/AnnotSV
curl -O https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.4.6.tar.gz
tar -xzf Annotations_Human_3.4.6.tar.gz
