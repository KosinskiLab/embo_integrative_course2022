# Anaconda 3
```
cd /beegfs/desy/group/school/software/tmp
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
#Enter /beegfs/desy/group/school/software/anaconda3 as installation location
#Do you wish the installer to initialize Anaconda3
#by running conda init? [yes|no]
#-> no !!!!!!
conda create --name ENV3 python=3.8
conda init bash #necessary?
conda activate ENV3
#Your shell prompt should display sth like: (ENV3) shool...
conda install -c salilab imp
#Test:
python -c "import IMP"
conda install -c salilab modeller
#Edit /beegfs/desy/group/school/software/anaconda3/envs/ENV3/lib/modeller-9.22/modlib/modeller/config.py
and replace XXXX with your Modeller license key (MODELIRANJE)
#test
python -c "import modeller"
conda install -c salilab imp-bayesianem
conda install scikit-learn
#Test:
python -c "from sklearn import datasets"
conda install scipy
conda install numpy
#Test:
python -c "import numpy"
python -c "import scipy"
conda install matplotlib
#Test
python -c "import matplotlib"
conda install -c conda-forge jupyterlab
#to test:
jupyter notebook
#and follow the displayed instructions
pip install calysto_bash
#to test:
jupyter notebook
# Download https://raw.githubusercontent.com/salilab/foxs_tutorial/master/foxs/nup133/FoXS.ipynb
# Open FoXS.ipynb
conda install -c salilab pyrmsd
conda install pandas msgpack-python
pip install pdb-tools
conda deactivate

#Add 
export PATH="/beegfs/desy/group/school/software/anaconda3/bin:$PATH"
to .bashrc for every user
```
