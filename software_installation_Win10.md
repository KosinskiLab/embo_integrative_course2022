```
#Create directory C:\Users\schulu16\Software (replace schulu16 with the proper username)
#Add Software and Desktop to Quick Lunch toolbar


#Pymol
#FastX3
#Sublime Text
#WinSCP
1. Copy the folders from INSTALLATION\Software-just_copy_content_to_software_dir to C:\Users\schulu16\Software (replace schulu16 with the proper username)
2. Make shortcuts on the Desktop
3. Configure connection for WinSCP

Host name: schulu1.desy.de
User name: Linux username
Password: Linux password
(save password and Save session, create shortcut)

Host name: max-display004.desy.de or max-display005.desy.de depending on the user
User name: Linux username
Password: Linux password
(save password and Save session, create shortcut)


4. Configure connection for FastX2
name and host: schulu1.desy.de
User: Linux user
Test the connection
5. Configure connection for FastX3
name and host: max-display004.desy.de or max-display005.desy.de depending on the user
User: Linux user
Test the connection

#VMD
#Execute and install to C:\Users\schulu16\Software\VMD
#Create a shortcut on the Desktop

#ChimeraX
#Install using ChimeraX-0.9.exe from INSTALLATION/ folder
#Installation path C:\Users\schulu16\Software\ChimeraX
#Create a shortcut on the Desktop

#UCSF Chimera 1.13
#Execute and install to C:\Users\schulu16\Software\Chimera 1.13.1
#Click yes to create a shortcut on the Desktop

#Install Anaconda3 by executing Anaconda3-2019.10-Windows-x86_64.exe from INSTALLATION/ folder
#Start Windows Command Line cmd.exe (doesn't work with PowerShell)
C:
cd C:\Users\schulu16\Software
conda.exe create python=3.7 --prefix ENV3
activate C:\Users\schulu16\Software\ENV3
#Your shell prompt should display sth like: (ENV3) shool...
conda install -c salilab modeller
#Edit C:\Users\schulu16\Software\ENV3\Library\modeller\modlib\modeller\config.py
and replace XXXX with your Modeller license key (MODELIRANJE)
#test
python -c "import modeller"


### DO NOT DO FOR NOW ###
conda install -c salilab imp
#Test:
python -c "import IMP"
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

```
