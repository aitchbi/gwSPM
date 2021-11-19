# gwSPM

This toolbox contains code implementing the method presented in the following [paper](https://bme.lth.se/fileadmin/biomedicalengineering/Personal_folders/Hamid_Behjat/HBehjat_NeuroImage2015.pdf): 

> H. Behjat, N. Leonardi, L. Sornmo, D. Van De Ville, "Anatomically-adapted graph wavelets for improved group-level fMRI activation mapping", NeuroImage, 123, pp. 185-199, 2015.

An overview of the toolbox is given in this [poster](https://bme.lth.se/fileadmin/biomedicalengineering/Personal_folders/Hamid_Behjat/HBehjat_OHBM2018a-poster.pdf), also shown below:

![gwSPM toolbox overview poster](figs/HBehjat_OHBM2018_poster.jpg?raw=true) 
 
The implementation is developed as a toolbox for the [SPM](https://www.fil.ion.ucl.ac.uk/spm/) software package, which is a free open-source software package. As such, the prerequisite for using this toolbox is to have SPM installed. In particular, only [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) is supported in the present release.  

### Usage 
Download the entire repository, place it in a folder named **gwSPM**, and then place this folder inside the **toolbox** directory of your spm12 package. In MATLAB, initiate SPM by typing `spm fmri` in the command window. The gwSPM toolbox should automatically be recognized by SPM and should show up in the **Toolboxes** drop down menu in the main SPM GUI. Click on gwSPM and off you go. 
