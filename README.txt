gwSPM : graph-based, wavelet-based Statistical Parametric Mapping Toolbox
(c) Hamid Behjat, 2016

contact: hamid.behjat@bme.lth.se

This toolbox contains code implementing the framework presented in the 
following paper: 

[1] H. Behjat, N. Leonardi, L. Sornmo, D. Van De Ville, 
"Anatomically-adapted graph wavelets for improved group-level fMRI 
activation mapping", NeuroImage, 123, pp. 185-199, 2015.

The implementation is developed as a toolbox for the SPM software package. 
As such, the prerequisite for using this toolbox is to have SPM installed. 
In particular, only SPM12 is supported in the present release. SPM8 
compatibility may also be considered in future updates, especially, 
upon request. 
 
Installation of the toolbox is very simple. Place the whole folder, i.e. 
'gwSPM', inside the 'toolbox' directory of spm12. If you now start SPM:

>> spm fmri

the toolbox should automatically be recognized by SPM and show up in the
'Toolboxes' drop down menu in the main SPM GUI. Click on gwSPM and off you 
go. 


I would be more than happy to answer questions, and to help out if you have 
any problems with analyzing your data. Also, make sure you keep an eye out 
for future demos, extensions and upgrades.

Good Luck!     

Hamid Behjat
June 2016


Download: 
The toolbox can be downloaded at https://github.com/hbehjat/gwSPM


License: 
The gwSPM toolbox is a Matlab library released under the GPL.

The gwSPM toolbox is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The gwSPM toolbox is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the gwSPM toolbox.  If not, see <http://www.gnu.org/licenses/>.


