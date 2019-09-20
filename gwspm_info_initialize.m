% This function is part of the toolbox:
%       gwSPM: Graph-based, Wavelet-based Statistical Parametric Mapping
%       (v1.00)
%
% 	Author: Hamid Behjat
% 
%   Biomedical Signal Processing Group, 
%   Dept. of Biomedical Engineering,
%   Lund University, Sweden
% 
%   June 2016
%
function info = gwspm_info_initialize(Im)

info.contrastDirName = 'gwspm_import_contrasts';
info.templateDirName = 'gwspm_templates';
info.templateDirRoot =  fileparts(deblank(Im(1,:)));
info.templateName = 'w1Template_6.nii,1';
info.maskName = 'GM';
