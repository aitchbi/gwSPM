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
%
function gwspm_postprocess_extractions(varargin)
data = varargin{1};

% Old Normalise: Estimate
%--------------------------------------------------------------------------
clear job
job.subj.source = {strcat(data.templateDirRoot,filesep,data.templateDirName,filesep,'w1Template_6.nii,1')};
job.subj.wtsrc = '';
job.subj.resample = {
    strcat(data.templateDirRoot,filesep,data.templateDirName,filesep,'GM_cerebellum.nii')
    strcat(data.templateDirRoot,filesep,data.templateDirName,filesep,'GM_cerebrum.nii')
    };
job.eoptions.template = {strcat(spm('Dir'),filesep,'tpm',filesep,'TPM.nii,1')};
job.eoptions.weight = '';
job.eoptions.smosrc = 8;
job.eoptions.smoref = 0;
job.eoptions.regtype = 'mni';
job.eoptions.cutoff = 25;
job.eoptions.nits = 16;
job.eoptions.reg = 1;
job.roptions.preserve = 0;
job.roptions.bb = [-78 -112 -70;78 76 85];
job.roptions.interp = 7;
job.roptions.wrap = [0 0 0];

for voxDim = [1, 2, 2.5, 3]
    job.roptions.vox = [voxDim, voxDim, voxDim];
    switch voxDim
        case 1, job.roptions.prefix = 'w1_';
        case 2, job.roptions.prefix = 'w2_';
        case 2.5, job.roptions.prefix = 'w2point5_';
        case 3, job.roptions.prefix= 'w3_';
    end
    
    spm_run_normalise(job);
    
    % Making sure we have a single connected cerebrum and a single connected cerebellum
    for region = [{'cerebrum'}, {'cerebellum'}]
        
        vol_info  = spm_vol(strcat(data.templateDirRoot,filesep,...
            data.templateDirName,filesep,job.roptions.prefix,...
            'GM_',region{:},'.nii'));
  
        gm = spm_read_vols(vol_info);
        
        ind = find(gm);
        numClean  = round(0.7*numel(ind)/2);
        gm_adj = gwspm_adjust_mask(gm,numClean,6,region{:});
        
        [pathstr,name,ext] = fileparts(vol_info.fname);
        vol_info.fname = strcat(pathstr,filesep,name,'_adjusted',ext);
        spm_write_vol(vol_info,gm_adj);
        
        % Threshold and binarize (to create downsampled mask)
        mask_adj = gm_adj;
        mask_adj(mask_adj<=0.5) = 0;
        
        numel(find(isnan(mask_adj))) 
        mask_adj(isnan(mask_adj)) = 0;
        
        indMask = find(mask_adj);
        mask_adj(indMask) = 1;

        vol_info.fname = strcat(pathstr,filesep,name,'_mask_adjusted',ext);
        spm_write_vol(vol_info,mask_adj);
        
    end
end

