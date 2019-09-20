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
function  gwspm_preprocess_contrasts(N_subjects,images,contrasts,meanVols,options)

%-Preliminaries 
%--------------------------------------------------------------------------

% Make gwspm directory for contrasts
clear job
for iS =1:N_subjects
    dummy = contrasts{iS}(1);
    pathstr = fileparts(dummy{:});

    job.parent = {pathstr};
    job.name = options.contrastDirName ;
    cfg_run_mkdir(job);
end

% File names of 'coregistered' & 'normalised' contrasts
c_contrasts = cell(size(contrasts));
w_contrasts = cell(size(contrasts));
for iS =1:N_subjects
    temp1 = cell(size(contrasts{iS}));
    temp2 = cell(size(contrasts{iS}));
    for i=1:size(contrasts{iS},1)
        [pathstr,name,ext] = fileparts(contrasts{iS}{i});
        temp1(i) = {strcat(pathstr,filesep,'c_',name,ext)};
        temp2(i) = {strcat(pathstr,filesep,'wc_',name,ext)};
    end
    c_contrasts(iS) = {temp1};
    w_contrasts(iS) = {temp2};
end


%-Coregister contrasts with strutural data (using the mean functional volume)
%--------------------------------------------------------------------------
% Even if the input data have been coregistered, its good to make sure. 
% Also a change of naming is required. 

clear job
job.eoptions.cost_fun = 'nmi';
job.eoptions.sep = [4 2];
job.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 ...
    0.001 0.01 0.01 0.01 0.001 0.001 0.001];
job.eoptions.fwhm = [7 7];
job.roptions.interp = 4;
job.roptions.wrap = [0 0 0];
job.roptions.mask = 0;
job.roptions.prefix = 'c_';

for iS =1:N_subjects
    job.ref = images(iS);
    job.source = meanVols(iS);
    job.other = contrasts{iS};
    spm_run_coreg(job);
end


%-DARTEL Normalise to MNI Space
%--------------------------------------------------------------------------
clear job
job.template =  {strcat(options.templateDirRoot,...
    filesep,options.templateDirName,filesep,...
    'Template_6.nii')};

for iS=1:N_subjects
    [pathstr,name,ext] = fileparts(images{iS});
    
    job.data.subj(iS).flowfield = {strcat(pathstr,...
        filesep,'u_rc1',name,'_Template',ext)};
    
    job.data.subj(iS).images = c_contrasts{iS}; 
end

job.bb = [-78 -112 -70; 78 76 85];
job.preserve = 0;
job.fwhm = [0 0 0];

for voxDim = [2, 2.5, 3]
    job.vox = [voxDim, voxDim, voxDim];
    
    spm_dartel_norm_fun(job); 
    
    for iS=1:N_subjects
       
        clear job2
        job2.files = w_contrasts{iS};
        pathstr = fileparts(contrasts{iS}{1});
        job2.action.moveren.moveto = {strcat(...
            pathstr,filesep,options.contrastDirName)};
        job2.action.moveren.patrep.pattern = 'w';
        switch voxDim
            case 2
                job2.action.moveren.patrep.repl = 'w2';
            case 2.5
                job2.action.moveren.patrep.repl = 'w2point5';
            case 3
                job2.action.moveren.patrep.repl = 'w3';
        end
        job2.action.moveren.unique = false;
        
        cfg_run_file_move(job2);
        
    end
    
end





