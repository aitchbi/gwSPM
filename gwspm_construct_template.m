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
function [options,t] = gwspm_construct_template(images,options)

spmVer = spm('Ver','',1);
spmDir = spm('Dir');
subjOneAnatDir= options.templateDirRoot;

if strcmp(spmVer,'SPM12') == 0
    error('The current code is only compatible with SPM12.')
end

t = zeros(3,1);

%-Segmentation & Dartel import format file construction
%-----------------------------------------
clear job
job.channel.vols = images;
job.channel.biasreg = 0.001;
job.channel.biasfwhm = 60;
job.channel.write = [1 1];
job.tissue(1).tpm = {strcat(spmDir,filesep,'tpm',...
    filesep,'TPM.nii,1')};
job.tissue(1).ngaus = 1;
job.tissue(1).native = [1 1];
job.tissue(1).warped = [0 0];
job.tissue(2).tpm = {strcat(spmDir,filesep,'tpm',...
    filesep,'TPM.nii,2')};
job.tissue(2).ngaus = 1;
job.tissue(2).native = [1 1];
job.tissue(2).warped = [0 0];
job.tissue(3).tpm = {strcat(spmDir,filesep,'tpm',...
    filesep,'TPM.nii,3')};
job.tissue(3).ngaus = 2;
job.tissue(3).native = [1 1];
job.tissue(3).warped = [0 0];
job.tissue(4).tpm = {strcat(spmDir,filesep,'tpm',...
    filesep,'TPM.nii,4')};
job.tissue(4).ngaus = 3;
job.tissue(4).native = [1 1];
job.tissue(4).warped = [0 0];
job.tissue(5).tpm = {strcat(spmDir,filesep,'tpm',...
    filesep,'TPM.nii,5')};
job.tissue(5).ngaus = 4;
job.tissue(5).native = [1 1];
job.tissue(5).warped = [0 0];
job.tissue(6).tpm = {strcat(spmDir,filesep,'tpm',...
    filesep,'TPM.nii,6')};
job.tissue(6).ngaus = 2;
job.tissue(6).native = [0 0];
job.tissue(6).warped = [0 0];
job.warp.mrf = 1;
job.warp.cleanup = 1;
job.warp.reg = [0 0.001 0.5 0.05 0.2];
job.warp.affreg = 'mni';
job.warp.fwhm = 0;
job.warp.samp = 3;
job.warp.write = [1 1];
job.warp.vox = 1;

tic
spm_preproc_run(job,'run');
t(1) =toc;


%-Dartel template construction
%-----------------------------------------
images_rc1 = cell(size(images));
images_rc2 = cell(size(images));

for i=1:size(images,1)
    [pathstr,name,ext] = fileparts(images{i});
    images_rc1(i) = {strcat(pathstr,filesep,'rc1',name,ext)};
    images_rc2(i) = {strcat(pathstr,filesep,'rc2',name,ext)};
end

clear job
job.images = {images_rc1, images_rc2};
job.settings.template = 'Template';
job.settings.rform = 0;
job.settings.param(1).its = 3;
job.settings.param(1).rparam = [4 2 1e-06];
job.settings.param(1).K = 0;
job.settings.param(1).slam = 16;
job.settings.param(2).its = 3;
job.settings.param(2).rparam = [2 1 1e-06];
job.settings.param(2).K = 0;
job.settings.param(2).slam = 8;
job.settings.param(3).its = 3;
job.settings.param(3).rparam = [1 0.5 1e-06];
job.settings.param(3).K = 1;
job.settings.param(3).slam = 4;
job.settings.param(4).its = 3;
job.settings.param(4).rparam = [0.5 0.25 1e-06];
job.settings.param(4).K = 2;
job.settings.param(4).slam = 2;
job.settings.param(5).its = 3;
job.settings.param(5).rparam = [0.25 0.125 1e-06];
job.settings.param(5).K = 4;
job.settings.param(5).slam = 1;
job.settings.param(6).its = 3;
job.settings.param(6).rparam = [0.25 0.125 1e-06];
job.settings.param(6).K = 6;
job.settings.param(6).slam = 0.5;
job.settings.optim.lmreg = 0.01;
job.settings.optim.cyc = 3;
job.settings.optim.its = 3;

tic
spm_dartel_template(job);
t(2) = toc;


%-Move constructed templates
%--------------------------------------------
clear job
job.parent = {subjOneAnatDir};
job.name = options.templateDirName;
cfg_run_mkdir(job);

clear job
job.files = {
    strcat(subjOneAnatDir,filesep,'Template_0.nii')
    strcat(subjOneAnatDir,filesep,'Template_1.nii')
    strcat(subjOneAnatDir,filesep,'Template_2.nii')
    strcat(subjOneAnatDir,filesep,'Template_3.nii')
    strcat(subjOneAnatDir,filesep,'Template_4.nii')
    strcat(subjOneAnatDir,filesep,'Template_5.nii')
    strcat(subjOneAnatDir,filesep,'Template_6.nii')
    };
job.action.moveto = {strcat(subjOneAnatDir,...
    filesep,options.templateDirName)};
cfg_run_file_move(job);


%-Normalise: Estimate and Write
%--------------------------------------------
clear job
job.subj.vol = {strcat(subjOneAnatDir,filesep,...
    options.templateDirName,filesep,'Template_6.nii,1')};
job.subj.resample = job.subj.vol;
job.eoptions.biasreg = 0.0001;
job.eoptions.biasfwhm = 60;
job.eoptions.tpm = {strcat(spmDir,filesep,'tpm',filesep,'TPM.nii')};
job.eoptions.affreg = 'mni';
job.eoptions.reg = [0 0.001 0.5 0.05 0.2];
job.eoptions.fwhm = 0;
job.eoptions.samp = 3;
job.woptions.bb = [-78 -112 -70; 78 76 85];
job.woptions.vox = [1 1 1];
job.woptions.interp = 7;
job.woptions.prefix = 'w1';

tic
spm_run_norm(job);
t(3) =toc;

