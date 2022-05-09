function gwspm_postprocess_extractions(varargin)
data = varargin{1};

% Old Normalise: Estimate--------------------------------------------------
clear job
job.subj.source = {fullfile(data.templateDirRoot,data.templateDirName,'w1Template_6.nii,1')};
job.subj.wtsrc = '';
job.subj.resample = {
    fullfile(data.templateDirRoot,data.templateDirName,'GM_cerebellum.nii')
    fullfile(data.templateDirRoot,data.templateDirName,'GM_cerebrum.nii')
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
    job.roptions.vox = [voxDim,voxDim,voxDim];
    switch voxDim
        case 1
            job.roptions.prefix = 'w1_';
        case 2
            job.roptions.prefix = 'w2_';
        case 2.5
            job.roptions.prefix = 'w2point5_';
        case 3
            job.roptions.prefix= 'w3_';
    end
    
    spm_run_normalise(job);
    
    % ensure cerebrum & cerebellum maska are connected
    for region = [{'cerebrum'}, {'cerebellum'}]
        
        h = spm_vol(...
            fullfile(...
            data.templateDirRoot,...
            data.templateDirName,...
            [job.roptions.prefix,'GM_',region{:},'.nii']...
            )...
            );
  
        gm = spm_read_vols(h);
        
        numClean = round(0.7*numel(find(gm))/2);
        gm = gwspm_adjust_mask(gm,numClean,6,region{:});
        
        [p,n,e] = fileparts(h.fname);
        h.fname = fullfile(p,[n,'_adjusted',e]);
        spm_write_vol(h,gm);
        
        % threshold & binarize (to create downsampled mask)
        mask = gm;
        mask(mask<=0.5) = 0;
        
        nnz(isnan(mask))
        mask(isnan(mask)) = 0;
        
        mask(find(mask)) = 1; %#ok<FNDSB>

        h.fname = fullfile(p,[n,'_mask_adjusted',e]);
        spm_write_vol(h,mask);
    end
end

