function  gwspm_preprocess_contrasts(Nsubj,images,contrasts,meanVols,opts)

%-Preliminaries------------------------------------------------------------ 
% make gwspm directory for contrasts
for iS =1:Nsubj
    d = contrasts{iS}(1);
    job.parent = {fileparts(d{:})};
    job.name = opts.contrastDirName ;
    cfg_run_mkdir(job);
end

% file names of 'coregistered' & 'normalised' contrasts
c_contrasts = cell(size(contrasts));
w_contrasts = cell(size(contrasts));
for iS =1:Nsubj
    d1 = cell(size(contrasts{iS}));
    d2 = cell(size(contrasts{iS}));
    for iC=1:size(contrasts{iS},1)
        [p,n,e] = fileparts(contrasts{iS}{iC});
        d1(iC) = {fullfile(p,['c_',n,e])};
        d2(iC) = {fullfile(p,['wc_',n,e])};
    end
    c_contrasts(iS) = {d1};
    w_contrasts(iS) = {d2};
end

%-Coregister---------------------------------------------------------------
% coreg contrasts with strutural data (using the mean functional volume)
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

for iS =1:Nsubj
    job.ref = images(iS);
    job.source = meanVols(iS);
    job.other = contrasts{iS};
    spm_run_coreg(job);
end

%-DARTEL Normalise to MNI Space--------------------------------------------
clear job
job.template = {...
    fullfile(...
    opts.templateDirRoot,...
    opts.templateDirName,...
    'Template_6.nii')...
    };

for iS=1:Nsubj
    [p,n,e] = fileparts(images{iS});
    job.data.subj(iS).flowfield = {...
        fullfile(p,'u_rc1',[n,'_Template',e])...
        };
    job.data.subj(iS).images = c_contrasts{iS};
end

job.preserve = 0;
job.fwhm     = [0 0 0];
job.bb       = [-78 -112 -70; 78 76 85];

for voxdim = [2, 2.5, 3]
    job.vox = [voxdim,voxdim,voxdim];
    spm_dartel_norm_fun(job); 
    for iS=1:Nsubj
        clear job2
        job2.files = w_contrasts{iS};
        d = fileparts(contrasts{iS}{1});
        job2.action.moveren.moveto = {...
            fullfile(d,opts.contrastDirName)...
            };
        job2.action.moveren.patrep.pattern = 'w';
        switch voxdim
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





