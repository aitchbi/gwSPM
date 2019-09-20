%
%     gwSPM: Graph-based, Wavelet-based Statistical Parametric Mapping
%       (v1.00)
%
% 	Author: 
%	Hamid Behjat 
% 
%   Biomedical Signal Processing Group, 
%   Dept. of Biomedical Engineering,
%   Lund University, Sweden
%   
%   gwSPM is a toolbox designed as an ad-on for SPM. For a 
%   thorough description of the implemented methodology, 
%   please refer to [1].
%   
%
%   gwSPM is an extension of the previously developed 
%   Wavelet-based, SPM (WSPM) [2] toolbox, developed by: 
% 
%   Dimitri Van De Ville
%   Medical Image Processing Lab
%   Institute of Bioengineering, EPFL, Lausanne, Switzerland.
%   Department of Radiology, University of Geneva, Switzerland.
%
%   As such, a significant part of the code for the statistical 
%   analysis and wavelet processing (step 5.) has been adopted 
%   from that work, but adapted to the graph setting.
%   
%
%   The toolbox can be downloaded at https://github.com/hbehjat/gwSPM
%
%
% 	References
%
% 	[1] H. Behjat, N. Leonardi, L. S?rnmo, D. Van De Ville
%           "Anatomically-adapted graph wavelets for improved
%            group-level fMRI activation mapping", 
%            NeuroImage, 123, pp. 185-199, Dec 2015.
%          
% 	[2] D. Van De Ville, T. Blu, M. Unser,
%           "Integrated Wavelet Processing and Spatial Statistical
%            Testing of fMRI Data", NeuroImage, 23(4), pp 1472-1485,
%            Dec 2004.
%

spmVer=spm('Ver');
gwspmVer=1.00;

SPMid=spm('FnBanner','gwSPM',sprintf('%3.2f',gwspmVer)); 
[Finter,Fgraph,CmdLine] = spm('FnUIsetup',strcat('gwSPM',...
    ' (v',sprintf('%3.2f',gwspmVer),')'));
addpath(fileparts(mfilename('fullpath')));
addpath(genpath(strcat(fileparts(mfilename('fullpath')),...
    filesep,'external')));

choice = spm_input('Select action',1,'m',...
    'step 1. Initial inputs & gray matter (GM) template construction |step 2. cerebrum & cerebellum GM extraction & data normalisation |step 3. classical SPM second-level analysis |step 4. design graph transform, transform data & estimate GLM |step 5. integrated wavelet processing & spatial statistical testing |step 6. display results (detected parameter maps) ',[1 2 3 4 5 6]);

switch choice
    case 1, action = 'templateDesign';
    case 2, action = 'extract';
    case 3, action = 'secondLevelSPM';
    case 4, action = 'design&estimate';
    case 5, action = 'reconstruct';
    case 6, action = 'results';
end

if strcmp(action,'templateDesign')
    
    % Input Data: 
    % 1. a strutural scan
    % 2. a mean functional volume 
    % 3. a set of first-level contrast maps 
    % for a set of subjects 
    
    N_subjects = spm_input('Number of Subjects','+1','e',[],1);
    
    message = sprintf('For each of the %d subjects, you will be prompted to specify:\n\n(1) a T1 structural scan\n\n(2) the mean functional volume (created when realigning functional data)\n\n(3) first-level contrast maps (one or several).',N_subjects);
    uiwait(msgbox(message));
    
    message = sprintf('First, you will be promted to select the structural scans for the %d subjetcs.\n\n------------------------------------------------------------\n\nBefore you proceed, note that for each of the subjetcs'' T1 scan, the Anterior Commissure (AC) should have been set to the origin (0 0 0).\n\nIf you are not sure whether this has been done, please quit and make sure this has been done.\n\nThis can be done using spm ''Display'' button. If the origin is not around AC, it can be easily fixed as follows:\n\n1) load the T1 volume\n2) place the crosshair in good proximity of the AC\n3) click ''Set Origin''\n4) click ''Reorient...''\n5) click ''Done''\n6) click ''No''.',N_subjects);
    uiwait(msgbox(message,'Info'));
    
    if ~spm_input('Proceed or exit to check AC s?','+1','b','proceed|Exit',[1 0]);
        return
    end

    [Im,sts] = spm_select(N_subjects,'image',sprintf('select the structural (T1) scans of the  %d subjetcs',N_subjects));
    if sts==0,
        message = sprintf('There was an error loading your files.\n\n');
        uiwait(msgbox(message, 'Error','error'));
        return
    end
    
    images = cell(N_subjects,1);
    for iS=1: N_subjects
        images(iS) = {deblank(Im(iS,:))};
    end
    
    gwspm_info = gwspm_info_initialize(Im);
    
    message = sprintf('Now, you will be promted to specify the location of the mean functional volume for each subject, separately.\n\nThis file results from the ''realign'' phase when preprocessing the functional volumes. Typically, it should be located in the functional data folder (in the first session folder if multiple sessions exists), and is prefixed with the tag: ''mean''\n\nThe file name of the structural scan of each subject will be specified.\n\nPlease be patient, and make sure you correctly specifty the mean functional volume for the ''specified'' subject.');
    uiwait(msgbox(message, 'Instructions'));
    
    meanVols = cell(N_subjects,1);
    for iS=1: N_subjects
        message = sprintf(strcat('Select the mean functional volume for the following subject:\n\n',images{iS}));
        uiwait(msgbox(message, 'Instructions'));
        
        [Im,sts] = spm_select(1,'^mean','select mean functional volume for the specified subject');
        meanVols(iS) = {Im};
    end
    
    message = sprintf('Now, you will be promted to specify the location of the contrast map(s) for each subject, separately.\n\nNote that contrast map(s) are prefixed with the tag: ''con''\n\nThe file name of the structural scan of each subject will be specified.\n\nPlease be patient, and make sure you correctly specifty the contrasts for the ''specified'' subject.');
    uiwait(msgbox(message, 'Instructions'));
    
    contrasts = cell(1,N_subjects);
    for iS=1: N_subjects
        message = sprintf(strcat('Select the first-level contrast map(s) for the following subject:\n\n',images{iS}));
        uiwait(msgbox(message, 'Instructions'));
        
        [Im,sts] = spm_select([1 Inf],'^con','select first-level contrast(s) for the specified subject');
        nn = size(Im,1);
        temp = cell(nn,1);
        for iC=1:size(Im,1)
            temp(iC) = {Im(iC,:)};
        end
        contrasts(iS) = {temp};
    end
    
    % construct group template (~ 3-5 hours)
    [gwspm_info,t] = gwspm_construct_template(images,gwspm_info);

    % normalise the contrast maps to the teplate spcae 
    gwspm_preprocess_contrasts(N_subjects,images,contrasts,meanVols,gwspm_info);
    
    % The 'gwspm_info' file will be loaded in the next steps to
    % locate the constructed templates and normalised data
    try
        if ~exist(strcat(fileparts(mfilename('fullpath')),filesep,'temp'),'dir')
            mkdir(strcat(fileparts(mfilename('fullpath')),filesep,'temp'))
        end
        save(strcat(fileparts(mfilename('fullpath')),...
            filesep,'temp',filesep,'gwspm_info.mat'),'gwspm_info')
    catch
        
    end

    save(strcat(gwspm_info.templateDirRoot,...
        filesep,gwspm_info.templateDirName,...
        filesep,'gwspm_info.mat'),'gwspm_info')
    
    
    clc
    fprintf('-------------------------------------------------------------------------\n')
    fprintf('A file named ''gwspm_info.mat'' was save in: %s \n', strcat(gwspm_info.templateDirRoot,filesep,gwspm_info.templateDirName))
    fprintf('\n')
    fprintf('Remember this location, as you will be asked to upload this file for the next step.\n')
    %fprintf('*****\n')
    %fprintf('A backup of the file is also save in: %s \n',strcat(fileparts(mfilename('fullpath')),filesep,'temp',filesep))
    %fprintf('But note that this file will be overwriten if you run step 1. for another dataset.\n')
    fprintf('-------------------------------------------------------------------------\n')
    
elseif strcmp(action,'extract')
    
    
    message = sprintf('Now, you will be promted to specify the location of the gwspm_info.mat file.\n\nThis file was created in step 1.\n\nThe file is located in a folder named: \n\ngwspm_templates \n\nwhich should be located in the same directory as the directory where the T1 strutural scan of the first subject you uploded.');
    uiwait(msgbox(message, 'Instructions'));
    
    [dummy,sts] = spm_select(1,'^gwspm','select the gwspm_info.mat file');
    if sts
        load(dummy)
    else
        error('gwspm_info.mat file not loaded.')
    end
    
    % extract cerebrum and cerebellum templates
    hObject = gwspm_extract(gwspm_info);
    
    waitfor(hObject)
    
    % downsample to MNI 2, 2.5 and 3 mm cubic resolutions
    gwspm_postprocess_extractions(gwspm_info)
    
    
elseif strcmp(action,'secondLevelSPM')
    
    message = sprintf('Now, you will be promted to specify the location of the gwspm_info.mat file.\n\nThis file was created in step 1.\n\nThe file is located in a folder named: \n\ngwspm_templates \n\nwhich should be located in the same directory as the directory where the T1 strutural scan of the first subject you uploded.');
    uiwait(msgbox(message, 'Instructions'));
    
    [dummy,sts] = spm_select(1,'^gwspm','select the gwspm_info.mat file');
    if sts
        load(dummy)
    else
        error('gwspm_info.mat file not loaded.')
    end
    
    dummy = strcat(sprintf('NOTE: This information is also diplayed in your MATLAB command line.'),...
        sprintf('\n\n----------------------------------------'));
    message = strcat(...
        sprintf('\n\nFor this step, you need to perform second-level analysis with SPM.'),...
        sprintf('\n\nIn the directories where the first-level analysis contrasts of each subject are located, a folder named:'),...
        sprintf('\n\n%s',gwspm_info.contrastDirName),...
        sprintf('\n\nhas been created. Inside these folders you will find the first-level contrasts which have been normalised'),...
        sprintf(' to the constructed template space (which is also MNI), at 3 different resolutions: 2, 2.5 and 3 mm cubic.'),...
        sprintf('\n\nPlease proceed as follows:'),...
        sprintf('\n\n\nStep 1. Decide which resolution you want to consider. (recommended: 3 mm cubic, as otherwise the computational cost can be too high.)'),...
        sprintf('\n\nStep 2. Use SPM to smooth the contrats that you want to work with, for instance the ones starting with ''w3c_'' if you have considered working at 3 mm cubic resoluton'),...
        sprintf('\n\nIMPORTANT note: please do not change the default spm prefix for smoothing volumes.. i.e. ''s'' '),...
        sprintf('\n\nStep 3. Perform second level analysis with SPM. Note that only T-contrast are supported in the current implementation. For example, you can run one-sample T-test on the first-level contrats.'),...
        sprintf('\n\nWhen you are done with your second level analysis, return to the toolbox and run the next step in the pipeline, i.e. Step 4.'));
    
    clc
    fprintf('-------------------------------------------------------------------------\n')
    fprintf(strcat(message,'\n'))
    fprintf('-------------------------------------------------------------------------\n')
    
    uiwait(msgbox(strcat(dummy,message)));
    
    
elseif ~strcmp(action, 'results')
    
    currDir = pwd;
    
    CONSV_THRESH_EST_N=40;
    
    CONST_DO_ANALYSIS=1;
    
    message = sprintf('Now, you will be promted to specify the location of the gwspm_info.mat file.\n\nThis file was created in step 1.\n\nThe file is located in a folder named: \n\ngwspm_templates \n\nwhich should be located in the same directory as the directory where the T1 strutural scan of the first subject you uploded.');
    uiwait(msgbox(message, 'Instructions'));
    
    [dummy,sts] = spm_select(1,'^gwspm','select the gwspm_info.mat file');
    if sts
        load(dummy)
    else
        error('gwspm_info.mat file not loaded.')
    end
    
    correctSPMs = 0;
    while ~correctSPMs
        
        if strcmp(spmVer,'SPM2') || strcmp(spmVer,'SPM5')
            error(strcat('Maybe its time to update your SPM -- your current SPM version: ', spmVer))
        elseif strcmp(spmVer,'SPM8') || strcmp(spmVer,'SPM12')
            [SPMs,sts]=spm_select([1,Inf],'^SPM\.mat$','Select SPM.mat files');
            if ~sts, SPM = []; xSPM = []; return; end
            %swd = spm_str_manip(currSPM,'H');
        else
            fprintf('! Error: unknown SPM version\n');
            return
        end
        
        % Check whether the files associated to the selected SPM s are in the same MNI resolution space
        if size(SPMs,1)>1
            for i=1:size(SPMs,1)
                load(deblank(SPMs(i,:)));
                if iscell(SPM.xY.P)
                    [pathstr,name,ext] = fileparts(deblank(SPM.xY.P{1}));
                else
                    [pathstr,name,ext] = fileparts(deblank(SPM.xY.P(1,:)));
                end
                switch name(2:4)
                    case 'w2c'
                        check1 = 2;
                    case 'w2p'
                        check1 = 2.5;
                    case 'w3c'
                        check1 = 3;
                    otherwise
                        
                        correctSPMs = 0;
                        
                        message = sprintf('The following selected SPM is not compatible for gwSPM:\n\n %s \n\nThe SPM files that you select should be those that result form ''second-level analysis'' on a set of smoothed, preprocessed contrasts that result from step 2 (these files are located in the folders named: ''gwspm_import_contrasts'').\n\n You will be prompted to either:\n\n i) re-select a set of compatible SPM files \n\n or \n\n ii) exit the program',deblank(SPMs(i,:)));
                        uiwait(msgbox(message, 'Instructions'));
                        
                        reSelect=spm_input('','+1','b','Re-select SPM s|Exit',[1 0]);
                        
                        if reSelect
                            break
                        else
                            return
                        end
                end
                
                if i==1
                    check2 = check1;
                else
                    if check1~=check2
                        
                        correctSPMs = 0;
                        
                        message = sprintf('The files associated to the %d chosen SPM s,\n are not in the same MNI resolution space.\n\n You will be prompted to either:\n\n i) re-select a set of compatible SPM files \n\n or \n\n ii) exit the program',size(SPMs,1));
                        uiwait(msgbox(message, 'Instructions'));
                        
                        reSelect=spm_input('','+1','b','Re-select SPM s|Exit',[1 0]);
                        
                        if reSelect
                            break
                        else
                            return
                        end
                    else
                        correctSPMs = 1;
                    end
                end
                clear SPM
            end
        else
            correctSPMs = 1;
        end
        disp('dummy')
    end
    
    defaults = spm_get_defaults;
    
    switch action
        
        case 'design&estimate'
            
            for iSPM=1:size(SPMs,1)
                
                clear SPM
                
                currSPM = deblank(SPMs(iSPM,:));
                load(currSPM);
                
                try
                    cd(SPM.swd)
                catch
                    error('The loaded SPM file has been moved from its original directory.')
                end
                
                try
                    WS=SPM.gWavelet;
                    WS.num=WS.num+1;
                    mkdir(strcat(SPM.swd,filesep,'gwspm',filesep,sprintf('gwspm_analysis_%02d',WS.num)))
                catch
                    clear WS;
                    WS.num=1;
                    mkdir(strcat(SPM.swd,filesep,'gwspm'))
                    mkdir(strcat(SPM.swd,filesep,'gwspm',filesep,sprintf('gwspm_analysis_%02d',WS.num)))
                    
                    WS.files = cell(SPM.nscan,1);
                    
                    for i = 1:SPM.nscan
                        if iscell(SPM.xY.P)
                            [pathstr,name,ext] = fileparts(deblank(SPM.xY.P{i}));
                        else
                            [pathstr,name,ext] = fileparts(deblank(SPM.xY.P(i,:)));
                        end
                        
                        WS.files{i} = strcat(pathstr,filesep,name(2:end),ext);
                        if i==1
                            switch name(2:4)
                                case 'w2c', WS.spmResolution = 2;
                                case 'w2p', WS.spmResolution = 2.5;
                                case 'w3c',  WS.spmResolution = 3;
                            end
                        end
                    end
                    
                end
                
                gwspm_info.resolution = WS.spmResolution;
                
                try
                    gwspm_info.atomsDir = strcat(fileparts(deblank(SPMs(1,:))),...
                        filesep,'gwspm',filesep,'atoms');
                catch
                    
                end
                
                if iSPM==1 % The selected set of SPM s will use the same graph transform design
                    sts=0;
                    while ~sts
                        [frameInfo,graphInfo,resol,parCheck,sts] = gwspm_design_transform(gwspm_info);
                        
                        if parCheck
                            delete(gcp('nocreate'))
                        end
                        
                        if ~sts
                            message = sprintf('It seems that you had not confirmed the settings in the GUI for designing the transform.\n\n You will be prompted to either:\n\n i) go back to the GUI and re-check the design. Make sure you ''check'' the box before closing the GUI.\n\n or\n\n ii) exit the program\n\n');
                            uiwait(msgbox(message, 'Instructions'));
                            
                            reSelect=spm_input('','+1','b','Go back to GUI|Exit',[1 0]);
                            if reSelect==0
                                return
                            end
                            continue
                        end
                        
                        if gwspm_info.resolution~=resol
                            sts = 0;
                            
                            message = sprintf(strcat('Incompatible tranform design and data.\n\n',...
                                ' The files associate to the selected SPM(s), are in %s cubic mm MNI',...
                                ' resolution.\n\n The selected design is in %s cubic mm MNI resolution.',...
                                '\n\n You will be prompted to either:\n\n i) go back to the GUI and',...
                                ' choose a design in %s mm resolution.\n\n or\n\n ii) exit the program,',...
                                ' and if you wish, do second level analysis with the contrast files that are',...
                                ' in %s mm resolution and then return, for performing anaysis at %s',...
                                ' resolution.\n\n',gwspm_info.resolution,resol,gwspm_info.resolution,resol,resol));
                            uiwait(msgbox(message, 'Instructions'));
                            
                            reSelect=spm_input('','+1','b','Go back to GUI|Exit',[1 0]);
                            if reSelect==0
                                return
                            end
                            continue
                            
                        end
                    end
                end
                
                % Check whether the chosen graph wavelet transform already exists.
                matchFound = 0;
                for iter=1:WS.num-1,
                    if WS.trans(iter).wav_scales==frameInfo.wav_scales...
                            && WS.trans(iter).chebyOrder== frameInfo.chebyOrder...
                            && WS.trans(iter).shift== frameInfo.shift
                        
                        matchFound=iter;
                        break;
                    end
                end
                
                if matchFound
                    message = sprintf(strcat('Found a matching graph decomposition in SPM structure.\n\n',...
                        ' You will be prompted to either:\n\n 1) ',...
                        ' exit the program\n\n or\n\n 2)',...
                        ' overwrite the existing associated files'));
                    uiwait(msgbox(message, 'Instructions'));
                    
                    switch spm_input('','+1','b','Exit|Overwrite',[1 2]);
                        case 1
                            return
                        case 2
                            WS.num = matchFound;
                    end
                end
                
                
                % Construct transformation info structure
                %----------------------------------------------------------------------
                clear trans

                trans.wav_scales = frameInfo.wav_scales;
                trans.shift = frameInfo.shift;
                trans.chebyOrder = frameInfo.chebyOrder;
                
                trans.cbr.g = frameInfo.g_cbr;
                trans.cbr.arange = frameInfo.arange_cbr;
                trans.cbr.c = frameInfo.c_cbr;
                trans.cbr.L = graphInfo.L_cbr;
                trans.cbr.indices = graphInfo.indCbr;
                
                trans.cbl.g = frameInfo.g_cbl;
                trans.cbl.arange = frameInfo.arange_cbl;
                trans.cbl.c = frameInfo.c_cbl;
                trans.cbl.L = graphInfo.L_cbl;
                trans.cbl.indices = graphInfo.indCbl;
                                
                trans.wav_dim = SPM.xY.VY(1).dim(1:3);

                trans.descrip = strcat('SGWT_scales',num2str(frameInfo.wav_scales),'_tunningFactor',num2str(frameInfo.shift),'_chebyOrder',num2str(frameInfo.chebyOrder));
                trans.descrip_short = sprintf('sgwt_sc%d_sh%d_or%d',frameInfo.wav_scales,frameInfo.shift,frameInfo.chebyOrder);
                
                
                % Compute Wavelet Transform for each image/volume
                %----------------------------------------------------------------------
                
                % Initialize dimensions
                rsz=trans.wav_dim;%SPM.xY.VY(1).dim(1:3); HB june 2016
                sz1=rsz(1);
                sz2=rsz(2);
                sz3=rsz(3);
                sz4=sum(SPM.nscan);
                psz3=rsz(3)*(trans.wav_scales+1); % number of slices of wavelet coefficients for each volume
                Y1=zeros(rsz(1:3));
                Y2=zeros(rsz(1:3));
                
                spm_progress_bar('Init',100,'Computing Graph wavelet Transforms','');
                
                % 24 June !
                clear files;
                files = cell(SPM.nscan,1);
                
                for iter=1:sum(SPM.nscan),
                    
                    if iscell(WS.files), % new multi-session convention
                        fname=char(WS.files{iter});
                    else % back-ward compatibility: old mono-session convention
                        fname=WS.files(iter,:);
                    end
                    
                    tmp=[ length(fname) strfind(fname,',') ];
                    fname22=fname(1:tmp(end)-1); % construct filename without volume separator
                    if strcmp(spmVer,'SPM8') || strcmp(spmVer,'SPM12')
                        if exist(fname22,'file')==0, % file not found: strip directory path (perhaps has moved)
                            tmp=[ 0 strfind(fname,filesep)];
                            fname=fname(tmp(end)+1:end);
                        end
                    end
                    
                    Vi=spm_vol(fname);
                    
                    [pathstr,name,ext] = fileparts(fname);
                    decomposed_fname = strcat(pathstr,filesep,name,'_',trans.descrip_short,ext(1:strfind(ext,',')-1));
                    
                    fprintf('%-40s: %s/%s\n','Computing graph wavelet transform',num2str(iter),num2str(sum(SPM.nscan)));
                    
                    [Y,XYZ]=spm_read_vols(Vi,0);
                    Y1=Y1+Y;
                    Y2=Y2+Y.^2;
                    
                    Y(isnan(Y))=0; % second-level analysis: masked areas are put to NaN
                    
                    % graph wavelet decomposition (cerebrum subgraph)
                    f_cbr=Y(trans.cbr.indices);
                    YYt_cbr=sgwt_cheby_op(f_cbr(:),trans.cbr.L,trans.cbr.c,trans.cbr.arange);
                    
                    % graph wavelet decomposition (cerebellum subgraph)
                    f_cbl=Y(trans.cbl.indices);
                    YYt_cbl=sgwt_cheby_op(f_cbl(:),trans.cbl.L,trans.cbl.c,trans.cbl.arange);
                    
                    Yt=zeros([rsz(1),rsz(2),psz3]);
                    for i=1:trans.wav_scales+1
                        Yt(trans.cbr.indices+(i-1)*numel(Y))=YYt_cbr{i};
                        Yt(trans.cbl.indices+(i-1)*numel(Y))=YYt_cbl{i};
                    end
                    
                    V=struct('fname',decomposed_fname,'dim',[sz1,sz2,psz3],'dt',[16 spm_platform('bigend')],'mat',Vi(1).mat,'descrip',strcat('WAV Y',num2str(iter),'-',trans.descrip));
                    
                    spm_write_vol(V,Yt);
                    
                    % 24 June !
                    %files(iter,:)=decomposed_fname;
                    files{iter} = decomposed_fname;
                    
                    spm_progress_bar('Set',100*iter/sum(SPM.nscan));
                
                end
                
                spm_progress_bar('Clear');
                
                trans.VY=files;
                Y1=Y1/sum(SPM.nscan); % mean contrast volume
                Y2=Y2/sum(SPM.nscan); % mean squared contrast volume
                Ys=Y2-Y1.^2;
                clear Y1 Y2;
                
                
                %======================================================================
                % Estimate a design using wavelets
                %======================================================================
                
                % Loading WT volumes
                %----------------------------------------------------------------------
                fprintf('%-40s: ','Loading graph wavelet images');
                Vi=spm_vol(trans.VY);
                
                spm_progress_bar('Init',100,'Running Wavelet Domain GLM','');
                
                clear Y;
                B=NaN*ones(size(SPM.xX.X,2),sz1*sz2,psz3);
                
                X = SPM.xX.nKX; % from SPM: design matrix K*X       
                A = SPM.xX.pKX; % from SPM: pseudo-inverse K*W*X    
                W = SPM.xX.W;   % from SPM: whitening matrix W      
                erdf = SPM.xX.erdf; % from SPM: effective degrees of freedom   
                Qs=zeros(sz1*sz2*psz3,1);                         
                
                
                % Each decomposition volume has sz3 slices.
                % Thus, there are in total, psz3 slices of decomposition coefficients
                % A volume (across subjects) is constructed for each slice
                % GLM is then fitted for all the graph wavelt coefficents that lie in that slice
                
                checkk=0;
                for p=1:psz3
                    pp=mod(p,sz3);
                    if pp==0,
                        pp=sz3;
                        fprintf('%-40s: %s/%s\n','Round',num2str(p),num2str(psz3));
                    end
                    
                    if p==1
                        disp('')
                        disp('GLM on scaling function coefficients')
                    elseif pp==1 && p~=1
                        scGlm=floor(p/sz3);
                        fprintf('%-40s: scale %s/%s\n','GLM on wavelet coefficients',...
                            num2str(scGlm),num2str(trans.wav_scales));
                    end
                    
                    fprintf('%-40s: %s/%s\n','Round',num2str(p),num2str(psz3));
                    
                    % find graph indices in the current slice, at the current scale
                    dummy1 = trans.cbr.indices>=(1+(pp-1)*(sz1*sz2));
                    dummy2 = trans.cbr.indices<=(pp*(sz1*sz2));
                    dummy3 = trans.cbr.indices(dummy1&dummy2);
                    
                    dummy1 = trans.cbl.indices>=(1+(pp-1)*(sz1*sz2));
                    dummy2 = trans.cbl.indices<=(pp*(sz1*sz2));
                    dummy4 = trans.cbl.indices(dummy1&dummy2);
                    
                    indiceRef = [dummy3;dummy4];
                    checkk = checkk+numel(indiceRef);
                    
                    % Adjust indiceref to appropriate scale
                    indiceP = indiceRef+(floor(p/sz3)*(sz1*sz2*sz3));
                    [xP yP zP] = ind2sub([sz1,sz2,psz3],indiceP);
                    xyzP = [xP,yP,zP]';
                    
                    % Construct data set to run GLM for current slice (p)
                    Y = zeros(sz4,numel(indiceP));
                    for i=1:sz4,
                        Y(i,:) = spm_get_data(Vi(i),xyzP);
                    end
                    Y = Y.*repmat(SPM.xGX.gSF,1,numel(indiceP)); % global scaling ???
                    
                    
                    % GLM for each graph wavelet coefficient
                    %----------------------------------------------------------------------
                    KWY = spm_filter(SPM.xX.K,W*Y);
                    clear Y;
                    B(:,indiceRef-(pp-1)*(sz1*sz2),p) = A*KWY;
                    KWY2 = zeros(sz4,sz1*sz2);
                    KWY2(:,indiceRef-(pp-1)*(sz1*sz2)) = KWY;
                    res = spm_sp('r',SPM.xX.xKXs,KWY2);
                    Qs(1+(p-1)*sz1*sz2:p*sz1*sz2) = squeeze(sum(res.^2))/erdf;
                    clear res;
                    
                    spm_progress_bar('Set',100*p/psz3);
                end
                
                Qs = reshape(Qs,sz1,sz2,psz3);
                                
                % Write volumes with B (model estimate parameters) and Qs (errors)
                for iter = 1:size(SPM.xX.X,2),
                    
                    fname = strcat(SPM.swd,filesep,'gwspm',filesep,...
                        sprintf('gwspm_analysis_%02d',WS.num),filesep,...
                        'B',num2str(iter),'_sc',num2str(trans.wav_scales),...
                        '_sh',num2str(trans.shift),'_or',...
                        num2str(trans.chebyOrder),'.nii');
                    
                    V = struct('fname',fname,'dim',[sz1,sz2,psz3],'dt',...
                        [16 spm_platform('bigend')],'mat',Vi{1}.mat,...
                        'descrip',strcat('WAV B',num2str(iter),' - ',...
                        trans.descrip));
                    
                    tmp = reshape(squeeze(B(iter,:,:)),sz1,sz2,psz3); 
                    trans.Vbeta(iter) = spm_write_vol(V,tmp);
                end
                
                % write VResMS to disk
                fname = strcat(SPM.swd,filesep,'gwspm',filesep,...
                    sprintf('gwspm_analysis_%02d',WS.num),filesep,...
                    'V','_sc',num2str(trans.wav_scales),'_sh',...
                    num2str(trans.shift),'_or',num2str(trans.chebyOrder),...
                    '.nii'); 
                
                V = struct('fname',fname,'dim',[sz1,sz2,psz3],'dt',...
                    [16 spm_platform('bigend')],'mat',Vi{1}.mat,...
                    'descrip',strcat('WAV VResMS ',trans.descrip));
                
                trans.VResMS = spm_write_vol(V,Qs);
                
                % write VResMS2 to disk
                fname = strcat(SPM.swd,filesep,'gwspm',filesep,...
                    sprintf('gwspm_analysis_%02d',WS.num),filesep,...
                    'Ys','_sc',num2str(trans.wav_scales),'_sh',...
                    num2str(trans.shift),'_or',num2str(trans.chebyOrder),...
                    '.nii');
                
                V = struct('fname',fname,'dim',rsz(1:3),'dt',...
                    [16 spm_platform('bigend')],'mat',Vi{1}.mat,...
                    'descrip',strcat(' VResMS2',trans.descrip));
                
                trans.VResMS2 = spm_write_vol(V,Ys);
                
                spm_progress_bar('Clear');
                
                WS.trans(WS.num) = trans;
                SPM.gWavelet=WS;
                
                save(currSPM,'SPM'); 
                
                clear trans
                
                clear files
                
            end
            
        case 'reconstruct'
            %======================================================================
            % Compute results for a contrast vector (over all selected SPM s)
            %======================================================================
            
            for iSPM=1:size(SPMs,1)
                
                clear SPM
                
                currSPM = deblank(SPMs(iSPM,:));
                load(currSPM);
                
                try
                    cd(SPM.swd)
                catch
                    error('The loaded SPM file has been moved from its original directory.')
                end
                
                try
                    WS=SPM.gWavelet;
                catch
                    message = sprintf(strcat('No graph wavelet transforms are available',...
                        ' for this experiment.\n\nMake sure you first run steps 1 to 4,',...
                        ' before running this step (i.e. step 5).\n\n By clicking OK,',...
                        ' the algorithm will be terminated.\n\n'));
                    uiwait(msgbox(message, 'Instructions'));
                    
                    fprintf('%s','No graph wavelet transforms are available for this experiment!');
                    return
                end
                
                
                % Select graph wavelet transform set
                %----------------------------------------------------------------------
                clear sets;
                for iter=1:numel(WS.trans)
                    sets(iter)={WS.trans(iter).descrip_short};
                end
                
                if iSPM==1 && size(SPMs,1)>1 && numel(WS.trans)>1
                    message = sprintf(strcat('You will be prompted to select one of the graph ',...
                        ' wavelet transform sets that were constructed in step 4.\n\nThe same ',...
                        ' transform should exist for all the %d SPM files that you have loaded,',...
                        'and they will be automatically uploaded.\n\n'),size(SPMs,1));
                    uiwait(msgbox(message, 'Instructions'));
                end
                
                if iSPM==1
                    select=spm_input('Select a graph wavelet transform set',1,'m',sets);
                    trans=WS.trans(select);
                    transNumm=select;
                    tempData.transform = WS.trans(select).descrip_short;
                else
                    dummy = zeros(numel(WS.trans),1);
                    for iter=1:numel(WS.trans)
                        dummy(iter) = strcmp(WS.trans(iter).descrip_short,tempData.transform);
                    end
                    
                    if ~isempty(find(dummy,1))
                        trans=WS.trans(find(dummy,1));
                    else
                        error(strcat(sprintf('gwspm error: \n\n The graph wavelet transform: \n\n ''%s'' ',...
                            tempData.transform),...
                            sprintf(' \n\n does not exist for the data associated to the following selected SPM analysis'),...
                            sprintf(': \n\n %s ',currSPM),...
                            strcat(sprintf('\n\n Please make sure all selected SPM s have been decomposed'),...
                            ' with the same graph wavelet transform in the previous step of the pipeline.')))
                    end
                    transNumm=find(dummy);
                end
                tempData.transNum(iSPM) = transNumm; % for creating directories and to load Qs
                
                % Select statistical test (todo: F-test)
                %-----------------------------------------
                clear sets sets_info;
                sets_idx=1;
                for iter=1:length(SPM.xCon),
                    if SPM.xCon(iter).STAT=='T'
                        if ~isfield(WS,'wCon') || ( isfield(WS,'wCon') && (length(WS.wCon)<iter || isempty(WS.wCon(iter).siglevel)) ), % short-circuit version to prevent error in case wCon does not exist
                            sets(sets_idx)={ sprintf('{%s}: %s',SPM.xCon(iter).STAT,SPM.xCon(iter).name) };
                            sets_info(sets_idx).num=iter;
                            sets_idx=sets_idx+1;
                        end
                    end
                end
                
                if sets_idx==1,
                    message = strcat(sprintf('The current Implementation is only suited for T-tests.'),...
                        sprintf('\n\nNo T-contrasts were defined for the SPM analysis:'),...
                        sprintf('\n\n %s \n\n',currSPM));
                    uiwait(msgbox(message, 'Instructions'));
                    error('No T-contrasts defined!'); 
                end
                
                
                % As their names might be different, its not possible to 
                % select it for first SPM and then extend to other SPMs.
                % Thus, will be prompted, for each SPM. However, if there 
                % is only a single contrats, there will be no prompt, 
                % and the contrast will be autoselected. 
                select=spm_input('Select test','+1','m',sets);
                name=char(SPM.xCon(sets_info(select).num).name);
                xCon=SPM.xCon(sets_info(select).num);
                
                wCon.mask=[];
                
                % Input significance level
                if iSPM==1
                    wCon.typeI=spm_input('Type I error control','+1','m',...
                        { 'Strong (corrected Bonferroni)',...
                        'Weak (corrected FDR)','Uncorrected'});
                    
                    if wCon.typeI<3,  % corrected
                        wCon.siglevel=spm_input('p value (significance level)','+1','r',0.05,1,[0,500]);
                    else % uncorrected
                        wCon.siglevel=spm_input('p value (significance level)','+1','r',0.001,1,[0,0.5]);
                    end
                    tempData.wConTypeI = wCon.typeI;
                    tempData.wConSiglevel = wCon.siglevel;
                else
                    wCon.typeI = tempData.wConTypeI;
                    wCon.siglevel = tempData.wConSiglevel;
                end
                
                % Check whether analysis already exists
                dummy = gwspm_check_xCon(SPM,name,wCon,trans.descrip_short);
                
                if dummy>0, % match found
                    if spm_input('Contrast with these parameters already exists!',...
                            1,'b',{'Stop','Recompute'},[1 0]),
                        CONST_DO_ANALYSIS=0;
                    end
                    fprintf('Found same analysis in contrast library: replacing contrast %d\n',checkContrast);
                    num=dummy;
                else
                    num=length(SPM.xCon)+1;
                end
                
                if CONST_DO_ANALYSIS==1, 
                    
                    % Initialize dimensions
                    %------------------------------------------------------------------
                    rsz=SPM.xY.VY(1).dim(1:3);
                    sz1=trans.wav_dim(1);
                    sz2=trans.wav_dim(2);
                    sz3=trans.wav_dim(3);
                    sz4=sum(SPM.nscan);
                    
                    % Loading V (model residual errors)
                    %------------------------------------------------------------------
                    Qs=spm_read_vols(trans.VResMS);
                    if length(trans.VResMS2)>0,
                        Ys=spm_read_vols(trans.VResMS2);
                    end
                    
                    psz3=sz3*(trans.wav_scales+1);
                    
                    c=xCon.c;     % contarst vector
                    X=SPM.xX.nKX; % design matrix
                    A=SPM.xX.pKX; 
                    erdf=SPM.xX.erdf; % effective degrees of freedom
                    
                    if size(c,2)==1,
                        sc=c;
                    else 
                        sc=sum(xCon.c')';
                    end
                    
                    % LOADING B (model parameters estimates)
                    %------------------------------------------------------------------
                    %if prod([size(SPM.xX.X,2) sz1 sz2 psz3 8])/(1024^3)>0.1, % too much memory required?
                    if 1, % always memory safe (more and more high resolution studies...)
                        Bm=zeros(sz1,sz2,psz3);
                        for iter=1:size(SPM.xX.X,2), 
                            if c(iter)~=0,
                                B=spm_read_vols(trans.Vbeta(iter),0); 
                                Bm=Bm+c(iter)*B; % Bm: contrast map
                            end
                        end
                    else
                        B=spm_read_vols(trans.Vbeta,0);
                        B=reshape(B,sz1*sz2*psz3,size(SPM.xX.X,2))';
                        Bm=reshape(c'*B,sz1,sz2,psz3);
                    end
                    clear B;
                    
                    % Compute alpha
                    %------------------------------------------------------------------
                    if wCon.typeI<3, % Strong or weak typeI error control (corrected) 
                        wCon.alpha=1-wCon.siglevel/(numel(trans.cbr.indices)+(numel(trans.cbl.indices)));   % Adjusting for  multiple comparison problem
                    else   % uncorrected
                        wCon.alpha=1-wCon.siglevel;
                    end
                    
                    switch(xCon.STAT)
                        case 'T', % t-test
                            Qs=Qs*(c'*SPM.xX.Bcov*c); % old: Qs*(c'*inv(X'*X)*c) % H- this is our estimate of variance at each voxel (i.e. sigma^=sqrt(Qs)) (Ref: slide 25)
                            Q=Bm./sqrt(Qs);           
                            Qm=Bm;
                    end
                    
                    wCon.t_high=max(abs(Q(:))); 
                    
                    
                    % Build map with t-values
                    %------------------------------------------------------------------
                    a=1-wCon.alpha; 
                    
                    if wCon.typeI<3,     % corrected Bonferroni or corrected FDR
                        CONST_TYPEI=wCon.typeI;
                    else
                        CONST_TYPEI=1; % uncorrected
                    end
                    
                    % Bonferroni corrected significance level
                    alphaB = a; % initialize threshold for FDR  % H- this is the alphaB which we use in equation 27, Ref: 2004
                    
                    if iSPM==1
                        
                        if sum(SPM.nscan)<CONSV_THRESH_EST_N; % good estimation of noise variance?
                            CONST_SIGMA_UNKNOWN=1;
                            fprintf('Less than %d scans: switching to more conservative threshold estimation\n',CONSV_THRESH_EST_N);
                            % H- i.e. switching to 'General case - true sigma are unknown' (Ref: page 1478, year:2004)
                        else
                            CONST_SIGMA_UNKNOWN=0;
                        end
                        
                        % Initialize threshold values
                        if CONST_TYPEI==1,  % H- either uncorrected or Bonferroni corrected
                            sprintf('Computing wavelet domain and spatial doamin thresholds\n')
                            if CONST_SIGMA_UNKNOWN==1,
                                % sigma assumed unknown
                                [TW,TS]=threshold_search2d(min(150,sz4),min(150,erdf),a)
                                %TW=8.21; TS=0.8257;
                            else
                                % sigma assumed known
                                TW=sqrt(real(-LambertW(-1,-2*pi*a^2)))  % H- equation 27 in Ref:2004
                                TS=1/TW
                            end
                        elseif CONST_TYPEI==2,   %H- for FDR corrected (weak corrected)
                            
                            %TW=min(sqrt(real(-LambertW(-1,-2*pi*a^2))),wCon.t_high);%Dimitri
                            TW=sqrt(real(-LambertW(-1,-2*pi*(a*...
                                (numel(trans.cbr.indices)+numel(trans.cbl.indices))...
                                )^2)));
                            
                            TS=1/TW;
                            prev_dec=0;
                        end
                    end
                    
                    
                    % temporarily save the data, to load them back in,
                    % after having computed the time consuming absolute wavelet reconstructions 
                    dummy = strcat(SPM.swd,filesep,'gwspm',filesep,...
                        sprintf('gwspm_analysis_%02d',tempData.transNum(iSPM)),...
                        filesep,'temp_data');
                    
                    if exist(dummy,'dir')~=7
                        mkdir(dummy)
                    end
                    
                    save(strcat(dummy,filesep,'pre_rs1_data.mat'),...
                        'Qs','Qm','Q','TW','TS','CONST_TYPEI','alphaB',...
                        'CONST_SIGMA_UNKNOWN','erdf','sz4','rsz',...
                        'wCon','xCon','name','transNumm','num',...
                        'WS','trans','Bm')
                    
                    clear Qs Qm Q CONST_TYPEI alphaB ...
                        erdf sz4 rsz ...
                        wCon xCon name transNumm num ...
                        WS trans Bm
                    
                    % TW and TS are only computed for iSPM==1, and are not deleted so that 
                    % their computation can be skipped for the next SPM s . This is because 
                    % their value should bethe same for all the SPM s, since for all of them, 
                    % the factors that affect their value, i.e. the sixe of the graph graph, 
                    % the dimensions of the data, and the chosen error control and significance 
                    % level ''should'' be the same. 
                    
                end
            end
            

            % Absolute-value wavelet synthesis for all size(SPMs,1) datasets
            %------------------------------------------------------------------

            % load one of the trans files
            dummy = load(strcat(fileparts(deblank(SPMs(iSPM,:))),filesep,...
                'gwspm',filesep,sprintf('gwspm_analysis_%02d',...
                tempData.transNum(1)),filesep,'temp_data',...
                filesep,'pre_rs1_data.mat')); %fixed BUG: trans.mat [gwspm v1.00] >> pre_rs1_data.mat [gwspm v1.01] -- [6 Aug 2017].
            
            trans = dummy.trans;
            
            atomOptions.atomsDir = strcat(fileparts(deblank(SPMs(1,:))),...
                filesep,'gwspm',filesep,'private',filesep,'atoms_chunks',...
                filesep,'atoms_',trans.descrip_short);
            
            atomOptions.cbr_atomsDir = strcat(atomOptions.atomsDir,...
                filesep,'cerebrum');
            
            atomOptions.cbl_atomsDir = strcat(atomOptions.atomsDir,...
                filesep,'cerebellum');
            
            atomOptions.chunkTag = 'atoms_chunk_';
            atomOptions.loadAtoms = 0; 
            atomOptions.saveAtoms = 0; % change to '1' if you would like to save the wavelet constructions >> but requires upto 25 gigabytes of free memory on disk.. watch out.  
      
            szChunks =200;
            
            Qss = cell(size(SPMs,1),1);
            for  iSPM=1:size(SPMs,1)
                dummy = load(strcat(fileparts(deblank(SPMs(iSPM,:))),...
                    filesep,'gwspm',filesep,sprintf('gwspm_analysis_%02d',...
                    tempData.transNum(iSPM)),filesep,'temp_data',...
                    filesep,'pre_rs1_data.mat'),'Qs');
                Qss{iSPM} = dummy.Qs;
            end
           
            % Estimate compute time for the absolute value wavelet phase 
            [tMx,tMn,option] = gwspm_estimate_computeTime(trans,[],[],[]);
            
            
            message = strcat(...
                sprintf('An estimation of the approximate time required for running the remaining, final part of the ''gwspm analysis pipeline'' was performed.'),...
                sprintf('\n\n Based on the estimation, the computation time will be:'),...
                sprintf('\n\n %d to %d hours',round(tMn),round(tMx)),...
                sprintf('\n\nAs the required time is long, it is best if this part is perfomed over night, and if more time is needed'),...
                sprintf(' maybe running this step (step 4.) over a weekend would be best.'),...
                sprintf('\n\nOf course, if your computer is free, and you don''t need to use MATLAB meanwhile, well then please let the analysis proceed :-)'),...
                sprintf('\n\nYou will be prompted to make a choice.'));
            uiwait(msgbox(message, 'Instructions'));
            
            if spm_input('','+1','b','Proceed|Exit',[0 1])
                return
            end
            
            
            %------------------------------------------------------------------
            % The absolute-value reconstruction for all input SPM s
            % ***This follwing line is time consumin part of the algorithm.
            % (computes the denominator in Eq.(10) in Behjat, et. al, NeuroImage, 2015)
            Rs1s = gwspm_wrapper_synthesis_abs(trans,Qss,szChunks,atomOptions);
            %------------------------------------------------------------------
            
                       
            for  iSPM=1:size(SPMs,1)
                clear SPM
                currSPM = deblank(SPMs(iSPM,:));
                load(currSPM);
                Rs1 = Rs1s{iSPM};
                save(strcat(SPM.swd,filesep,'gwspm',filesep,...
                    sprintf('gwspm_analysis_%02d',tempData.transNum(iSPM)),...
                    filesep,'temp_data',filesep,'Rs1.mat'),'Rs1')
                clear Rs1
            end
            
            
            
            for  iSPM=1:size(SPMs,1)
                clear SPM
                currSPM = deblank(SPMs(iSPM,:));
                load(currSPM);
                
                dummy = strcat(SPM.swd,filesep,'gwspm',filesep,...
                        sprintf('gwspm_analysis_%02d',tempData.transNum(iSPM)),...
                        filesep,'temp_data');
                    
                load(strcat(dummy,filesep,'pre_rs1_data.mat'))
                load(strcat(dummy,filesep,'Rs1.mat'))
            
                % Normal wavelet synthesis
                Rm0 = gwspm_wrapper_synthesis(Qm,trans);
                
                % Start reconstruction (single iteration if Bonferroni; loop if FDR)
                CONST_ITERATE=1;
                while CONST_ITERATE>0,
                    tmp = Qm;                
                    idx = find(abs(Q)<=TW);  
                    tmp(idx) =0;            
                    
                    Rm=gwspm_wrapper_synthesis(tmp,trans);
                    Rs = Rs1;
                    
                    % spatial bias reduction
                    Rm = min(Rm,Rm0);
                    
                    % secon-level parameter map 
                    S = Rm./Rs+TW;
                    
                    S_det = S;
                    indNan = find(isnan(S));
                    S_det(indNan) = 0;
                    S_det(S<(TW+TS)) = 0;
                    
                    dec=numel(find(S>=TW+TS));

                    wCon.TW = TW;
                    wCon.TS = TS;
                    wCon.active = dec;
                 
                    try
                        dummy = pwd;
                        cd(strcat(SPM.swd,filesep,'gwspm',filesep,sprintf('gwspm_analysis_%02d',tempData.transNum(iSPM)),filesep,'temp_data'))
                        save Rs Rs
                        save S S
                        save S_det S_det
                        save dec dec
                    catch
                        error('something is fishy..')
                    end
                    
                    
                    if CONST_TYPEI==1,  
                        
                        fprintf('BON: detections=%d, thresholds [TW,TS]=%f,%f\n',dec,TW,TS);
                        CONST_ITERATE=0; 
                        
                    else
                        % False Discovery Rate
                        
                        fprintf('FDR: detections=%d, thresholds [TW,TS]=%f,%f\n',dec,TW,TS);
                        
                        % 1) count detections <- dec
                        
                        % 2) adapt threshold values
                        a=alphaB*max(1,dec);
                        if CONST_SIGMA_UNKNOWN==1,
                            % sigma assumed unknown
                            [TW,TS]=threshold_search2d(min(150,sz4),min(150,erdf),a);
                        else
                            % sigma assumed known
                            TW=sqrt(real(-LambertW(-1,-2*pi*a^2)));
                            TS=1/TW;
                        end
                        
                        % 3) stop?
                        if dec==prev_dec || CONST_ITERATE==15,
                            CONST_ITERATE=0;
                        else
                            CONST_ITERATE=CONST_ITERATE+1;
                            prev_dec=dec;
                        end
                        
                    end 
                    
                end 
                
                if wCon.typeI==1,
                    tmp='SIG';
                elseif wCon.typeI==2,
                    tmp='SIG-FDR';
                else
                    tmp='SIG-UNC';
                end
                
                if wCon.siglevel>=0.01,
                    tmp=[tmp sprintf(' (%.2f)',wCon.siglevel)];
                else
                    tmp=[tmp sprintf(' (%.0e)',wCon.siglevel)];
                end
                
                descrip=strcat(sprintf('%s %s(%4.2f) %s WAV',name,xCon.STAT,TW+TS,tmp),'(',trans.descrip_short,')'); 
                
                dummy = strcat(SPM.swd,filesep,'gwspm',filesep,sprintf('gwspm_analysis_%02d',WS.num));
                fname1 = strcat(dummy,filesep,sprintf('con_%04d.nii',num));
                fname2 = strcat(dummy,filesep,sprintf('dncon_%04d.nii',num));
                fname3 = strcat(dummy,filesep,sprintf('wspmT_%04d.nii',num));
                fname4 = strcat(dummy,filesep,sprintf('detects_%04d.nii',num));
                fname5 = strcat(dummy,filesep,sprintf('wavecon_%04d.nii',num));
                fname6 = strcat(dummy,filesep,sprintf('wavevar_%04d.nii',num));

                V=struct('fname',' ','dim',' ','dt',[16 spm_platform('bigend')],'mat',xCon.Vcon.mat,'descrip',' ','pinfo',[1 0 0]');

                clear Y
                Y=Rm0(1:rsz(1),1:rsz(2),1:rsz(3));
                if 0,  Y(find(M<0.5))=NaN; end   
               
                V.fname = fname1;
                V.descrip = strcat(xCon.STAT,'-contrast:',descrip);
                V.dim = rsz(1:3);
                wCon.Vcon = spm_write_vol(V,Y);
                
                clear Y
                Y=Rm(1:rsz(1),1:rsz(2),1:rsz(3));
                if 0,  Y(find(M<0.5))=NaN; end   
               
                V.fname = fname2;
                V.descrip = strcat(xCon.STAT,'-denoised contrast:',descrip);
                V.dim = rsz(1:3);
                wCon.Vdncon=spm_write_vol(V,Y);
                
                clear Y
                Y=S(1:rsz(1),1:rsz(2),1:rsz(3));
                if 0,  Y(find(M<0.5))=NaN; end   
                
                V.fname=fname3;
                V.descrip=strcat(xCon.STAT,'- wspmT:',descrip); 
                V.dim = rsz(1:3);
                wCon.Vspm = spm_write_vol(V,S);
                
                clear Y
                Y=S_det(1:rsz(1),1:rsz(2),1:rsz(3));
                if 0,  Y(find(M<0.5))=NaN; end   
                
                V.fname=fname4;
                V.descrip=strcat(xCon.STAT,'- detects:',descrip); 
                V.dim = rsz(1:3);
                wCon.detects = spm_write_vol(V,S_det);
                
                sQm=size(Qm); 
                if length(sQm)==2, sQm(3)=1; end
                
                V.fname = fname5;
                V.descrip = strcat(xCon.STAT,'-wavecontrast:',descrip);
                V.dim = sQm;
                wCon.Wcon = spm_write_vol(V,Qm);
                
                V.fname=fname6;
                V.descrip=strcat(xCon.STAT,'-wavevar:',descrip);
                V.dim = sQm;
                wCon.Wvar = spm_write_vol(V,Qs);
                
                xCon.name=descrip;
                xCon.Vcon=wCon.Vcon; 
                xCon.Vspm=wCon.Vspm;
                
               
                % adapt wCon structure if necessary
                clear wCon2;
                if isfield(WS,'wCon'), 
                    FN=fieldnames(WS.wCon(1));
                    for iter=1:length(FN),
                        if isfield(wCon,FN{iter})==0,
                            wCon.(FN{iter})=[];
                        end
                    end
                    FN=fieldnames(wCon);
                    for iter=1:length(FN),
                        if isfield(WS.wCon(1),FN{iter})==0,
                            WS.wCon(1).(FN{iter})=[];
                        end
                    end
                    FN=fieldnames(WS.wCon(1));
                    for iter=1:length(FN),
                        wCon2.(FN{iter})=wCon.(FN{iter});
                    end
                    WS.wCon(num)=wCon2;
                else
                    WS.wCon(num)=wCon;
                end
                
                SPM.gWavelet=WS;
                SPM.xCon(num)=xCon;
                
               clear Qs Qm Q TW TS CONST_TYPEI alphaB ...
                        CONST_SIGMA_UNKNOWN erdf sz4 rsz ...
                        prev_dec wCon xCon name transNumm...
                        num WS trans tComp
                    % and maybe more..
                    
                    if exist('CONST_DO_ANALYSIS','var') && CONST_DO_ANALYSIS==0,
                        % do not save SPM - no analysis done
                    elseif ~strcmp(action,'design&estimate') % HB since already saved
                        fprintf('%-40s: ','Saving SPM configuration')
                        save(currSPM,'SPM'); % use currSPM instead
                        fprintf('%30s\n','..SPM.mat saved')
                    end
                    
            end
            
            cd(currDir)
    end 
    
else
    
    message = sprintf('Now, you will be promted to specify the location of the gwspm_info.mat file.\n\nThis file was created in step 1.\n\nThe file is located in a folder named: \n\ngwspm_templates \n\nwhich should be located in the same directory as the directory where the T1 strutural scan of the first subject you uploded.');
    uiwait(msgbox(message, 'Instructions'));
    
    [dummy,sts] = spm_select(1,'^gwspm','select the gwspm_info.mat file');
    if sts
        load(dummy)
    else
        error('gwspm_info.mat file not loaded.')
    end
    
    message = strcat(...
        sprintf('The detected parameter map will be displayed, overlayed on the constructed group template.'),...
        sprintf('\n\nBetter visulaizations and more display options will be provided in a near future upgrade of gwSPM.'));
    uiwait(msgbox(message, 'Instructions'));
    
    dummy = '^detects';
    imgs{1} = spm_select(1,dummy,'select the detected parameter map -- starts with tag ''detects'' ');
    imgs{2} = strcat(gwspm_info.templateDirRoot,filesep,gwspm_info.templateDirName,filesep,'w1_GM_cerebrum_adjusted.nii');
    imgs{3} = strcat(gwspm_info.templateDirRoot,filesep,gwspm_info.templateDirName,filesep,'w1_GM_cerebellum_adjusted.nii');
    gwspm_plot_overlay(imgs,[],[],[],'results')
    
end

%-End: Cleanup GUI
%==========================================================================
spm_clf(Finter)
spm('FigName','Wavelets: finished',Finter,CmdLine);
spm('Pointer','Arrow')
fprintf('\n\n')
return

