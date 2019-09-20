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
function [V,rs1_cbr,rs1_cbl,tComp] = gwspm_wrapper_synthesis_abs(trans,Qss,szChunks,option)

if ~isempty(option) && isstruct(option)
    if option.saveAtoms || option.loadAtoms
        saveLoadAtoms = 1;
        if ~gwspm_check_dirs(option);
            error('Incompatible directories.')
        end
    else
        saveLoadAtoms = 0;
    end
else
    saveLoadAtoms = 0;
    option.saveAtoms = 0;
    option.loadAtoms = 0;
end

% considering several Qs volumes in one go 
sz = numel(Qss);
tmps =  cell(sz,1);
for j=1:numel(Qss)
    tmps{j}=sqrt(Qss{j});
end

if ~option.loadAtoms
    % use parpool if available
    if gwspm_check_par([])
        optionPar = 'parallel';
    else
        optionPar = 'sequential';
    end
end

% compute separately for cerebrum and cerebellum subgraphs
for subG = [{'cbr'}, {'cbl'}]
    
    tic
    switch subG{:}
        case 'cbr'
            indice = trans.cbr.indices;
            if saveLoadAtoms
                partialPath = strcat(option.cbr_atomsDir,filesep,option.chunkTag);
            end
            printTag = ' Cerebrum';
            
        case 'cbl'
            indice = trans.cbl.indices;
            if saveLoadAtoms
                partialPath = strcat(option.cbl_atomsDir,filesep,option.chunkTag);
            end
            printTag = ' Cerebellum';
            
    end
    
    gSize = numel(indice);
    
    iLast = gSize*(trans.wav_scales+1);
    
    indiceT=[];
    for i=1:trans.wav_scales+1,
        indiceT=[indiceT;indice+(i-1)*prod(trans.wav_dim)]; %#ok<AGROW>
    end
    
    fprintf(strcat('Absolute value wavelet reconstruction on c',printTag(3:end),'...\n'))
    spm_progress_bar('Init',100,strcat('Absolute Value Wavelet Reconstruction - ',printTag),'');
    
    
    % initialize rs1 vectors
    rs1 = cell(sz,1);
    for j = 1:sz
        rs1{j} = zeros(size(indice));
    end
    
    N_chunks = ceil(iLast/szChunks);
    nChunk =1;
    
    for iChunk=1:szChunks:iLast-mod(iLast,szChunks)
        
        if option.loadAtoms
            load(strcat(partialPath,num2str(nChunk),'.mat'))
        else
            atoms = gwspm_construct_atoms(trans,iChunk,iChunk+szChunks-1,subG{:},optionPar);
            if option.saveAtoms
                save(strcat(partialPath,num2str(nChunk),'.mat'),'atoms')
            end
        end
        
        count = 1;
        for iter = iChunk:iChunk+szChunks-1
            i=indiceT(iter);
            for j = 1:sz
                rs1{j}=rs1{j}+tmps{j}(i)*abs(atoms(:,count));
            end
            count = count +1;
        end
        clc
        clear atoms
        nChunk = nChunk+1;
        spm_progress_bar('Set',100*nChunk/N_chunks);
    end
    
    % last chunk -- probably has less than szChunk atoms
    if option.loadAtoms
        load(strcat(partialPath,num2str(nChunk),'.mat'))
    else
        atoms = gwspm_construct_atoms(trans,iLast-mod(iLast,szChunks)+1,iLast,subG{:},optionPar);
        if option.saveAtoms
            save(strcat(partialPath,num2str(nChunk),'.mat'),'atoms')
        end
    end
    
    count =1;
    for iter= iLast-mod(iLast,szChunks)+1:iLast
        i=indiceT(iter);
        for j = 1:sz
            rs1{j}=rs1{j}+tmps{j}(i)*abs(atoms(:,count));
        end
        count = count +1;
    end
    clc
    clear atoms
    
    spm_progress_bar('Clear');
    fprintf(strcat('Done. \n'))
    
    switch subG{:}
        case 'cbr'
            rs1_cbr = rs1;
            tComp.cbr = toc;
        case 'cbl'
            rs1_cbl = rs1;
            tComp.cbl = toc;
    end
    
    clear rs1
    
end

if ~option.loadAtoms
    if strcmp(optionPar,'parallel')
        delete(gcp('nocreate'))
    end
end

V= cell(sz,1);
for j=1:sz
    V{j} = zeros(trans.wav_dim);
    V{j}(trans.cbr.indices) = rs1_cbr{j};
    V{j}(trans.cbl.indices) = rs1_cbl{j};
end






















%         case 'parallel'
%
%
%             for j=1:sz
%                 eval(['dummy_',num2str(j),'= zeros(numel(indice),numel(indice)*wav_scales+1);'])
%             end
%
%             parfor iter = iFirst:iLast
%
%                 i=indiceT(iter);
%                 f = zeros(numel(indice), wav_scales+1);
%                 dummy1=ceil(iter/gSize);
%                 dummy2 = iter-(dummy1-1)*gSize;
%                 f(dummy2,dummy1)=1;
%                 wave=sgwt_inverse(f,L,c,arange);
%                 %rs1=rs1+tmps(i)*abs(wave);
%
%
%                 %for j = 1:sz
%                 %dummy(:,iter)=tmps{j}(i)*abs(wave);
%                 %end
%
%                 % I am aware that the following code could be simply written as above :-)
%                 % This is just a requirement for 'parfor' to work better.
%                 if sz==1
%                     dummy_1(:,iter)=tmps_1(i)*abs(wave);
%                 else
%                     dummy_2(:,iter)=tmps_2(i)*abs(wave);
%                     if sz>2
%                         dummy_3(:,iter)=tmps_3(i)*abs(wave);
%                         if sz>3
%                             dummy_4(:,iter)=tmps_4(i)*abs(wave);
%                             if sz>4
%                                 dummy_5(:,iter)=tmps_5(i)*abs(wave);
%                                 if sz>5
%                                     dummy_6(:,iter)=tmps_6(i)*abs(wave);
%                                     if sz>6
%                                         dummy_7(:,iter)=tmps_7(i)*abs(wave);
%                                         if sz>7
%                                             dummy_8(:,iter)=tmps_8(i)*abs(wave);
%                                             if sz>8
%                                                 dummy_9(:,iter)=tmps_9(i)*abs(wave);
%                                                 if sz>9
%                                                     dummy_10(:,iter)=tmps_10(i)*abs(wave);
%                                                     if sz>10
%                                                         dummy_11(:,iter)=tmps_11(i)*abs(wave);
%                                                         if sz>11
%                                                             dummy_12(:,iter)=tmps_12(i)*abs(wave);
%                                                             if sz>12
%                                                                 dummy_13(:,iter)=tmps_13(i)*abs(wave);
%                                                                 if sz>13
%                                                                     dummy_14(:,iter)=tmps_14(i)*abs(wave);
%                                                                     if sz>14
%                                                                         dummy_15(:,iter)=tmps_15(i)*abs(wave);
%                                                                         if sz>15
%                                                                             dummy_16(:,iter)=tmps_16(i)*abs(wave);
%                                                                             if sz>16
%                                                                                 dummy_17(:,iter)=tmps_17(i)*abs(wave);
%                                                                                 if sz>17
%                                                                                     dummy_18(:,iter)=tmps_18(i)*abs(wave);
%                                                                                     if sz>18
%                                                                                         dummy_19(:,iter)=tmps_19(i)*abs(wave);
%                                                                                         if sz>19
%                                                                                             dummy_20(:,iter)=tmps_20(i)*abs(wave);
%                                                                                             if sz>20
%                                                                                                 dummy_21(:,iter)=tmps_21(i)*abs(wave);
%                                                                                             else
%                                                                                                 [pathstr,name] = fileparts(mfilename('fullpath'));
%                                                                                                 message = sprintf(strcat('A maximum of 21 SPMs was considered when writing the code for this implementation,',...
%                                                                                                     ' which we thought should be sufficient; i.e. studying a maximum 21 different contrasts.\n\n',...
%                                                                                                     'If you have bumped into this error message, well, you have proved me wrong. To consider more contrats,',...
%                                                                                                     ' and also use the parallel processing option, please easily edit the code in the following function:\n\n',...
%                                                                                                     pathstr,filesep,name,'.m','\n\n As you will see, adjusting the code will be a piece of cake ;-)\n\nThe',...
%                                                                                                     ' location of the file will be also printed in the command line.\n\nSimply type:',...
%                                                                                                     ' ''edit gwspm_wrapper_synthesis_abs'' in the command line to reach the file and edit it.\n\nYou will have to re-run this step of analysis.\n\n'));
%                                                                                                 sprintf(strcat(pathstr,filesep,name,'.m'))
%                                                                                                 uiwait(msgbox(message, 'Error'));
%                                                                                             end
%                                                                                         end
%                                                                                     end
%                                                                                 end
%                                                                             end
%                                                                         end
%                                                                     end
%                                                                 end
%                                                             end
%                                                         end
%                                                     end
%                                                 end
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%
%             end
%             t = toc;
%             rs1 = cell(sz,1);
%             for j=1:sz
%                 eval(['rs1{j} = sum(dummy_',num2str(j),',2);'])
%             end
%
%             tComp.parallel = t/numel(iFirst:iLast);
%             clear t
%             delete(gcp('nocreate'))
%
% %         case 'parallel_easy'
% %             dummy = cell(sz,1);
% %             for j=1:sz
% %                 dummy{j} = zeros(numel(indice),numel(indice)*wav_scales+1);
% %             end
% %
% %             tic
% %             parfor iter = iFirst:iLast
% %
% %                 i=indiceT(iter);
% %                 f = zeros(numel(indice), wav_scales+1);
% %                 dummy1=ceil(iter/gSize);
% %                 dummy2 = iter-(dummy1-1)*gSize;
% %                 f(dummy2,dummy1)=1;
% %                 wave=sgwt_inverse(f,L,c,arange);
% %                 %rs1=rs1+tmps(i)*abs(wave);
% %
% %
% %                 %for j = 1:sz
% %                 %dummy(:,iter)=tmps{j}(i)*abs(wave);
% %                 %end
% %
% %                 for j = 1:sz
% %                     dummy{j}(:,iter)=tmps{j}(i)*abs(wave);
% %                 end
% %             end
% %             t = toc
% %
% %             rs1 = cell(sz,1);
% %             for j=1:sz
% %                 rs1{j} = sum(dummy{j},2);
% %             end
% %
% %             tComp.parallel = t/numel(iFirst:iLast);
% %             clear t
% %
%         case 1
%             for crapFor=1:1
%                 tic
%                 iter=iFirst-1;
%                 for i=indiceT(iFirst:iLast)',
%                     iter=iter+1;
%                     %if rem(iter,10)==0, iter, end
%
%                     fIndex=ceil(iter/gSize);
%                     f{fIndex}(iter-(fIndex-1)*gSize)=1;
%                     wave=sgwt_inverse(f,trans.L,trans.c,trans.arange); % <<<
%                     rs1=rs1+tmps(i)*abs(wave);
%                     f{fIndex}(iter-(fIndex-1)*gSize)=0;
%                 end
%                 toc
%             end
%         case 2
%             for crapFor=1:1
%                 tic
%                 iter=iFirst-1;
%                 for i=indiceT(iFirst:iLast)',
%                     iter=iter+1;
%                     %if rem(iter,10)==0, iter, end
%
%                     for iF=1:wav_scales+1,
%                         f{iF} = zeros(size(indice));
%                     end
%                     dummy1=ceil(iter/gSize);
%                     f{dummy1}(iter-(dummy1-1)*gSize)=1;
%                     wave=sgwt_inverse(f,trans.L,trans.c,trans.arange);
%                     rs1=rs1+tmps(i)*abs(wave);
%                 end
%                 toc
%             end
%     end
%
%     switch subG{:}
%         case 'cbr'
%             rs1_cbr = rs1;
%         case 'cbl'
%             rs1_cbl = rs1;
%     end
%
% end
%
% V = cell(sz,1);
% for j=1:sz
%     V{j} = zeros(trans.wav_dim);
%     V{j}(trans.cbr.indices) = rs1_cbr{j};
%     V{j}(trans.cbl.indices) = rs1_cbl{j};
% end
%
% %
% % V=zeros(trans.wav_dim);
% % V(trans.cbr.indices)=rs1_cbr;
% % V(trans.cbl.indices)=rs1_cbl;
%
% %%
% if 0
%     tmps=sqrt(Qs);
%     trans=SPM.Wavelet.trans(1);
%     indice=trans.indice;
%     gSize=numel(indice);
%     wav_dim=trans.wav_dim;
%     wav_scales=trans.wav_scales;
%     sizeVol=prod(wav_dim);
%
%     indiceT=[];
%     for i=1:trans.wav_scales+1,
%         indiceT=[indiceT;indice+(i-1)*sizeVol]; %#ok<AGROW>
%     end
%
%     ttt=zeros([wav_dim(1:2),wav_dim(3)*(wav_scales+1)]);
%     rs1=zeros(size(indice)); %Rs1: output
%
%     for j=1:trans.wav_scales+1,
%         f{j}=ttt(indice+(j-1)*sizeVol); %#ok<SAGROW>
%     end
%
%     % Here is where the the loop takes place..
%     %(needs to be repeated 100k and above)
%     switch option
%         case 'parallel'
%             dummy = zeros(numel(indice),numel(indice)*wav_scales+1);
%
%             tic
%             parfor iter = iFirst:iLast,
%
%                 i=indiceT(iter);
%                 f = zeros(numel(indice), trans.wav_scales+1);
%                 dummy1=ceil(iter/gSize);
%                 dummy2 = iter-(dummy1-1)*gSize;
%                 f(dummy2,dummy1)=1;
%                 wave=sgwt_inverse(f,trans.L,trans.c,trans.arange);
%                 %rs1=rs1+tmps(i)*abs(wave);
%                 dummy(:,iter)=tmps(i)*abs(wave);
%             end
%             %toc
%             rs1 = sum(dummy,2);
%
%         case 'sequential'
%             dummy = zeros(numel(indice),numel(indice)*wav_scales+1);
%             tic
%             for iter = iFirst:iLast,
%
%                 i=indiceT(iter);
%                 f = zeros(numel(indice), trans.wav_scales+1);
%                 dummy1=ceil(iter/gSize);
%                 dummy2 = iter-(dummy1-1)*gSize;
%                 f(dummy2,dummy1)=1;
%                 wave=sgwt_inverse(f,trans.L,trans.c,trans.arange);
%                 %rs1=rs1+tmps(i)*abs(wave);
%                 dummy(:,iter)=tmps(i)*abs(wave);
%             end
%             toc
%             rs1 = sum(dummy,2);
%
%         case 1
%             for crapFor=1:1
%                 tic
%                 iter=iFirst-1;
%                 for i=indiceT(iFirst:iLast)',
%                     iter=iter+1;
%                     %if rem(iter,10)==0, iter, end
%
%                     fIndex=ceil(iter/gSize);
%                     f{fIndex}(iter-(fIndex-1)*gSize)=1;
%                     wave=sgwt_inverse(f,trans.L,trans.c,trans.arange); % <<<
%                     rs1=rs1+tmps(i)*abs(wave);
%                     f{fIndex}(iter-(fIndex-1)*gSize)=0;
%                 end
%                 toc
%             end
%         case 2
%             for crapFor=1:1
%                 tic
%                 iter=iFirst-1;
%                 for i=indiceT(iFirst:iLast)',
%                     iter=iter+1;
%                     %if rem(iter,10)==0, iter, end
%
%                     for iF=1:trans.wav_scales+1,
%                         f{iF} = zeros(size(indice));
%                     end
%                     dummy1=ceil(iter/gSize);
%                     f{dummy1}(iter-(dummy1-1)*gSize)=1;
%                     wave=sgwt_inverse(f,trans.L,trans.c,trans.arange);
%                     rs1=rs1+tmps(i)*abs(wave);
%                 end
%                 toc
%             end
%
%     end
% end
