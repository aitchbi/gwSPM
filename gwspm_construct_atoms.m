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
function [atoms,tComp,option] = gwspm_construct_atoms(trans,iFirst,iLast,subG,option)

if isempty(option)
    if gwspm_check_par([])
        option = 'parallel';
        switchPoolOff = 1;
    else
        option = 'sequential';
        switchPoolOff =0;
    end
else
    switchPoolOff =0;
end

wav_dim=trans.wav_dim;
wav_scales=trans.wav_scales;
sizeVol=prod(wav_dim);

switch subG
    case 'cbr'
        indice = trans.cbr.indices;
        L = trans.cbr.L;
        c = trans.cbr.c;
        arange = trans.cbr.arange;
    case 'cbl'
        indice = trans.cbl.indices;
        L = trans.cbl.L;
        c = trans.cbl.c;
        arange = trans.cbl.arange;
end
gSize = numel(indice);

if isempty(iFirst)
    iFirst = 1;
end
if isempty(iLast)
    iLast = gSize*(wav_scales+1);
end

indiceT=[];
for i=1:wav_scales+1,
    indiceT=[indiceT;indice+(i-1)*sizeVol]; %#ok<AGROW>
end

szChunk = iLast-iFirst+1;
atoms = zeros(numel(indice),szChunk);

switch option
    case 'sequential'
        
        tic
        for index = 1:szChunk
            iter = (iFirst-1)+index;
            f = zeros(numel(indice), wav_scales+1);
            dummy1=ceil(iter/gSize);
            dummy2 = iter-(dummy1-1)*gSize;
            f(dummy2,dummy1)=1;
            atoms(:,index)=sgwt_inverse(f,L,c,arange);
        end
        
    case 'parallel'
        
        tic
        parfor index = 1:szChunk
            iter = (iFirst-1)+index;
            f = zeros(numel(indice), wav_scales+1);
            dummy1=ceil(iter/gSize);
            dummy2 = iter-(dummy1-1)*gSize;
            f(dummy2,dummy1)=1;
            atoms(:,index)=sgwt_inverse(f,L,c,arange);
        end
end

tComp.total = toc;
tComp.single = tComp.total /numel(iFirst:iLast);

if switchPoolOff
    delete(gcp('nocreate'))
end
