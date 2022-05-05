function [atoms,tComp,RunType] = gwspm_construct_atoms(trans,iFirst,iLast,subG,RunType)

if isempty(RunType)
    if gwspm_check_par([])
        RunType = 'parallel';
        SwitchPoolOff = 1;
    else
        RunType = 'sequential';
        SwitchPoolOff = 0;
    end
else
    SwitchPoolOff = 0;
end

wav_dim = trans.wav_dim;
wav_scales = trans.wav_scales;
Nv = prod(wav_dim);

switch subG
    case 'cbr'
        indices = trans.cbr.indices;
        L = trans.cbr.L;
        c = trans.cbr.c;
        arange = trans.cbr.arange;
    case 'cbl'
        indices = trans.cbl.indices;
        L = trans.cbl.L;
        c = trans.cbl.c;
        arange = trans.cbl.arange;
end
Ng = length(indices);

if isempty(iFirst)
    iFirst = 1;
end

if isempty(iLast)
    iLast = Ng*(wav_scales+1);
end

indiceT=[];
for i = 1:wav_scales+1
    indiceT = [
        indiceT
        indices+(i-1)*Nv
        ]; %#ok<AGROW>
end

szChunk = iLast-iFirst+1;

atoms = zeros(Ng,szChunk);

switch RunType
    case 'sequential'
        
        tic
        for index = 1:szChunk
            iter = (iFirst-1)+index;
            f = zeros(Ng, wav_scales+1);
            d1 = ceil(iter/Ng);
            d2 = iter-(d1-1)*Ng;
            f(d2,d1) = 1;
            atoms(:,index) = sgwt_inverse(f,L,c,arange);
        end
        
    case 'parallel'
        
        tic
        parfor index = 1:szChunk
            iter = (iFirst-1)+index;
            f = zeros(Ng,wav_scales+1);
            d1 = ceil(iter/Ng);
            d2 = iter-(d1-1)*Ng;
            f(d2,d1) = 1;
            atoms(:,index) = sgwt_inverse(f,L,c,arange);
        end
end

tComp.total = toc;
tComp.single = tComp.total /length(iFirst:iLast);

if SwitchPoolOff
    delete(gcp('nocreate'))
end
