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
function V=gwspm_wrapper_synthesis(Y,trans)

V=zeros(trans.wav_dim);

for subG = [{'cbr'}, {'cbl'}]
    
    switch subG{:}
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
    
    for i=1:trans.wav_scales+1
        f{i}=Y(indice+(i-1)*prod(trans.wav_dim)); %#ok<AGROW>
    end
    
    V(indice) = sgwt_inverse(f,L,c,arange);
        
end
