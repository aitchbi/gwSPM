function V = gwspm_wrapper_synthesis(Y,trans)
V = zeros(trans.wav_dim);
for subG = [{'cbr'}, {'cbl'}]
    switch subG{:}
        case 'cbr'
            inds = trans.cbr.indices;
            L = trans.cbr.L;
            c = trans.cbr.c;
            arange = trans.cbr.arange;
        case 'cbl'
            inds = trans.cbl.indices;
            L = trans.cbl.L;
            c = trans.cbl.c;
            arange = trans.cbl.arange;
    end
    for k = 1:trans.wav_scales+1
        f{k} = Y(inds+(k-1)*prod(trans.wav_dim)); %#ok<AGROW>
    end
    V(inds) = sgwt_inverse(f,L,c,arange);
end
