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
function check = gwspm_check_dirs(option)

check = 1;
if option.saveAtoms || option.loadAtoms
    dummy = pwd;
    try
        cd(option.atomsDir)
    catch
        if option.loadAtoms
            check = 0;
            return
        end
        if option.saveAtoms
            mkdir(option.atomsDir)
        end
    end
    try
        cd(option.cbr_atomsDir)
    catch
        if option.loadAtoms
            check = 0;
            return
        end
        if option.saveAtoms
            mkdir(option.cbr_atomsDir)
        end
    end
    try
        cd(option.cbl_atomsDir)
    catch
        if option.loadAtoms
            check = 0;
            return
        end
        if option.saveAtoms
            mkdir(option.cbl_atomsDir)
        end

    end
    cd(dummy)
end
