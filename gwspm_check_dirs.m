function chk = gwspm_check_dirs(opts)
chk = 1;
if opts.saveAtoms || opts.loadAtoms
    cwd = pwd;
    try
        cd(opts.atomsDir)
    catch
        if opts.loadAtoms
            chk = 0;
            return
        end
        if opts.saveAtoms
            mkdir(opts.atomsDir)
        end
    end
    try
        cd(opts.cbr_atomsDir)
    catch
        if opts.loadAtoms
            chk = 0;
            return
        end
        if opts.saveAtoms
            mkdir(opts.cbr_atomsDir)
        end
    end
    try
        cd(opts.cbl_atomsDir)
    catch
        if opts.loadAtoms
            chk = 0;
            return
        end
        if opts.saveAtoms
            mkdir(opts.cbl_atomsDir)
        end

    end
    cd(cwd);
end
