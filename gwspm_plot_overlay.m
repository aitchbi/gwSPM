function gwspm_plot_overlay(imgs,curr_orient,curr_slicenum,hh,opts)

if ischar(imgs)
    imgs = cellstr(imgs);
end

obj = slover;

cscale = [];
obj.cbar = [];

switch opts
    case 'atom'
        %-Atom.
        obj.img(1).vol = spm_vol(imgs{1});
        obj.img(1).type = 'truecolour';
        dcmap = 'flow.lut';
        [mx,mn] = slover('volmaxmin', obj.img(1).vol);
        drange = [mn mx];
        cscale = [cscale 1];
        obj.cbar = [obj.cbar 1];
        obj.img(1).cmap = slover('getcmap',dcmap); 
        dummy = max(abs(drange));
        obj.img(1).range = [-dummy; dummy]; 
        obj.img(1).prop = 1;

        %-GM Template.
        obj.img(2).vol = spm_vol(imgs{2});
        obj.img(2).type = 'truecolour';
        obj.img(2).cmap = gray;
        [mx,mn] = slover('volmaxmin', obj.img(2).vol);
        obj.img(2).range = [mn mx];
        cscale = [cscale 2];
        obj.img(2).prop = 0.5;

    case 'results'
        obj.img(1).vol = spm_vol(imgs{1});
        obj.img(1).type = 'truecolour';
        obj.img(1).cmap = gray;
        [mx,mn] = slover('volmaxmin', obj.img(1).vol);
        obj.img(1).range = [mn mx];
        cscale = [cscale 2];
        obj.img(1).prop = 1;
        
        %-Cerebrum Template.
        obj.img(2).vol = spm_vol(imgs{2});
        obj.img(2).type = 'truecolour';
        obj.img(2).cmap = gray;
        [mx,mn] = slover('volmaxmin', obj.img(2).vol);
        obj.img(2).range = [mn mx];
        cscale = [cscale 2];
        obj.img(2).prop = 0.25;
        
        %-Cerebellum Template.
        obj.img(3).vol = spm_vol(imgs{3});
        obj.img(3).type = 'truecolour';
        obj.img(3).cmap = gray;
        [mx,mn] = slover('volmaxmin', obj.img(3).vol);
        obj.img(3).range = [mn mx];
        cscale = [cscale 2];
        obj.img(3).prop = 0.25;
end

if ~isempty(curr_orient)
    switch curr_orient
        case 1
            obj.transform = 'sagital';
        case 2
            obj.transform = 'coronal';
        case 3
            obj.transform = 'axial';
    end
end

if isempty(hh)
    obj.figure = spm_figure('GetWin', 'Graphics');
else
    obj.figure = hh;
end

obj = fill_defaults(obj);
sls = obj.slices;

switch opts
    case 'atom'
        obj.slices = sls(curr_slicenum);
    case 'results'
        obj.slices = sls;
end

paint(obj)

