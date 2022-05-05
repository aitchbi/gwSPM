function [f,x] = gwspm_view_design(g,arange,varargin)
% An extended version of sgwt_filter_design.m from the SGWT toolbox.

p = inputParser;
addParameter(p,'showLegend','yes');
addParameter(p,'Graph',0);
addParameter(p,'eigsToInterp',[]);
addParameter(p,'subSampleWarping',0);
addParameter(p,'warping','none');
addParameter(p,'plotLineWidth',1);
addParameter(p,'lambda',[]);
addParameter(p,'guiHandle',[]);
addParameter(p,'chebyOrder',[]);
addParameter(p,'onlyChebyApprox','no');
parse(p,varargin{:});
opts = p.Results;

lw = opts.plotLineWidth;
gh = opts.guiHandle;
cOrd = opts.chebyOrder;
eTI  = opts.eigsToInterp;
Graph = opts.Graph;

if isstruct(Graph)
    d = Graph.E;
    d = d(2:end)-d(1:end-1);
    I_ui = [1; find(logical(d)) + 1]; % indices of unidentical eigs
    
    if strcmp(opts.warping,'none')
        if isempty(eTI)
            f = Graph.E(I_ui);
            x = f;
        else 
            f = eTI;
            x = f;
        end
    else
        if isempty(eTI)
            f = Graph.E(I_ui);
            x = opts.warping(I_ui)*Graph.lmax;
        else
            if opts.subSampleWarping
                d1 = Graph.E(I_ui);
                interp_x = d1(1:opts.subSampleWarping:end-1);
                interp_x(end+1) = d1(end);
                
                d2 = opts.warping(I_ui)*Graph.lmax;
                interp_y = d2(1:opts.subSampleWarping:end-1);
                interp_y(end+1) = d2(end);
                
                f = eTI(eTI <= interp_x(end));
                x = abs(gsp_mono_cubic_warp_fn(interp_x,interp_y,f)); 
            else
                f = eTI;
                interp_x = Graph.E(I_ui);
                interp_y = opts.warping(I_ui)*Graph.lmax;
                x = abs(gsp_mono_cubic_warp_fn(interp_x,interp_y,f));
            end
        end
    end
else 
    x = linspace(arange(1),arange(2),1e3);
    f = x;
end


J = length(g)-1;
co = get(gca,'ColorOrder'); 
co = [co;co;co];
G = 0*x;
G_cheby = 0*x;

if ~isempty(cOrd)
    c = cell(1,length(g));
    for k = 1:length(g)
        c{k} = sgwt_cheby_coeff(g{k},cOrd,cOrd+1,arange);
    end
end

for n=0:J
    if isempty(gh)
        if isempty(cOrd)
            plot(...
                f,...
                g{1+n}(x),...
                'Color',co(1+n,:),...
                'LineWidth',lw...
                );
        elseif strcmp(opts.onlyChebyApprox,'yes')
            plot(...
                f,...
                sgwt_cheby_eval(x,c{1+n},arange),...
                'Color',co(1+n,:),...
                'LineWidth',lw...
                );
        else
            plot(...
                f,...
                g{1+n}(x),...
                'Color',co(1+n,:),...
                'LineWidth',lw...
                );
            plot(...
                f,...
                sgwt_cheby_eval(x,c{1+n},arange),...
                '-.',...
                'Color',co(1+n,:),...
                'LineWidth',lw...
                );
        end
    else
        if isempty(cOrd)
            plot(...
                f,...
                g{1+n}(x),...
                'Color',co(1+n,:),...
                'LineWidth',lw,...
                'Parent',gh...
                );
        elseif strcmp(opts.onlyChebyApprox,'yes')
            plot(...
                f,...
                sgwt_cheby_eval(x,c{1+n},arange),...
                'Color',co(1+n,:),...
                'LineWidth',lw,...
                'Parent',gh...
                );
        else
            plot(...
                f,...
                g{1+n}(x),...
                'Color',co(1+n,:),...
                'LineWidth',lw,...
                'Parent',gh...
                );
            plot(...
                f,...
                sgwt_cheby_eval(x,c{1+n},arange),...
                '-.',...
                'Color',co(1+n,:),...
                'LineWidth',lw,...
                'Parent',gh...
                );
        end
    end
    
    if n==0
        if isempty(gh)
            hold on
        else
            hold(gh,'on')
        end
    end
    
    if isempty(cOrd)
        G = G+g{1+n}(x).^2;
    else
        G = G+g{1+n}(x).^2;
        G_cheby = G_cheby + sgwt_cheby_eval(x,c{1+n},arange).^2;
    end
end

if isempty(gh)
    if isempty(cOrd)
        plot(f,G,'k:','LineWidth',lw);
    elseif strcmp(opts.onlyChebyApprox,'yes')
        plot(f,G_cheby,'k:','LineWidth',lw);
    else
        plot(f,G,'k:','LineWidth',lw);
        plot(f,G_cheby,'k-.','LineWidth',lw);
    end
else
    if isempty(cOrd)
        plot(f,G,'k:','LineWidth',lw,'Parent',gh);
    elseif strcmp(opts.onlyChebyApprox,'yes')
        plot(f,G_cheby,'k:','LineWidth',lw,'Parent',gh);
    else
        plot(f,G,'k:','LineWidth',lw,'Parent',gh);
        plot(f,G_cheby,'k-.','LineWidth',lw,'Parent',gh);
    end
end

leglabels{1} = 'h';
for j=1:J
    leglabels{1+j} = sprintf('g_{%d}',j);
end
leglabels{J+2} = 'G';

if ~isempty(opts.lambda)
    y = -.5*ones(size(opts.lambda,1),1);
    stem(opts.lambda,y, 'k');
    leglabels{J+3} = 'Lambda';
end

switch opts.showLegend
    case 'yes'
        legend(leglabels)
        title([...
            'Scaling function kernel h(x), ',...
            'Wavelet kernels g(t_j x), ',...
            'Sum of Squares G, ',...
            'and Frame Bounds'
            ]);
    case 'no'
end

if isempty(gh)
    hold off
else
    hold(gh,'off')
end
end