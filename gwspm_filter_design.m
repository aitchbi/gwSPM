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
function [g,gp,t,gb,gb2,gb3]=gwspm_filter_design(lmax,Nscales,hamShift,varargin)
% This function is a modified version of the function: 
% nsgwt_filter_design.m 
% by
% Nora Leonardi, Dimitri Van De Ville (2011)
% which itself is an extension based on the function
% sgwt_filter_design.m 
% from the the SGWT toolbox.

control_params={'designtype','default','lpfactor',20,...
    'a',2,'b',2,'t1',1,'t2',2,'wav_type','meyer','dmax',lmax/2,...
    };

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);
g=cell(Nscales+1,1);
gp=cell(Nscales+1,1);

switch designtype
    case 'default'
        
        lmin=lmax/lpfactor;
        t=sgwt_setscales(lmin,lmax,Nscales);
        
        gl = @(x) exp(-x.^4);
        glp = @(x) -4*x.^3 .*exp(-x.^4);
        
        gb= @(x) sgwt_kernel(x,'a',a,'b',b,'t1',t1,'t2',t2);
        gbp = @(x) sgwt_kernel_derivative(x,'a',a,'b',b,'t1',t1,'t2',t2);
        
        for j=1:Nscales
            g{j+1}=@(x) gb(x*t(j));
            gp{j+1}=@(x) gbp(t(j)*x)*t(j);
        end
        
        f=@(x) -gb(x);
        xstar=fminbnd(f,1,2);
        gamma_l=-f(xstar);
        lminfac=.6*lmin;
        gb2=@(x) gamma_l*gl(x); 
        
        g{1}=@(x) gamma_l*gl(x/lminfac);
        gp{1} = @(x) gamma_l*glp(x/lminfac)/lminfac;
        
    case 'mh'
        lmin=lmax/lpfactor;
        t=sgwt_setscales(lmin,lmax,Nscales);
        gb=@(x) sgwt_kernel(x,'gtype','mh');
        gl = @(x) exp(-x.^4);
        for j=1:Nscales
            g{j+1}=@(x) gb(t(j)*x);
        end
        lminfac=.4*lmin;  
        gb2=@(x) 1.2*exp(-1)*gl(x); 
        g{1}=@(x) 1.2*exp(-1)*gl(x/lminfac);

    case 'mey'
        switch wav_type
            case 'meyer'
                gb=@(x) sgwt_meyer(x,'t1',t1);
                gb2=@(x) sgwt_mey_h(x,'t1',t1);
                gb3=@(x) sgwt_meyer_end(x,'t1',t1);
            case 'simonc'
                gb=@(x) sgwt_simonc(x,'t1',t1);
                gb2=@(x) sgwt_simonc_h(x,'t1',t1);
                gb3=@(x) sgwt_simonc_end(x,'t1',t1);
            case 'shannon'
                gb=@(x) sgwt_shannon(x,'t1',t1);
                gb2=@(x) sgwt_shannon_h(x,'t1',t1);
            case 'papadakis'
                gb=@(x) sgwt_papadakis(x,'t1',t1);
                gb2=@(x) sgwt_papadakis_h(x,'t1',t1);
            case 'held'
                gb=@(x) sgwt_held(x,'t1',t1);
                gb2=@(x) sgwt_held_h(x,'t1',t1);
            otherwise
                error('Unknown wavelet type');
        end
        t=sgwt_setscales_mey(lmax,Nscales,hamShift,'t1',t1);
        
        switch wav_type
            case {'meyer', 'simonc'}
                for j=1:Nscales-1
                    g{j+1}=@(x) gb(t(j)*x);
                    gp{j+1}=@(x) 0;
                end
                g{Nscales+1}=@(x) gb3(t(end)*x); 
            otherwise
                for j=1:Nscales
                    g{j+1}=@(x) gb(t(j)*x);
                    gp{j+1}=@(x) 0;
                end
        end
        g{1}=@(x) gb2(t(1)*x);
        gp{1}=@(x) 0;
        
    otherwise
        error('Unknown design type');
end
