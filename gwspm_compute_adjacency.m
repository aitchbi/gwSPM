function A = gwspm_compute_adjacency(mask,conn,varargin)
% Initial version of code by Elena Najdenovska & Nora Leonardi, 2011.

[w,pow,s,dist,sw]=process_options(...
    varargin,...
    'Weight','no',...
    'pow',[],...
    'S',[],...
    'Dist','Eucl',...
    'SW',1,...
    );

dim = size(mask);
N = numel(mask);

%the indices of the non-zero mask elements
indices=find(mask);
[alli, allj, allk]=ind2sub(dim,indices);

switch conn
    case 6
        nN = 3;
    case 10
        nN = 5;
    case 18
        nN = 9;
    case 26
        nN = 13;
    otherwise
        error('undefined 3D connectivity neighbourhood..!')
end

% coordinates of the center points
ci=repmat(alli,nN,1);
cj=repmat(allj,nN,1);
ck=repmat(allk,nN,1);

switch conn
    case 6
        ni=[alli  ; alli+1; alli  ];
        nj=[allj+1; allj  ; allj  ];
        nk=[allk  ; allk  ; allk+1];
        
    case 10 
        ni=[alli  ; alli+1; alli+1; alli+1; alli  ];
        nj=[allj+1; allj  ; allj+1; allj-1; allj  ];
        nk=[allk  ; allk  ; allk  ; allk  ; allk+1];
    case 18
        ni=[alli  ;alli+1;alli+1;alli+1;alli  ;alli  ;alli  ;alli+1;alli+1];
        nj=[allj+1;allj  ;allj-1;allj+1;allj  ;allj+1;allj+1;allj  ;allj  ];
        nk=[allk  ;allk  ;allk  ;allk  ;allk+1;allk-1;allk+1;allk-1;allk+1];
    case 26
        ni=[alli;alli+1;alli+1;alli+1;alli;alli;alli;alli+1;alli+1;alli+1;alli+1;alli+1;alli+1];
        nj=[allj+1;allj;allj-1;allj+1;allj;allj+1;allj+1;allj;allj;allj-1;allj+1;allj+1;allj-1];
        nk=[allk;allk;allk;allk;allk+1;allk-1;allk+1;allk-1;allk+1;allk-1;allk-1;allk+1;allk+1];
end

maskZ=cat(2,mask,zeros(dim(1),1,dim(3)));
maskZ=cat(1,maskZ,zeros(1,dim(2)+1,dim(3)));
maskZ=cat(3,maskZ,zeros(dim(1)+1,dim(2)+1,1));

valid=(ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2) & nk>=1 & nk<=dim(3));
ni=ni(valid);
nj=nj(valid);
nk=nk(valid);
ci=ci(valid);
cj=cj(valid);
ck=ck(valid);

tt=sub2ind(size(maskZ),ni,nj,nk);
ee=maskZ(tt);
valid=logical(ee);
ni=ni(valid);
nj=nj(valid);
nk=nk(valid);
ci=ci(valid);
cj=cj(valid);
ck=ck(valid);

cInd=sub2ind(dim,ci,cj,ck); 
nInd=sub2ind(dim,ni,nj,nk); 

switch w
    case 'yesNorm_additive'
         
        pro=1.5;
        p=(mask(cInd).*mask(nInd)*pro).^pow;
        A=sparse([cInd,nInd],[nInd,cInd],[p,p],N,N);
                
    case 'yesNorm_subtractive'
        
        pro=1;
        p=(mask(cInd).*mask(nInd)*pro).^(-pow);
        A=sparse([cInd,nInd],[nInd,cInd],[p,p],N,N);
        
    case 'yesNora'
        d=mask(cInd)-mask(nInd); 
        if isempty(s) 
            s=mean(abs(d))*sw;
        end
        
        if strcmp(dist,'Gaus')
            sim=exp(-d.^2/(2*s^2)); 
        else
            sim=1./(1+d.^2/(2*s^2)); 
        end
        A=sparse([cInd,nInd],[nInd,cInd],[sim,sim],N,N);
        
    case 'no'
        A=sparse([cInd,nInd],[nInd,cInd],ones(1,2*numel(ni)),N,N);
    case 'yesDistance'
        
        [aha1, aha2, aha3] = ind2sub(dim, cInd);
        [nababa1, nababa2, nababa3] = ind2sub(dim, nInd);
        d = sqrt(sum(([aha1-nababa1,aha2-nababa2,aha3-nababa3]).^2,2));
        
        if isempty(s) 
            s=mean(abs(d))*sw;
        end
        
        if strcmp(dist,'Gaus')     
            sim=exp(-d.^2/(2*s^2));    % gaussian similarity
        elseif strcmp(dist,'adjEucl')
            sim=1./(1+d.^2/(2*s^2));   % adjusted Euclidean similarity
        elseif strcmp(dist,'Eucl')
            sim=1./d;                  % Euclidean similarity
        end
        A=sparse([cInd,nInd],[nInd,cInd],[sim,sim],N,N);
        
end

c=find(~sum(A,1));
c=c(~ismember(c,indices));

A(:,c)=[];
A(c,:)=[];
end


