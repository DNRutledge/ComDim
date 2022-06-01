function [res_calib]=comdim_ICA_2020(col, Options);
%comdim			- Finding common dimensions in multitable data (saisir format)
%
% With data compression methods INSIDE jade_DNR_2014
%
%
% USAGE :
%---------------
% [res_calib]=comdim_ICA_2020(col, Options);
%
% INPUT :
%---------------
% col : vector of saisir files (the numbers "nrow" of rows in each table must be equal) .
%
% Options for ComDim:
% Options.ndim : Number of common dimensions
% Options.normalise : 0== no; 1== yes (default)
% Options.threshold : If the "difference of fit" < threshold
% then break the iterative loop
%
% % % % To accelerate the ComDim calculations with big matrices
% % % % Options.CompMethod = 'Normal', 'Kernel', 'PCT', 'Tall_Wide', 'Tall' or 'Wide'
% % % % Options.Partitions = number of partition if X 'Tall' or 'Wide'
% % % % Options.rows_par = number of row partition if X 'Tall' or 'Tall_Wide'
% % % % Options.cols_par = number of columnpartition if X 'Tall_Wide' or 'Wide'
%
% Options for ICA:
% To accelerate the ICA calculations with big matrices
% Options.Data= 'Data', 'Loadings'
%
% Compression is done INSIDE jade_DNR_2014
% For compatability 
% Options.Method=Options.CompMethod; 
%
% Options.Method = 'Normal' , 'Tall', 'Wide', 'Kernel', 'Tall_Wide'
% Options.Partitions = number of partition if X 'Tall' or 'Wide'
% Options.rows_par = number of row partition if X 'Tall' or 'Tall_Wide'
% Options.cols_par = number of columnpartition if X 'Tall_Wide' or 'Wide'
%
% Options.SortICs =0; % Don't sort ICs by variance
%
% Options.loquace : 1== displays calculation times; 0== no (default)
%
% Options.Output : To limit the number of extra parameters calculated
% To save time and memory
% 'Q' (Global Scores) and Statistice will always be output
% if Output=='A' or 'TPL' , all extra parameters are calculated
% if Output=='[]', no extra parameters are calculated
% 'T' : Local Scores
% 'P' : Scaled Loadings
% 'L' : Unscaled Loadings

% Options.FileName : To save results for successive CCs to files
% if FileName=='Toto', then results saved to 'Toto_01', Toto_02', etc
%
%
% OUTPUT :
%-----------------
% res_calib with fields:
% Q : Global scores_calib (nrow x ndim)
% P : Scaled Global loadings calculated from Q and tables (nvars x ndim)
% P_Loc : Structure with ntable of Scaled Local loadings in matrices of (local vars x ndim)
% T : Structure with ntable of Local scores_calib calculated from local Loadings and tables (nrow x ndim)
% Lx : Unscaled Global loadings for X
% Lx_Loc : Unscaled Local loadings for X
% saliences : weight of the original tables in each dimension (ntable x ndim)
% Sum_saliences_Tab - Sum of Saliences for each Tab in a Dimension
% Sum_saliences_Dim - Sum of Saliences for each Dim for a Table
% explained : percentage explanation given by each dimension (1 x ndim)
%
% NORM : ntable norms of the tables (1 x 1)
% MEAN : ntable means of tables (1 x nvars)
%
%
% SingVal : vector of singular values (1 x ndim)
%
%
% CALLS:
%-----------------
% mlr_DB
% Normalise_DB
% Compress_Data_2020
% ColumnsPartition
% PCA_TW_DNR
% PCA_Tall_PCT_DNR
%   which calls
%   RowsPartition
%
% jade_DNR_2014
%   which calls
%   jadeR_2005
%   IC_sort
%
% REFERENCES :
%-----------------
% Method published by:
% E.M. Qannari, I. Wakeling, P. Courcoux and H. J. H. MacFie
% in Food quality and Preference 11 (2000) 151-154
%
% D. Juan-Rimbaud Bouveresse, D. N. Rutledge
% Independent Components Analysis with the Jade Algorithm
% TrAC., (2013), 50, 22–32
% 
% D. N. Rutledge and D. Juan-Rimbaud Bouveresse
% Corrigendum to "Independent Components Analysis with the Jade Algorithm"
% TrAC., (2015), 67, 220

%
% Example application
% Chemometric Tools to Highlight Possible Migration of Compounds
% from Packaging to Sunflower Oils
% J. Maalouly, N. Hayeck, A. Kassouf,D. N. Rutledge and V. Ducruet
% J. Agric. Food Chem. 2013, 61, 10565?10573


% EXAMPLE :
%-----------------
% (suppose 3 SAISIR matrices "spectra1","spectra2","spectra3")
% collection(1)=spectra1; collection(2)=spectra2; collection(3)=spectra3
% [res_calib]=comdim_ICA(collection, Options);

tic
ntable=size(col,2);
[nR,nC]=size(col(1).d);

if(nargin<2)
    ndim=ntable;
    threshold=1E-10;
    FileName=[];
end

if exist('Options','var')
    if isfield(Options,'ndim')
        ndim= Options.ndim;
    else
        ndim=ntable;
    end

    if isfield(Options,'normalise')
        normalise= Options.normalise;
    else
        Options.normalise=1;
        normalise= Options.normalise;
    end

    if isfield(Options,'threshold')
        threshold= Options.threshold;
        if size(threshold)==0
            threshold=1E-10;
            Options.threshold=threshold;
        end
    else
        threshold=1E-10;
        Options.threshold=threshold;
    end

    if isfield(Options,'Output')
        Output= Options.Output;
    else
        Options.Output= 'TPL';
        Output= Options.Output;
    end
    if Output== 'A'
        Options.Output= 'TPL';
        Output= Options.Output;
    end

    if isfield(Options,'FileName')
        FileName= Options.FileName;
    else
        Options.FileName= [];
        FileName= Options.FileName;
    end

    if isfield(Options,'loquace')
        loquace= Options.loquace;
    else
        Options.loquace= 0;
        loquace= Options.loquace;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(Options,'Test')
        Test= Options.Test;
    else
        Options.Test=0;
        Test= Options.Test;
    end
    
% For compatability with ICA
% Options.Method=Options.CompMethod; 
    if isfield(Options,'CompMethod')
        Options.Method=Options.CompMethod;
    end
end

DimLabels=[repmat('D',ndim,1),num2str([1:ndim]'),repmat(' ',ndim,6)];
DimLabels=DimLabels(:,1:6);

TableLabels=[repmat('t',ntable,1),num2str([1:ntable]'),repmat(' ',ntable,6)];
TableLabels=TableLabels(:,1:6);

isT=((strfind(Output,'T')));
isL=((strfind(Output,'L')));
isP=((strfind(Output,'P')));
isLP=[isP,isL];
isTLP=[isT,isP,isL];

NoFile=isempty(FileName);

nrow=size(col(1,1).d,1);

duree=toc;
if loquace==1
    disp(['Initialisation finished after ',num2str(duree)]);
end

%%
% Do Normalisation of blocks
% s{i} --->>> s_n{i}

tic
X_mat=[];
Xnorm_mat=[];

for i=1:ntable
    if normalise==1
        % Normalise original blocks
        [X_Normed, Norm.d(i), Mean(i).d]=Normalise_DB(col(i).d); % MS
    
        res_calib.MEAN(i).d= Mean(i).d';
        res_calib.MEAN(i).i= col(i).v';
        res_calib.MEAN(i).v= TableLabels(i,:);
        
        temp_tabCalib{i}.d=X_Normed;
        s_n{i}=X_Normed;

    else
        % Do nothing
        Norm.d(i)= 1;
        
        res_calib.MEAN(i).d= zeros(size(col(i).d,2),1);
        res_calib.MEAN(i).i= col(i).v';
        res_calib.MEAN(i).v= TableLabels(i,:);
        
        temp_tabCalib{i}.d=col(i).d;
        s_n{i}=col(i).d;
    end
    
    % Difficult to preallocate (Total number of variables)
    X_mat=[X_mat,col(i).d]; % X blocs concaténés
    Xnorm_mat=[Xnorm_mat,temp_tabCalib{i}.d]; % tab blocs concaténés
    
end

res_calib.NORM.d= Norm.d;
res_calib.NORM.i= 'Norm';
res_calib.NORM.v= TableLabels;

clear Norm;

[nR,nC]=size(Xnorm_mat);

duree=toc;
if loquace==1
    disp(['Normalisation finished after ',num2str(duree)]);
end
%%%%%%%%%%%

%% Do NO Compression
% Compression is done INSIDE jade_DNR_2014
% For compatability 
% Options.Method=Options.CompMethod; 

% [s_r]=Compress_Data_2020(s_n, Options);
[s_r]=s_n;

%% Do ComDim-ICA
h1 = waitbar(0,'Doing ComDim...');

% main loop
unexplained=ntable; % Warning!: true only if the tables were set to norm=1

% %%%%%%%%%%%%%%%%%%%%%%%
ICs=ndim;
% Start with ndim ICs
% and decrease by 1 each time
% or use ndim ICs each time
% %%%%%%%%%%%%%%%%%%%%%%%

% Do ComDim
tic
for dim=1:ndim
    waitbar(dim/ndim,h1,['CC : ', sprintf('%0.1u',dim) '/', sprintf('%0.1u',ndim)]);

    previousfit=10000;
    deltafit=1000000;
    lambda=ones(ntable,1);

%     qini=Xs{1,1}(:,1);
    qini=s_r{1,1}(:,1);
    qini=qini/sqrt(qini'*qini); % random initialization of t + set to unit length
    qold=100;

    iters=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    q_Mat{dim}=qini;
    while and(norm(qold-qini)>threshold, iters<100)
        iters=iters+1;
        qold=qini;
        q=0;

        W=[]; % according to your code
        for j=1:ntable
            W = [W, sqrt(lambda(j)).*s_r{j}];
        end

        %%%%%% ICA
%         ICs=ndim;
        [Signal,Scores] =  jade_DNR_2014(W,ICs, Options);  % Compression is done INSIDE jade_DNR_2014

        sv=sqrt(Scores(:,1)'*Scores(:,1));
        qtmp=Scores(:,1)/sv;
    
        % takes each table into account for lambda after ICA
        for j=1:ntable
%             lambda(j)=qtmp'*(s_r{j}*s_r{j}')*qtmp;
            % SLOWER but takes up less memory
            toto=qtmp'*s_r{j};
            titi=s_r{j}'*qtmp;
            lambda(j)=toto*titi;
            
            q=q+(lambda(j)*qtmp); 
        end
            
        q=q./sqrt(q'*q); % standardizes t
        if abs(min(q))>abs(max(q))
            q=-q;
        end
        qini=q; % updates the initialization of t and so on
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Test==1
            q_Mat{dim}=[q_Mat{dim},qini];
        end
        
    end %deltafit>threshold

    saliences.d(:,dim)=lambda;
    Q.d(:,dim)=q;

    %%%%%
    % 'SingVal' is square of the SVD singular value
    % because here it is based on X, not X.X'
    % Does NOT work for Kernel & Wide & Tall
    res_calib.SingVal.d(1,dim)=sv(1,1)*sv(1,1); % from SVD on W
    res_calib.SingVal.i='Singular value';
    res_calib.SingVal.v=DimLabels; % Dimensions

    varexp(:,dim)=res_calib.SingVal.d(1,dim)*res_calib.SingVal.d(1,dim);

    % Deflate blocks
    nrow=size(q,1);
    aux=eye(nrow,nrow)-q*q';
    for j=1:ntable
        s_r{j}=aux*s_r{j};
    end
    
% %%%%%%%%%%%%%%%%%%%%%%%
	ICs=ICs-1;
% Start with ndim ICs
% and decrease by 1 each time
% %%%%%%%%%%%%%%%%%%%%%%%

end

Q.i=col(1).i;
Q.v=DimLabels;
res_calib.Q=Q;

%%
duree=toc;
if loquace==1
    disp(['Scores finished after ',num2str(duree)]);
end

close(h1);

clear s_n;

%%
% explained.i='% explained';
% explained.v=DimLabels; % Dimensions
% res_calib.explained=explained;
% clear explained;

%%
explained.d=varexp/sum(varexp)*100;
explained.i='% explained';
explained.v=DimLabels; % Dimensions
res_calib.explained=explained;
clear varexp explained;

% Calculate Sums of saliences for each Dim
Sum_saliences_Dim.d=sum(saliences.d);
Sum_saliences_Dim.i='Sum Dim Saliences';
Sum_saliences_Dim.v=DimLabels; % Dimensions

res_calib.Sum_saliences_Dim=Sum_saliences_Dim;
clear Sum_saliences_Dim;

% Calculate Sums of saliences for each Table
Sum_saliences_Tab.d=sum(saliences.d,2)';
Sum_saliences_Tab.i='Sum Tab Saliences';
Sum_saliences_Tab.v=TableLabels;% tables

res_calib.Sum_saliences_Tab=Sum_saliences_Tab;
clear Sum_saliences_Tab;

saliences.i=TableLabels;% tables
saliences.v=DimLabels; % Dimensions

res_calib.saliences=saliences;

%% Calculate Normalised concatenated Xs ('Calib') from col
% Calculate concatenated CD loadings
% Reorganise Loadings - 1 matrix / LV

%%%%%% If Output==[], nothing else is calculated
tic;
if ~isempty(Output)
    clear L_CD_Vec Calib;
    L_CD_Vec=[];
    L_X_Vec=[];
    Calib=[];
 
    nCalib= size(col(1,1).d,1);
end

%%
isT=((strfind(Output,'T')));
isL=((strfind(Output,'L')));
isP=((strfind(Output,'P')));
isTLP=[isT,isP,isL];

% Calculate concatenated CD loadings
% Reorganise Loadings - 1 matrix / LV
for j=1:ndim
    T_mat=[];

    for i=1:ntable
        % Q.d are orthonormal in ICA & PCA
        temp=temp_tabCalib{i}.d'*Q.d(:,j); % Scaled CD Loadings "locaux" calculés

        if  ~isempty(isP)
            L_CD{i}(:,j)=temp;
        end

        if  ~isempty(isL)
            % Q.d are orthonormal in ICA & PCA
            L_X{i}(:,j)=col(i).d'*Q.d(:,j); % Unscaled X Loadings "locaux" calculés;
        end
        
%         T_mat=[T_mat,temp_tabCalib{i}.d*temp/(temp'*temp)]; % Scores "locaux"];
        T_Loc=temp_tabCalib{i}.d*temp/(temp'*temp);
        T_mat=[T_mat,T_Loc]; % Scores "locaux"];

        if  ~isempty(isT)
%             T{i}(:,j)=temp_tabCalib{i}.d*temp/(temp'*temp); % Scores "locaux"
            T{i}(:,j)=T_Loc; % Scores "locaux"
        end

        % Deflate each temp_tabCalib
        temp_tabCalib{i}.d=temp_tabCalib{i}.d-Q.d(:,j)*temp';
end

    %     % For each CC
    %     % MLR b-coefficients between Local and Global Scores
    [b0,b(:,j)]=mlr_DB(T_mat,Q.d(:,j),0);
    
%     for i=1:size(T_mat,2)
%         [b0,b(i,j)]=mlr_DB(T_mat(:,i),Q.d(:,j),0);
%         i
%         j
%         b(i,j)
%     end

end

%%%%%% If Output==[], nothing else is calculated
if ~isempty(Output)
    
    % Calculate Global Loadings
    if  ~isempty(isP)
        L_CD_Vec=Xnorm_mat'*Q.d; % Scaled CD Loadings "globaux"
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Test==1
            P_Mat={};
            for dim=1:ndim
                P_Mat{dim}=Xnorm_mat'*q_Mat{dim};
            end
        end
        
    end

    if  ~isempty(isL)
        L_X_Vec=X_mat'*Q.d; % Unscaled X Loadings "globaux"
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if Test==1
            Lx_Mat={};
            for dim=1:ndim
                Lx_Mat{dim}=X_mat'*q_Mat{dim};
            end
        end
        
    end
end

%%%%%% If Output==[], nothing else is calculated
duree=toc;
if loquace==1
    disp(['Loadings finished after ',num2str(duree)]);
end

clear X_mat temp_tabCalib temp;


%% Output
tic;

res_calib.b.i=TableLabels;% tables
res_calib.b.v=DimLabels; % Dimensions
res_calib.b.d=b;

if  ~isempty(isT)
    res_calib.T.i=res_calib.Q.i;% samples
    res_calib.T.v=DimLabels; % Dimensions
    res_calib.T.d=T;
end


% T is no longer needed 
clear T T_mat;

if  ~isempty(isP)
    res_calib.P.i=[1:nC]';% numbers of all variables;
    res_calib.P.v=DimLabels; % Dimensions
    res_calib.P.d= L_CD_Vec;

    res_calib.P_Loc.i=TableLabels;
    res_calib.P_Loc.v=DimLabels; % Dimensions
    res_calib.P_Loc.d= L_CD;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Test==1
        res_calib.P_Mat=P_Mat;
        clear P_Mat;
    end
    
end
clear L_CD_Vec;
clear L_CD;

if  ~isempty(isL)
    res_calib.Lx.i=[1:nC]';% numbers of all variables;
    res_calib.Lx.v=DimLabels; % Dimensions
    res_calib.Lx.d= L_X_Vec;

    res_calib.Lx_Loc.i=TableLabels; % Tables
    res_calib.Lx_Loc.v=DimLabels; % Dimensions
    res_calib.Lx_Loc.d=L_X;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Test==1
        res_calib.Lx_Mat=Lx_Mat;
        clear Lx_Mat;
    end
    
end
clear L_X_Vec;
clear L_X;

% Singular value calculated from b and saliences 
% Since :
%     b_T_Q(:,j)=saliences.d(:,j).*saliences.d(:,j)/SingVal.d(1,j);
for j=1:ndim
    SingVal.d(:,j)=b(:,j)'*(saliences.d(:,j).^2)/(b(:,j)'*b(:,j));
end
SingVal.i='Singular value';
SingVal.v=DimLabels;
res_calib.SingVal=SingVal; % DB, 13.03.2014
clear SingVal;

clear saliences;

%%%%%% If Output==[], nothing else is calculated
if  ~isempty(Output)
    if normalise==1
        res_calib.Xnorm_mat.i=res_calib.Q.i;% samples
        res_calib.Xnorm_mat.v=[1:nC]';% numbers of all variables;
        res_calib.Xnorm_mat.d=Xnorm_mat;
    end
end
clear Xnorm_mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Test==1
    res_calib.q_Mat=q_Mat;
    clear q_Mat;
end

%%%%%% If Output==[], nothing else is calculated

clear Q;


duree=toc;
if loquace==1
    disp(['Output finished after ',num2str(duree)]);
end

if NoFile==0
    tic;

    s=['SaveFile= ','''',FileName,'_D.mat',''';'];
    eval(s);
    save(SaveFile, 'res_calib', '-mat');

    duree=toc;
    if loquace==1
        disp(['File saved after ',num2str(duree)]);
    end
end

