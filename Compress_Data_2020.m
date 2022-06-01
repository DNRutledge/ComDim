function [s_r]=Compress_Data_2020(s_n, Options);
% Do Compression

tic

h1 = waitbar(0,'Doing Compression...');

ntable=size(s_n,2);
[nR,nC]=size(s_n(1));

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
    
end

%%%%
if exist('Options','var')
    if isfield(Options,'CompMethod')
        switch Options.CompMethod
            case 'Tall'
                if isfield(Options,'Partitions')
                    partitions=Options.Partitions;
                else
                    Options.Partitions=1;
                    partitions=Options.Partitions;
                end
                for i=1:ntable
                    MaxRank=min(size(s_n{i}));
                    %% Do Tall segmented PCT
                    
                    s_r{i} = PCA_Tall_PCT_DNR(s_n{i}, MaxRank, partitions);
                end
                
            case 'Wide'
                if isfield(Options,'Partitions')
                    partitions=Options.Partitions;
                else
                    Options.Partitions=1;
                    partitions=Options.Partitions;
                end
                
                for i=1:ntable
                    %% Concatenate Scores of Segments
                    Xi = ColumnsPartition(s_n{i}, partitions);
                    T_all=[];
                    for p = 1:partitions
                        [U_in,S_in,V_in]= svd(Xi{p}, 'econ');
                        T_in=U_in*S_in;
                        T_all=[T_all T_in];
                    end
                    s_r{i}=T_all;
                end
                
            case 'Kernel'
                for i=1:ntable
                    Kernel=(s_n{i}*s_n{i}')/nR;
                    
                    [U_in,S_in,V_in]= svd(Kernel,'econ')	;
                    s_r{i}=U_in*S_in;
                end
                
            case 'PCT'
                for i=1:ntable
                    [U_in,S_in,V_in]= svd(s_n{i},'econ')	;
                    s_r{i}=U_in*S_in;
                end
                
            case 'Normal'
                for i=1:ntable
                    s_r{i} =s_n{i};
                end
        end
    else
        Options.CompMethod='Normal';
        Options.Method='Normal';
        for i=1:ntable
            s_r{i} =s_n{i};
        end
    end
    
else
    Options.CompMethod='Normal';
    Options.Method='Normal';
    for i=1:ntable
        s_r{i} =s_n{i};
    end
end

duree=toc;
if loquace==1
    disp(['Compression finished after ',num2str(duree)]);
end

close(h1);

