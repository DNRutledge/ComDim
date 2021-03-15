function [T Vq] = PCA_TW_DNR(X, PCs, rows_par, cols_par, fullRank)
% 28/1/2010
% INPUT :
%      X        : the xdata
%      PCs      : number of PCs to extract
%      rows_par : number of partitions along rows
%      cols_par : number of partitions along columns
%      fullRank : if=0 -> full Rank; if>0 -> number of PCs in PCT

% Segment X into blocks, Xi
Xi = ColumnsPartition(X, cols_par);
Tx = cell(1, cols_par);
Px = cell(1, cols_par);

[rows cols] = size( Xi{1} );
if ( fullRank == 0 )
    fullRank = min(rows,cols);
end

% for i = 1:cols_par
%     %         disp(['Segment : ' num2str(i) ]);
%             [Tx{i}] = PCA_Tall_PCT_DNR( Xi{i}, fullRank, rows_par);
% end

for i = 1:cols_par
    Yi = ColumnsPartition(Xi{i}', rows_par);
    T_all=[];
    for p = 1:rows_par
        [U_in,S_in,V_in]= svd(Yi{p}', 'econ');
        T_in=Xi{i}*V_in;
        T_all=[T_all T_in];
    end
    %% Concatenated Scores of Segments
    [u,sv,v]= svd(T_all', 'econ');
    
%     Tx{i}=u(:,1:fullRank);
    Tx{i}=u;

end


Q = Tx{1}*Tx{1}';
for i = 2:cols_par
    Q = Q + Tx{i}*Tx{i}';
end

%     [Uq Sq Vq] = svds(Q,PCs); % next point to optimize
%     T = Uq * sqrt(Sq);
T = NIPALS_D( Q, PCs);

invT = inv(T'*T);
for i = 1:cols_par
    V{i}=(invT*(T'*Xi{i}))';
end

Vq=[];
for i = 1:cols_par
    Vq=[Vq; V{i}];
end

