function [T] = PCA_Tall_PCT_DNR(X, PCs, partitions)
% Segments the tall X matrix into blocks and does PCA on each block
%
% INPUT :
% X = Tall data matrix
% PCs = number of PCs to extract from X
% partitions= Number of segments
%
% OUTPUT :
% T = PCA scores
%
% CALLS :
% RowsPartition
%

% NaNs are  produced in X by normalisation !
% So ...
% [X]=NaN2Median(X')';
% [X]=NaN2Median(X);

W = RowsPartition(X, partitions);
P=[];
T=[];
cols = size(X,2);
% I = eye( cols, cols );

D=W{1}'*W{1};

for par = 2:partitions
    D = D + (W{par}'*W{par});
end % partitions

[Vx Sx]=eig(D);
% Vx=flipud(Vx);
Vx=fliplr(Vx);

T = [];
for par = 1:partitions
    t = W{par}*Vx;
    T = [T ; t(:,1:PCs)];
end


