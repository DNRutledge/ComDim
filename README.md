# ComDim
ComDim
ComDim is an unsupervised multi-block method that aims to simultaneously consider multiple data tables
to find the components that are common to all the tables and those that are specific to each data table,
with the contribution of each of the tables to each of these components.
ComDim determines a common space describing the dispersion of the samples in all the blocks,
each block having a specific weight (salience, λ) associated with each dimension in this common space.
Significant differences in the saliences for a given dimension reflect the fact
that the dimension contains different amounts of information coming from each block.
In addition to the saliences, Local loadings for each analyzed block and two different sets of scores are obtained.
The first set corresponds to the Local scores for each analyzed block
while the second set is composed of the Global scores, common to all the blocks.

%
% USAGE :
%---------------
% [res_calib]=comdim_PCA_2020(col, Options);
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
% 'Tall_Wide' not yet available
% % % % Options.Partitions = number of partition if X 'Tall' or 'Wide'
% % % % Options.rows_par = number of row partition if X 'Tall' or 'Tall_Wide'
% % % % Options.cols_par = number of columnpartition if X 'Tall_Wide' or 'Wide'
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
% ColumnsPartition
% PCA_TW_DNR
% PCA_Tall_PCT_DNR
%   which calls
%   RowsPartition
%
% REFERENCES :
%-----------------
% Method published by:
% E.M. Qannari, I. Wakeling, P. Courcoux and H. J. H. MacFie
% in Food quality and Preference 11 (2000) 151-154
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
% [res_calib]=comdim_PCA_2020(collection, Options);

