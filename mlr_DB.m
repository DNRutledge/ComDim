function [b0,b]=mlr_DB(X,Y,intercept)
%#  function [b0,b]=mlr_DB(X,Y,intercept)				
%#								
%#  AIM:	Computation of the coefficients of a MLR model.	
%#								
%#  PRINCIPLE:  Builds a regression model, based on minimising	
%#		the sum of squared residuals of the calibration	
%#		objects.					
%#								
%#  INPUT:						
%#	X	: (row x col) Calibration matrix of col 	
%#		  independent variables				
%#	Y	:  (row x 1) Dependent variable			
%#	intercept: (optional) Set this parameter equal to 1 if	
%#		   you want to compute an intercept, and to 0	
%#		   if you do not want to compute an intercept.	
%#		   BY DEFAULT, INTERCEPT = 1 			
%#								
%#  OUTPUT:							
%#	b0	: Value of the intercept of the MLR model	
%#		  (set to 0 if no intercept was required)	
%#	b	: Vector of regression coefficients		
%#								
%#  AUTHOR:     Delphine Jouan-Rimbaud				
%#								


if nargin==2, intercept=1; end	%default value of intercept

if intercept==1,		% compute an intercept
   row=size(X,1);
   X=[ones(row,1),X];		% add a column of 1 to X,
				% the corresponding b is the intercept
   b=pinv(X'*X)*X'*Y;
   b0=b(1,:); b(1,:)=[];	% computation of b and b0
else				% no intercept
   b=pinv(X'*X)*X'*Y;
   b0=0;			% computation of b and b0
end
