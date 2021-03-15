function [X_Normed, Norm_X, Mean_X]=Normalise_DB(X);

    Mean_X=mean(X);
    X_mean=X-ones(size(X,1),1)*Mean_X;

    XX=X_mean.*X_mean;
    Norm_X=sqrt(sum(XX(:)));

    X_Normed=(X_mean)/Norm_X;
    
end