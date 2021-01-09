function [X2,verror] = fun_SSA_filling_a(X,id, MM, KK)

X2 = X;
ind = id < 4;
X = X(ind);
X = X(:);

ind_nan = id == 3;

AMP = max([nanmedian( abs(X) ), 1D-14]); % rescale the values to be around 1
X = X / AMP;

[X_F, LAMBDA, RC, htest, EOF, PC] = ssa_missing_iterative(X, MM, KK);

X2(ind_nan) = X_F * AMP;

ind_nan = isnan(X);
verror =  std(X(~ind_nan)-sum(RC(~ind_nan,htest==1),2))*AMP;

end