function [X_F, LAMBDA, RC, htest, EOF, PC] = ssa_missing_iterative(X, MM, KK)
% {{SSA; series}}

% X: time series, gap filled by NaN
% MM: window length
% KK: number of RC, KK <= MM

if nargin == 0
    fun_ex1;
    return;
end

ithreshold = 1D-2; % when the iteration stops

% tt = 1:numel(X);

if numel(KK)>1
    error('K should be a single integer');
end

if numel(MM) > 1
    error('M should be a single integer');
end

if KK > MM
    error('K = %d should < MM = %d\n',KK,MM);
end

X = X(:);
N = numel(X);

ind_nan = isnan(X);
N1 = N - MM+1;

if sum(ind_nan) > 0
    max_iteration = 10000; % max number of iteration
    iK = 1;
else
    max_iteration = 1;
    iK = KK;
end

X1 = zeros(size(X));

% Y = NaN(N1,MM);
EOF = NaN(MM, MM);
PC = NaN(N1, MM);
RC = NaN(N,iK);

% colors = mycolors(MM);

% figure('position',[1,41,1280,683]);
for iter = 1:max_iteration
    
    v_diff = rms(X(ind_nan) - X1(ind_nan))/rms(X(ind_nan));
    
    %     fprintf('iter = %d, K = %d: rms(diff) = %.5e\n',iter,iK,v_diff);
%     X_save = X;
    X(ind_nan) = X1(ind_nan);
    
    if v_diff < ithreshold
        if iK == KK
            break;
        else
            iK = iK+1;
        end
    end
    
%     Y_save = Y;
    Y = zeros(N1,MM);
    for m = 1:MM
        Y(:,m) = X((1:N1) + m-1);
%         Y_nan(:,m) = ind_nan((1:N1)+m-1);
    end
    
%     EOF_save =EOF;
%     PC_save = PC;
    
    [U,S,V] = svd(Y);
    LAMBDA = diag(S.^2);
    EOF = V;
    PC = U*S;
    
    %% Calculate reconstructed components RC
%     RC_save = RC;
    RC = zeros(N,iK);
    for m = 1:iK
        
        buf = PC(:,m)*EOF(:,m)'; % invert projection
        buf = flipud(buf);
        for n = 1:N % anti-diagonal averaging
            RC(n,m)=mean( diag(buf,-N1+n) );
        end
    end
    
    X1 = sum(RC(:,1:iK),2);
    
end

X_F = X1(ind_nan);

htest(KK) = 0;
for ii = 1:KK
    [b,~] = periodogram( EOF(:,ii) );
    cum_pdf = cumsum(b)/max(cumsum(b));
    htest(ii) = cum_pdf(round(length(b)/2))>0.9;
end

end

%%
function fun_ex1()

M = 40;    % window length = embedding dimension
K = 10;


% Create time series X
% First of all, we generate a time series, a sine function of length N with
% observational white noise

% rng('default');

load('csr06_n20M40K20_200204-202008_ssa-filling-a.mat');
X = vCS2(1,:);
t = 1:numel(X);

ind_nan = isnan(X);

[X_F, lambda, RC, htest ] = ssa_missing_iterative(X, M, K);

% figure;
subplot(2,2,1)
a = lambda(1:min([20,numel(lambda)])); a = a/sum(a);
x = 1:numel(a);
yyaxis left;
plot(x,a,'o-');
yyaxis right;
plot(x,cumsum(a),'o-');
title('eigenvalues LAMBDA');

subplot(2,2,2)
for ii = 1:K
    plot(t,RC(:,ii));
    hold on
end
hold off;
legend;
title('RC')

subplot(2,2,3)
plot(t,X,'o-','color',[1,1,1]*0.6);
hold on;
plot(t,sum(RC(:,htest(1:K) == 1),2),'k-');
plot(t(ind_nan),X_F,'rx');


legend('Original time series','Reconstruction with RCs','Filled value');

hold off;

end