% clear;
% SH = load_SH('csr06');
% [vCS,LIST] = SH.reformat('vCS');
% tt = [SH.tt];
% 
% ser = [tt',vCS(7,:)'];
% save('example_C30.mat','ser');

%% input of time series
% Only the time series of C30 is provided. The time spans from 2002.4 to 2020.8
clear;
load('example_C30.mat'); 
% ser(:,1): time in decimal years
% ser(:,2): value of C(3,0)

%% generate uniformly spaced time series
[tt1,X1] = uniform_time(ser(:,1),ser(:,2), [2002,4,2020,8]);
% tt1: equal-spaced time, from April 2002 to August 2020 (defined by the third input)
% X1: rearranged C(3,0), NaN is assigned to gaps

%% SSA-filling-a
ind_nan = isnan(X1);
id = zeros(size(tt1)); % classify observations and gaps by id
id(tt1<2017.5 & ~ind_nan) = 1; % 1: GRACE
id(tt1>2017.5 & ~ind_nan) = 2; % 2: GFO
id(tt1<2017.5 & ind_nan) = 3;  % 3: gaps within GRACE
id(tt1>2017.5 & ind_nan) = 4;  % 4: the 11-month gap & a gap within GFO

MM = 72; % Window size
KK = 10; % Maximum number of RCs to be used
[X2,verror1] = fun_SSA_filling_a(X1,id, MM, KK);
% X2: results after SSA-filling-a gaps (id = 3) are filled.
% verror1: error estimation, based on fitting residuals

%% SSA-filling-b
Mlist = 60:12:96; 
Klist = [1,2:2:18]; 
% The following code traverses Mlist & Klist to implement the cross validation to find 
% the optimal parameter set.
% If both Mlist and Klist consist of only one element, the value will be used
% directly.
[X3,verror2,opt_MK] = fun_SSA_filling_b(tt1,X2,Mlist,Klist);
% X3: final output, all gaps are filled
% verror2: error esimation, based on the cross validation (if implemented, 
%          otherwise based on fitting residuals).

%% plot
figure('position',[1,1,1028,303]);
plot(tt1,X3,'o-','color',[1,1,1]/2);
hold on; 
errorbar(tt1(id == 3),X2 (id==3), ones(sum(id==3),1)*verror1,'ro','markerfacecolor','r');
errorbar(tt1(id == 4),X3 (id==4), ones(sum(id==4),1)*verror2,'bo','markerfacecolor','b');
hold off;
legend('Final series','SSA-filling-a','SSA-filling-b','location','best');
title(sprintf('Optimal parameter: M=%d, K=%d',opt_MK));
