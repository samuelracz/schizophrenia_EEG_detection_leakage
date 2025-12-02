function [ output ] = PermEn( ts, dim, lag )
% algorithm for calculating the permutation entropy (PermEn) of time series
% at diven dimension (d) and time lag (lag)
% Reference: Bandt C, Pompe B, 2002.


if max(dim) < 3
    error('dim has to be larger than 2')
end

if size(ts,1) < size(ts,2)
    ts = ts';
end

N = length(ts);
tstart = (dim-1)*lag + 1;
tmax = N - (dim-1)*lag;
x = zeros(tmax,dim);
ps0 = factorial(dim);
ps1 = factorial(dim-1);

for d = 1:max(dim)                                                         % embedding
    x(:,d) = ts(tstart-(d-1)*lag:end-(d-1)*lag,1);
end

[~, perms0] = sort(x, 2, 'ascend');
[~, perms1] = sort(x(:,1:end-1), 2, 'ascend');



output = struct(...
    'Hd',[],...
    'Hd_norm',[],...
    'hd',[],...
    'SortEn',[],...
    'Hd_1',[],...
    'Hd_1norm',[],...
    'hd_1',[],...
    'AR0',[],...
    'AR1',[]);

[~,uniques0,symb0] = unique(perms0,'rows');                           % find appearind permutations (attach symbolic value - symb)
pi0 = histcounts(symb0,length(uniques0));                                   % count the appearence of each permutation
pi0 = pi0/length(symb0);
output.Hd = -sum( (pi0./sum(pi0)) .* log2(pi0./sum(pi0)) );
output.Hd_norm = (output.Hd)/log2(ps0);
output.hd = (output.Hd)/(max(dim)-1);

[~,uniques1,symb1] = unique(perms1,'rows');
pi1 = histcounts(symb1,length(uniques1));
pi1 = pi1/length(symb1);
output.Hd_1 = -sum( (pi1./sum(pi1)) .* log2(pi1./sum(pi1)) );
output.Hd_1norm = (output.Hd_1)/log2(ps1);
output.hd_1 = (output.Hd_1)/(max(dim)-2);

output.SortEn = output.Hd - output.Hd_1;






end

