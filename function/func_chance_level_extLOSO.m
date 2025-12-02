function [ p ] = func_chance_level_extLOSO(acc,CVpartitions,K)

acc_true = mean(acc);
acc_perm = zeros(K,1);
Ncv = length(CVpartitions);
LOSO_lib = [0,0;0,1;1,0;1,1];
for k = 1:K
    acc_cv = zeros(Ncv,1);
    for cv = 1:Ncv
        y_true = [0,1];
        rand_cv = randperm(4);
        y_pred = LOSO_lib(rand_cv(1),:);
        acc_cv(cv) = mean(y_true == y_pred);
    end
    acc_perm(k) = mean(acc_cv);
end

p = (sum(acc_perm >= acc_true)+1)/(K+1);

end