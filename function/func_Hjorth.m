function [C, M, A] = func_Hjorth(eeg)

if size(eeg,2) > size(eeg,1)
    eeg = eeg';
end

eeg_d1 = diff(eeg,[],1);
eeg_d2 = diff(eeg_d1,[],1);

m0 = var(eeg,[],1);
m2 = var(eeg_d1,[],1);
m4 = var(eeg_d2,[],1);

A = m0;
M = sqrt(m2./m0);
C = sqrt(m4./m2)./sqrt(m2./m0);


end