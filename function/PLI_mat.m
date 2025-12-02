function [ R_mat, PLI_mat, c_mat ] = PLI_mat( data )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% data should be organized as a matrix of size data points * channels!!!

[L, nch] = size(data);

z = hilbert(data);
x = real(z);
x_hat = imag(z);

A = sqrt(x.^2 + x_hat.^2);
fi = atan(x_hat./x);

d_Fi = reshape(fi',nch,1,L) - reshape(fi',1,nch,L);

R_mat = abs(mean(exp(1i*d_Fi),3));
R_mat(eye(nch)==1) = 0;
PLI_mat = abs(mean(sign(d_Fi),3));

num_mat = pagemtimes(reshape(A',nch,1,L),reshape(A',1,nch,L)).*exp(1i*d_Fi);
denom_mat = sqrt(mean(A.^2,1)'*mean(A.^2,1));
c_mat = abs(imag(mean(num_mat,3)./denom_mat));

end

