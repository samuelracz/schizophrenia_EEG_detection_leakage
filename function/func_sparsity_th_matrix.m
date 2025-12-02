function [ mat_th ] = func_sparsity_th_matrix( mat, th )
% mat - matrix
% th - number of edges remaining

mat_th = zeros(size(mat));
matu = triu(mat);
matl = tril(mat);

nc = ceil(th*size(mat,1)*(size(mat,1)-1)/2);

edge_vecu = matu(:);
[Wu,Iu] = sort(edge_vecu,'descend');
mat_th(Iu(1:nc)) = Wu(1:nc);

edge_vecl = matl(:);
[Wl,Il] = sort(edge_vecl,'descend');
mat_th(Il(1:nc)) = Wl(1:nc);

end

