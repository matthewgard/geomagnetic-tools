function ind = legendre_ind(n,m)
% Indices for degree n and order m in the legendre
% SEE ALSO: legendre_schmidt
ind = (n*(n+1)/2)+m+1;
return
