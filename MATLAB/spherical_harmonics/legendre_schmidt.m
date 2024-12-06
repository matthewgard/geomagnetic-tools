function [Pnm,dPnmdt] = legendre_schmidt(N,theta)
% Compute Pnm and dPndmdt for given theta and degree N
% SEE ALSO: legendre_ind
% Edited from design_SHA_matlab by Nils Olsen (DSRI), March 2021

% Initialise
ndeg_terms   = N*(N+3)/2+1;
Pnm   = zeros(1,ndeg_terms);
dPnmdt = zeros(1,ndeg_terms);

% Create the first few terms manually so that you can recursively solve for
% the rest n=0 to 1
Pnm(1:3) = [1 cos(theta) sin(theta)];
dPnmdt(1:3) = [0 -sin(theta) cos(theta)];

% Computed schmidt semi-normalized associated legendre polynomials and
% their derivatives
for n=2:N
    % Solve the m = n terms first
    ind_nn = legendre_ind(n,n);
    ind_nn1 = legendre_ind(n-1,n-1);
    Pnm(ind_nn) = sqrt((2*n-1)/(2*n)) * sin(theta) * Pnm(ind_nn1);
    dPnmdt(ind_nn) = sqrt((2*n-1)/(2*n)) * (sin(theta) * dPnmdt(ind_nn1) + cos(theta) * Pnm(ind_nn1));

    % Solve for m = 0 to (n-1)
    for m = 0:1:(n-1)
        ind_nm = legendre_ind(n,m);
        ind_n1m = legendre_ind(n-1,m);
        ind_n2m = legendre_ind(n-2,m);
        Pnm(ind_nm) = (((2*n)-1) * cos(theta) * Pnm(ind_n1m) - sqrt((n-1+m)*(n-1-m))*Pnm(ind_n2m))/sqrt((n+m)*(n-m));
        dPnmdt(ind_nm) = (((2*n)-1) * (cos(theta) * dPnmdt(ind_n1m) ...
            - sin(theta) * Pnm(ind_n1m)) - sqrt((n-1+m)*(n-1-m)) * dPnmdt(ind_n2m))/sqrt((n+m)*(n-m));
    end
end
return
