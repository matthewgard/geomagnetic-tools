function [Pnm,dPnmdt] = legendre_schmidt(N,theta)
% Compute Pnm and dPndmdt for given theta and degree N
% SEE ALSO: legendre_ind

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















% 
% % Computed schmidt semi-normalized associated legendre polynomials and
% % their derivatives
% for n=0:N
%     % Pnm - MATLAB already has a function, might as well use it
%     % WELL nevermind it has a lot of overhead. Write my own for this?
%     % Very slow
%     %temp_Pnm = legendre(n, cos(theta), 'sch')';
%     temp_Pnm = legendre_mine(n,theta);
%     % DO THE BELOW, BUT DO IT ALL AT ONCE. Generate pxyz using pxyznm and
%     % pnmdt using below? Or just use iaga school method but remove 
% 
% 
%     % Pnmdt - derivative wrt theta
%     % Edited from design_SHA_matlab by Nils Olsen (DSRI), March 2021
%     if n == 0
%         tempd_Pnmdt = zeros(size(temp_Pnm));
%     else
%         tempd_Pnmdt = NaN(size(temp_Pnm));
%         % m = 0
%         tempd_Pnmdt(:,1) = -temp_Pnm(:,2).*sqrt(n*(n+1)/2);
%         % m = n
%         if n == 1
%             tempd_Pnmdt(:,2) = sqrt(2)*temp_Pnm(:,n).*sqrt(n/2);
%         else
%             tempd_Pnmdt(:,n+1) =  temp_Pnm(:,n).*sqrt(n/2);
%         end
%         % m = 1, n > 1
%         if n > 1
%             tempd_Pnmdt(:,2) = (temp_Pnm(:,1).*sqrt(2*n*(n+1))-temp_Pnm(:,3).*sqrt((n+2)*(n-1)))/2;
%         end
%         % 2 < m < (n-1)
%         for m = 2:(n-1)
%             tempd_Pnmdt(:,m+1)=(temp_Pnm(:,m).*sqrt((n+m)*(n-m+1))-temp_Pnm(:,m+2).*sqrt((n+m+1)*(n-m)))/2;
%         end
%     end
% 
%     % Append to legendre vectors
%     if n == 0
%         Pnm(1,1) = temp_Pnm;
%         dPnmdt(1,1) = tempd_Pnmdt;
%     else
%         ind_start = legendre_ind(n,0);
%         ind_end = legendre_ind(n+1,0)-1;
%         Pnm(1,ind_start:ind_end) = temp_Pnm;
%         dPnmdt(1,ind_start:ind_end) = tempd_Pnmdt;
%     end
% 
% end
% return


