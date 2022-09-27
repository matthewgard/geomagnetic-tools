function sshc = data2sh(theta,phi,r,a,N,Br,Bt,Bp,varargin)
% data2sh(theta,phi,r,a,N,Br,Bt,Bp) computes the schmidt semi-normalized
% coefficients given a set of colatitude, longitude, radius and geomagnetic
% field values using least squares (lsqr)
%
% sshc = data2sh(theta,phi,r,a,N,Br,Bt,Bp)
%
%   Inputs:
%     -theta: Co-latitude of points to use in calculation of sshc (degrees)
%     -phi:   Longitude of points to use in calculation of sshc (degrees)
%     -r:     Radius of points i.e. altitude + radius of the reference
%             sphere (in km)
%     -a:     Radius of reference sphere e.g. radius of Earth (in km)
%     -N:     1x2 vector containing the minimum and maximum spherical 
%             harmonic degree and order to use e.g. [16 90] or [0 13]
%     -Br:    Radial component of the magnetic field (nT)
%     -Bt:    N component of the magnetic field (nT)
%     -Bp:    E component of the magnetic field (nT)
%   Optional inputs:
%     -tol:   Tolerance for least squares solution (default 1e-6)
%             Must be a positive number
%     -iter:  Maximum iterations for least squares solution (default 1000)
%             Must be a positive integer
%   theta, phi, Br, Bt and Bp should be vector arrays of numeric data of
%   equal size (single row/column)
% 
%   Outputs:
%       -sshc:      Schmidt semi-normalised spherical harmonic coefficients
%                   Single column of g/h pairs
% 
% SEE ALSO: data2sh, lsqr
%
% Dr. Matthew Gard, 2022

% Pre-define sshc in case of early return due to an error
sshc = NaN;
max_N_permitted = 300;
tol = 1e-6;
iter = 1000;


% Check inputs
%--------------
if any(size(theta)~=size(phi)) || any(size(phi)~=size(Br)) || ...
        any(size(Br)~=size(Bt)) || any(size(Bt)~=size(Bp))
    fprintf('Theta, phi, Br, Bt and Bp must all be the same size\n')
    return
end
[rowN,colN]=size(N);
if (rowN ~= 1 && colN ~= 1) || (max([rowN,colN]) > 2)
    fprintf('N must be a 1x2 or 1x1 vector of positive integers\n')
    return
end
% If N is a single number, append with a 0 for minimum N
if max(size(N)) == 1
    N = [0 N];
end
% Check valid degree and order - N=300 arbitrarily picked
if (min(N) < 0) || (max(N) >= max_N_permitted) || any(mod(N,1) ~= 0) || max(size(N)) > 2
    error('Invalid maximum degree N: Must be an integer between 1 and %d',max_N_permitted)
end
% Check theta and phi values
if any((theta < 0) | (theta > 180)) || any(phi < -180)
    error('Theta/phi value error; 0 <= theta <= 180, -180 <= phi')
end
if any(phi > 360)
    warning('Phi was found to have values > 180. Wrapping values to -180 180')
    phi = wrapTo180(phi);
end

% Optional inputs for tol and iter
if length(varargin) > 2
    warning('Only two optional arguments can be taken, but more than two have been specified')
elseif length(varargin)==2
    if isnumeric(varargin{1})
        % Iter must be positive
        tol = abs(varargin{1});
    else
        warning('Optional argument tol must be a number.')
    end
    if isnumeric(varargin{2})
        % Iter must be a positive integer, round it down
        iter = floor(abs(varargin{2}));
    else
        warning('Optional argument iter must be a number.')
    end
elseif length(varargin)==1
    if isnumeric(varargin{1})
        tol = abs(varargin{1});
    else
        warning('Optional argument tol must be a number.')
    end
end




% Processing
%------------
% Convert degrees to radians
% I could just ask for this as the input, but most users are probably
% looking/wanting for degrees
theta = deg2rad(theta);
phi = deg2rad(phi);
maxN = max(N);
minN = min(N);

% Preallocate requested coefficients vector and terms
% Coefficients in the form:
% n = 0   1       2           3 
%     [g] [g g h] [g g h g h] [g g h g h g h]
coeffs_length = (maxN+1)^2 - (minN)^2;

% Preallocate the rest of the 'terms' in the SH equation
brterms = zeros(1,coeffs_length);
btterms = zeros(1,coeffs_length);
bpterms = zeros(1,coeffs_length);

% Merge mag data into single vector for mldivide, sets of Br, Bt, Bp next
% to each other i.e.
% Br(1), Bt(1), Bp(1), Br(2), Bt(2), Bp(2)
mag_data = [Br(:) Bt(:) Bp(:)]';
mag_data = mag_data(:);

% Preallocate the combined terms vector
terms = zeros([length(mag_data),coeffs_length]);

% Create the (a/r) terms
a_r = (a/r).^((minN:1:maxN)+2);

% Loop through each data point
ind_terms = 1;
progress_updates = floor((0.1:0.1:1)*length(theta));
fprintf('data2sh - Looping through location data:\n')
for i = 1:length(theta)
    % Get the schmidt semi-normalized associated legendre polynomials and
    % their derivatives
    % This runs quickly, haven't bothered to subset from minN to maxN
    % Do it at the end with indices instead
    [Pnm,dPnmdt] = legendre_schmidt(maxN,theta(i));

    % Loop through degree
    for n = minN:1:maxN
        % Solve for m = 0
        ind_legendre = legendre_ind(n,0);
        ind_sshc = n^2 + 1 - minN^2; % h terms not required for m = 0
        brterms(ind_sshc) = (n+1) * a_r(n+1-min(N)) * Pnm(ind_legendre);
        btterms(ind_sshc) = -a_r(n+1-min(N)) * dPnmdt(ind_legendre);
        bpterms(ind_sshc) = 0;
        
        % Loop through order m = 1:N
        for m = 1:n
            ind_legendre = legendre_ind(n,m);
            ind_sshc = n^2+2*m - minN^2;

            % Compute g and h pairs
            brterms(ind_sshc) = (n+1) * a_r(n+1-min(N)) * cos(m*phi(i)) * Pnm(ind_legendre);
            brterms(ind_sshc+1) = (n+1) * a_r(n+1-min(N)) * sin(m*phi(i)) * Pnm(ind_legendre);

            btterms(ind_sshc) = -a_r(n+1-min(N)) * cos(m*phi(i)) * dPnmdt(ind_legendre);
            btterms(ind_sshc+1) = -a_r(n+1-min(N)) * sin(m*phi(i)) * dPnmdt(ind_legendre);

            if (theta(i)) == (0) || (theta(i)) == (pi())
                bpterms(ind_sshc) = -cos(theta(i)) * a_r(n+1-min(N)) * (-sin(m*phi(i))) * dPnmdt(ind_legendre);
                bpterms(ind_sshc+1) = -cos(theta(i)) * a_r(n+1-min(N)) * (cos(m*phi(i))) * dPnmdt(ind_legendre);
            else
                bpterms(ind_sshc) = -(m/sin(theta(i))) * a_r(n+1-min(N)) * (-sin(m*phi(i))) * Pnm(ind_legendre);
                bpterms(ind_sshc+1) = -(m/sin(theta(i))) * a_r(n+1-min(N)) * (cos(m*phi(i))) * Pnm(ind_legendre);
            end
        end
    end
    terms(ind_terms:(ind_terms+2),:) = [brterms;btterms;bpterms];
    ind_terms = ind_terms + 3;

    if any(i == progress_updates)
        fprintf('%s%%...',num2str(floor(i/length(theta)*100)))
    end
end
fprintf('\n')

% mldivide is real slow for higher sizes
% sshc = terms\mag_data;
% Decided to use lsqr with customisable tol and iter settings 
[sshc,flag,relres,iter,resvec,lsvec] = lsqr(terms,mag_data,tol,iter);
fprintf('lsqr converged at iteration %d to a solution with relative residual %e and least-squares residual %e\n',[iter,relres,lsvec(end)])

% Pad result with zeros if subset was selected
sshc = [zeros([(minN)^2 1]);sshc];

return