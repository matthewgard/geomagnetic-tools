function [Br,Bt,Bp] = sh2data(sshc,theta,phi,r,a,N)
% sh2data(sshc,theta,phi,r,a,N) computes the magnetic field values at
% specified co-latitude, longitude and radius given a set of spherical
% harmonic coefficients (Schmidt semi-normalized) in column format
%
% [Br,Bt,Bp] = sh2data(sshc,theta,phi,r,a,N)
%
%   Inputs:
%       -sshc:      Schmidt semi-normalised spherical harmonic coefficients
%                   Single column of g/h pairs (see reformatSHC if required)
%       -theta:     Co-latitude of points to calculate the magnetic field (degrees)
%       -phi:       Longitude of points to calculate the magnetic field (degrees)
%       -r:         Radius of points to calculate the magnetic field at i.e.
%                   altitude + radius of the reference sphere (in km)
%       -a:         Radius of reference sphere e.g. radius of Earth (in km)
%       -N:         1x2 vector containing the minimum and maximum spherical
%                   harmonic degree and order to use e.g. [16 90] or [0 13]
%   theta and phi should be vector arrays of numeric data of equal size
% 
%   Outputs:
%       -Br:   Radial component of the magnetic field (nT)
%       -Bt:   N component of the magnetic field (nT)
%       -Bp:   E component of the magnetic field (nT)
%   Br, Bt and Bp will be the same size as the input arrays theta/phi/r
% 
% SEE ALSO: data2sh, reformatSHC, legendre_schmidt
%
% Dr. Matthew Gard, 2022

% Arbitrary limitation on N just to limit processing
max_N_permitted = 300;

% Preallocate outputs
%--------------------
Br = NaN(size(theta));
Bt = NaN(size(theta));
Bp = NaN(size(theta));

% Check inputs
%--------------
% Check dimensions theta/phi
if any(~(size(theta)==size(phi)))
    error('Unequal dimensions of lat and long')
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
% Check size of sshc is larger than N selected at a minimum and also a
% single column vector
if size(sshc,2)~=1
    error('SH coefficients are not a single column vector of g/h')
elseif ((max(N)+1)^2) > length(sshc)
    error('Maximum degree (N) selected exceeds length of SH coefficients provided')
end
% Check that the size matches sshc with a monopole, if it doesnt check if
% adding it in satisfies it
sshc_sizes = ((0:1:max_N_permitted)+1).^2;
ind = any(length(sshc)==sshc_sizes);
if ~ind
    ind2 = any((length(sshc)+1)==sshc_sizes);
    if ind2
        warning('SH coefficients provided do not include the monopole term. Adding a 0 monopole term for processing.')
        sshc = [0;sshc];
    else
        error('SH coefficients size mismatch with expected size input.')
    end
end
clear sshc_sizes


% Processing
%------------
% (a/r)^(N+2) terms for all n = minN:maxN
a_r = (a/r).^((min(N):1:max(N))+2);

% Convert degrees to radians
% I could just ask for this as the input, but most users are probably
% looking/wanting for degrees
theta = deg2rad(theta);
phi = deg2rad(phi);

% Loop through each location (r,theta(i),phi(i))
progress_updates = floor((0.1:0.1:1)*length(theta));
fprintf('sh2data - Looping through location data\n')
for i = 1:length(theta)
    % Compute the Schmidt-semi normalized associated legendre polynomials
    % and their derivatives
    % This runs quickly, haven't bothered to subset from minN to maxN
    % Do it at the end with indices instead
    [Pnm,dPnmdt] = legendre_schmidt(max(N),theta(i));

    % Calculate the rest of the terms for theta(i) by looping through
    % degrees and order n,m
    br_terms = zeros(max(N)*(max(N)+3)/2+1 - (min(N)*(min(N)+3)/2),1);
    bt_terms = zeros(max(N)*(max(N)+3)/2+1 - (min(N)*(min(N)+3)/2),1);
    bp_terms = zeros(max(N)*(max(N)+3)/2+1 - (min(N)*(min(N)+3)/2),1);

    ind_bterms = 0;
    ind_sshc_nm = min(N)^2;
    for n = min(N):1:max(N)
        ind_bterms = ind_bterms + 1;
        ind_sshc_nm = ind_sshc_nm + 1;
        % m = 0
        br_terms(ind_bterms)=(n+1)*sshc(ind_sshc_nm)*a_r(n+1-min(N));
        bt_terms(ind_bterms)=-sshc(ind_sshc_nm)*a_r(n+1-min(N));
        %bp_terms = 0 at m=0 and doesnt require changing
        % 0 < m <= n
        for m = 1:1:n
            ind_sshc_nm = ind_sshc_nm + 2;
            ind_bterms = ind_bterms + 1;

            br_terms(ind_bterms) = (n+1) * a_r(n+1-min(N)) * ((sshc(ind_sshc_nm-1)*cos(m*phi(i)) + sshc(ind_sshc_nm)*sin(m*phi(i))));
            bt_terms(ind_bterms) = -a_r(n+1-min(N)) * (sshc(ind_sshc_nm-1)*cos(m*phi(i)) + sshc(ind_sshc_nm)*sin(m*phi(i)));
            % If at the poles, dividing by sin(theta) is undefined
            if (theta(i)) == (0) || (theta(i)) == (pi())
                bp_terms(ind_bterms) = -cos(theta(i)) * a_r(n+1-min(N)) * (-sshc(ind_sshc_nm-1)*sin(m*phi(i)) + sshc(ind_sshc_nm)*cos(m*phi(i)));
            else
                bp_terms(ind_bterms) = -(m/sin(theta(i))) * a_r(n+1-min(N)) * (-sshc(ind_sshc_nm-1)*sin(m*phi(i)) + sshc(ind_sshc_nm)*cos(m*phi(i)));
            end
        end
    end
    
    % Legendre and its derivatives computed for all n=0:N. Need to subset
    % them
    ind_min = legendre_ind(min(N),0);
    ind_max = legendre_ind(max(N),max(N));


    Br(i) = Pnm(ind_min:ind_max) * br_terms;
    Bt(i) = dPnmdt(ind_min:ind_max) * bt_terms;
    % If at the poles, dividing by sin(theta) is undefined
    if (theta(i)) == 0 || (theta(i)) == pi()
        Bp(i) = dPnmdt(ind_min:ind_max) * bp_terms;
    else
        Bp(i) = Pnm(ind_min:ind_max) * bp_terms;
    end

    if any(i == progress_updates)
        fprintf('%s%%...',num2str(floor(i/length(theta)*100)))
    end
end


fprintf('\n')

return