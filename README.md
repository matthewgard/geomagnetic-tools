# geomagnetic-tools
A collection of scripts relating to processing geomagnetic data, spherical harmonics etc.

![Capture](https://user-images.githubusercontent.com/48192690/192556396-0c946764-e318-410e-a09e-eaf78519e280.PNG)

MATLAB/spherical_harmonics/
  - data2sh.m
      - Computes a solution of schmidt semi-normalized spherical harmonic coefficients given a set of geomagnetic field values using least squares
  - sh2data.m
      - Computes geomagnetic field values at specified locations given a set of schmidt semi-normalized spherical harmonic coefficients
  - legendre_schmidt.m
      - Calculate the set of schmidt semi-normalized associated legendre polynomials and their derivatives for a given degree n
  - legendre_ind.m
      - Get index within Pnm and dPnmdt calculated using legendre_schmidt for the associated degree (n) and order (m)

MATLAB/general_scripts/
  - declination.m
      - Compute declination from xyz data
  - inclination.m
      - Compute inclination from xyz data
  - readIAGA.m
      - Parse IAGA2002 format geomagnetic data into a MATLAB structure array for easier processing
  - xyz2hdz.m
      - Convert xyz geomagnetic data to hdz
  - hdz2xyz.m
      - Convert hdz geomagnetic data to xyz

MATLAB/plotting/
  - bluewhitered.m
      - Custom colormap for plotting geomagnetic data (blue->white->red)
