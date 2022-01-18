%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The pseudo inclination method for magnetic reduction to the pole at low latitudes
% Author: Lianghui Guo (guolh@cugb.edu.cn)
% Organization: China University of Geosciences (Beijing), School of Geophysics and Information Technology
% Compiled version: MATLAB R2017b
% Reference: 
%       Guo L H, Shi L, Meng X H. The antisymmetric factor method for magnetic reduction
%       to the pole at low latitudes. Journal of Applied Geophysics, 2013, 92, 103-109.
% Description of the input parameters: 
%       infile: original magnetic anomaly data file
%       inc: inclination, unit: ¡ã
%       dec: declination, unit: ¡ã
%       inc2: pseudo inclination, unit: ¡ã
% Description of the output parameters: £º
%       outfile: result file
% Description of primary identifiers£º
%       x, y: x, y verctor
%       nx, ny: number of points in x and y directions
%       dx, dy: spacing in x and y directions
%       g: anomaly before rtp
%       npts: extension points
%       grtp: anomaly after rtp
% Description of subroutine function£º
%       readgrd.m: read surfer text grd file
%       taper2d.m: cosine attenuation edge extension
%       rtp_pi.m: pseudo inclination method for magnetic reduction to the pole (rtp)
%       savegrd.m: save surfer text grd file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
infile = 'anomaly.grd'; 
%%%%%%%%%%%% I/O parameters %%%%%%%%%%%%%
inc = 0; 
dec = 0; 
inc2 = 30; 
outfile = 'anomaly_rtp_pi.grd'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[g,x,y,nx,ny,dx,dy] = readgrd(infile);
[X,Y] = meshgrid(x,y);
nmax = max([nx ny]);
npts = 2^nextpow2(nmax);
% The pseudo inclination method for magnetic reduction to the pole (rtp)
grtp = rtp_pi(g,npts,ny,nx,dec,inc,inc2);
% Save
savegrd(grtp,x,y,nx,ny,outfile);