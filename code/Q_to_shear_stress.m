% Alex plan
% 1. load Q data
% 2. Fourier decompose it into Qn
% 3. Estimate shear stress from Qn




% 1. load Q data
data = load('./test/data.txt');

% need to look into the original Excel to 
% gets which column is which

t = data(:,0);
Q = data(:,3); % not really

% just for the test let's take some data from
% the readme_pulsatile.m


clf; clear all; clc;

r=0.001; grid=64; timestep=28; 
ru=1060; freq=1.0; mu=0.0035; 

%p0=0.0; pn=[7.99 4.43 0.99 0.74]; phi=[3.1219 -1.6675 0.4821 0.3938]; 
q=[-7.7183,-8.2383,-8.6444,-8.8797,-9.6337,-10.5957,-11.8705,-10.0942,-6.2839,-1.1857,2.6043,4.4323,6.1785,7.8211,9.1311,9.9138,10.3447,10.4011,10.2807,9.8951,8.0597,5.6717,2.5232,1.3301,1.4405,1.9094,1.8145,0.8738,0.7055,0.7343,0.7788,0.7495,0.6711,-0.4796,-1.6541,-2.8643,-3.4902,-4.1714,-5.6581,-6.8024];
[q0,qn,phi] = FourierSeries(q,0.025,4,1/freq);


% 2. Fourier decompose Q

