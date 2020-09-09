% Alex plan
% 1. load Q data
% 2. Fourier decompose it into KQ
% 3. Estimate shear stress from KQ


%%

% 1. load Q data
% data = load('./test/data.txt');
% metadata = load('./test/info.txt');

% run the reader 
clear all
read_data_xlsx

%%
%% Allocate imported array to column variable names
% Normalisedtime1 = data(:,1);
% PA_Surgery_90_Aortic1 = data(:,2);
% PA_Surgery_90_Branch1 = data(:,3);

t = Normalisedtime1;
% Q = PA_Surgery_90_Aortic1;


Q = PA_Surgery_90_Branch1;
Q=Q/1000/60; % idit divided in 60000 to converge L/min to m^3/s

% cleaning
ind = ~isnan(t);
t = t(ind);
Q = Q(ind);


figure
plot(t,Q,'-o');
title('original')


% 90 RPM 1.5 Hz
% 115RPM 1.9167 Hz
% 135RPM 2.25 Hz
% T = 1/1.5; % we use normalized time and take BPM values (in Hz)
t = t/1.5; % we trust t for being 0-1 sec and normalize by BPM to get physical time
T = max(t);
figure
plot(t,Q,'-o');
xlabel('t (sec)');
title('scaled to physical units')

nf = 10; % we can take n frequencies


% aortic diameter [mm]:	42
% branch diameter [mm]	18

a = 0.042/2; % (m) - aorta radius
a = 0.018/2; % (m) - branch radius

% we need to resample at fixed frequency

dt = diff(t(1:2));

% remove nans (from excel import)




% just for the test let's take some data from
% the readme_pulsatile.m

%%
% r=0.001; grid=64; timestep=28; 
%  
% 


% 2. Fourier decompose Q

% [KQ0,KQ,pshi] = FourierSeries(Q,nf,t);
fourier_decomposition


%%
% 3. Q to \tau

% flow properties: density, dynamic and kinematic viscosity
rho=1060; 
mu=0.0035; 
nu = mu/rho; 

nt = length(t);
temp2 = ones(1,nt);  % temp2 = temp2(1,:);


w0 = 2*pi/T;			    % fundamental radian frequency
alpha0 = a*sqrt(w0/nu);	    % Womersley parameter, fundamental frequency
w = w0*[1:nf];			    % array of radian frequencies
alpha = a*sqrt(w/nu);	    % array of Womersley parameters
K = 0*KQ;					% initialize pressure gradient array
Q = KQ0*temp2;				% initialize flow rate array with DC flow component
Pp = -8*nu*Q/(pi*a^4);	% initialize pressure gradient array with DC component
Tau_steady = -a*Pp/2;	% initialize wall shear stress array with DC component
Tau = Tau_steady;
j = sqrt(-1);
j32 = j^1.5;


ny = 101;		% # points in profile (over diameter)

mean_u = mean(Q)*100 /(6*(pi*a*a));% calculate mean velocity [cm/s]
Re = mean_u*(2*a)/nu;% calculate Reynolds number

y = [-1: 2/(ny-1): 1]';				% column vector of y's
u = (2*KQ0/(pi*a^2))*((1 - y.*y)*temp2);

for n = 1:nf		% Sum over all frequencies
   K(n) = KQ(n)*j*w(n)/(pi*a^2*(1 - 2*besselj(1,j32*alpha(n))/(j32*alpha(n)*besselj(0,j32*alpha(n)))));
   Ktau(n) = rho*K(n)*a*besselj(1,j32*alpha(n))/(j32*alpha(n)*besselj(0,j32*alpha(n)));
   for k = 1:nt		% Calculate for each point in time
      Q(k) = Q(k)+real(KQ(n)*exp(j*w(n)*t(k)));
      Pp(k) = Pp(k)+real(K(n)*exp(j*w(n)*t(k)));
      Tau(k) = Tau(k)+real(Ktau(n)*exp(j*w(n)*t(k)));
      for m = 1:ny		% Calculate for each spatial location y
         u(m,k) = u(m,k) + real((K(n)*a^2/(nu*alpha(n)^2*i))*(1-besselj(0,j32*alpha(n)*y(m))/besselj(0,j32*alpha(n)))*exp(i*w(n)*t(k)));
      end	% y loop
   end	% t loop
end	% f loop

%%
nExp = 1;

figure;
% Plot time waveforms of flow, pressure gradient, wall shear stress, and centerline velocity
subplot(2,2,1); plot(t,Q); xlabel('time (s)'); ylabel('Flow Rate (ml/s)'); title(sprintf('set %g',nExp));
subplot(2,2,2); plot(t,Tau,t,Tau_steady); xlabel('time (s)'); ylabel('Wall Shear Stress (dyn/cm^2)');title(sprintf('set %g',nExp));
subplot(2,2,3); plot(t,Pp); xlabel('time (s)'); ylabel('Pressure Gradient (dyn/cm^3)');title(sprintf('set %g',nExp));
subplot(2,2,4); plot(t,u((ny+1)/2,:)); xlabel('time (s)'); ylabel('Centerline Velocity (cm/s)');title(sprintf('set %g',nExp));


figure
plot(t,Q/max(Q),'r-',t,Tau/max(Tau),'b--',t,Pp/max(Pp),'g-.');title('Observe phase')