% Alex plan
% 1. load Q data
% 2. Fourier decompose it into Qn
% 3. Estimate shear stress from Qn


%%

% 1. load Q data
% data = load('./test/data.txt');
% metadata = load('./test/info.txt');

sample,PH_DA,PH_IA,normalisedtime,PA_AorticFlowLmin,PA_BranchFlowmLmin = importdata_excel('test/Phase_Average.xlsx');

%%
t = normalisedtime;
Q = PA_AorticFlowLmin; % not really

% we need to resample at fixed frequency

dt = diff(t(1:2));

% remove nans (from excel import)
 
ind = ~isnan(t);

t = t(ind);
Q = Q(ind);

plot(t,Q)
xlabel('t^*')

plot(t,Q,'-o')

% just for the test let's take some data from
% the readme_pulsatile.m

%%
% r=0.001; grid=64; timestep=28; 
% ru=1060; freq=1.0; mu=0.0035; 
% 


% 2. Fourier decompose Q
nfreq = 10; % we can take n frequencies
T = 1; % we use normalized time

[Q0,Qn,phi] = FourierSeries(Q,dt,nfreq,T);

% 3. Q to \tau

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


k=[1 13 23 38 51 63 76];
col=['b';'k';'r';'g';'y';'c';'m'];

subplot(2,4,nExp);
for j=1:7
subplot(2,4,j);
for i=1:7
    b = fh(t1(i),y,T,nu,a,KQ,KQ0,nf,rho);
    plot(y,b,'LineStyle','--','Color',sprintf('%c',col(i)));hold on;
    temp=d(nExp).(sprintf('t%g',i-1));
    plot(temp(:,1)/max(temp(:,1)),100*temp(:,2),'LineStyle','none','Marker','o','Color',sprintf('%c',col(i)));
    hold on;
   
end

xlabel('r/R');
ylabel('u [cm/s]');
title(sprintf('set %g alpha=%g R=%g',nExp,alpha0,a));
grid on;
axis([-1 1 -2 15]);

end
umin = min(min(u));  umax = max(max(u));



%Make an animated velocity profile plot
%{
clear M;
for n = 1:nt
   h1 = figure(3);
   plot(u(:,n),y); axis([1.1*umin 1.1*umax -1 1]);
   xlabel(['U(y,t=' num2str((n-1)*T/nt) ')']); ylabel('y'); title('Velocity Profiles');
   M(n) = getframe(h1);
end

moviereps = 30;	% number of repetions of the movie
speedratio = 0.5;	% ratio of movie speed to actual speed
%movie(M,moviereps,speedratio*nt/T);
movie2avi(M,'velocity_profile_femoral.avi','quality',100) % record the animation as AVI files
%}

%subplot(1,2,2);
%Plot selected velocity profiles on a static plot
%t = [0: T/nt : T];
%hold on;
%for i = 1:10:100
    %plot(y,fh(t(i),y,T,nu,a,KQ,KQ0,nf,rho));
%end
%legend('0%','10%','20%','30%','40%','50%','60%','70%','80%','90%');
%title(sprintf('Velocity Profiles set%g, alpha=%g',nExp,alpha0)); xlabel('U(y;percent of period) [cm/s]'); ylabel('r/R');
%grid on;
%axis([-1 1 -2 20]);

%hold off
%figure;
%Plot time waveforms of flow, pressure gradient, wall shear stress, and centerline velocity
%subplot(2,2,1); plot(t,Q); xlabel('time (s)'); ylabel('Flow Rate (ml/s)');title(sprintf('set %g',nExp));
%subplot(2,2,2); plot(t,Tau,t,Tau_steady); xlabel('time (s)'); ylabel('Wall Shear Stress (dyn/cm^2)');title(sprintf('set %g',nExp));
%subplot(2,2,3); plot(t,Pp); xlabel('time (s)'); ylabel('Pressure Gradient (dyn/cm^3)');title(sprintf('set %g',nExp));
%subplot(2,2,4); plot(t,u((ny+1)/2,:)); xlabel('time (s)'); ylabel('Centerline Velocity (cm/s)');title(sprintf('set %g',nExp));



% save model_time_stamp.mat ts