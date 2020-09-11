% Alex plan
% 1. load Q data
% 2. Fourier decompose it into KQ
% 3. Estimate shear stress from KQ

% run the reader
clc, clear all, close all
%%
% general, global data
nf = 10; % we can take n frequencies
% flow properties: density, dynamic and kinematic viscosity
rho=1060;
mu=0.0035;
nu = mu/rho;

ny = 101;		% # points in profile (over diameter)

% read_data_xlsx

%% Import the data
[~, ~, raw] = xlsread('test\data.xlsx','Sheet1','A2:AE249');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[4,8,12,16,20]);
raw = raw(:,[1,2,3,5,6,7,9,10,11,13,14,15,17,18,19,21,22,23,24,25,26,27,28,29,30,31]);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
data = reshape([raw{:}],size(raw));

% %% Allocate imported array to column variable names
% Normalisedtime1 = data(:,1);
% PA_Surgery_90_Aortic1 = data(:,2);
% PA_Surgery_90_Branch1 = data(:,3);
%
%
% Normalisedtime2 = data(:,4);
% PA_Surgery_135_Aortic = data(:,5);
% PA_Surgery_135_Branch = data(:,6);
%
% Normalisedtime3 = data(:,7);
% PA_Hybrid_90_Aortic = data(:,8);
% PA_Hybrid_90_Branch = data(:,9);
%
%
% Normalisedtime4 = data(:,10);
% PA_Hybrid_115_Aortic = data(:,11);
% PA_Hybrid_115_Branch = data(:,12);
%
% Normalisedtime5 = data(:,13);
% PA_Hybrid_135_Aortic = data(:,14);
% PA_Hybrid_135_Branch = data(:,15);
%
% Normalisedtime6 = data(:,16);
% PA_Chimney_90_Aortic = data(:,17);
% PA_Chimney_90_Branch = data(:,18);
%
% Normalisedtime7 = data(:,20);
% PA_Chimney_115_Aortic = data(:,21);
% PA_Chimney_115_Branch = data(:,22);
%
% Normalisedtime8 = data(:,24);
% PA_Chimney_135_Aortic = data(:,25);
% PA_Chimney_135_Branch = data(:,26);

%% Clear temporary variables
clearvars raw cellVectors R;

%% we have 8 datasets

bpm = [90,135,90,115,135,90,115,135];
index=0; %idit
counter=0; %idit

for col = [1]%,4,7,10,13,16,20,24] % columns in the Excel file
    index=index+1; % idit added
    % 8 datasets, see above
    for data_type = 1%:2
        counter=counter+1;
        % 1 is Aorta, 2 is branch
        if data_type == 1
            % aortic diameter [mm]:	42
            % branch diameter [mm]	8
            a = 0.042/2; % (m) - aorta radius
            Qcoef=0.001;
        elseif data_type == 2
            a = 0.008/2; % (m) - branch radius
            Qcoef=1e-6;
        end
        
        t = data(:,col); % the NormalizedTime
        Q = data(:,col + data_type); % or Aorta or Branch
        
        
        Q=Q*Qcoef/60; % idit divided in 60000 to converge L/min to m^3/s
        
        
        
        % cleaning
        ind = ~isnan(t);
        t = t(ind);
        Q = Q(ind);
        
        
        % figure
        % plot(t,Q,'-o');
        % title('original')
        
        
        % 90 RPM 1.5 Hz
        % 115RPM 1.9167 Hz
        % 135RPM 2.25 Hz
        % T = 1/1.5; % we use normalized time and take BPM values (in Hz)
        
        freq = bpm(col)/60;
        
        
        
        t = t/freq; % we trust t for being 0-1 sec and normalize by BPM to get physical time
        T = max(t);
        
        % debugging
        % figure
        % plot(t,Q,'-o');
        % xlabel('t (sec)');
        % title('scaled to physical units')
        
        
        
        %% 2. Fourier decompose Q
        Qf = fft(Q);
        nt = length(t);
        KQ0 = real(Qf(1))/nt;
        KQ = zeros(nf,1);
        for n = 1:nf
            KQ(n) = 2*Qf(n+1)/nt;
            % w(n)	= 2*pi*n/T;
        end
        %% 3. Q to \tau
        
        
        
        % nt = length(t);
        temp2 = ones(1,nt);  % temp2 = temp2(1,:);
        
        
        w0 = 2*pi/T;			    % fundamental radian frequency
        alpha0 = a*sqrt(w0/nu);	    % Womersley parameter, fundamental frequency
        w = w0*(1:nf);			    % array of radian frequencies
        alpha = a*sqrt(w/nu);	    % array of Womersley parameters
        
        K = 0*KQ;					% initialize pressure gradient array
        Qr = KQ0*temp2;				% initialize flow rate array with DC flow component
        Pp = -8*nu*Qr/(pi*a^4);	% initialize pressure gradient array with DC component
        Tau_steady = -a*Pp/2;	% initialize wall shear stress array with DC component
        Tau = Tau_steady;
        j = sqrt(-1);
        j32 = j^1.5;
        
        
        
        
        mean_u = mean(Qr)*100 /(6*(pi*a*a));% calculate mean velocity [cm/s]
        mean_u = mean(Qr) /((pi*a*a))%;% idit changed to m/s  and deleted the 6
        Re = mean_u*(2*a)/nu%;% calculate Reynolds number
        
        y = (-1: 2/(ny-1): 1)';				% column vector of y's
        u = (2*KQ0/(pi*a^2))*((1 - y.*y)*temp2);
        
        Ktau = zeros(nf,1);
        
        for n = 1:nf		% Sum over all frequencies
            K(n) = KQ(n)*j*w(n)/(pi*a^2*(1 - 2*besselj(1,j32*alpha(n))/(j32*alpha(n)*besselj(0,j32*alpha(n)))));
            Ktau(n) = rho*K(n)*a*besselj(1,j32*alpha(n))/(j32*alpha(n)*besselj(0,j32*alpha(n)));
            for k = 1:nt		% Calculate for each point in time
                Qr(k) = Qr(k)+real(KQ(n)*exp(j*w(n)*t(k)));
                Pp(k) = Pp(k)+real(K(n)*exp(j*w(n)*t(k)));
                Tau(k) = Tau(k)+real(Ktau(n)*exp(j*w(n)*t(k)));
                for m = 1:ny		% Calculate for each spatial location y
                    u(m,k) = u(m,k) + real((K(n)*a^2/(nu*alpha(n)^2*1i))*(1-besselj(0,j32*alpha(n)*y(m))/besselj(0,j32*alpha(n)))*exp(1i*w(n)*t(k)));
                end	% y loop
            end	% t loop
        end	% f loop
        
        %%
        % nExp = 1;
        
        
        %idit
        %         figure;
        %         % Plot time waveforms of flow, pressure gradient, wall shear stress, and centerline velocity
        %         subplot(2,2,1); plot(t,Q,t,Qr); xlabel('time (s)'); ylabel('Flow Rate (ml/s)'); % title(sprintf('set %g',nExp));
        %         subplot(2,2,2); plot(t,Tau,t,Tau_steady); xlabel('time (s)'); ylabel('Wall Shear Stress (dyn/cm^2)');% title(sprintf('set %g',nExp));
        %         subplot(2,2,3); plot(t,Pp); xlabel('time (s)'); ylabel('Pressure Gradient (dyn/cm^3)');% title(sprintf('set %g',nExp));
        %         subplot(2,2,4); plot(t,u((ny+1)/2,:)); xlabel('time (s)'); ylabel('Centerline Velocity (cm/s)');% title(sprintf('set %g',nExp));
        %
        %
        %         figure
        %         plot(t,Qr/max(Qr),'r-',t,Tau/max(Tau),'b--',t,Pp/max(Pp),'g-.');title('Observe phase')
        %
        %         figure
        %         hold on
        %         title(sprintf('Velocity profiles for \alpha = %f',alpha0));
        %
        %         for j=1:size(u,2)
        %             plot(u(:,j))
        %         end
        
        
        PI= (max(Q)-min(Q))/mean(Q);
        minmax=abs(min(Tau)/max(Tau));
        Ppr = length(t(Tau < 0))/length(t)*100; % percentage of negative shear stress
        OSI=0.5*(1-(trapz(t,Tau)/trapz(t,abs(Tau))));
        
        
        summ(counter,:)=[freq*60, alpha0, Re, PI, minmax, Ppr, OSI]; %idit
        
        
        
        
    end
end
summ  %idit
