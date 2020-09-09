% Fourier decomposition

test = false;

if test 
    Q =[-7.7183,-8.2383,-8.6444,-8.8797,-9.6337,-10.5957,-11.8705,-10.0942,-6.2839,-1.1857,2.6043,4.4323,6.1785,7.8211,9.1311,9.9138,10.3447,10.4011,10.2807,9.8951,8.0597,5.6717,2.5232,1.3301,1.4405,1.9094,1.8145,0.8738,0.7055,0.7343,0.7788,0.7495,0.6711,-0.4796,-1.6541,-2.8643,-3.4902,-4.1714,-5.6581,-6.8024];
end

% T = 1;
% nf = 10;
    
Qf = fft(Q);

nt = length(t);

dt = linspace(0,T,length(Q));
dt1=dt(1:end);

KQ0 = real(Qf(1))/nt;

j = sqrt(-1);

for n = 1:nf
   KQ(n) = 2*Qf(n+1)/nt;
   w(n)	= 2*pi*n/T;
end
                                                               
KQmag = [KQ0, abs(KQ)];
KQphase = [0.0, angle(KQ)];

t1 =[0:0.001*T:(T-(dt1(2)-dt1(1)))]';
nt1 = length(t1);
% max(nt1);

Qr = KQ0*ones(1,nt1)';

for n = 1:nf
   for k = 1:nt1
      Qr(k) = Qr(k) + real(KQ(n)*exp(j*w(n)*t1(k)));
   end
end

%fprintf(1,['Mean volume flow rate = ' num2str(mean(Q)) '(ml/s) \n']);
% t2=T*[0 0.125 0.25 0.375 0.5 0.625 0.75];

% disp(length(dt1))
% disp(length(Qraw))
% disp(length(t1))
% disp(length(Q))

figure;
plot(dt1,Q,t1,Qr);
hold on;
% apparently q(nExp, :) was the experimental data
% and EX3 was the reconstructed Q spline

% plot(t2,q(nExp,:),'LineStyle','none','Marker','o'); 
xlabel('Time (s)'); 
ylabel('Q (ml/s)'); 
legend('Original','Fourier Series');