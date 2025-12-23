close all;
clearvars;
clc;

t_start = 0;
t_end = 3;
N = 5e4;
t = t_start : (t_end - t_start)/N : t_end - (t_end - t_start)/N;

A0 = 10;
A1 = 5; A2 = 3; A3 = 1;
W1 = 3*pi; W2 = pi; W3 = (pi/4);
m_t = A0 + A1 * sin(W1*t) + A2 * (cos(W2*t)).^3 + A3 * sin(W3*t);
Len_m = length(m_t);

figure(1);
plot(t,m_t);
title('m(t)');
xlabel('time(s)');

fs = 500;
Ts = 1 / fs;

samp_m = zeros(1,Len_m/fs);

for i = 1: 100
    samp_m(i) = m_t(fs*i);
end

%samp_m = m_t(1:fs:length(m_t));
len_samp = length(samp_m);

t_s2 = (t_end - t_start)/len_samp;
t_samp = t_start : t_s2 : t_end-t_s2;

figure(2)

subplot(2,1,1)
plot(t_samp,samp_m,'p');
title('sampled m(t)');

subplot(2,1,2)
plot(t_samp,samp_m,'b');


%%%% JIM ) Quantization

q_n = 32;
dist = max(m_t)-min(m_t);
dist2 = dist/(q_n-1);
Q_lvl = zeros(1,q_n);
for i = 1:q_n
    Q_lvl(i) = min(m_t)+(i-1)*dist2;
end

Q_m_t = zeros(1,len_samp);
for i = 1:length(samp_m)
    Q_m_t(i) = round((samp_m(i) - min(m_t))/(dist2)) * dist2 + min(m_t);
end


figure(3)
subplot(2,1,1);
plot(t_samp,Q_m_t,'x');
ylim([0 20]);
title("Quantized m_t");
subplot(2,1,2);
plot(t_samp,Q_m_t,'r');



%%%% DAL ) 

fs_tri = 1000;
ts_tri = (1 - 0) / fs_tri;
N_tri = (1-0)/ts_tri;
t3_tri= 0 :ts_tri: 1 - ts_tri;

data = load('p.mat');
tri = cell2mat(struct2cell(data));
figure(4)
plot(t3_tri , tri);
title("Pulse Of Amplitude modulation");

tri_ENG = sum(tri.*tri);

tri_lvl = -31:2:31;

N_P = fs_tri*len_samp;
P_m_t = zeros(1,N_P);
for i = 1:length(Q_m_t)
    for j = 1:length(Q_lvl)
        if ( Q_m_t(i) == Q_lvl(j) )
            P_m_t(1000*(i-1)+1:1000*i) = tri_lvl(j)*tri;
        end
    end
end

P_power = sum(P_m_t.*P_m_t)/length(P_m_t);

t_start_P = 0; t_end_P = 100; ts_P = 1/fs_tri;
t_P = t_start_P:ts_P:t_end_P - ts_P;

figure(5)
plot(t_P , P_m_t)
title("AMP modulation");
xlabel("time");



% GrayCode: 
bin_to_gray = [0, 1, 3, 2, 6, 7, 5, 4, 12, 13, 15, 14, 10, 11, 9, 8,   24, 25, 27, 26, 30, 31, 29, 28, 20, 21, 23, 22, 18, 19, 17, 16];
gray = dec2bin(bin_to_gray , 5);

G_m_t = [];
for i = 1:length(Q_m_t)
    for j = 1:length(Q_lvl)
        if ( Q_m_t(i) == Q_lvl(j) )
            k = find(bin_to_gray == j-1);
            G_m_t = [G_m_t ; gray(k,:)];
        end
    end
end



%%%% HE )

SNR = 2;
N0 = P_power/(10^(SNR/10));
noise = randn(1,N_P)*sqrt(N0);
P_m_with_noise = noise + P_m_t;


figure(6)
plot(t_P , P_m_with_noise, 'r')
title("digital sig with noise");
xlabel("time");






%%%% VAV and Z )

D_q_m_t = zeros(1,len_samp); J = [];
for i = 1:len_samp
    rec_sig =  sum(P_m_with_noise(fs_tri*(i-1)+1:fs_tri*(i)) .* tri)/tri_ENG;
    [~ , j] = min(abs(rec_sig - tri_lvl));
    D_q_m_t(i) = Q_lvl(j);
    J = [J;j];
end
D_m_t = spline(t_samp , D_q_m_t , t);
noise_ENG = tri_ENG/SNR;

figure(7)
plot(t_samp ,D_q_m_t, '*');
title("decoded sig");

figure(8)
subplot(2,1,1)
plot(t , m_t , 'g');
title("analog sig VS decoded dig sig");
subplot(2,1,2)
plot(t , D_m_t , 'r')

sum_err = sum(Q_m_t ~= D_q_m_t);
SIG_ERR = sum_err/len_samp;
ERROR = immse(m_t , D_m_t);









