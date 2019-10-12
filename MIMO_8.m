
clc; clear all;
%Task 1
%Jakes Fading Filter
 
N = 34;
M = (N/2 - 1)/2;
fdmax = (50*1.6*800*10^9)/(60*60*3*10^8);%Hz(from calculations)
  
hn = 0;
Ts = 10^-3/30;%(33us)
Had = [1 1;1 -1];
for t = 1:10^5
    hn1 = 0;hn2 = 0;hn3 = 0;hn4 = 0;
    for n = 1:M
        Beta = pi*n/(M+1);
        Theta = 2*pi*n/N;
        gamma = 2*pi*(2)*n/(M+1);
        hn1 = hn1+ power(exp(1),1i*Beta)*cos(2*pi*t*Ts*fdmax*cos(Theta)+gamma);
        gamma = 2*pi*(4)*n/(M+1);
        hn2 = hn2+ power(exp(1),1i*Beta)*cos(2*pi*t*Ts*fdmax*cos(Theta)+gamma);
        gamma = 2*pi*(6)*n/(M+1);
        hn3 = hn3+ power(exp(1),1i*Beta)*cos(2*pi*t*Ts*fdmax*cos(Theta)+gamma);
        gamma = 2*pi*(8)*n/(M+1);
        hn4 = hn4+ power(exp(1),1i*Beta)*cos(2*pi*t*Ts*fdmax*cos(Theta)+gamma);
    end
    h1(t) = Had(1,1)*hn1;h2(t) = Had(2,1)*hn2;
    h3(t) = Had(1,2)*hn3;h4(t) = Had(2,2)*hn4;
end
 
%Avg power of each sample to be normalised to 1
E1 = 0;E2 = 0;E3 = 0;E4 = 0;
for k = 1:10^5
    E1 = E1 + (abs(h1(k)))^2;
    E2 = E2 + (abs(h2(k)))^2;
    E3 = E3 + (abs(h3(k)))^2;
    E4 = E4 + (abs(h4(k)))^2;
end
E1 = E1/10^5;E2 = E2/10^5;E3 = E3/10^5;E4 = E4/10^5;
 
% normalised samples
for i = 1:10^5
    H1(i) = h1(i)/E1;H2(i) = h2(i)/E2;H3(i) = h3(i)/E3;H4(i) = h4(i)/E4;
end
%plots
for i = 1:500
x(i) = i;
y1(i) =10*log10(abs(H1(i)));
y2(i) =10*log10(abs(H2(i)));
y3(i) =10*log10(abs(H3(i)));
y4(i) =10*log10(abs(H4(i)));
end
figure;
plot(x,y1,x,y2,x,y3,x,y4);%
xlabel('Sample');
ylabel('amplitude in db');
title("Rayleigh Fading channel" );
 
%Task 2
Mt = 2;Mr = 1;
nSample = 10^5;
for l = 1:4
    c(:,:,l) = sqrt(2)*[exp(1i*(l-1)*pi/2) 0;0 exp(1i*(l-1)*pi/2)];
end
s1(:,:,1) = sqrt(Mt).*[1 0;0 1];
s2(:,:,1) = sqrt(Mt).*[1 0;0 1];
s3(:,:,1) = sqrt(Mt).*[1 0;0 1];
s4(:,:,1) = sqrt(Mt).*[1 0;0 1];
rho_dB = 0:15;
rho = (10.^(rho_dB./10));
pep1 = zeros(1,length(rho_dB));
pep2 = zeros(1,length(rho_dB));
for n = 2:4
    s1(:,:,n) = c(:,:,1)*s1(:,:,n-1)/sqrt(Mt);
    s2(:,:,n) = c(:,:,2)*s2(:,:,n-1)/sqrt(Mt);
    s3(:,:,n) = c(:,:,3)*s3(:,:,n-1)/sqrt(Mt);
    s4(:,:,n) = c(:,:,4)*s4(:,:,n-1)/sqrt(Mt);
end
%Hc = (randn(Mt,Mr)+1i*randn(Mt,Mr))/sqrt(2);
%nSample = 50000;
for itr = 1:nSample
    for n = 1:4
            Np(:,:,n) = (randn(Mt,Mr)+1i*randn(Mt,Mr))/sqrt(2);
            if(n == 1)
                NpHat(:,:,n) = Np(:,:,n);
            else
                
                NpHat(:,:,n) = Np(:,:,n)-sqrt(1/Mt)*c(:,:,2)*Np(:,:,n-1);
            end
    end
    H11 = [H1(itr);H2(itr)];%H22 =[H2(itr);H2(itr)];
        for nn = 1:length(rho_dB)
            %transmit c2 at T1 and T2
            y1_1 = sqrt(rho(nn)/Mt).*(s1(:,:,1)*H11)+Np(:,:,1);
            %y2_1 = sqrt(1/Mt).*(c(:,:,1)*y1_1)+NpHat(:,:,2);
            y1_2 = sqrt(rho(nn)/Mt).*(s2(:,:,1)*H11)+Np(:,:,1);
            y2_2 = sqrt(1/Mt).*(c(:,:,2)*y1_2)+NpHat(:,:,2);
            y1_3 = sqrt(rho(nn)/Mt).*(s3(:,:,1)*H11)+Np(:,:,1);
            %y2_3 = sqrt(1/Mt).*(c(:,:,3)*y1_3)+NpHat(:,:,2);
            y1_4 = sqrt(rho(nn)/Mt).*(s4(:,:,1)*H11)+Np(:,:,1);
            %y2_4 = sqrt(1/Mt).*(c(:,:,4)*y1_4)+NpHat(:,:,2);
            
            c1hat = (norm(y2_2 - sqrt(1/Mt).*(c(:,:,1)*y1_2),2))^2;
            c2hat = (norm(y2_2 - sqrt(1/Mt).*(c(:,:,2)*y1_2),2))^2;
            c3hat = (norm(y2_2 - sqrt(1/Mt).*(c(:,:,3)*y1_2),2))^2;
            c4hat = (norm(y2_2 - sqrt(1/Mt).*(c(:,:,4)*y1_2),2))^2;
            if((c1hat<c2hat)||(c3hat<c2hat)||(c4hat<c2hat))
                pep1(nn) = pep1(nn)+1;
            end
            end
end
pep1 =  pep1/nSample;
Mr = 2;
%Task 3
for itr = 1:nSample
    for n = 1:4
            Np1(:,:,n) = (randn(Mt,Mr)+1i*randn(Mt,Mr))/sqrt(2);
            if(n == 1)
                NpHat1(:,:,n) = Np1(:,:,n);
            else
                
                NpHat1(:,:,n) = Np1(:,:,n)-sqrt(1/Mt)*c(:,:,2)*Np1(:,:,n-1);
            end
    end
    H11 = [H1(itr) H2(itr);H3(itr) H4(itr)];%H22 =[H2(itr);H2(itr)];
        for nn = 1:length(rho_dB)
            %transmit c2 at T1 and T2
            y1_1 = sqrt(rho(nn)/Mt).*(s1(:,:,1)*H11)+Np1(:,:,1);
            y2_1 = sqrt(1/Mt).*(c(:,:,1)*y1_1)+NpHat1(:,:,2);
            y1_2 = sqrt(rho(nn)/Mt).*(s2(:,:,1)*H11)+Np1(:,:,1);
            y2_2 = sqrt(1/Mt).*(c(:,:,2)*y1_2)+NpHat(:,:,2);
            y1_3 = sqrt(rho(nn)/Mt).*(s3(:,:,1)*H11)+Np1(:,:,1);
           % y2_3 = sqrt(1/Mt).*(c(:,:,3)*y1_3)+NpHat(:,:,2);
            y1_4 = sqrt(rho(nn)/Mt).*(s4(:,:,1)*H11)+Np1(:,:,1);
           % y2_4 = sqrt(1/Mt).*(c(:,:,4)*y1_4)+NpHat(:,:,2);
 
            c1hat = (norm(y2_2 - sqrt(1/Mt).*(c(:,:,1)*y1_2),2))^2;
            c2hat = (norm(y2_2 - sqrt(1/Mt).*(c(:,:,2)*y1_2),2))^2;
            c3hat = (norm(y2_2 - sqrt(1/Mt).*(c(:,:,3)*y1_2),2))^2;
            c4hat = (norm(y2_2 - sqrt(1/Mt).*(c(:,:,4)*y1_2),2))^2;
            if((c1hat<c2hat)||(c3hat<c2hat)||(c4hat<c2hat))
                pep2(nn) = pep2(nn)+1;
            end
            end
end
pep2 =  pep2/nSample;
figure;
semilogy(rho_dB,pep1,'*-',rho_dB,pep2,'^-');
xlabel('SNR(dB)')
ylabel('BER');
legend('Mt=2,Mr=1','Mt=2,Mr=2');
title('Differential ST code SNR vs BER');
grid on;
                
