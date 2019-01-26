%% loading data
close all;
load proj18.mat 



%% Excercise D i) - using ElGeneina model on Kassala

close all; 

% Importing rain-data
nvdi_ka = Kassala.nvdi; 
nvdi_eg = ElGeneina.nvdi;

% Converting nvdi from [0,255] to [-1, 1]
nvdi_ka = uintconv(nvdi_ka); 
nvdi_eg = uintconv(nvdi_eg);

%% Plotting NVDI index for Kassala against ElGeneina 
figure
plot(nvdi_ka,'r') 
hold on
plot(nvdi_eg,'b')
title('NVDI index for Kassala (red) and ElGeneina (blue)')
xlabel('Year') 
ylabel('NVDI Index')
print -depsc D1_NVDI_ElG_vs_Kas
hold off

%% Plotting precipitation for Kassala against ElGeneina

figure
plot(Kassala.rain_org,'r') 
hold on
plot(ElGeneina.rain_org,'b')
title('Precipitation for Kassala (red) and ElGeneina (blue)')
xlabel('Time') 
ylabel('Precipitation')
print -depsc D3_rain_ElG_vs_Kas
hold off
hold off

sum(Kassala.rain_org)
sum(ElGeneina.rain_org)
sum(ElGeneina.rain_org)/sum(Kassala.rain_org)

%%

% Performning basic analysis on data 
basicanalysis(nvdi_ka, 50)
print -depsc D2_acfpacfnormplotoverview


% Splitting the data into modelling and validation sets 2/3 - model 1/3 -
% val

mean_nvdi_ka = mean(nvdi_ka); 
nvdi_ka = nvdi_ka - mean(nvdi_ka);
n_est = floor(2/3 * length(nvdi_ka)); 
nvdi_ka_m = nvdi_ka(1:n_est); 
nvdi_ka_v = nvdi_ka(n_est + 1:end); 




%% Using model for ElGeneina on Kassala
nvdi_data_ka = iddata(nvdi_ka_m')
e_ka1 = resid(nvdi_data_ka, model_init_eg_p); 

basicanalysis(e_ka1.y, 50)
print -depsc D4_acfpacfnormplotnoinput
figure
whitenessTest(e_ka1.y) 
print -depsc D5_whitenesstestnoinput


%% Reestimating parameters without input
close all
S=36;
A = [1 0 0];
C = [1 0 zeros(1, S-2) -0.9]; 
model_init_ka = idpoly(A, [], C); 
model_init_ka.Structure.c.Free = [0 1 zeros(1,S-2) 1];
model_init_ka_p = pem(nvdi_data_ka, model_init_ka); 
e_ka2 = resid(nvdi_data_ka, model_init_ka_p);

basicanalysis(e_ka2.y, 100)
print -depsc D_BAARMAreest

figure
whitenessTest(e_ka2.y) 
print -depsc D_WTARMAreest
%% Saving parameters

A_eg = model_init_eg_p.A;  
C_eg = model_init_eg_p.C; 

A_ka = model_init_ka_p.A;
C_ka = model_init_ka_p.C;

%% 1-step for ElGeneina model not reestimated
[Fk, Gk] = diophantine(C_eg,A_eg,1); 

yhat_k1_ka = filter(Gk, C_eg, nvdi_ka_v); 
yhat_k1_ka = yhat_k1_ka(max(length(Gk), length(C_eg):end)); 

figure
plot(nvdi_ka_v(length(nvdi_ka_v)-length(yhat_k1_ka)+1:end)+mean_nvdi_ka, 'b')
hold on
plot(yhat_k1_ka+mean_nvdi_ka, 'r') 
title('1-step prediction (red) versus real data (blue)')
legend('Original data', '1-step prediction', 'Location', 'North')
print -depsc D6_1stepnoinput

ehat1_ka = nvdi_ka_v(length(nvdi_ka_v)-length(yhat_k1_ka)+1:end)-yhat_k1_ka;  %y-yhat_k1
noise_var = var(ehat1_ka)


%% plotting K-step prediction for 4<=k<=10 ElGeneina model no reestimation

kstepvector = [4:10];
%P_error_ARMA = kstepb(model_init_ka_p,kstepvector,nvdi_ka_v,noise_var);
j = 1; 

for i = 1:length(kstepvector) 
    A_eg = model_init_eg_p.A;  
    C_eg = model_init_eg_p.C; 

    
    % Performing polynominal division 
    [Fk, Gk] = diophantine(C_eg,A_eg,kstepvector(i)); 
  
    % Performing k-step prediction 
    yhat = filter(Gk, C_eg, nvdi_ka_v); 
    
    % Removing some samples
    yhat = yhat(max(length(Gk), length(C_eg):end)); 
    
    % Making data equal length 
    tempData = nvdi_ka_v(length(nvdi_ka_v) - length(yhat) + 1:end); 
   
    % Plotting the k-step prediction 
    figure
    plot(tempData + mean_nvdi_ka, 'b')
    hold on
    plot(yhat + mean_nvdi_ka, 'r') 
    title(string(kstepvector(i))+ '-step prediction (red) versus real data (blue)')

    %Saving prediction error
    ehat_ka_BJ = nvdi_ka_v(length(nvdi_ka_v)-length(yhat)+1:end)-yhat;  %y-yhat_k1
    noise_var = var(ehat_ka_BJ);
    prediction_error(i) = noise_var;
end

prediction_error


%% 1-step for ElGeneina model reestimated
[Fk, Gk] = diophantine(C_ka,A_ka,1); 

yhat_k1_ka_reest = filter(Gk, C_ka, nvdi_ka_v); 
yhat_k1_ka_reest = yhat_k1_ka_reest(max(length(Gk), length(C_ka):end)); 

figure
plot(nvdi_ka_v(length(nvdi_ka_v)-length(yhat_k1_ka_reest)+1:end) + mean_nvdi_ka, 'b')
hold on
plot(yhat_k1_ka_reest + mean_nvdi_ka, 'r') 
title('1-step prediction (red) versus real data (blue)')
print -depsc D_1stepARMAreest

ehat1_ka_reest = nvdi_ka_v(length(nvdi_ka_v)-length(yhat_k1_ka_reest)+1:end)-yhat_k1_ka_reest;  %y-yhat_k1
noise_var_reest = var(ehat1_ka_reest)

%% plotting K-step prediction for 4<=k<=10 ElGeneina model reestimated

kstepvector = [4:10];
%P_error_ARMA = kstepb(model_init_ka_p,kstepvector,nvdi_ka_v,noise_var);
j = 1; 

for i = 1:length(kstepvector) 
    A_ka = model_init_ka_p.A;  
    C_ka = model_init_ka_p.C; 

    
    % Performing polynominal division 
    [Fk, Gk] = diophantine(C_ka,A_ka,kstepvector(i)); 
  
    % Performing k-step prediction 
    yhat = filter(Gk, C_ka, nvdi_ka_v); 
    
    % Removing some samples
    yhat = yhat(max(length(Gk), length(C_ka):end)); 
    
    % Making data equal length 
    tempData = nvdi_ka_v(length(nvdi_ka_v) - length(yhat) + 1:end); 
   
    % Plotting the k-step prediction 
    figure
    plot(tempData + mean_nvdi_ka, 'b')
    hold on
    plot(yhat + mean_nvdi_ka, 'r') 
    title(string(kstepvector(i))+ '-step prediction (red) versus real data (blue)')
    
    %Saving prediction error
    ehat_ka_reest = nvdi_ka_v(length(nvdi_ka_v)-length(yhat)+1:end)-yhat;  %y-yhat_k1
    noise_var_ka = var(ehat_ka_reest);
    prediction_error_reest(i) = noise_var_ka;
end

prediction_error_reest



%% Reconstruct Kassala rain with Kalman filter




%% % Filling missing samples with zeros
rain_org_mod = zeros(1, 3 * length(Kassala.rain_org));
k = 3;
i = 1; 
while i < length(Kassala.rain_org)
    rain_org_mod(k) = Kassala.rain_org(i); 
    k = k + 3; 
    i = i + 1; 
end

close all;
clear i;
clear k;

%% Initializing Kalman parameters

close all;
Ckalman = [1 1 1];
Re = [ 1 0 0; 0 0 0; 0 0 0]; %var e_rain
Rw = 0.001; % should be low, 
Akalman = [1 0 0; 1 0 0; 0 1 0];
z = kalmanfunction(Akalman,Re,Rw,Ckalman,rain_org_mod);



%% creating a vector with all the x:es
i=1;
k=3;
reconstructed_z_ka = zeros(1,length(z)/3);
while i < length(z)
    reconstructed_z_ka(i) = z(3,k);
    reconstructed_z_ka(i+1) = z(2,k);
    reconstructed_z_ka(i+2) = z(1,k);
    k=k+3;
    i=i+3;
end

%% Zero-mean and Creating model and validation sets for reconstructed rain
rain_ka = reconstructed_z_ka - mean(reconstructed_z_ka);
rain_ka = rain_ka(793:end);
rain_ka_m = rain_ka(1:length(nvdi_ka_m));
rain_ka_v = rain_ka(length(nvdi_ka_m)+1:end);

%% Using Box Jenkins model on Kassala 

z = iddata(nvdi_ka_m',rain_ka_m');
ehat_ka_BJ=resid(MboxJ,z);
basicanalysis(ehat_ka_BJ.y,30);
print -depsc D7_acfpacfnormplotBJ
present(MboxJ)
figure()
whitenessTest(ehat_ka_BJ.y)
print -depsc D8_whitenessTestBJ


%% Reestimating Box Jenkins model parameters

close all; 
p=1
q=0;
r=2;
s=1;
d=2;

A1 = [1 zeros(1,p)];
A2 = [1 zeros(1,r)];
B = [zeros(1,d) zeros(1,s)];
C = [1];
Mi = idpoly(1,B,C,A1,A2);
Mi.structure.b.Free = [zeros(1,d), ones(1,s)];
z = iddata(nvdi_ka_m',rain_ka_m');
MboxJ_reest = pem(z,Mi);
ehat=resid(MboxJ_reest,z);
basicanalysis(ehat.y,30);
print -depsc D_BAreestBJ 
crosscorrplot(30,rain_m,ehat.y);
print -depsc D_CCreestBJ
present(MboxJ_reest)
figure
whitenessTest(ehat.y)
print -depsc D_WTreestBJ


%% k-step prediction using Box Jenkins
close all;

A_rain = rain_model.A; 
C_rain = rain_model.C; 

A1 = MboxJ.D;
A2 = MboxJ.F; 
C1 = 1;
B = MboxJ.B;
A_tilde = conv(A1,A2);
C_tilde = conv(A2,C1);
B_tilde = conv(A1, B); 

A1_reest = MboxJ_reest.D;
A2_reest = MboxJ_reest.F;
B_reest = MboxJ_reest.B;
A_tilde_reest = conv(A1_reest,A2_reest);
C_tilde_reest = conv(A2_reest,C1);
B_tilde_reest = conv(A1_reest, B_reest); 


%% Performing 1-step prediction no reestimation

[Fy_1, Gy_1] = diophantine(C_tilde,A_tilde,1); 
 
BF_1 = conv(Fy_1,B_tilde); 
[Fu_1,Gu_1] = diophantine(BF_1,C_tilde,1);
uhat_k1 = filter(Gu_1,C_tilde,rain_ka_v); 
uhat_k1 = uhat_k1(max(length(Gu_1),length(C_tilde)):end);%remove samples
y1hat_k1 = filter(Gy_1,C_tilde,nvdi_ka_v); 
y1hat_k1 = y1hat_k1(max(length(Gy_1),length(C_tilde)):end); %remove samples
yhat_k1_BJ = y1hat_k1(1:end)+uhat_k1;

figure()
plot(yhat_k1_BJ + mean_nvdi_ka,'r')
hold on
plot(nvdi_ka_v(length(nvdi_ka_v) - length(yhat_k1_BJ) + 1:end) + mean_nvdi_ka,'b');
title('1-step prediction versus real data');
axis([0 length(yhat_k1_BJ) 0 0.5])
legend('1-step prediction', 'Original data', 'Location', 'North')
print -depsc D9_1stepBJ

%% Calculating the prediction error 

ehat1_BJ = nvdi_ka_v(length(nvdi_ka_v)-length(yhat_k1_BJ)+1:end)-yhat_k1_BJ;  %y-yhat_k1
noise_var_BJ = var(ehat1_BJ)

%% 4-10 step prediction no reestimation

pe_var_BJ = zeros(1,length(kstepvector));
for i = 1:length(kstepvector) 
    
    % Prediction of input 
    [F, G] = diophantine(C_rain, A_rain, kstepvector(i)); 
    xhat = filter(G, C_rain, rain_ka_v); 
    xhat = xhat(max(length(G), length(C_rain)):end); 
    
    % k-step prediction of output
    [Fy, Gy] = diophantine(C_tilde,A_tilde,kstepvector(i)); 
    BF = conv(Fy,B_tilde); 
    [Fu,Gu] = diophantine(BF,C_tilde,kstepvector(i));
    uhat = filter(Gu,C_tilde,rain_ka_v); 
    uhat = uhat(max(length(Gu),length(C_tilde)):end);%remove samples
    y1hat = filter(Gy,C_tilde,nvdi_ka_v); 
    y1hat = y1hat(max(length(Gy),length(C_tilde)):end); %remove samples
    
    
    yhat_p = y1hat(length(y1hat)-length(xhat)+1:end)+uhat(length(uhat)-length(xhat)+1:end) + filter(Fu, 1, xhat);

    % Plot against original data
    figure()
    plot(yhat_p + mean_nvdi_ka,'r')
    hold on
    plot(nvdi_ka_v(length(nvdi_ka_v) - length(yhat_p) + 1:end) + mean_nvdi_ka,'b');
    title(string(kstepvector(i)) + '-step prediction with predicted input')
    legend('k-step prediction', 'Original data', 'Location', 'North')
    axis([0 length(yhat_p) 0 0.5])
    
    pe_var_BJ(i) = var( nvdi_ka_v(length(nvdi_ka_v) - length(yhat_p) + 1:end) - yhat_p); 
end

pe_var_BJ


%% Performing 1-step prediction with reestimated parameters

[Fy_1, Gy_1] = diophantine(C_tilde_reest,A_tilde_reest,1); 
 
BF_1 = conv(Fy_1,B_tilde_reest); 
[Fu_1,Gu_1] = diophantine(BF_1,C_tilde_reest,1);
uhat_k1 = filter(Gu_1,C_tilde_reest,rain_ka_v); 
uhat_k1 = uhat_k1(max(length(Gu_1),length(C_tilde_reest)):end);%remove samples
y1hat_k1 = filter(Gy_1,C_tilde_reest,nvdi_ka_v); 
y1hat_k1 = y1hat_k1(max(length(Gy_1),length(C_tilde_reest)):end); %remove samples
yhat_k1_BJ = y1hat_k1(1:end)+uhat_k1;

figure()
plot(yhat_k1_BJ + mean_nvdi_ka,'r')
hold on
plot(nvdi_ka_v(length(nvdi_ka_v) - length(yhat_k1_BJ) + 1:end) + mean_nvdi_ka,'b');
title('1-step prediction versus real data');
axis([0 length(yhat_k1_BJ) 0 0.5])
legend('1-step prediction', 'Original data', 'Location', 'North')
print -depsc D9_1stepBJreest

%% Calculating the prediction error 

ehat1_BJ_reest = nvdi_ka_v(length(nvdi_ka_v)-length(yhat_k1_BJ)+1:end)-yhat_k1_BJ;  %y-yhat_k1
noise_var_BJ_reest = var(ehat1_BJ_reest)

%% 4-10 step prediction reestimation

pe_var_BJ_reest = zeros(1,length(kstepvector));
for i = 1:length(kstepvector) 
    
    % Prediction of input 
    [F, G] = diophantine(C_rain, A_rain, kstepvector(i)); 
    xhat = filter(G, C_rain, rain_ka_v); 
    xhat = xhat(max(length(G), length(C_rain)):end); 
    
    % k-step prediction of output
    [Fy, Gy] = diophantine(C_tilde_reest,A_tilde_reest,kstepvector(i)); 
    BF = conv(Fy,B_tilde_reest); 
    [Fu,Gu] = diophantine(BF,C_tilde_reest,kstepvector(i));
    uhat = filter(Gu,C_tilde_reest,rain_ka_v); 
    uhat = uhat(max(length(Gu),length(C_tilde_reest)):end);%remove samples
    y1hat = filter(Gy,C_tilde_reest,nvdi_ka_v); 
    y1hat = y1hat(max(length(Gy),length(C_tilde_reest)):end); %remove samples
    
    
    yhat_p = y1hat(length(y1hat)-length(xhat)+1:end)+uhat(length(uhat)-length(xhat)+1:end) + filter(Fu, 1, xhat);

    % Plot against original data
    figure()
    plot(yhat_p + mean_nvdi_ka,'r')
    hold on
    plot(nvdi_ka_v(length(nvdi_ka_v) - length(yhat_p) + 1:end) + mean_nvdi_ka,'b');
    title(string(kstepvector(i)) + '-step prediction with predicted input')
    legend('k-step prediction', 'Original data', 'Location', 'North')
    axis([0 length(yhat_p) 0 0.5])
    
    
    pe_var_BJ_reest(i) = var( nvdi_ka_v(length(nvdi_ka_v) - length(yhat_p) + 1:end) - yhat_p); 
end

pe_var_BJ_reest




