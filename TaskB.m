% Loading data 
close all;
load proj18.mat 
rain_org = ElGeneina.rain_org;


%% Excercise B i) - non-recursive fitting without input

close all; 

% Importing rain-data
nvdi_eg = ElGeneina.nvdi; 

% Converting nvdi from [0,255] to [-1, 1]
nvdi_eg = uintconv(nvdi_eg); 

% Plotting the data 
T = linspace(datenum(1982,1,1),datenum(1999,12,31),length(nvdi_eg));
figure
plot(T, nvdi_eg) 
axis([0 length(nvdi_eg) 0 0.42]);
title('Overview of original NVDI-index data')
datetick('x');
xlabel('Year')
ylabel('NVDI')
print -depsc 0_NVDIoverview

% Performning basic analysis on data 
basicanalysis(nvdi_eg, 50)
print -depsc 0_acfpacfnormNVDIoverview
% Splitting the data into modelling and validation sets 2/3 - model 1/3 -
% val
mean_nvdi = mean(nvdi_eg); 
nvdi_eg = nvdi_eg - mean(nvdi_eg);

%% 
figure
bcNormPlot(nvdi_eg) 
%% 
n_est = floor(2/3 * length(nvdi_eg)); 
nvdi_m = nvdi_eg(1:n_est); 
nvdi_v = nvdi_eg(n_est + 1:end); 

%% Fitting ARMA-polynominals to data
close all
S=36;
nvdi_data = iddata(nvdi_m'); 
A = [1 0 0];
C = [1 0 zeros(1, S-2) -0.9]; 
model_init_eg = idpoly(A, [], C); 
model_init_eg.Structure.c.Free = [0 1 zeros(1,S-2) 1];
model_init_eg_p = pem(nvdi_data, model_init_eg); 
e = resid(nvdi_data, model_init_eg_p); 

basicanalysis(e.y, 50)
print -depsc B1_acfpacfnormnoinput

figure
whitenessTest(e.y) 
print -depsc B2_whitenessTestnoinput
% The data have strong jumps at a seasonality of 36. By applying a season
% filter of 36 and then modelling ARMA-process will first force an
% 36-seasonality dependancy into the data which we then try to model with
% ARMA-polynominals which seems to work fine.

%% Performing 1-step prediction 

A = model_init_eg_p.A;  
C = model_init_eg_p.C; 

[Fk, Gk] = diophantine(C,A,1); 

yhat_k1 = filter(Gk, C, nvdi_v); 
yhat_k1 = yhat_k1(max(length(Gk), length(C):end)); 


%% plot against validation set
figure
plot(nvdi_v(length(nvdi_v)-length(yhat_k1)+1:end) + mean_nvdi, 'b')
hold on
plot(yhat_k1 + mean_nvdi, 'r') 
title('1-step prediction versus given data')
legend({'Given data', '1-step prediction'}, 'Location', 'north')
ylabel('NVDI')
print -depsc B3_1stepnoinput

ehat1 = nvdi_v(length(nvdi_v)-length(yhat_k1)+1:end)-yhat_k1;  %y-yhat_k1
noise_var = var(ehat1);


%% plotting K-step prediction for 4<=k<=10
graph = 4; 
kstepvector = 4:10;
% P_error_ARMA = kstepb(model_init_eg_p,kstepvector,nvdi_v,noise_var, nvdi_mean);
pe_var_arma = zeros(1,length(kstepvector)); 
for i = 1:length(kstepvector) 
    % Performing polynominal division 
    [Fk, Gk] = diophantine(C,A,kstepvector(i)); 
    
    % Performing k-step prediction 
    yhat = filter(Gk, C, nvdi_v); 
    
    % Removing some samples
    yhat = yhat(max(length(Gk), length(C):end));
    
    % Making data equal length 
    tempData = nvdi_v(length(nvdi_v) - length(yhat) + 1:end); 
   
    % Plotting the k-step prediction 
    figure
    plot(tempData + mean_nvdi, 'b')
    hold on
    plot(yhat + mean_nvdi, 'r') 
    title(string(kstepvector(i))+ '-step prediction versus given data')
    ylabel('NVDI')
    legend({'Given data', 'k-step prediction'}, 'Location', 'north')
    %Saving prediction error
    ehat = nvdi_v(length(nvdi_v)-length(yhat)+1:end)-yhat;  %y-yhat_k1
    basicanalysis(ehat, 50); 
    pe_var_arma(i) = var(ehat);
end

%% B ii) using rain as input building a box jenkins model

close all
rain = reconstructed_z(793:end);

T = linspace(datenum(1982,1,1),datenum(1999,12,31),length(rain));
figure
plot(T, rain) 
axis([0 length(rain) 0 100]);
title('Overview of Kalman-interpolated rain data')
datetick('x');
xlabel('Year')
ylabel('Accumulated rain')
print -depsc B5_rainoverview

mean_rain = mean(rain); 
rain = rain - mean(rain); 
 

%% Creating model and validation set 

rain_m = rain(1:length(nvdi_m));
rain_v = rain(length(nvdi_m)+1:end);

%% Checking the plot of Kalman interpolation vs linear interpolation.

T = linspace(datenum(1982,1,1),datenum(1999,12,31),length(reconstructed_z(793:end)));
plot(T,reconstructed_z(793:end))
hold on
plot(T,ElGeneina.rain(793:end))
datetick('x');
xlabel('Year') 
ylabel('Accumulated rain')
title('Kalman interpolated rain versus linear interpolated rain')
legend('Reconstructed rain', 'Linear interpolated rain', 'Location' ,'NorthWest')
print -depsc AX_rainvskalmanrain80-99

%% Basic analysis of the rain 
basicanalysis(rain,100)
print -depsc B6_basicanalysisrain
%% Trying to fit an ARMA-polynominal to pre-whithen input

rain_data = iddata(rain_m');
Arain = [1 0 0];
Crain = [1 zeros(1,35) -1];
rain_poly = idpoly(Arain, [], Crain);
rain_poly.structure.c.Free = [0 0 1 1 1 zeros(1,31) 1];
%rain_poly.structure.a.Free = [0 1 0 0 1];
rain_model = pem(rain_data,rain_poly);
rain_pw = resid(rain_model,rain_data);

%% Basic analysis and whiteness test of data

basicanalysis(rain_pw.y,50)
print -depsc B7_basicanalysisrainpw
figure
whitenessTest(rain_pw.y)
print -depsc B8_whitenesstestrainpw 
%% Pre-whiten nvdi data with the model found in prev. section, check crosscorr.

nvdi_m_data = iddata(nvdi_m');
nvdi_pw = resid(rain_model,nvdi_m_data); 

%%

basicanalysis(nvdi_pw.y,50) % Should not be white? 
crosscorrplot(40,rain_pw.y,nvdi_pw.y)
title('Cross-correlation between pre whitened input and pre whitened output')
print -depsc B9_crosscorrplotinputpw-outputpw

%% based on crosscorrplot, determine d, s, r
d=2;
s=1;
r=2;

A2 = [1 zeros(1,r)];
B = [zeros(1,d) ones(1,s) ];
%B = [0 1];
Mi = idpoly([1],[B],[],[],[A2]);
Mi.Structure.b.Free = [zeros(1,d) ones(1,s)];
z_pw = iddata(nvdi_pw.y,rain_pw.y);
Mba2 = pem(z_pw,Mi)
vhat=resid(Mba2,z_pw);
crosscorrplot(40,rain_pw.y,vhat.y)
title("Cross-correlation prewhitened rain and model residual")
print -depsc B10_crosscorrrainpwBA2resid

%% Forming the residual x = output - B(z)/A_2(z) * input

x = nvdi_m - filter(Mba2.b,Mba2.f,rain_m);
basicanalysis(x,30);
crosscorrplot(60,rain_m,x)
title('Cross correlation between external rain input and residual')
print -depsc B11_CCrainandresidual
% there were significant crosscorrelation even in lab using the hair dryer
% data so it is not unlikely that we have significant correlation between
% our input and the residual

%% Estimating ARMA process for the residual
x_data = iddata(x');
Cresid = 1;
Aresid = [1 0];
ARMA_poly=idpoly(Aresid,[],Cresid);
ARMA_model=pem(x_data,ARMA_poly);
x_res = resid(ARMA_model,x_data);
basicanalysis(x_res.y,30)
print -depsc B12_BAarmaresidx
%%
whitenessTest(x_res.y)
print -depsc B13_wtARMARESIDX
%% 
p=1;
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
z = iddata(nvdi_m',rain_m');
MboxJ = pem(z,Mi);
ehat=resid(MboxJ,z);
basicanalysis(ehat.y,50);
print -depsc B14_BABJresid
crosscorrplot(30,rain_m,ehat.y);
title("Cross correlation between model input and Box Jenkins residuals")
print -depsc B15_CCrainBJresid
present(MboxJ)
figure
whitenessTest(ehat.y)
print -depsc B16_WTBJresid
 
%% k-step prediction using ARMAX
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

%% Performing 1-step prediction 

[Fy_1, Gy_1] = diophantine(C_tilde,A_tilde,1); 
 
BF_1 = conv(Fy_1,B_tilde); 
[Fu_1,Gu_1] = diophantine(BF_1,C_tilde,1);
uhat_k1 = filter(Gu_1,C_tilde,rain_v); 
uhat_k1 = uhat_k1(max(length(Gu_1),length(C_tilde)):end);%remove samples
y1hat_k1 = filter(Gy_1,C_tilde,nvdi_v); 
y1hat_k1 = y1hat_k1(max(length(Gy_1),length(C_tilde)):end); %remove samples
yhat_k1_BJ = y1hat_k1(1:end)+uhat_k1;

figure()
plot(yhat_k1_BJ + mean_nvdi,'r')
hold on
plot(nvdi_v(length(nvdi_v) - length(yhat_k1_BJ) + 1:end) + mean_nvdi,'b');
title("1-step prediction Box Jenkins model")
legend('1-step prediction', 'Given data', 'Location', 'North')
axis([0 length(yhat_k1_BJ) 0 0.5])
print -depsc B17_1stepBJ

%% Calculating the prediction error 

ehat1_BJ = nvdi_v(length(nvdi_v)-length(yhat_k1_BJ)+1:end)-yhat_k1_BJ;  %y-yhat_k1
noise_var_BJ = var(ehat1_BJ);

%% 4-10 step prediction

pe_var_armax = zeros(1,length(kstepvector));
for i = 1:length(kstepvector) 
    
    % Prediction of input 
    [F, G] = diophantine(C_rain, A_rain, kstepvector(i)); 
    xhat = filter(G, C_rain, rain_v); 
    xhat = xhat(max(length(G), length(C_rain)):end); 
    
    % k-step prediction of output
    [Fy, Gy] = diophantine(C_tilde,A_tilde,kstepvector(i)); 
    BF = conv(Fy,B_tilde); 
    [Fu,Gu] = diophantine(BF,C_tilde,kstepvector(i));
    uhat = filter(Gu,C_tilde,rain_v); 
    uhat = uhat(max(length(Gu),length(C_tilde)):end);%remove samples
    y1hat = filter(Gy,C_tilde,nvdi_v); 
    y1hat = y1hat(max(length(Gy),length(C_tilde)):end); %remove samples
    
    yhat_np = y1hat(length(y1hat)-length(xhat)+1:end)+uhat(length(uhat)-length(xhat)+1:end);
    yhat_p = y1hat(length(y1hat)-length(xhat)+1:end)+uhat(length(uhat)-length(xhat)+1:end) + filter(Fu, 1, xhat);

    % Plot against original data
    figure()
    plot(yhat_p + mean_nvdi,'r')
    hold on
    plot(nvdi_v(length(nvdi_v) - length(yhat_p) + 1:end) + mean_nvdi,'b');
    title(string(kstepvector(i)) + '-step prediction with predicted input')
    axis([0 length(yhat_p) 0 0.5]) 
    legend('k-step prediction', 'Given data', 'Location', 'North')
    
%     figure()
%     plot(yhat_np + mean_nvdi,'r')
%     hold on
%     plot(nvdi_v(length(nvdi_v) - length(yhat_np) + 1:end) + mean_nvdi,'b');
%     title(string(kstepvector(i)) + '-step prediction without predicted input')
    
    
    pe_var_armax(i) = var( nvdi_v(length(nvdi_v) - length(yhat_p) + 1:end) - yhat_p); 
end

%% Testing prediction with a naive model 

% Creating AR(1)-process for nvdi_m 

A_naive = [1 0]; 
data_naive = iddata(nvdi_m'); 
poly_naive = idpoly(A_naive, [], 1); 
param_naive = pem(data_naive, poly_naive); 
resid_naive = resid(data_naive, param_naive); 

figure
whitenessTest(resid_naive.y) 
basicanalysis(resid_naive.y, 50) 

%% Predicting with naive model 

% 1-step prediction 

[F_naive, G_naive] = diophantine(1, param_naive.A, 1); 
y_naive_1 = filter(G_naive, 1, nvdi_v); 
y_naive_1 = y_naive_1(length(G_naive):end); 

plot(y_naive_1 + mean_nvdi, 'b')
hold on
plot(nvdi_v(length(nvdi_v) - length(y_naive_1) + 1:end) + mean_nvdi, 'r')
hold off 
title('1-step prediction using naive AR(1) process') 

naive_err_1 = nvdi_v(length(nvdi_v) - length(y_naive_1)) - y_naive_1; 
naive_var_1 = var(naive_err_1); 

%% 4-10 step prediction using naive AR(1)-process 

kstepvector = 4:10;
% P_error_ARMA = kstepb(model_init_eg_p,kstepvector,nvdi_v,noise_var, nvdi_mean);
pe_var_naive = zeros(1,length(kstepvector)); 
for i = 1:length(kstepvector) 
    % Performing polynominal division 
    [Fk, Gk] = diophantine(param_naive.C,param_naive.A,kstepvector(i)); 
    
    % Performing k-step prediction 
    yhat = filter(Gk, param_naive.C, nvdi_v); 
    
    % Removing some samples
    yhat = yhat(max(length(Gk), length(param_naive.C):end));
    
    % Making data equal length 
    tempData = nvdi_v(length(nvdi_v) - length(yhat) + 1:end); 
   
    % Plotting the k-step prediction 
    figure
    plot(tempData + mean_nvdi, 'b')
    hold on
    plot(yhat + mean_nvdi, 'r') 
    title(string(kstepvector(i))+ '-step naive prediction (red) versus real data (blue)')
    
    %Saving prediction error
    ehat = nvdi_v(length(nvdi_v)-length(yhat)+1:end)-yhat;  %y-yhat_k1
    %basicanalysis(ehat, 50); 
    pe_var_naive(i) = var(ehat);
end

