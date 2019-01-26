%% 3.1
clear all
n=500;
A1 = [1 -0.65]; 
A2 = [1 0.90 0.78];
C = 1; 
B = [0 0 0 0 0.4];
e = sqrt(1.5)*randn(n + 100, 1);
w = sqrt(2)*randn(n + 200, 1);
A3 = [1 .5]; 
C3 = [1 -.3 .2];
u = filter(C3,A3,w); 
u = u(101:end);
y = filter(C,A1,e) + filter(B,A2,u);
u = u(101:end);
y = y(101:end);
clear A1; clear A2; clear C; clear B; clear e; clear w; clear A3; clear C3;

%% task1: determine A3 & C3
basicanalysis(u,30)
%%
u_data = iddata(u);
u_poly = idpoly([1 0], [], [1 0 0]);
u_poly.structure.a.Free = [0 1];
u_model = pem(u_data,u_poly)
u_pw = filter(u_model.a,u_model.c,u);
u_pw = u_pw(10:end);
basicanalysis(u_pw,30)
% Q1 Started of with AR(1),peak at 3 eliminated with MA(3)

%% task2
y_pw = filter(u_model.a,u_model.c,y);
y_pw = y_pw(10:end);
basicanalysis(y_pw,30)
crosscorrplot(40,u_pw,y_pw)
%%
%Q2 From the cross-correlation plot, we can guess that (d=4, s=8, r=2)
d=4;
s=1;
r=2;

A2 = [1 zeros(1,r)];
B = [zeros(1,d) ones(1,s) ];
Mi = idpoly([1],[B],[],[],[A2]);
Mi.Structure.b.Free = [zeros(1,d) ones(1,s)];
z_pw = iddata(y_pw,u_pw);
Mba2 = pem(z_pw,Mi)
vhat=resid(Mba2,z_pw);
basicanalysis(vhat.y,30);
crosscorrplot(40,u_pw,vhat.y)

%% task 3
x = y - filter(Mba2.b,Mba2.f,u);
basicanalysis(x,30);
crosscorrplot(30,u,x) %???? ask at lab
% u and x show high correlation
%%
x_data = iddata(x);
C = [1];
A = [1 0];
ARMA_poly=idpoly([A],[],[C]);
ARMA_model=pem(x_data,ARMA_poly);
%ARMA_model.structure.c.Free = [1 zeros(1,3) 1];
x_res = filter(ARMA_model.a,ARMA_model.c,x);
basicanalysis(x_res,30)

%Q3 started of with AR(2), added ARMA(4,4) to eliminate peak at 4
% dependence from u(t) not removed from x, tried various combinations

%% task 4
% (p,q,r,s,d = 1,0,2,1,4)

p=1;
q=0;
r=2;
s=1;
d=4;

A1 = [1 zeros(1,p)];
A2 = [1 zeros(1,r)];
B = [zeros(1,d) zeros(1,s)];
C = [1 zeros(1,q)];
Mi = idpoly(1,B,C,A1,A2);
Mi.structure.b.Free = [zeros(1,d), ones(1,s)];
z = iddata(y,u);
MboxJ = pem(z,Mi);
ehat=resid(MboxJ,z);
basicanalysis(ehat.y,30);
crosscorrplot(30,u,ehat.y);
present(MboxJ)

% The correlation is low, with a few peaks slightly outside
% the confint. some parameters not sign. diff from 0

%% 3.2
clear all
close all;
load tork.dat
tork = tork -repmat(mean(tork),length(tork),1);
y = tork(:,1); 
u = tork(:,2);
z = iddata(y,u);
figure()
plot(z(1:300))

u = u(1:300);
y = y(1:300);

%% task1
close all;
basicanalysis(u,30)
%%
u_data = iddata(u);
u_poly = idpoly([1 0],[],[]);
u_model = pem(u_data,u_poly);
u_pw = filter(u_model.a,u_model.c,u);
basicanalysis(u_pw,30)
%%
y_pw = filter(u_model.a,u_model.c,y);
basicanalysis(y_pw,40)
crosscorrplot(40,u_pw,y_pw)
%%
d=3;
s=2;
r=2;

A2 = [1 zeros(1,r)];
B = [zeros(1,d) ones(1,s+1) ];
Mi = idpoly([1],[B],[],[],[A2]);
Mi.Structure.b.Free = [zeros(1,d) ones(1,s+1)];
z_pw = iddata(y_pw,u_pw);
Mba2 = pem(z_pw,Mi)
vhat=resid(Mba2,z_pw);
basicanalysis(vhat.y,30);
crosscorrplot(40,u_pw,vhat.y)

whitenessTest(vhat.y)

%%
x = y - filter(Mba2.b,Mba2.f,u);
basicanalysis(x,30);
crosscorrplot(30,u,x)

%%
x_data = iddata(x);
C = 1;
A = [1 0];
ARMA_poly=idpoly([A],[],[C]);
ARMA_model=pem(x_data,ARMA_poly);
x_res = filter(ARMA_model.a,ARMA_model.c,x);
basicanalysis(x_res,30)
%% 
p=1;
q=0;
r=2;
s=2;
d=3;

A1 = [1 zeros(1,p)];
A2 = [1 zeros(1,r)];
B = [zeros(1,d) zeros(1,s+1)];
C = 1;
Mi = idpoly(1,B,C,A1,A2);
Mi.structure.b.Free = [zeros(1,d), ones(1,s+1)];
z = iddata(y,u);
MboxJ = pem(z,Mi);
ehat=resid(MboxJ,z);
basicanalysis(ehat.y,30);
crosscorrplot(30,u,ehat.y);
present(MboxJ)


%Q5
%Delay 0.08*d = 0.24s
%The residual is white, (not really )uncorrelated and parameters
%sign. different from zero


%% 3.3 ARMA prediction of temperature in Svedala
clear all
load svedala
y=svedala;
A = [1 -1.79 0.84];
C = [1 -0.18 -0.11];

%% calculating residual e1
% Q6 putting k=1 to get noise variance = 0.3751;
k1=1;
[Fk1,Gk1] = diophantine(C,A,k1);
yhat_k1 = filter(Gk1, C, y);
y1 = y(max(length(Gk1),length(C)):end);
yhat_k1 = yhat_k1(max(length(Gk1),length(C)):end);
ehat1 = y1-yhat_k1;
noise_var = var(ehat1);
%% Calculating yhat for k=3 & k=26

k3=3;
k26=26;

[Fk3,Gk3] = diophantine(C,A,k3);
yhat_k3 = filter(Gk3, C, y);
y3 = y(max(length(Gk3),length(C)):end);
yhat_k3 = yhat_k3(max(length(Gk3),length(C)):end);

[Fk26,Gk26] = diophantine(C,A,k26);
yhat_k26 = filter(Gk26, C, y);
y26 = y(max(length(Gk26),length(C)):end);
yhat_k26 = yhat_k26(max(length(Gk26),length(C)):end);

%% Q7

ehat3 = y3-yhat_k3;
ehat26 = y26-yhat_k26;

%1
m3 = mean(ehat3); %=-0.0071
m26 = mean(ehat26); %=-0.0435

%2 Expectation is zero

%% Q8
err_var_estk3 = var(ehat3); %=2.5937
err_var_estk26 = var(ehat26); %=10.6925
%estimated

var3 = noise_var*norm(Fk3)^2; %=2.7475
var26 = noise_var*norm(Fk26)^2; %=12.6179
%theoretical

%% Q9
ci_3high = m3 + sqrt(var3)*1.96; 
ci_3low = m3 - sqrt(var3)*1.96;


ci_26high = m26 + sqrt(var26)*1.96;
ci_26low = m26 - sqrt(var26)*1.96;


nbrError3 = sum(abs(ehat3)>ci_3high) + sum(abs(ehat3)<ci_3low);
p_Error3 = nbrError3/length(ehat3);

nbrError26 = sum(abs(ehat26)>ci_26high) + sum(abs(ehat26)<ci_26low);
p_Errork26 = nbrError26/length(ehat26);

% Q10 see workspace

%% Q11
close all;
figure()
plot(y3)
hold on
plot(yhat_k3)
hold off

figure()
plot(y26)
hold on
plot(yhat_k26)
hold off

% cov_est3 = covf(ehat3,30);
% figure()
% stem(cov_est3)
% hold on
% plot(0:30,ones(1,31)*ci_3high);
% plot(0:30,ones(1,31)*ci_3low);

basicanalysis(ehat3,30)

% Residuals don't look like MA(k-1) processes, acf peaking around 24

%% 3.4 ARMAX prediction using sturup as input

% Q12 Delay?

load sturup;
u = sturup;
A_u = [1 -1.49 0.57];
B_u = [0 0 0 0.28 -0.26];
C_u = [1];

[F_u3,G_u3] = diophantine(C_u,A_u,k3);
[F_u26,G_u26] = diophantine(C_u,A_u,k26);

BF3 = conv(B_u,F_u3);
BF26 = conv(B_u,F_u26);
[Fhat3,Ghat3] = diophantine(BF3,C_u,k3);
[Fhat26,Ghat26] = diophantine(BF26,C_u,k26);

%%
% Computing yhat for k=3;
uhat_k3 = filter(Ghat3,C_u,u); 
uhat_k3 = uhat_k3(max(length(Ghat3),length(C_u)):end);%remove samples
y1hat_k3 = filter(Gk3,C_u,y); 
y1hat_k3 = y1hat_k3(max(length(Ghat3),length(C_u)):end); %remove samples
yhat_k3 = y1hat_k3+uhat_k3;
%%
% Computing yhat for k=26;
uhat_k26 = filter(Ghat26,C_u,u); 
uhat_k26 = uhat_k26(max(length(Ghat26),length(C_u)):end); %remove samples
u_k26 = u(max(length(Ghat26),length(C_u)):end);
y1hat_k26 = filter(Gk26,C_u,y); 
y1hat_k26 = y1hat_k26(max(length(Ghat26),length(C_u)):end); %remove samples
yhat_k26 = y1hat_k26+uhat_k26;

%%
% Plot
close all;
figure()
plot(y((length(y)-length(yhat_k3)+1):end),'b')
hold on
plot(yhat_k3,'r')
hold off

figure()
plot(y((length(y)-length(yhat_k26)+1):end),'b')
hold on
plot(yhat_k26,'r')
hold off

%% Q13

ARMAX_ehat3 = y(length(y)-length(yhat_k3)+1:end)-yhat_k3; %y-yhat;
ARMAX_ehat26 = y(length(y)-length(yhat_k26)+1:end)-yhat_k26; %y-yhat

%1
ARMAX_m3 = mean(ARMAX_ehat3);
ARMAX_m26 = mean(ARMAX_ehat26);

ARMAX_var3 = var(ARMAX_ehat3);
ARMAX_var26 = var(ARMAX_ehat26);

basicanalysis(ARMAX_ehat3,30)

%% 3.5 SARIMA process
close all;
load svedala
plot(svedala)
S = 24;
AS = [1 zeros(1,S-1) -1];
ys = filter(AS,1,svedala);
ys = ys(S+1:length(ys)); % remove samples

basicanalysis(ys,30)

%adding AR(2) and c24 coefficient

ys_poly = idpoly([1 0 0], [], [1 zeros(1,23) -1],[],[]);
ys_poly.Structure.c.Free = [0 zeros(1,23) 1];

ysdata = iddata(ys);
model_init = pem(ysdata,ys_poly);
ys_res = resid(model_init,ysdata);
basicanalysis(ys_res.y,30)
%Q15 AR2 c24
%%
Astar = conv(AS,model_init.a); 
model2 = model_init; %creating a copy
model2.a = Astar; % replace A polynomial with SARIMA
svedala_data = iddata(svedala);
svedala_res = resid(model2,svedala_data);


A = model2.a;
C = model2.c;

[Fk3s,Gk3s] = diophantine(C,A,k3);
yhat_sk3 = filter(Gk3s,C,svedala);
yhat_sk3 = yhat_sk3(max(length(Gk3s),length(C)):end);


[Fk26s,Gk26s] = diophantine(C,A,k26);
yhat_sk26 = filter(Gk26s,C,svedala);
yhat_sk26 = yhat_sk26(max(length(Gk26s),length(C)):end);

close all;
figure()
plot(svedala(length(svedala)-length(yhat_sk3)+1:end),'b')
hold on
plot(yhat_sk3,'r')
hold off

close all;
figure()
plot(svedala(length(svedala)-length(yhat_sk26)+1:end),'b')
hold on
plot(yhat_sk26,'r')
hold off
%%
SARIMA_ehat3  = svedala(length(svedala)-length(yhat_sk3)+1:end) -yhat_sk3;
SARIMA_ehat26 = svedala(length(svedala)-length(yhat_sk26)+1:end) -yhat_sk26;

SARIMA_m3 = mean(SARIMA_ehat3);
SARIMA_m26 = mean(SARIMA_ehat26);

SARIMA_var3 = var(SARIMA_ehat3);
SARIMA_var26 = var(SARIMA_ehat26);
%%
basicanalysis(SARIMA_ehat3,30)
