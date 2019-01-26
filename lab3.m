%% Task 3.1 Recursive LS estmation
clear
load tar2.dat
load thx.dat
subplot(311)
plot(tar2)
subplot(312)
plot(thx)
na = 2;
nb = 0;
nk = 0;
model = [na]; %[na,nb,nk]
lambda = 0.95;
[Aest,yhat,covAest,yprev]= rarx(tar2,model,'ff',lambda);
subplot(313)
plot(Aest)

%%
n = 100;
lambda_line = linspace(0.85,1,n);
ls2 = zeros(n,1);
for i=1:length(lambda_line)
    [Aest,yhat,CovAest,trash] = rarx(tar2,[2],'ff',lambda_line(i));
    ls2(i) = sum((tar2-yhat).^2);
end
figure
plot(lambda_line,ls2)

% Q2- ls2 vector contains quadratic error for different values of lambda. minimum at 0.95

%% initial guess for Task 3.2

A = [1 0 0];
C = 1;
tar2init = tar2(1:20);
tar2initdata = iddata(tar2init);
model_init = idpoly(A,[],C);
model_param = pem(tar2initdata,model_init)

%% Task 3.2 Kalman filtering of time series

%   Simulate data
y=tar2;

% Length of process
N = length(y);

% Define the state space equations
A = [1 0; 0 1];
Re = 0.001*[1 0; 0 0]; % Hidden state noise covariance matrix
Rw = 1.25; % Observation variance, (what does it do?)

%usually C should be set here to
%but in this case C is a function of time.



% Set initial values

Raa_1 = 0.04*eye(2); % Initial variance
att_1 = [model_param.a(2) model_param.a(3)]'; % Initial state, 

% Vector to store values in
z=zeros(2,N);

% Kalman filter. Start from k=3, since we need old values of y.
for k=3:N
  % C is a function of time.
  C = [-y(k-1) -y(k-2)];
   
  % Update
  Ryy = C*Raa_1*C' + Rw; % scalar? 1 x 1
  Kt = (Raa_1*C')/Ryy; % 2 x 1
  att = att_1 + (Kt*(y(k) - C*att_1)); % 2 x 1
  Raa = (eye(2)-Kt*C)*Raa_1; % 2 x 2 
 
  % Save
 z(:,k) = att;
  
  %1-step prediction
  
  att_1 = A*att;
  Raa_1 = A*Raa*A' + Re;
  
  
end

%% plotting

subplot(311) %True
plot(thx)
title('True Values')

subplot(312) %RLS lambda = 0.95
plot(Aest)
title('RLS lambda = 0.95')

subplot(313) %Kalman
plot(z')
title('Kalman')

%% Task 3.3 Quality control of a process
%% generate u (input) as a markov chain

P = [7/8 1/8; 1/8 7/8];
chain_length = 1000;
u = markovsimul(P,chain_length);

%%


close all
b = 20;
sigma2_e = 1;
sigma2_w = 4;


e_t = 1*randn(500,1);
v_t = 2*randn(500,1);

%generating x
x = zeros(500,1);
x(1) = e_t(1);
for i=2:length(x)
    x(i) = x(i-1)+e_t(i);
end

%generating y
y = zeros(500,1);
for i=1:length(y)
    y(i) = x(i) + b*u(i) + v_t(i);
end


%% Kalman filter
%   Simulate data

% Length of process
N = length(y);

% Define the state space equations
A = [1 0; 0 1];
Re = [1 0; 0 0];
Rw = 1.25; 

% Set initial values

Rzz_1 = var(y(1:5))*eye(2); % Initial variance, why take variance of y? 
ztt_1 = [0 0]'; % Initial state, 

% Vector to store values in
z=zeros(2,N);

% Kalman filter. Start from k=3, since we need old values of y.
for k=2:N
  % C is a function of time.
  C = [1 u(k)];
   
  % Update
  Ryy = C*Rzz_1*C' + Rw;
  Kt = (Rzz_1*C')/Ryy;
  ztt = ztt_1 + (Kt*(y(k) - C*ztt_1));
  Rzz = (eye(2)-Kt*C)*Rzz_1;
 
  % Save
  z(:,k) = ztt;
  
  %1-step prediction
  
  ztt_1 = A*ztt;
  Rzz_1 = A*Rzz*A' + Re;
  
  
end


%%
% close all;
% Re = [sigma2_e 0; 0 0];
% Rw = sigma2_w;
% A = [1 0; 0 1];
% C = [1 u];
% z = kalmanfunction(A,Re,Rw,u,y);
figure()
plot(z(1, :),'b')
hold on
plot(z(2, :), 'r')
hold on
plot(x,'g')
%% 3.4 Recursive temperature modeling using RLS
close all;
load svedala94.mat
plot(svedala94)
Diff = [1 zeros(1,5) -1];
y_diff = filter(Diff,1,svedala94);
y_diff = y_diff(length(Diff)+1:end);
T = linspace(datenum(1994,1,1),datenum(1994,12,31),length(y_diff));
subplot(211)
plot(T,svedala94(length(svedala94)-length(y_diff)+1:end));
datetick('x');
subplot(212)
plot(T,y_diff);
datetick('x');

%% Q.2 
th = armax(y_diff,[2 2]);
th_winter = armax(y_diff(1:540-length(Diff)),[2 2]);
th_summer = armax(y_diff(907-length(Diff):1458),[2 2])

%% Q.3 

th0 = [th_winter.A(2:end) th_winter.C(2:end)]; 
[thr, yhat] = rarmax(y_diff, [2 2], 'ff', 0.99, th0); 

subplot(311) 
plot(T, svedala94(length(svedala94)-length(y_diff)+1:end)); 
datetick('x')

subplot(312) 
plot(thr(:, 1:2))
hold on 
plot(repmat(th_winter.A(2:end), [length(thr) 1]), 'b:'); 
plot(repmat(th_summer.A(2:end), [length(thr) 1]), 'r:'); 
axis tight
hold off 

subplot(313) 
plot(thr(: ,3:end))
hold on 
plot(repmat(th_winter.C(2:end), [length(thr) 1]), 'b') 
plot(repmat(th_summer.C(2:end), [length(thr) 1]), 'r:')
axis tight
hold off

% The recursive estimation 
% What is happening when the parameters switch? 
% Better coincide at summer 
%% 3.5 Recursive temperature modeling using varying means
clear
close 
load svedala94 
start = 850; 
finish = 1100; 
y = svedala94(start:finish) - mean(svedala94(start:finish)); 
t = (1:length(y))'; 
U = [sin(2*pi*t/6) cos(2*pi*t/6)]; % Why division by six and not by any other number? 
Z = iddata(y, U); 
model = [3 [1 1] 4 [0 0]]; %[na [nb1 nb2] nc [nk1 nk2]];

thx = armax(Z,model); 

plot(y, 'b')
hold on
plot(U*cell2mat(thx.b)', 'r')

% There seems to be a difference in phase between the sinosodals and y

%% 

U = [sin(2*pi*t/6) cos(2*pi*t/6) ones(size(t))]; 
Z = iddata(y,U); 
m0 = [thx.A(2:end) thx.B 0 thx.C(2:end)]; % Initial values
m0 = cell2mat(m0); 
Re = diag([zeros(1,5) 1 zeros(1,4)]); % how do I initialize these? 
model = [3 [1 1 1] 4 0 [0 0 0] [1 1 1]]; %na nb nc nk (k is delay)
[thr, yhat] = rpem(Z, model, 'kf', Re, m0); 

m = thr(:, 6); 
a = thr(end, 4); 
b = thr(end, 5); 
y_mean = m + a * U(:,1) + b*U(:,2); 
y_mean = [0; y_mean(1:end-1)]; 

plot(y, 'b')
hold on
plot(y_mean, 'r')

% y and y_mean are very similar, but it seems like y_mean is oscillating
% more

%% Q10 

y = svedala94; 
y = y - y(1); 
t = (1:length(y))'; 
U = [sin(2*pi*t/6) cos(2*pi*t/6)]; % Why division by six and not by any other number? 
Z = iddata(y, U); 
model = [3 [1 1] 4 [0 0]]; %[na [nb1 nb2] nc [nk1 nk2]];

thx = armax(Z,model); 

%% 

U = [sin(2*pi*t/6) cos(2*pi*t/6) ones(size(t))]; 
Z = iddata(y,U); 
m0 = [thx.A(2:end) thx.B 0 thx.C(2:end)]; % Initial values
m0 = cell2mat(m0); 
Re = diag([zeros(1,5) 1 zeros(1,4)]); % how do I initialize these? 
model = [3 [1 1 1] 4 0 [0 0 0] [1 1 1]]; %na nb nc nk (k is delay)
[thr, yhat] = rpem(Z, model, 'kf', Re, m0); 

m = thr(:, 6); 
a = thr(end, 4); 
b = thr(end, 5); 
y_mean = m + a * U(:,1) + b*U(:,2); 
y_mean = [0; y_mean(1:end-1)]; 

plot(y, 'b')
hold on
plot(y_mean, 'r')
