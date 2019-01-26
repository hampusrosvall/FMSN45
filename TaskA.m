
% Loading data 
close all;
load proj18.mat 
rain_org = ElGeneina.rain_org;



%% % Filling missing samples with zeros
rain_org_mod = zeros(1, 3 * length(rain_org));
k = 3;
i = 1; 
while i < length(rain_org)
    rain_org_mod(k) = rain_org(i); 
    k = k + 3; 
    i = i + 1; 
end

close all;
clear i;
clear k;

%% Initializing Kalman parameters

close all;
Ckalman = [1 1 1];
Re = [ 1 0 0; 0 0 0; 0 0 0]; 
Rw = 0.001; % should be low, 
Akalman = [1 0 0; 1 0 0; 0 1 0];
z = kalmanfunction(Akalman,Re,Rw,Ckalman,rain_org_mod);



%% creating a vector with all the x:es
i=1;
k=3;
reconstructed_z = zeros(1,length(z)/3);
while i < length(z)
    reconstructed_z(i) = z(3,k);
    reconstructed_z(i+1) = z(2,k);
    reconstructed_z(i+2) = z(1,k);
    k=k+3;
    i=i+3;
end




%% plotting against rain
close all;
figure()
stem(rain_org_mod,'r')
hold on
stem(reconstructed_z,'b')
title('Reconstructed rain plotted against original rain year 1987-1989')
hold off
axis([972 1044 0 250])
legend('Original rain data', 'Reconstructed rain data', 'Location', 'North')
print -depsc A2_rain_plot2

sum(reconstructed_z)
sum(rain_org)

figure()
plot(reconstructed_z,'b')
hold on
plot(ElGeneina.rain,'r')
title('Reconstructed rain(blue) plotted against linearly interpolated rain(red)')
hold off
print -depsc A1_rain_plot1

% %% Finding the optimal Re-parameter in AR(1)-process
% Re_test = linspace(0.001,1, 1000); % To ensure stable system
% err = zeros(1,length(Re_test)); % To store errors for every a
% sq_err = zeros(1,length(Re_test)); % Store squared errors for every a
% % Creating squared error for every rain_org point i.e
% % y_i - (x_i + x_i-1 + x_i-2)
% for i=1:length(Re_test)
%     Re = [Re_test(i) 0 0; 0 0 0; 0 0 0];
%     z_test = kalmanfunction(Akalman,Re,Rw,Ckalman,rain_org_mod);
%     j=1;
%     k=3;
%     while j < length(z_test)/3
%         err(j) = (rain_org_mod(k) - (z_test(3,k) + z_test(2,k) + z_test(1,k)))^2;
%         k=k+3;
%         j=j+1;
%     end
%     sq_err(i) = sum(err);
% end
% 
% figure
% stem(sq_err)
% 
% [val, index] = min(sq_err);
% 
% % %% Finding the optimal a-parameter in AR(1)-process
% A_test = linspace(-1,1, 100); % To ensure stable system
% err = zeros(1,length(A_test)); % To store errors for every a
% sq_err = zeros(1,length(A_test)); % Store squared errors for every a
% % Creating squared error for every rain_org point i.e
% % y_i - (x_i + x_i-1 + x_i-2)
% for i=1:length(A_test)
%     A_temp = [A_test(i) 0 0; 1 0 0; 0 1 0];
%     z_test = kalmanfunction(A_temp,Re,Rw,Ckalman,rain_org_mod);
%     j=1;
%     k=3;
%     while j < length(z_test)/3
%         err(j) = (rain_org_mod(k) - (z_test(3,k) + z_test(2,k) + z_test(1,k)))^2;
%         k=k+3;
%         j=j+1;
%     end
%     sq_err(i) = sum(err);
% end
% 
% figure
% stem(sq_err)
% 
% [val, index] = min(sq_err);
% 
% %% Finding the optimal sigmaw-parameter in AR(1)-process
% Rw_test = linspace(0.001,1, 1000); 
% err = zeros(1,length(Rw_test)); % To store errors 
% sq_err = zeros(1,length(Rw_test)); % Store squared errors
% % Creating squared error for every rain_org point i.e
% % y_i - (x_i + x_i-1 + x_i-2)
% for i=1:length(Rw_test)
%     Rw_temp = Rw_test(i);
%     z_test = kalmanfunction(Akalman,Re,Rw_temp,Ckalman,rain_org_mod);
%     j=1;
%     k=3;
%     while j < length(z_test)/3
%         err(j) = (rain_org_mod(k) - (z_test(3,k) + z_test(2,k) + z_test(1,k)))^2;
%         k=k+3;
%         j=j+1;
%     end
%     sq_err(i) = sum(err);
% end
% 
% figure
% stem(sq_err)
% 
% [val, index] = min(sq_err);
% 
