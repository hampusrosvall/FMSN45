function [ z ] = kalmanfunction( A, Re, Rw, C, y )
%KALMANFUNCTION Summary of this function goes here
%   Detailed explanation goes here

% Length of process
N = length(y);

% Set initial values
Rzz_1 = 10*eye(3); % Initial variance
ztt_1 = [0 0 0]'; % Initial state

% Vector to store values in
z=zeros(3,N);

% Kalman filter. Start from k=3, since we need old values of y.
for k=3:N
  
  
  if mod(k,3) > 0
      % Update with missing y-sample
      ztt = ztt_1; 
      Rzz = Rzz_1; 
  else
       % Update with available y-sample 
      Ryy = C*Rzz_1*C' + Rw;
      Kt = Ryy\(Rzz_1*C');
      ztt = ztt_1 + (Kt*(y(k) - C*ztt_1));
      Rzz = (eye(3)-Kt*C)*Rzz_1;
  end
 

     for i=1:3
         if(ztt(i)<0)
             ztt(i)=0;
         end
     end

 
  % Save
  z(:,k) = ztt;
  %1-step prediction
  
  ztt_1 = A*ztt;
  Rzz_1 = A*Rzz*A' + Re;
  
  
end

  
 
 

end

