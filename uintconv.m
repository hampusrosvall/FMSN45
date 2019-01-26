function [ sorted ] = uintconv( uint )
%UINTCONV Summary of this function goes here
%   Detailed explanation goes here

sorted = zeros(1, length(uint)); 
for i = 1:length(uint)
    sorted(i) = uint(i) / (255/2) - 1;
end

