clear; clc; close all

n = 1000;
w = logspace(-1, 2, n); 
R = 15;

% 1. Use a separate variable name for the vector
u_vec = 0.1:0.1:100; 

N = zeros(length(u_vec), n); % Pre-allocate matrix

for k = 1:length(u_vec)
    % 2. Extract scalar into u_rle from u_vec
    u_rle = u_vec(k); 
    
    U_W = u_rle * w;
    arg_clamped = min(1, (pi/2) * R ./ U_W);
    
    N_temp = (4 / pi) * (R ./ U_W) .* exp(-1i * acos(arg_clamped));
    N_temp(U_W <= R) = 1;
    
    N(k, :) = N_temp;
end