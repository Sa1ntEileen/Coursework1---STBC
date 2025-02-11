`% A MUD receiver, which recovers 4 bits for each of 6 users from the signal
% received using 6 antennas and 4 timeslots. The 4 users' STBC signals interfere
% with each other and are impaired by quasi-static fading channels and by noise.
% Y - this is a 4x6 matrix of complex received symbols, where the 4 rows correpsond
%    to the 4 timeslots and the 6 columns correspond to the 6 receive antennas
% H - this is a 16x6 matrix of complex channel coefficients, where the 16
%    rows correspond to the 4 transmit antennas of the 4 users, while the 6
%    columns correspond to the 6 receive antennas
% B_hat - this is a 6x4 matrix of decoded bits, where the 6 rows correspond
%    to the 6 users and the 4 columns correpond to the 4 bits of each user
% rx_title - this is a string that describes the operation of the receiver
function [B_hat, rx_title] = receiver(Y, H)

% Describe the operation of the receiver
rx_title = 'A QPSK Alamouti STBC signal is received using 1 antenna and 2 timeslots';

% Set the parameters that are assumed by the receiver
K = 6; % number of users
L = 4; % number of bits per user
T = 4; % number of timeslots
P = 4; % number of transmit antennas per user
Q = 6; % number of receive antennas

% Check that the dimensions of the inputs agree with the parameters above
if ~isequal(size(Y), [T, Q])
    error('Soton:argChk','~isequal(size(Y), [T, Q])');
end
if ~isequal(size(H), [P*K,Q])
    error('Soton:argChk','~isequal(size(H), [P*K,Q])');
end

% Build matrices of the form Y2 = H2*X2+N2
Y2 = [Y(1,1); conj(Y(2,1))];
H2 = [H(1,1), H(2,1); conj(H(2,1)), -conj(H(1,1))];

% Perform single-user detection
Z2 = H2'*Y2;

% Zero will be used for the decoded bits of most users
B_hat = zeros(K, L);

% QPSK demodulation for user 1
B_hat(1,1) = imag(Z2(1,1)) < 0;
B_hat(1,2) = real(Z2(1,1)) < 0;
B_hat(1,3) = imag(Z2(2,1)) < 0;
B_hat(1,4) = real(Z2(2,1)) < 0;
  
end