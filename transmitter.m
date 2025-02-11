% The kth user's transmitter, which converts 4 bits into an STBC signal, 
% which is transmitted using 4 antennas and 4 timeslots
% b_k - this is the kth user's vector of four bits
% k - this is the index of the user, in the range 1 to 4
% X_k - this is a 4x4 matrix of complex STBC symbols, where the 4 rows correpsond 
%    to the 4 timeslots and the 4 columns correspond to the 4 transmit antennas
% tx_title - this is a string that describes the operation of the transmitter
function [X_k, tx_title] = transmitter(b_k, k)

% Describe the operation of the transmitter
tx_title = 'A QPSK Alamouti STBC signal is transmitted by 1 user using 2 antennas and 2 timeslots';

% Set the parameters that are assumed by the transmitter
L = 4; % number of bits per user
T = 4; % number of timeslots
P = 4; % number of transmit antennas per user

% Check that the dimensions of the input agree with the parameters above
if length(b_k) ~= L
    error('Soton:argChk','length(b_k) ~= L');
end

% Nothing will be transmitted by most users on most antennas in most timeslots
X_k = zeros(T, P);

% Only user 1 transmits 
if k == 1
    % Determine the transmissions in the first timeslot using BPSK 
    X_k(1,1) = (1-2*b_k(1,1)); 
    X_k(1,2) = (1-2*b_k(1,2));
    X_k(1,3) = (1-2*b_k(1,3));
    X_k(1,4) = (1-2*b_k(1,4));
    
    % Determine the transmissions in the second timeslot using quasi-orthogonal STBC
    X_k(2,1) = -conj(X_k(1,2)); 
    X_k(2,2) = conj(X_k(1,1));
    X_k(2,3) = -conj(X_k(1,4));
    X_k(2,4) = conj(X_k(1,3));
    
    % third time slot in quasi-orthoghonal 
    X_k(3,1) = -conj(X_k(1,3)); 
    X_k(3,2) = -conj(X_k(1,4));
    X_k(3,3) = conj(-X_k(1,1));
    X_k(3,4) = conj(X_k(1,2));

    % fourth time slot in quasi-orthoghonal 
    X_k(4,1) = (X_k(1,4)); 
    X_k(4,2) = -(X_k(1,3));
    X_k(4,3) = -(X_k(1,2));
    X_k(4,4) = (X_k(1,1));
    
end
end



