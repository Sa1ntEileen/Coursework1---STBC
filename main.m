clear all

K = 6; % number of users
L = 4; % number of bits per user
T = 4; % number of timeslots
P = 4; % number of transmit antennas per user
Q = 6; % number of receive antennas

inv_N_0_dB = -10:1:25; % Choose the x axis for our BER plot
min_BER_to_plot = 1e-3; % Choose the y axis for our BER plot
min_bit_errors_to_plot = 100; % Adjust this to control how smooth your BER plot is and how long your simulation takes to run

% Create a figure to plot the results.
figure
axes1 = axes('YScale','log');
ylabel('BER');
xlabel('1/N_0 (in dB)');
xlim([inv_N_0_dB(1), inv_N_0_dB(end)]);
ylim([min_BER_to_plot,1]);
hold on
plots = zeros(size(K));
key = cell(K,1);
cmap = colormap('Lines');
for k = 1:K
    plots(k) = plot(nan,'Parent',axes1,'Color',cmap(k,:));
    key{k} = ['User ',num2str(k)'];
end
leg1 = legend(key,'Location','SouthWest');

% Counters to store the number of bits and errors simulated so far
bit_counters=zeros(length(inv_N_0_dB),1);
bit_error_counters=zeros(length(inv_N_0_dB),K);

% Variables used to check if the average transmission energy per bit of each user is too high
E_b = zeros(K,1);
E_b_too_high = 0;

% Loop over the Eb/N0 values
for inv_N_0_dB_index = 1:length(inv_N_0_dB)
    
    % Calculate the noise power spectral density
    N_0 = 10^-(inv_N_0_dB(inv_N_0_dB_index)/10);
    
    % Continue the simulation until enough bit errors have been observed
    while min(bit_error_counters(inv_N_0_dB_index, :)) < min_bit_errors_to_plot
        
        % Increment the bit counter for this N_0 value
        bit_counters(inv_N_0_dB_index) = bit_counters(inv_N_0_dB_index) + L;
        
        % Generate 4 random bits for each of the 4 users
        B = round(rand(K, L));
        
        % Build the X matrix of the 4 user's STBC signals
        X = zeros(T, P*K);
        for k = 1:K
            
            % Call transmitter.m to obtain the STBC signal transmitted by each user
            X_k = transmitter(B(k, :),k);
            
            % Check that transmitter.m produces an output having the correct dimensions
            if ~isequal(size(X_k),[T,P])
                error('Soton:argChk','~isequal(size(X_k),[timeslots,tx_antennas_per_user])');
            end
            
            % Accumulate the transmission energy of the STBC signal
            E_b(k) = E_b(k) + sum(sum(abs(X_k.^2)));
            
            % Put each user's STBC signal into the X matrix
            start = (k-1)*P+1;
            stop = k*P;
            X(:,start:stop) = X_k;
        end
        
        % Generate Rayleigh-distributed complex channel coefficients
        H = sqrt(1/2)*(randn(P*K,Q)+1i*randn(P*K,Q));
        
        % Generate AWGN
        N = sqrt(N_0/2)*(randn(T,Q)+1i*randn(T,Q));
        
        % Obtain the received signal
        Y = X*H+N;
        
        % Call receiver.m to demodulate the received signal
        B_hat = receiver(Y, H);
        
        % Check that receiver.m produces an output having the correct dimensions
        if ~isequal(size(B_hat),size(B))
            error('Soton:argChk','~isequal(size(decoded_bits),size(bits))');
        end
        
        % Measure the number of bit errors in each user's recovered bit vector
        for k = 1:K
            bit_error_counters(inv_N_0_dB_index,k) = bit_error_counters(inv_N_0_dB_index,k) + sum(B(k,:) ~= B_hat(k,:));
        end
        
    end
    
    % Plot the BER vs 1/N0 results
    if E_b_too_high == 1
        delete(ann1);
        E_b_too_high = 0;
    end
    for k = 1:K
        set(plots(k),'XData',inv_N_0_dB);
        set(plots(k),'YData',bit_error_counters(:,k)./bit_counters);
        key{k} = ['User ',num2str(k),'   E_b^{(',num2str(k),')}  = ', num2str(E_b(k)/sum(bit_counters))];
        set(leg1, 'String', key);
        
        if E_b(k)/sum(bit_counters) > 1.05
            E_b_too_high = 1;
        end
    end
    if E_b_too_high == 1
        ann1 = annotation('textbox', [0.2, 0.4, 0.6, 0.2],'String','The average E_b of one or more of your users is too high!','FontSize',24, 'Color','red');
    end
    drawnow
    
    % Stop the simulation early if any of the BER curves disappear off the bottom of the BER plot
    if min(bit_error_counters(inv_N_0_dB_index,:))/bit_counters(inv_N_0_dB_index) < min_BER_to_plot
        break;
    end
    
end
