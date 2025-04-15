% STBC_3x4_simulation.m
clear; clf
Lfr = 130; % Frame length (fixed for now)
N_packet = input('Enter the number of packets (e.g., 4000): '); % Number of packets
NT = input('Enter the number of transmit antennas (e.g., 3): '); % Number of transmit antennas
NR = input('Enter the number of receive antennas (e.g., 4): '); % Number of receive antennas
b = input('Enter modulation order (bits per symbol, e.g., 2 for QPSK): '); % Modulation order
M = 2^b; % Modulation size (QAM or PSK)
SNRdBs = input('Enter the range of SNR in dB (e.g., [0:2:30]): '); % SNR range in dB
sq_NT = sqrt(NT); 
sq2 = sqrt(2);

% Initialize BER array for plotting
BER = zeros(1, length(SNRdBs));

for i_SNR = 1:length(SNRdBs)
    SNRdB = SNRdBs(i_SNR);
    sigma = sqrt(0.5 / (10^(SNRdB / 10))); % Noise power calculation
    
    % Initialize error tracking
    noeb_p = zeros(1, N_packet);

    % Simulation over N_packet
    for i_packet = 1:N_packet
        msg_symbol = randi([0 M-1], Lfr * b, 1); % Generate message symbols
        tx_bits = msg_symbol.'; 
        temp = []; temp1 = [];

        % Modulate message bits
        for i = 1:4
            [temp1, sym_tab, P] = modulator(tx_bits(i,:), b); 
            temp = [temp; temp1];
        end

        % Transmitter processing
        X = temp.'; % Transmit signal

        % Block coding for G3 STBC (Space Time Block Code)
        X1 = X(:, 1:3); X5 = conj(X1); 
        X2 = [-X(:, 2) X(:, 1) -X(:, 4)]; X6 = conj(X2);
        X3 = [-X(:, 3) X(:, 4) X(:, 1)]; X7 = conj(X3);
        X4 = [-X(:, 4) -X(:, 3) X(:, 2)]; X8 = conj(X4);
        
        % Generate random channel matrices for each antenna
        for n = 1:NT
            Hr(n, :, :) = (randn(Lfr, NT) + 1j * randn(Lfr, NT)) / sq2;
        end

        % Channel processing for each antenna
        for n = 1:NT
            H = reshape(Hr(n, :, :), Lfr, NT); 
            Hc = conj(H);
            Habs(:, n) = sum(abs(H).^2, 2);
            
            % Received signals for each STBC encoded symbol
            R1n = sum(H .* X1, 2) / sq_NT + sigma * (randn(Lfr, 1) + 1j * randn(Lfr, 1));
            R2n = sum(H .* X2, 2) / sq_NT + sigma * (randn(Lfr, 1) + 1j * randn(Lfr, 1));
            R3n = sum(H .* X3, 2) / sq_NT + sigma * (randn(Lfr, 1) + 1j * randn(Lfr, 1));
            R4n = sum(H .* X4, 2) / sq_NT + sigma * (randn(Lfr, 1) + 1j * randn(Lfr, 1));
            R5n = sum(H .* X5, 2) / sq_NT + sigma * (randn(Lfr, 1) + 1j * randn(Lfr, 1));
            R6n = sum(H .* X6, 2) / sq_NT + sigma * (randn(Lfr, 1) + 1j * randn(Lfr, 1));
            R7n = sum(H .* X7, 2) / sq_NT + sigma * (randn(Lfr, 1) + 1j * randn(Lfr, 1));
            R8n = sum(H .* X8, 2) / sq_NT + sigma * (randn(Lfr, 1) + 1j * randn(Lfr, 1));

            % STBC processing and combining received signals
            Z1_1 = R1n .* Hc(:, 1) + R2n .* Hc(:, 2) + R3n .* Hc(:, 3);
            Z1_2 = conj(R5n) .* H(:, 1) + conj(R6n) .* H(:, 2) + conj(R7n) .* H(:, 3);
            Z(:, n, 1) = Z1_1 + Z1_2;
            Z2_1 = R1n .* Hc(:, 2) - R2n .* Hc(:, 1) + R4n .* Hc(:, 3);
            Z2_2 = conj(R5n) .* H(:, 2) - conj(R6n) .* H(:, 1) + conj(R8n) .* H(:, 3);
            Z(:, n, 2) = Z2_1 + Z2_2;
            Z3_1 = R1n .* Hc(:, 3) - R3n .* Hc(:, 1) - R4n .* Hc(:, 2);
            Z3_2 = conj(R5n) .* H(:, 3) - conj(R7n) .* H(:, 1) - conj(R8n) .* H(:, 2);
            Z(:, n, 3) = Z3_1 + Z3_2;
            Z4_1 = -R2n .* Hc(:, 3) + R3n .* Hc(:, 2) - R4n .* Hc(:, 1);
            Z4_2 = -conj(R6n) .* H(:, 3) + conj(R7n) .* H(:, 2) - conj(R8n) .* H(:, 1);
            Z(:, n, 4) = Z4_1 + Z4_2;
        end

        % Symbol decision and calculation of bit errors
        for m = 1:P
            tmp = (-1 + sum(Habs, 2)) * abs(sym_tab(m))^2;
            for i = 1:4
                d(:, m, i) = abs(sum(Z(:, :, i), 2) - sym_tab(m)).^2 + tmp;
            end
        end

        % Decode the transmitted symbols
        Xd = [];
        for n = 1:4
            [yn, in] = min(d(:, :, n), [], 2);
            Xd = [Xd sym_tab(in).'];
        end

        % Calculate bit errors
        temp1 = X > 0; 
        temp2 = Xd > 0;
        noeb_p(i_packet) = sum(sum(temp1 ~= temp2));
    end
    
    % Compute BER for the current SNR value
    BER(i_SNR) = sum(noeb_p) / (N_packet * Lfr * b);
end

% Plot the results
semilogy(SNRdBs, BER, '-o', 'LineWidth', 2);
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('STBC 3x4 Simulation: BER vs SNR');
grid on;
axis([min(SNRdBs) max(SNRdBs) 1e-6 1e0]);
legend('STBC 3x4');

