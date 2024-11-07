clear all;close all;clc; 

%% MRC,EGC,SC 

nSym = 10e5; %Number of symbols to transmit 

EbN0dBs = -15:2:20; %Eb/N0 range in dB for simulation 

M = 2; %M-ary PSK modulation 

MOD_TYPE = 'PSK'; %Modulation type 

k=log2(M);EsN0dBs = 10*log10(k)+EbN0dBs; %EsN0dB calculation 

N = [2,4]; %number of diversity branches 

figure; 

for nRx = N %simulate for each # of received branchs 

    ser_MRC = zeros(1,numel(EsN0dBs));%simulated symbol error rates 

    ser_EGC = zeros(1,numel(EsN0dBs)); 

    ser_SC = zeros(1,numel(EsN0dBs));q=1; 

    %---------- Transmitter ----------- 

    d = ceil(M*rand(1,nSym));%uniform random symbols 1:M 

    s = modulate(MOD_TYPE,M,d);%Modulation 

    s_diversity = kron(ones(nRx,1),s);%nRx signal branches 

    for EsN0dB = EsN0dBs 

        h = sqrt(1/2)*(randn(nRx,nSym)+1j*randn(nRx,nSym));%Rayleigh flat-fading 

        signal = h.*s_diversity; %effect of channel on the modulated signal 

 

        gamma = 10.^(EsN0dB/10);%for AWGN noise addition 

        P = sum(abs(signal).^2,2)./nSym; %calculate power in each branch of signal 

        N0 = P/gamma; %required noise spectral density for each branch 

        %Scaled each row of noise with the calculated noise spectral density 

        noise = (randn(size(signal))+1j*randn(size(signal))).*sqrt(N0/2); 

 

        r = signal+noise;%received signal branches 

 

        %MRC processing assuming perfect channel estimates 

        s_MRC = sum(conj(h).*r,1)./sum(abs(h).^2,1); %detector decisions 

        d_MRC = demodulate(MOD_TYPE,M,s_MRC); %demodulation decisions 

 

        %EGC processing assuming perfect channel estimates 

        h_phases = exp(-1j*angle(h));%estimated channel phases 

        s_EGC = sum(h_phases.*r,1)./sum(abs(h),1); %detector decisions 

        d_EGC = demodulate(MOD_TYPE,M,s_EGC); %demodulation decisions 

 

        %SC processing assuming perfect channel estimates 

        [h_max,idx] = max(abs(h),[],1); %max |h| values along all branches 

        h_best = h(sub2ind(size(h),idx,1:size(h,2)));%best path's h estimate 

        y = r(sub2ind(size(r), idx, 1:size(r,2)));%selected output 

 

        s_SC = y.*conj(h_best)./abs(h_best).^2; %detector decisions 

        d_SC = demodulate(MOD_TYPE,M,s_SC); %demodulation decisions 

 

        ser_MRC(q) = sum(d_MRC ~= d)/nSym;%Error rates computation 

        ser_EGC(q) = sum(d_EGC ~= d)/nSym; 

        ser_SC(q) = sum(d_SC ~= d)/nSym;q=q+1; 

    end 

    semilogy(EbN0dBs,ser_MRC,'b-','lineWidth',1.5);hold on; 

    semilogy(EbN0dBs,ser_EGC,'k--','lineWidth',1.5); 

    semilogy(EbN0dBs,ser_SC,'r','lineWidth',1.5); 

end 

xlim([-20,36]);ylim([0.00001,1.1]); 

xlabel('Eb/N0(dB)');ylabel('Symbol Error Rate (P_s)') 

title('Receive diversity schemes in Rayleigh flat-fading') 

grid on 

 

%% Alamouti[Tx-2,Rx-1] 

 

nSym = 10e5; %Number of symbols to transmit 

EbN0dBs = -15:2:20; %Eb/N0 range in dB for simulation 

M = 2; %M-ary PSK modulation 

MOD_TYPE = 'PSK'; %Modulation type 

k=log2(M);EsN0dBs = 10*log10(k)+EbN0dBs; %EsN0dB calculation 

ser_sim = zeros(1,numel(EsN0dBs));q=1; %simulated symbol error rates 

for EsN0dB = EsN0dBs %simulate for each # of received branches 

    %---------- Transmitter ----------- 

    d = ceil(M*rand(nSym,1));%uniform random symbols 1:M 

    s = modulate(MOD_TYPE,M,d);%Modulation 

    ss = kron(reshape(s,2,[]),ones(1,2));%shape as 2xnSym vector 

     

    h = sqrt(1/2)*(randn(2,nSym/2)+1j*randn(2,nSym/2));%channel coeffs 

    H = kron(h,ones(1,2)); %shape as 2xnSym vector 

    H(:,2:2:end) = conj(flipud(H(:,2:2:end))); 

    H(2,2:2:end) = -1*H(2,2:2:end);%Alamouti coded channel coeffs 

     

    signal = sum(H.*ss,1);%effect of channel on the modulated signal 

     

    gamma = 10.^(EsN0dB/10);%for AWGN noise addition 

    P = sum(abs(signal).^2,2)./nSym; %calculate power in each branch of signal 

    N0 = P/gamma; %required noise spectral density for each branch 

    %Scale each row of noise with the calculated noise spectral density 

    noise = (randn(size(signal))+1j*randn(size(signal))).*sqrt(N0/2); 

     

    r = signal+noise; %received signal 

     

    %Receiver processing 

    rVec = kron(reshape(r,2,[]),ones(1,2)); %2xnSym format 

     

    Hest = H; %perfect channel estimation 

    Hest(1,1:2:end) = conj(Hest(1,1:2:end)); 

    Hest(2,2:2:end) = conj(Hest(2,2:2:end));%Hermitian transposed Hest matrix 

     

    y = sum(Hest.*rVec,1); %combiner output 

    sCap = y./sum(conj(Hest).*Hest,1);%decision vector for demod 

    dCap = demodulate(MOD_TYPE,M,sCap).';%demodulation 

    ser_sim(q) = sum(dCap ~= d)/nSym;q=q+1;%error rate calculation 

end 

semilogy(EbN0dBs,ser_sim,'m--','lineWidth',1.5);hold on;%plot simulated error rates 

title('2x1 Transmit diversity - Alamouti coding');grid on; 

xlabel('EbN0 (dB)');ylabel('Symbol error rate (Ps)'); 

 

%% Alamouti[Tx-2,Rx-2] 

 

N = 10e5; % number of bits or symbols 

Eb_N0_dB = [-15:2:20]; % multiple Eb/N0 values 

nRx = 2; 

for ii = 1:length(Eb_N0_dB) 

 

    % Transmitter 

    ip = rand(1,N)>0.5; % generating 0,1 with equal probability 

    s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 0 

 

    % Alamouti STBC 

    sCode = 1/sqrt(2)*kron(reshape(s,2,N/2),ones(1,2)) ; 

 

    % channel 

    h = 1/sqrt(2)*[randn(nRx,N) + j*randn(nRx,N)]; % Rayleigh channel 

    n = 1/sqrt(2)*[randn(nRx,N) + j*randn(nRx,N)]; % white gaussian noise, 0dB variance 

 

    y = zeros(nRx,N); 

    yMod = zeros(nRx*2,N); 

    hMod = zeros(nRx*2,N); 

    for kk = 1:nRx 

 

        hMod = kron(reshape(h(kk,:),2,N/2),ones(1,2)); % repeating the same channel for two symbols 

        hMod = kron(reshape(h(kk,:),2,N/2),ones(1,2)); 

        temp = hMod; 

        hMod(1,[2:2:end]) = conj(temp(2,[2:2:end])); 

        hMod(2,[2:2:end]) = -conj(temp(1,[2:2:end])); 

 

        % Channel and noise Noise addition 

        y(kk,:) = sum(hMod.*sCode,1) + 10^(-Eb_N0_dB(ii)/20)*n(kk,:); 

 

        % Receiver 

        yMod([2*kk-1:2*kk],:) = kron(reshape(y(kk,:),2,N/2),ones(1,2)); 

 

        % forming the equalization matrix 

        hEq([2*kk-1:2*kk],:) = hMod; 

        hEq(2*kk-1,[1:2:end]) = conj(hEq(2*kk-1,[1:2:end])); 

        hEq(2*kk,  [2:2:end]) = conj(hEq(2*kk,  [2:2:end])); 

 

    end 

 

    % equalization 

    hEqPower = sum(hEq.*conj(hEq),1); 

    yHat = sum(hEq.*yMod,1)./hEqPower; % [h1*y1 + h2y2*, h2*y1 -h1y2*, ... ] 

    yHat(2:2:end) = conj(yHat(2:2:end)); 

 

    % receiver - hard decision decoding 

    ipHat = real(yHat)>0; 

 

    % counting the errors 

    nErr(ii) = size(find([ip- ipHat]),2); 

 

end 

 

simBer = nErr/N; % simulated ber 

EbN0Lin = 10.^(Eb_N0_dB/10); 

theoryBer_nRx1 = 0.5.*(1-1*(1+1./EbN0Lin).^(-0.5)); 

 

p = 1/2 - 1/2*(1+1./EbN0Lin).^(-1/2); 

theoryBerMRC_nRx2 = p.^2.*(1+2*(1-p)); 

 

pAlamouti = 1/2 - 1/2*(1+2./EbN0Lin).^(-1/2); 

theoryBerAlamouti_nTx2_nRx1 = pAlamouti.^2.*(1+2*(1-pAlamouti)); 

 

 

% semilogy(Eb_N0_dB,theoryBer_nRx1,'bp-','LineWidth',2); 

hold on 

% semilogy(Eb_N0_dB,theoryBerMRC_nRx2,'kd-','LineWidth',2); 

% semilogy(Eb_N0_dB,theoryBerAlamouti_nTx2_nRx1,'c+-','LineWidth',2); 

semilogy(Eb_N0_dB,simBer,'g','LineWidth',2); 

semilogy(Eb_N0_dB,simBer,'ko','LineWidth',1.5); 

% axis([0 25 10^-5 0.5]) 

grid on 

% legend('theory (nTx=1,nRx=1)', 'theory (nTx=1,nRx=2, MRC)', 'theory (nTx=2, nRx=1, Alamouti)', 'sim (nTx=2, nRx=2, Alamouti)'); 

% xlabel('Eb/No, dB'); 

% ylabel('Bit Error Rate'); 

%% 

title('Performance Analysis using BPSK modulation in Rayleigh Flat Fading Channel'); 

axis tight 

xlabel('SNR(dB)'); 

ylabel('Symbol Error Rate(SER)'); 

%% 

legend('MRC 1x2','EGC 1x2','SC 1x2','MRC 1x4','EGC 1x4','SC 1x4','Alamouti 2x1','Alamouti 2x2') 

legend('Location','southwest') 
