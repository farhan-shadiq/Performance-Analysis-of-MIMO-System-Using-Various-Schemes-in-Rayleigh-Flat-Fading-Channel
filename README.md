# Performance-Analysis-of-MIMO-System-Using-Various-Schemes-in-Rayleigh-Flat-Fading-Channel
## This project analyzes the performance of a MIMO System using MRC, EGC, SC and Alamouti Scheme with BPSK Modulation in Rayleigh Flat Fading Channel.

### Introduction

Transmitted Message Signal undergoes distortion, fading, shadowing etc. due to the nature of the wireless channel. The received signal SNR is fluctuating in nature due to channel fading. Fading causes performance degradation in wireless communication severely. Diversity techniques are employed to make a wireless communication system robust and reliable even over varying channel conditions. Diversity techniques combat fading and interference by presenting the receiver with multiple uncorrelated copies of the same transmitted signal. Combining techniques are employed at the receiver to exploit multipath propagation characteristics of a channel. MIMO (Multiple Input Multiple Output) systems provide spatial diversity which can be used to resolve the issue of fading. Various techniques such as Maximal Ratio Combining (MRC), Equal Gain Combining (EGC), Selection Combining (SC), Alamouti Scheme etc. are used in MIMO Systems. In this project, we simulated these various techniques in MATLAB and compared there performance under Rayleigh Flat Fading Channel with AWGN noise using BPSK Modulation.

Diversity techniques are employed to make a communication system robust and reliable even over varying channel conditions. Diversity techniques exploit the channel variations rather than mitigating it. Diversity techniques combat fading and interference by presenting the receiver with multiple uncorrelated copies of the same information bearing signal. Essentially, diversity techniques are aimed at creating uncorrelated random channels, thereby creating uncorrelated copies of the same signal at the receiver front end. Combining techniques are employed at the receiver to exploit multipath propagation characteristics of a channel. Broadly speaking, the term diversity is categorized as: Time diversity, Frequency Diversity, Spatial diversity (antenna diversity), Polarization diversity (antenna diversity) etc. In Spatial diversity Aimed at creating uncorrelated propagation paths for a signal, spatial diversity is effected by the usage of multiple antennas in the transmitter and/or the receiver. Employing multiple antennas at the transmitter is called transmit diversity and multiple antennas at the receiver is called receive diversity. Diversity combining techniques like selection combining (SC), feedback or scanning combining, maximum ratio combining (MRC) can be employed by the receiver to exploit the multipath effects. Spatial diversity techniques can also be used to increase the data rates (spatial multiplexing) rather than improving the reliability of the channel. MIMO and space-time block coding (STBC) can be mentioned as their example.

The MIMO systems can be modeled as follows:

![MIMO System Model](https://github.com/farhan-shadiq/Performance-Analysis-of-MIMO-System-Using-Various-Schemes-in-Rayleigh-Flat-Fading-Channel/blob/main/MIMO%20System%20Model.jpg)

The system configuration typically contains M antennas at the transmitter and N antennas at the receiver front end as illustrated in the figure above. This is a multiple input multiple output (MIMO) channel model. Each receiver antenna receives not only the direct signal intended for it, but also receives a fraction of signal from other propagation paths. Thus, the channel response is expressed as a transmission matrix H. The direct path formed between antenna 1 at the transmitter and the antenna 1 at the receiver is represented by the channel response h11. The channel response of the path formed between antenna 1 in the transmitter and antenna 2 in the receiver is expressed as h21 and so on. Thus, the channel matrix is of dimension N ×M. The received vector y is expressed in terms of the channel transmission matrix H, the input vector x and noise vector n as, y = Hx+n; where the various symbols are,

![Symbols for Channel Transmission Matrix](https://github.com/farhan-shadiq/Performance-Analysis-of-MIMO-System-Using-Various-Schemes-in-Rayleigh-Flat-Fading-Channel/blob/main/Symbols%20for%20Channel%20Transmission%20Matrix.jpg)

For asymmetrical antenna configuration (M 6 = N), the number of data streams (or the number of uncoupled equivalent channels) in the system is always less than or equal to the minimum of the number of transmitter and receiver antennas min(M,N).

In mobile wireless systems, the envelope of the received signal is composed of superposition of slow and fast fading components. Slow fading is caused by the coarse variations in the terrain between the mobile and base stations. Fast fading is caused by the interference between the multipath waves scattered from the surrounding scatterers. The effects of rapid and deep fading are combated using spatial diversity techniques. In spatial diversity techniques, same information is sent across independent uncorrelated channels to com bat fading. When multiple copies of the same data are sent across independently fading channels, the amount of fade suffered by each copy of the data will be different. This guarantees that at-least one of the copies will suffer less fading compared to rest of the copies. At the receiver, the signal from the independent propagation paths are combined and hence the chance of properly receiving the transmitted data increases. In effect, this improves the reliability of the entire system. MRC, EGC, SC, Alamouti Scheme can be used to provide spatial diversity in MIMO systems.

MRC combines all the signals in a co-phased and weighted manner so as to have the highest achievable SNR at the receiver at all times. The following figure shows the working procedure of MRC:

![Working Procedure of MRC](https://github.com/farhan-shadiq/Performance-Analysis-of-MIMO-System-Using-Various-Schemes-in-Rayleigh-Flat-Fading-Channel/blob/main/Working%20Procedure%20of%20MRC.jpg)

Although, MRC provides the best average SNR at the output, the computational complexity of its implementation is very high.

EGC mitigates this problem by combining all the signals in a co-phased manner with unity weights for all signal levels so as to have the highest achievable SNR at the receiver at all times. This reduces the computational complexity with the trade-off in slightly lower SNR at the receiver side. The following figure shows the working procedure of EGC:

![Working Procedure of EGC](https://github.com/farhan-shadiq/Performance-Analysis-of-MIMO-System-Using-Various-Schemes-in-Rayleigh-Flat-Fading-Channel/blob/main/Working%20Procedure%20of%20EGC.jpg)

In SC technique, the best signal among all the signals received from different branches is selected at the receiver. This is a very simple technique to implement. But the average received SNR is poorer in this case compared to MRC and EGC which becomes significant enough with higher diversity branches.

Alamouti Scheme is a coding technique that does not require any feedback from the receiver nor does it require any bandwidth expansion. The computational complexity is similar to MRC receive diversity technique with one transmit and two receive antennas. The following figure shows its working procedure:

![Alamouti Scheme](https://github.com/farhan-shadiq/Performance-Analysis-of-MIMO-System-Using-Various-Schemes-in-Rayleigh-Flat-Fading-Channel/blob/main/Working%20Procedure%20of%20Alamouti%20Scheme.jpg)

Rayleigh Flat Fading Channel model is used to describe the statistical time varying nature of wireless signal with flat fading when Line-of-Sight (LOS) does not exist. In this channel model, the received voltage envelope provided by the Rayleigh distribution is:

![Rayleigh Distribution](https://github.com/farhan-shadiq/Performance-Analysis-of-MIMO-System-Using-Various-Schemes-in-Rayleigh-Flat-Fading-Channel/blob/main/Rayleigh%20Distribution.jpg)

Rayleigh Flat Fading Channel with BPSK Modulation was used to simulate MIMO system in this project. Then the performance of MRC, EGC, SC with Tx-1, Rx=2 (L=2) & Tx-1, Rx=4 (L=4), Alamouti 2x1 (L=2) & Alamouti 2x2 (L=4) were analysed by plotting the BER vs SNR graph.

### System

We performed the simulation using MATLAB. Firstly, 106 bits were generated randomly which would be used as the transmit symbols. SNR range was taken between -15dB to +20dB with an increment of 2dB. The bits were modulated with BPSK and then sent through a Rayleigh Flat Fading Channel. Signal energy was assumed to be equal for each transmission bit. For each SNR, Rayleigh Fading gain was applied and AWGN noise was added to Tx signal to generate the Received Signal. In the receiver side, perfect channel estimation was assumed. We performed the simulation according to the figures shown previously. After bit detection, Bit Error Rate (BER) was calculated by comparing the received bits with the sent ones. Then all the BER’s were plotted against the SNR values. The MATLAB code for this project can be found in [`main.m`](https://github.com/farhan-shadiq/Performance-Analysis-of-MIMO-System-Using-Various-Schemes-in-Rayleigh-Flat-Fading-Channel/blob/main/main.m).

### Result

BER vs. SNR plots for different simulations is shown below: 

![BER vs. SNR plot](https://github.com/farhan-shadiq/Performance-Analysis-of-MIMO-System-Using-Various-Schemes-in-Rayleigh-Flat-Fading-Channel/blob/main/BER%20vs.%20SNR%20plot.jpg)

### Conclusion 

We can see from the BER vs. SNR plot that, the better BER performance for all the systems can be noted down in the following order: 

MRC-1x4 > EGC-1x4 > Alamouti-2x2 ≈ SC-1x4 > MRC-1x2 > EGC-1x2 > SC-1x2 

MRC, EGC-1x4 outperforms Alamouti-2x2 although they have the same diversity L=4 because, same energy for each bit was assumed. The Alamouti-2x2 uses 2 transmit antennas which eventually divides the energy for each signal copy by a factor of 2. As a result, the received SNR for the Alamouti is less compared to MRC and EGC. Same thing happens to MRC, EGC-1x2 and Alamouti-2x1 although they have same diversity of L=2. Alamouti-2x2 and SC-1x4 shows equivalent BER vs SNR performance initially, but as SNR increases, Alamouti outperforms SC. 
