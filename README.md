# Subband averaging kurtogram with dual-tree complex wavelet packet
This repository provides the code for the paper Subband averaging kurtogram with dual-tree complex wavelet
packet transform for rotating machinery fault diagnosis
https://www.sciencedirect.com/science/article/pii/S0888327020301412

Matlab Implementation of SAK. This paper presents a method called subband 
averaging kurtogram (SAK), incorporating with dual-tree complex wavelet packet 
transform (DTCWPT), to improve performance of the fast kurtogram (FK) for rotating 
machinery fault diagnosis. The proposed method first segments a signal into 
M sub-signals by a sliding window, then computes the kurtosis of subbands 
obtained by DTCWPT of each sub-signal. Finally, average kurtosis of 
corresponding subbands are calculated to obtain the SAK, which indicates 
the optimal frequency band for the envelope analysis. The FK is easily 
misled by non-Gaussian noise (e.g., sporadic impulse interferences) 
whereas the SAK can overcome this problem. Moreover, the DTCWPT 
simultaneously subdivides bands at high and low frequencies, offers the 
desirable property of approximate shift-invariance and meanwhile remains 
less computationally expensive. When the original DTCWPT iterates filter 
banks on the high-pass channel, the obtained subbands of a signal are not 
arranged in monotone order of the center frequency. This problem can be 
resolved by exchanging the inverted filter banks based on their band-pass 
properties. The proposed method provides improved performance compared 
to FK, in particular, for extracting periodic transients from noisy signals containing
a variety of interferences.
