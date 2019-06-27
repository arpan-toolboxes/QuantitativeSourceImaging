function S=SNR_Gol(Stim, Base)
% SNR- Signal/Noise ratio
%  Usage
%    value=SNR(sig1,sig2)
%  Inputs   
%    sig1	Original reference signal
%    sig2	Restored or noisy signal
%  Outputs
%    value	Signal/Noise ratio.
%
ab_B_sigma=((max(Stim)-min(Stim))./std(Base)).^2;
% size(ab_B_sigma)
S=10*log10(mean(ab_B_sigma));
% size(S)
