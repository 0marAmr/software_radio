clc;
clear;

%% initialization
N = 11;         % number of bits in a realization
M = 5;          % number of realizations in the ensamble

PW  = 70;       % pulse width is 70ms
SR  = 10;       % DAC sampling rate
N_S = PW/SR;    % number of samples 
A=4;            % amplitude
ensamble = [];

%% polar NRZ
data = randi([0,1],M,N);            % Generate random Data
for i = 1 : M
    Tx = data(i,:);         % POLAR NRZ: maping for 0 to be –A, 1 to be A
    Tx2=repmat(Tx,N_S,1);
    Tx2 = polar_nrz(Tx2,A);
    Tx_out=reshape(Tx2,1, size(Tx2,1)* size(Tx2,2));
    ensamble = vertcat(ensamble,Tx_out);
end
realization1 = array
