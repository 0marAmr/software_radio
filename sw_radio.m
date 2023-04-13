clc;
clear;

%% initialization
N = 100;         % number of bits in a realization, an extra 1 added for delay
M = 500;              % number of realizations in the ensamble
PW  = 0.07;           % pulse width is 70ms
SR  = 0.01;           % DAC sampling rate
N_S = PW/SR;        % number of samples 
A=4;                % amplitude

%%ask the user about the desired line coding he wants to use
lineCode_var=input('please write polar or uni-polar : \n ','s');
polar_var = strcmp(lineCode_var,'polar');
uni_polar_var = strcmp(lineCode_var,'uni-polar');

while (~polar_var && ~uni_polar_var) %we will not accept any inputs except polar or uni -polar
lineCode_var=input('invalid word ,please write polar or uni-polar : \n','s');
polar_var = strcmp(lineCode_var,'polar');
uni_polar_var = strcmp(lineCode_var,'uni-polar');
end

%now if the user choose polar, we will ask him either polar-RZ or
%nonPolar-NRZ

if polar_var
TypePolar_var=input('please write RZ or NRZ: \n','s');
polar_RZ=strcmp(TypePolar_var,'RZ');
polar_NRZ=strcmp(TypePolar_var,'NRZ');
  
   while(~polar_RZ && ~polar_NRZ)
     TypePolar_var=input('invalid word, please write RZ or NRZ: \n','s');
     polar_RZ=strcmp(TypePolar_var,'RZ');
     polar_NRZ=strcmp(TypePolar_var,'NRZ');
   end
end

%%creation of ensemble
Tx_CompleteEnsemble=zeros(500,707);   % ensemble is a matrix of 500*700
 %707--> (700 SAMPLE +7 samples delay in each realization )   %500--->no.of realizations in ensemble
 for i=1:500  % i---> realization number iterator
    data=randi([0 1],1,101);
    if polar_var   %in case user choose polar line coding
        Tx= (2*data-1)*A;

        if polar_RZ   %in case user choose polar return to zero line coding
          Tx1=repmat(Tx,4,1);
          Tx1=[Tx1;zeros(3,101)];

        elseif polar_NRZ   %in case user choose polar non-return to zero line coding
           Tx1=repmat(Tx,7,1);
        end 

    elseif uni_polar_var
         Tx=data*A;
         Tx1=repmat(Tx,7,1);
    end 
 %now Tx1 is a matrix of 7*100
 Tx2=reshape(Tx1,1,size(Tx1,1)*size(Tx1,2)); %Tx2 is a materix of 1*700
 Tx_CompleteEnsemble(i,:)=Tx2;
 end
  
%% generation of time shift
time_shift=randi([0 6],500,1);
Tx_CompleteEnsemble_shifted=zeros(500,700);
 for i=1:500 % i---> realization number iterator
   if(time_shift(i)~=0) %ask if there is a time delay or not
     %in case the presence of a delay in a realization
     Tx_CompleteEnsemble_shifted(i,:)=Tx_CompleteEnsemble(i,[time_shift(i)+1:700+time_shift(i)]);
   else
       Tx_CompleteEnsemble_shifted(i,:)= Tx_CompleteEnsemble(i,[1:700]);
   end
 end
 
%%the realizations plotting

%1st realization
figure(1)
subplot(3,1,1,'LineWidth',2)
stairs(Tx_CompleteEnsemble_shifted(1,:),'B','LineWidth',2);
axis([1 701 -5 5])
title('1st realization');
xlabel('time');
ylabel('data');

%2nd realization
figure(1)
subplot(3,1,2,'LineWidth',5)
stairs(Tx_CompleteEnsemble_shifted(2,:),'R','LineWidth',2);
axis([1 701 -5 5])
title('2nd realization');
xlabel('time');
ylabel('data');

%3rd realization
figure(1)
subplot(3,1,3,'LineWidth',5)
stairs(Tx_CompleteEnsemble_shifted(3,:),'y','LineWidth',2);
axis([1 701 -5 5])
title('3rd realization');
xlabel('time');
ylabel('data');

%%statistcal mean
statistical_mean=zeros(1,700);
for k=1:700
  summation=0;
    for i=1:500
       summation= Tx_CompleteEnsemble_shifted(i,k)+summation;
    end
    statistical_mean(k)=summation/500;
end
fprintf('the mean of statistical mean = %f \n ',mean(statistical_mean));

%plotting
time=1:700; 
figure(2)
subplot(2,1,1,'LineWidth',3) 
plot(time,statistical_mean(1:700)); 
axis([1 700 -5 5]); 
title('Statistical Mean ');
xlabel('time'); 
ylabel('Statistical Mean');

%%time mean
time_mean=zeros(500,1);
for k =1:500
    summation=0;
    for i=1:700
        summation=Tx_CompleteEnsemble_shifted(k,i)+summation;
    end
    time_mean(k)=summation/700;
end
fprintf('the mean of time mean = %f \n ',mean(time_mean));

%plotting
time=1:500; 
figure(2)
subplot(2,1,2,'LineWidth',3)
plot(time,time_mean(1:500));
axis([1 500 -5 5]); 
title('Time Mean '); 
xlabel('time'); 
ylabel('Time Mean');
%% statistical autocorrelation
stat_autocorrelation = zeros(1, 700);

%Implementation of statstical autocorrelation function
for column = 1:700
    stat_autocorrelation(column) =(1/500) * sum(Tx_CompleteEnsemble_shifted(:, 1) .* Tx_CompleteEnsemble_shifted(:, column));
end

R_tau = [fliplr(stat_autocorrelation(2:end)) stat_autocorrelation]; %Flipping R_tau to the -ve quad

%Plotting statistical autocorrelation
figure(3)
x_axis = [1:1:1399];
subplot(2,1,1,'LineWidth',5)
plot(x_axis - 700,R_tau(x_axis),'r','LineWidth',2)
axis([-350 350 -3 18])
title('Statistical autocorrelation');
xlabel('tau');
ylabel('R(tau)');

%% time autocorrelation
autocorr_time = zeros(1, 700);
time_autocorrelation = zeros(1, 700);

for tau = 0:699
    accumulate = 0;
    
    %Summation over a constant tau across one realization    
    for count = 1:(700-tau)
        accumulate = accumulate + (Tx_CompleteEnsemble_shifted(1, count) * Tx_CompleteEnsemble_shifted(1, (count + tau))); 
    end
    %Divide by their number
    time_autocorrelation((tau+1)) = (1/(700-tau)) * accumulate;
end
R_tu_time = [fliplr(time_autocorrelation(2:end)) time_autocorrelation];

%Plotting time autocorrelation
figure(3)
subplot(2,1,2,'LineWidth',5)
plot(x_axis - 700,R_tu_time(x_axis),'b','LineWidth',2)
axis([-350 350 -3 18])
title('time autocorrelation');
xlabel('tau');
ylabel('R(tau)');

%%Bandwidth
L = 1400;             % Length of signal
n = 2^nextpow2(L);
T = 0.01;             % Sampling period
Fs = 1/T;             % Sampling frequency                      

t = (0:L-1)*T;        % Time vector

%Compute the Fourier transform of the signal.
PSD = fft(R_tau,n);
PSD = fftshift(PSD);

% f = Fs*(0:(L/2))/L;
f = Fs * (-n/2:n/2-1)/n;

%Compute the two-sided spectrum P2
P2 = abs(PSD/n).^2;
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);

Smoothed_PSD = smooth(P2); %Smooth the curve from the noise, to be plotted
figure (4)
subplot(1,1,1,'LineWidth',5)
plot(f,Smoothed_PSD) 
title("PSD")
xlabel("f (Hz)")
ylabel("S(f)")


