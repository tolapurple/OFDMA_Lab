%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  TRANSMISSION CODES  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OFDM Transmitter,

clc;
close all
clear all



% Training Sequence Generation

% Generate m-squences

tr1 = (mseq(2,6)).';

% Baseband signal generation (OFDM)

frame               = [];
frame_symbols = [];
num_OFDM_symbols    = 25+21;     %%Added 21 because of the pilot situation
N                   = 64;         % Number of Subcarriers/fftSize
N_data              = 48;          % Active subcarriers
CP_size             = 8;
subSpacing          = 15e3;        % Subcarrier Spacing
OOBE_reduction  = 0;  % Enables OOBE reduction
winLen = 0;  % Windowing length

message = ['We hold these truths 2 be self-evident, that all men are created equal, that they are endowed by their Creator with certain unalienable Rights, that among these are Life, Liberty and the pursuit of Happiness'];

decimal_vector_std_01 = dec2bin(double(char(message)));
binary_string_vector_std_01 = reshape(decimal_vector_std_01,1,size(decimal_vector_std_01,1)*size(decimal_vector_std_01,2));
binary_vector_std_01 = [];

%%% converts to string
for i=1:length(binary_string_vector_std_01)
    binary_vector_std_01(i) = str2num(binary_string_vector_std_01(i));
end

% Place the pilot bits, % please note this is only specific to this data size
buff = reshape(binary_vector_std_01(2:end),2,[]);  %%% reshape data to insert pilots
buff = [repmat([1 0],1,length(buff)/2); buff];
symbols = 2*[binary_vector_std_01(1) buff(:).']-1;
symbols = [symbols 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...        
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];      %%zero padding
                                     

used_subc_indices       = [N-(N_data/2:-1:1)+1, 1:N_data/2];  % Make sure to note that this is the order you should receive the data.


for ii=1:num_OFDM_symbols%(25+21=46)   %%Added 21
    
    ssc                     = zeros(1,N);
    ssc(used_subc_indices)  = symbols((1:N_data) + (ii-1)*N_data);
    time                    = sqrt(N)*ifft(ssc);
    
    time_cp                 = [time(end-CP_size+1:end) time];
    frame                   = [frame time_cp];
    
end

tx_signal   = [tr1 tr1 frame tr1];      % This is the main variable
% containing the basband signal
% that will be transmitted.


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
                    %%%%%%%%%SECTION 2 %%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % Plotting the power spectrum
figure(1)
 clf
 [Pxx f] = pwelch(tx_signal,[],[],[],subSpacing,'centered');
 plot(10*log10(Pxx));
 title('power spectrum of the ofdm signal with 48 carriers')
 
%%% Plotting the ccdf
figure (2)
[f]= cumsum(hist(abs(tx_signal)));
myccdf = 1-f;
plot(myccdf)
title('CCDF of the ofdm signal with 48 carriers')
xlabel ('X')
ylabel ('CCDF')

% % % %%%% Applying Windowing as an OOBE technique

L= length(tx_signal);
wind = window(@tukeywin,L,1);
New_Signal = tx_signal .* wind'; %(Scalar Multiplication)
Freq = fft(New_Signal);
figure(3)
plot(abs(New_Signal));

% % % %%%% Applying Clipping as a PAPR Tecnique
avg=0.5;
clipped_signal=tx_signal;
for i=1:length(clipped_signal)
	if clipped_signal(i) > avg
		clipped_signal(i) = avg;
    end
    if clipped_signal(i) < -avg
		clipped_signal(i) = -avg;
    end
end
dd=(abs(clipped_signal));
hist(dd)

[g]= cumsum(hist(abs(clipped_signal)));
myccdf = 1-g;
figure(4)
plot(myccdf)
title('CCDF of the clipped ofdm signal')
xlabel ('X')
ylabel ('CCDF')


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
                    %%%%%%%%%RECEIVER DESIGN%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%fft(abs(tx_signal))
seqlength = length(tr1);
seqfrmseq = seqlength+length(frame)+seqlength;

% % % % Coarse Time synchronization
%R=(zeros);
for i = 1:(length(tx_signal)-seqfrmseq)  
a =abs(tx_signal(i+1:i+seqlength));  %estimated training sequence
b =abs(tx_signal(i+(seqlength+length(frame)+1):i+seqfrmseq));
R(i,1)= myAutocorr(a,b);
end 
[l,k]= max(R);


tx1_signal = tx_signal(k+1 : k+seqfrmseq);
frame2 = tx1_signal(k+1 : k+length(frame));

for i = 1 : (length(frame2))
%frame_ofdmSym1 = frame2(i : (N+CP_size)+i);
 frame_ofdmSym1= reshape(frame2,[(N+CP_size),(num_OFDM_symbols)]);
end
frame_ofdmSym1(1:8,:)=[]; %%% Removing the CP

fft_frame_ofdmSym1= fft(frame_ofdmSym1,64);  %%% performing fft

fft_frame_ofdmSym1(25:40,:)=[];

fft_frame_ofdmSym2 = fft_frame_ofdmSym1/8;
%%%%%%%%%%%%%%%%
% Splitting in order to aid flipping the data 41-64,1-24
tmefrme1 = fft_frame_ofdmSym2(25:48,:); %%%25-48
tmefrme2 = fft_frame_ofdmSym2(1:24,:);  %%%%1-24
Frme = [tmefrme1; tmefrme2];       %%%%% Joining them back           
Frme2= reshape(Frme,1,[]);             %%%% converting to a row vector
% % % % % % % % % % % % % % removing pilots% % % % % % % % % % % % 
m = Frme2(1);
Frme2v2 = Frme2(2:(length(Frme2)));
Frme2v3 = Frme2v2(1:length(Frme2v2)-35);
Frme2v4=reshape(Frme2v3,[3,724]);

% Remove pilots
Frme2v4(1,:)=[];  %%Remove the first row             
Frme2v5= reshape(Frme2v4,1,[]);   %% Joining them back
 Frme2v6 = [m Frme2v5];
% % % % % % % % % % % % % % Decoding % % % % % % % % % % % % % % 
Rx_Bits=zeros;
for i= 1:length(Frme2v6)
if (Frme2v6(i) >= 0)
    Rx_Bits(i,1) = 1; 
end
if (Frme2v6(i) < 0)
    Rx_Bits(i,1) = 0; 
end
end
%%%Conversion to text
 %%%Selecting the data from the entire frame
Rx_text = bin2text(Rx_Bits);


 % % Time signal
 figure(5)
 plot(10*log10(abs(Frme2v6)),'b-')
 title('Time signal of the received OFDM symbols')

 figure(6)
 [Pxx f] = pwelch(Frme2v6,[],[],[],subSpacing,'centered');
 plot(10*log10(Pxx));
 title('power spectrum of the received ofdm symbols')
 
  %Constellation diagram
 figure(7)
 scatterplot(Frme2v6)



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
                    %%%%%%%%%Section 4%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
load('OFDM_CableReference.mat', 'Y')
for i = 1:(length(Y)-(2*length(tr1)))  
    p = (Y(i+1:i+length(tr1)));
    q = (tr1)';
    Z(i,1)= myAutocorr(p,q);
end 
[v,w]= max(Z);

for i = 1:(length(Y)-(2*length(tr1)))   
 s = (Y(i+1:i+length(tr1)));
 u = (Y(i+length(tr1)+1:i+length(tr1)+length(tr1)));
RW(i,1)= myAutocorr(s,u);
end 
[l,k]= max(RW);

 figure(8)
 plot(abs(RW))

function c = myAutocorr(m,n)
xc = m.*n;
c = sum(xc);
end
