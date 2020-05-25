%{
Communication Theory and Systems
CIE 337 - Spring 20 20
By: Ahmed Ashraf Hussein
Email: s-ahmedashraf@zewailcity.edu.eg
Final Assessment
PART1_MATCHED_FILTER
%}
%{
Prompt:
It's required to design an optimum ISI-free system
Channel impulse response c(t) = delta(t)
Noise w(t) is AWGN ~ G(0, No/2)
%}
clear;
clc;
T_b = 1; %assuming bit dirutaion is 1 second
time_step = 0.001;
t = 0:time_step:T_b; %time
%definging the channel impulse response c(t)
c = dirac(t);
idx = c == Inf; % find Inf
c(idx) = 1;
%defining the pulse shape g(t)
g = zeros(1,length(t));
g(t >= 0 & t <= 0.5*T_b) = 0.5;
g(t >= 0.5*T_b & t <= T_b) = -0.5;
%defining the filter impulse response h(t)
%it's the y-axis mirrored and shifted version of g(t)
h = circshift(fliplr(g),T_b*(1/time_step));
A = [1,1,-1,1,-1];
stream = [];
for i = 1:length(A)
    stream = [stream A(i)*g];
end
rx = conv(h, stream);
k = (1/(2*time_step));
time_samples = [];
samples_rx = [];
A_rx = [];
for i = 1:length(A)
    time_samples = [time_samples i];
    samples_rx = [samples_rx rx(i*(1/time_step)+1)];
    if rx(i*(1/time_step)+1) > 0
        A_rx = [A_rx -1]; 
    elseif rx(i*(1/time_step)+1) < 0
        A_rx = [A_rx 1];        
    end
end
h1 = zeros(1,(length(t)));
h1(t >= 0 & t <= 0.5*T_b) = 1;

rx_h1 = conv(h1, stream);
samples_rx_h1 = [];
A_rx_h1 = [];
for i = 1:length(A)
    samples_rx_h1 = [samples_rx_h1 rx_h1(i*(1/time_step)+1)];
    if rx_h1(i*(1/time_step)+1) > 0
        A_rx_h1 = [A_rx_h1 1]; 
    elseif rx_h1(i*(1/time_step)+1) < 0
        A_rx_h1 = [A_rx_h1 -1];        
    end
end
%% Run this section in the presence of AWGN only
%defining the noise w(t)
N_o = 200;
w = sqrt(N_o/2).*randn(1,length(stream));
Noisy_stream = stream + w;

rx_noisy = conv(h, Noisy_stream);
k = (1/(2*time_step));
samples_rx_noisy = [];
A_rx_noisy = [];
for i = 1:length(A)
    samples_rx_noisy = [samples_rx_noisy rx_noisy(i*(1/time_step)+1)];
    if rx_noisy(i*(1/time_step)+1) > 0
        A_rx_noisy = [A_rx_noisy -1];
    elseif rx_noisy(i*(1/time_step)+1) < 0
        A_rx_noisy = [A_rx_noisy 1];    
    end
end
h1 = zeros(1,(length(t)));
h1(t >= 0 & t <= 0.5*T_b) = 1;

rx_h1_noisy = conv(h1, Noisy_stream);
samples_rx_h1_noisy = [];
A_rx_h1_noisy = [];
for i = 1:length(A)
    samples_rx_h1_noisy = [samples_rx_h1_noisy rx_h1_noisy(i*(1/time_step)+1)];
    if rx_h1_noisy(i*(1/time_step)+1) > 0
        A_rx_h1_noisy = [A_rx_h1_noisy 1];
    elseif rx_h1_noisy(i*(1/time_step)+1) < 0
        A_rx_h1_noisy = [A_rx_h1_noisy -1];        
    end
end
%%
plot(t, g)
axis([0 1 -0.55 0.55])
xticks([0 0.5*T_b T_b])
hline = refline(0, 0);
hline.Color = 'k';
xticklabels({'0', 'T_b/2', 'T_b'})
xlabel('time (Seconds)')
ylabel('Amplitude')
title('Manchester signaling pulse-shaping filter')
%%
plot(t, h)
axis([0 1 -0.55 0.55])
xticks([0 0.5*T_b T_b])
hline = refline(0, 0);
hline.Color = 'k';
xticklabels({'0', 'T_b/2', 'T_b'})
xlabel('time (Seconds)')
ylabel('Amplitude (Multiplied by an arbitrary constant K)')
title('Optimum Receiver Impulse Response')
%%
t_output = 0:time_step:length(rx)*time_step-time_step;
plot(t_output, rx/k)
hold on;
stem(time_samples, samples_rx/k)
axis([0 5.5 -0.55 0.55])
xticks([0 0.5*T_b T_b 1.5*T_b 2*T_b 2.5*T_b 3*T_b 3.5*T_b 4*T_b 4.5*T_b 5*T_b])
hline = refline(0, 0);
hline.Color = 'k';
xticklabels({'0', '0.5T_b', 'T_b', '1.5T_b', '2T_b', '2.5T_b', '3T_b', '3.5T_b', '4T_b', '4.5T_b', '5T_b'})
xlabel('time (Seconds)')
ylabel('Normalized Amplitude (*K)')
title('Output of the Optimum Receiver')
legend('Output','Sampled Output')
%%
plot(t, h1)
axis([0 1 0 1.05])
xticks([0 0.5*T_b T_b])
hline = refline(0, 0);
hline.Color = 'k';
xticklabels({'0', 'T_b/2', 'T_b'})
xlabel('time (Seconds)')
ylabel('Amplitude')
title('Receiver Impulse Response h1(t)')
%%
t_output = 0:time_step:length(rx_h1)*time_step-time_step;
plot(t_output, rx_h1/k)
hold on;
stem(time_samples, samples_rx_h1/k)
axis([0 5.5 -0.55 0.55])
xticks([0 0.5*T_b T_b 1.5*T_b 2*T_b 2.5*T_b 3*T_b 3.5*T_b 4*T_b 4.5*T_b 5*T_b])
hline = refline(0, 0);
hline.Color = 'k';
xticklabels({'0', '0.5T_b', 'T_b', '1.5T_b', '2T_b', '2.5T_b', '3T_b', '3.5T_b', '4T_b', '4.5T_b', '5T_b'})
xlabel('time (Seconds)')
ylabel('Normalized Amplitude (*K)')
title('Output of the Optimum Receiver')
legend('Output','Sampled Output')
%%
t_output = 0:time_step:length(rx_noisy)*time_step-time_step;
subplot(2,1,1)
plot(t_output, rx_noisy/k)
hold on;
stem(time_samples, samples_rx_noisy/k)
axis([0 5.5 -0.55 0.55])
xticks([0 0.5*T_b T_b 1.5*T_b 2*T_b 2.5*T_b 3*T_b 3.5*T_b 4*T_b 4.5*T_b 5*T_b])
hline = refline(0, 0);
hline.Color = 'k';
xticklabels({'0', '0.5T_b', 'T_b', '1.5T_b', '2T_b', '2.5T_b', '3T_b', '3.5T_b', '4T_b', '4.5T_b', '5T_b'})
xlabel('time (Seconds)')
ylabel('Normalized Amplitude (*K)')
title('Output of the Optimum Receiver')
legend('Output','Sampled Output')

t_output = 0:time_step:length(rx_h1_noisy)*time_step-time_step;
subplot(2,1,2)
plot(t_output, rx_h1_noisy/k)
hold on;
stem(time_samples, samples_rx_h1_noisy/k)
axis([0 5.5 -0.65 0.65])
xticks([0 0.5*T_b T_b 1.5*T_b 2*T_b 2.5*T_b 3*T_b 3.5*T_b 4*T_b 4.5*T_b 5*T_b])
hline = refline(0, 0);
hline.Color = 'k';
xticklabels({'0', '0.5T_b', 'T_b', '1.5T_b', '2T_b', '2.5T_b', '3T_b', '3.5T_b', '4T_b', '4.5T_b', '5T_b'})
xlabel('time (Seconds)')
ylabel('Normalized Amplitude (*K)')
title('Output of the Optimum Receiver')
legend('Output','Sampled Output')