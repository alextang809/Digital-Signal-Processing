diary Tutorial3

%Question 1a
% Define the time vector
t   = 0 : 0.001 : 0.999;

%Define the signal
y = 2 + 5*sin(4*pi*t) + 10*cos(8*pi*t) + 15*sin(10*pi*t + pi/4) + 10*cos(14*pi*t);

%Plot the signal
plot(t, y);
xlabel('Time(sec)');
ylabel('Amplitude');
title('signal y(t)');

%Question 1b
% Define the time vector
%%t = 0:0.001:0.999;

% Define the frequency vector
%%f = 0:9;

% Initialize the sinusoidal and cosine basis function matrices
%%B_sin = zeros(length(t), length(f));
%%B_cos = zeros(length(t), length(f));

% Project
%%p_i = y * B_sin;
%%q_i = y * B_cos;

% Compute square root
%%m_i = zeros(size(f));
%%for i = 1:length(f)
%%    m_i(i) = sqrt(p_i(i)^2 + q_i(i)^2);
%%end

% Define the time vector
t = 0:0.001:0.999;

% Define the frequency vector
f = 0:9;

% Define the signal
y = 2 + 5*sin(4*pi*t) + 10*cos(8*pi*t) + 15*sin(10*pi*t + pi/4) + 10*cos(14*pi*t);

% Initialize the sinusoidal and cosine basis function matrices
B_sin = zeros(length(t), length(f));
B_cos = zeros(length(t), length(f));

% Generate each basis function and store it in the corresponding matrix
for i = 1:length(f)
    B_sin(:,i) = sin(2*pi*f(i)*t);
    B_cos(:,i) = cos(2*pi*f(i)*t);
end

% Project the signal onto the basis functions
p_i = y * B_sin;
q_i = y * B_cos;
m_i = sqrt(p_i.^2 + q_i.^2);

% Create a figure with four subplots
figure;
subplot(2, 2, 1);
plot(t, y);
title('Signal y(t)');
xlabel('Time (sec)');
ylabel('Amplitude');

subplot(2, 2, 2);
stem(f, p_i);
title('Projection coefficients p(i)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot(2, 2, 3);
stem(f, q_i);
title('Projection coefficients q(i)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

subplot(2, 2, 4);
stem(f, m_i);
title('Magnitude of projection coefficients m(i)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%Question 1c
%The values of m_i indicate the magnitude of the projection coefficients for each basis function i,
%The magnitude of the projection coefficient indicates the amplitude of the corresponding frequency component in the signal.

%The largest value of m_i occurs at i = 3, corresponding to the frequency component of the sinusoidal function with a frequency of f = 6 Hz. This indicates that the frequency component with a frequency of 6 Hz has the largest amplitude in the signal y(t).
%The second largest value of m_i occurs at i=4, corresponding to the frequency component of the sine function with a frequency of f=8 Hz. This indicates that the frequency component with 8 Hz frequency has the second largest amplitude in the signal y(t).
%Smaller values of m_i for other frequencies indicate that the corresponding frequency component has a smaller amplitude in the signal y(t).

%Question 2

Y = fft(y, 1000);

%help
Y_shift = fftshift(fftfreq(1000, 0.001));
figure;
plot(Y_shift, abs(Y));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude spectrum of DFT of y');

% Compute the real and imaginary parts of Y
Y_real = real(Y);
Y_imag = imag(Y);

% Compute the magnitude of Y
Y_mag = abs(Y);

% Plot the real, imaginary, and magnitude parts of Y using stem and subplot
figure;
subplot(3,1,1);
stem(Y_real);
title('Real part of Y');
subplot(3,1,2);
stem(Y_imag);
title('Imaginary part of Y');
subplot(3,1,3);
stem(Y_mag);
title('Magnitude of Y');
xlabel('Frequency (Hz)');

%2c

%The real part Y_real of the DFT is equivalent to the cosine projection q_i of the signal on the basis function, while the imaginary part Y_imag is equivalent to the sine projection p_i on the basis function.
%The amplitude spectrum of the DFT is equivalent to the projection length m(i) calculated from the sine and cosine projections.
%The magnitude spectrum plot of the DFT has the same shape and peak as the magnitude Y_mag plotted in 2.2. The frequency values in both plots are also the same.
%The real and imaginary parts of the DFT show the phase relationship between the sine and cosine components. In the amplitude spectrogram in Question 1(b), the phase information is not shown.
%The graph drawn with the stem command used in 2.2 makes it easier to see the individual frequency components of the DFT, while the graph in Question 1(b) shows the overall shape of the amplitude spectrum.

%Question 3
% Define the signal yq3
yq3 = [1 zeros(1,7)];

% Compute the FFT of yq3 and compute the magnitude spectrum
Yq3 = fft(yq3);
Yq3_mag = abs(Yq3);

% Plot yq3 and its magnitude spectrum using subplot and stem
figure;
subplot(2,1,1);
stem(yq3);
title('yq3');
ylabel('Amplitude');
subplot(2,1,2);
stem(Yq3_mag);
title('Magnitude Spectrum of yq3');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%Question 4
% Define the signal yq4
yq4 = ones(1,8);

% Compute the FFT of yq4 and compute the magnitude spectrum
Yq4 = fft(yq4);
Yq4_mag = abs(Yq4);

% Plot yq4 and its magnitude spectrum using subplot and stem
figure;
subplot(2,1,1);
stem(yq4);
title('yq4');
ylabel('Amplitude');
subplot(2,1,2);
stem(Yq4_mag);
title('Magnitude Spectrum of yq4');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%Q5
% Define the time vector
t = 0:0.0001:0.9999;

% Define the signal yq5
f = 300; % Hz
yq5 = cos(2*pi*f*t);

% Compute the FFT of yq5 and compute the magnitude spectrum
Yq5 = fft(yq5);
Yq5_mag = abs(Yq5);

% Plot yq5 and its magnitude spectrum using subplot and stem
figure;
subplot(2,1,1);
stem(yq5);
title('yq5');
ylabel('Amplitude');
subplot(2,1,2);
stem(Yq5_mag);
title('Magnitude Spectrum of yq5');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Find the frequency with the highest magnitude
[max_mag, max_idx] = max(Yq5_mag);
freq_resolution = 1 / (t(2) - t(1));
freq_max = (max_idx - 1) * freq_resolution;
fprintf('The frequency with the highest magnitude is %f Hz.\n', freq_max);
fprintf('Its index in the FFT array is %d.\n', max_idx);

%Q6

% Define the time vector
t = 0:0.0001:0.9999;

% Define the signal yq6
f = 400; % Hz
yq6 = sin(2*pi*f*t);

% Compute the FFT of yq6 and compute the magnitude spectrum
Yq6 = fft(yq6);
Yq6_mag = abs(Yq6);

% Plot yq6 and its magnitude spectrum using subplot and stem
figure;
subplot(2,1,1);
stem(yq6);
title('yq6');
ylabel('Amplitude');
subplot(2,1,2);
stem(Yq6_mag);
title('Magnitude Spectrum of yq6');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Find the frequency with the highest magnitude
[max_mag, max_idx] = max(Yq6_mag);
freq_resolution = 1 / (t(2) - t(1));
freq_max = (max_idx - 1) * freq_resolution;
fprintf('The frequency with the highest magnitude is %f Hz.\n', freq_max);
fprintf('Its index in the FFT array is %d.\n', max_idx);


diary off




%Function definition
function freqs = fftfreq(n, d)
% Compute the FFT frequency values for a signal of length n with sampling interval d
freqs = (0:n-1) / (n*d);
if mod(n, 2) == 0
    freqs(n/2+1) = -1 / (2*d);
    freqs(n/2+2:end) = flipud(freqs(2:n/2));
else
    freqs((n+1)/2:end) = flipud(freqs((n+1)/2:n))';
end
end



