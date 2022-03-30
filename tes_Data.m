clear all
close all
clc
% cd ('D:\Perkuliahan\Tugas Akhir (Yok pasti bisa)\Data')
file_out=('E:\Kuliah\Semester 8\Tugas Akhir\Pengolahan data\Data RF TA Kepulauan Mentawai\Sipura Island\*Signal*SOBS*.mat');
file_out2=dir(file_out);

flow=0.1;
fhigh=2;
for n=1:length(file_out2)
    filename= file_out2.name;
    A = importdata(filename);
    A.Vert.stats.sampling_rate
    
    fs= A.Vert.stats.sampling_rate;
end

[B1,A1] = butter(2,[flow/(fs/2) fhigh/(fs/2)]);
% freqz(B1,A1)
     alph=4;
     wl=.001;

if fs > 50
    L= 17000
else
    L= 8500
end

for n=1:length(file_out2)
    filename1=file_out2(n).name
    A=importdata(filename1);
    
    komponen = A.Vert.stats.channel
    samp_rate = A.Vert.stats.sampling_rate
    magnitude = A.Header(7)
    ArcDisatnace = A.Header(8)
   
    A.vert=A.Vert.data;
    A.vert_ori=A.vert;
    A.East=A.East.data;
    A.North=A.North.data;
    
    A.vert=filtfilt(B1,A1,detrend(A.vert));
    A.East=filtfilt(B1,A1,detrend(A.East));
    A.North=filtfilt(B1,A1,detrend(A.North));
    
    figure(1)
    plot(A.vert);title(filename1);legend(komponen);xlabel(samp_rate); hold on
    
    A.vert= fft(A.vert)
    P2 = abs(A.vert/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;
    figure(2)
    plot(f,P1) 
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    
    if fs > 50
        A.vert =A.Vert.data(100:2000);
    else
        A.vert =A.Vert.data(50:1000);
    end
    noise =filtfilt(B1,A1,detrend(A.vert));
    noise_power = noise.^2;
    noise_average_power = mean(noise_power);
    noise_power_db = 10 * log10(noise_average_power)
    
    if fs > 50
        A.vert =A.Vert.data(2000:17000);
    else
        A.vert =A.Vert.data(1000:8500);
    end
    signal =filtfilt(B1,A1,detrend(A.vert));
    signal_power = signal.^2;
    signal_average_power = mean(signal_power);
    signal_power_db = 10 * log10(signal_average_power)

    SNR_db = signal_power_db - noise_power_db
    SNRp = 10.^(SNR_db/10)
    
end

help detrend
