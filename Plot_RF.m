close all;
clear all;
clc


file_out=('E:\Kuliah\Semester 8\Tugas Akhir\Pengolahan data\Data RF TA Kepulauan Mentawai\Siberut Island\*Signal*Z05S*.mat');
file_out2=dir(file_out);
for n=1:length(file_out2)
    filename= file_out2.name;
    A = importdata(filename);
    fs= A.Vert.stats.sampling_rate;
    A.vert=A.Vert.data;
    taxis = 1:1:8500;
    maxi = 8501;
    taxis2 = (0:(maxi-1))/fs;
   
%     b = A.Vert.stats.sampling_rate;
    
    figure(1)
    plot(taxis2(1:8500),A.vert(1:8500));title(filename); hold on
    
end


