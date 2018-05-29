%used for processing gait video 
clear all
close all
clc
obj=VideoReader('lab.dv');%read in gait video
k=4;          
lambda=0.2;   
N=50;           
%Isotropy
K2=0.45;          
N2=10;         
%Variable
nframes = get(obj, 'NumberOfFrames');
for i=1:nframes
    Frame1 = read(obj, i);
    Frame2 = read(obj, i+1);
    Frame3 = read(obj, i+2);   
    %consecutive three frames
    Frame1 = sobel_ani(Frame1,k,lambda,N);
    Frame2 = sobel_ani(Frame2,k,lambda,N);
    Frame3 = sobel_ani(Frame3,k,lambda,N);

    [I,HFO] = isotropy(Frame1,Frame2,Frame3,K2,N2);  
    %Non-Maxima_Suppression
    [HFO_max] = Maxima_Suppression(HFO);
    [HFO_Hys] = Hysteresis_Thresholding(uint8(HFO_max));
    OutputDir = '/Users/Chestnuts/Dropbox/ACV/P3-Moving-edge detection by heat flow/S/test/out';%store all processed frames
    imwrite(uint8(HFO_Hys),[OutputDir,int2str(i),'.png']);
end
 

    
 