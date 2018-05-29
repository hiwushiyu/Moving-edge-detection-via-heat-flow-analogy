%used for transforming processed consecutive frames to a video 

vedio = VideoWriter('non-max.avi'); %initialize a ".avi" file
vedio.FrameRate = 30;%FrameRate
open(vedio);
for i=1:181  %number of frames
    fname=strcat('out',num2str(i,'%d'),'.png');
    frame = imread(fname);
    writeVideo(vedio,frame);
end
close(vedio);