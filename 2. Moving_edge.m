clear all;  
close all;  
clc;  
 
%Anisotropy
k=4;           %导热系数,控制平滑  
lambda=0.2;    %控制平滑  
N=50;           %迭代次数  

%Isotropy
K2=0.45;           %导热系数,控制平滑  
N2=10;           %迭代次数  

I1=imread('11.jpg');  
I1= double(rgb2gray(I1));

I2=imread('22.jpg');  
I2= double(rgb2gray(I2));

I3=imread('33.jpg');  
I3= double(rgb2gray(I3));

% figure;  
% imshow(uint8(I1));
% title('1');
% figure;  
% imshow(uint8(I2));
% title('2');
% figure;  
% imshow(uint8(I3));
% title('3');

I1 = anisotropy(I1,k,lambda,N);
I2 = anisotropy(I2,k,lambda,N);
I3 = anisotropy(I3,k,lambda,N);

% figure;  
% imshow(uint8(I1));
% title('1 anisotropy');
% figure;  
% imshow(uint8(I2));
% title('2 anisotropy');
% figure;  
% imshow(uint8(I3));
% title('3 anisotropy');

I1 = sobel(I1);
I2 = sobel(I2);
I3 = sobel(I3);

figure;  
imshow(uint8(I1));
title('1 sobel');
figure;  
imshow(uint8(I2));
title('2 sobel');
figure;  
imshow(uint8(I3));
title('3 sobel');

[I,HFO] = isotropy(I1,I2,I3,K2,N2);

figure;  
imshow(uint8(I));
title('Total Heat Flow');
figure;  
imshow(uint8(HFO));
title('HFO');

[HFO_max] = Maxima_Suppression(HFO);
figure;  
imshow(uint8(HFO_max));
title('Maxima Suppression');


[HFO_Hys] = Hysteresis_Thresholding(uint8(HFO_max));
figure;  
HFO_Hys_bin=im2bw(HFO_Hys,0.5);
imshow(HFO_Hys_bin);
title('Hysteresis Thresholding');


function [img] = sobel(img) 
P1 = [1 0 -1; 2 0 -2; 1 0 -1];
P2 = [1 2 1; 0 0 0; -1 -2 -1];

I1 = imfilter(img,P1);
I2 = imfilter(img,P2);

img = sqrt(I1.^2 + I2.^2);
end

function [img] = anisotropy(img,k,lambda,N) 
[m n]=size(img);  
imgn=zeros(m,n);  
for i=1:N  
  
    for p=2:m-1  
        for q=2:n-1  
            %当前像素的散度，对四个方向分别求偏导，局部不同方向上的变化量，  
            %如果变化较多，就证明是边界，想方法保留边界  
            NI=img(p-1,q)-img(p,q);  
            SI=img(p+1,q)-img(p,q);  
            EI=img(p,q-1)-img(p,q);  
            WI=img(p,q+1)-img(p,q);  
              
            %四个方向上的导热系数，该方向变化越大，求得的值越小，从而达到保留边界的目的  
            cN=exp(-NI^2/(k*k));  
            cS=exp(-SI^2/(k*k));  
            cE=exp(-EI^2/(k*k));  
            cW=exp(-WI^2/(k*k));  
              
            imgn(p,q)=img(p,q)+lambda*(cN*NI+cS*SI+cE*EI+cW*WI);  %扩散后的新值        
        end  
    end  
      
    img=imgn;       %整个图像扩散完毕，用已扩散图像的重新扩散。  
end 
end

function [I,HFO] = isotropy(I1,I2,I3,K2,N2) 

Cur = I2;
Total = 0;

Diff_Tem = 0;
Total_HFO =0;

for i=1:N2  
     Diff = I1 + I3 - 2*Cur;
     Cur = Cur + K2*Diff;
     Total = Total + abs(Diff);
     
     Diff_Tem =  Diff;
     Diff_Tem(Diff_Tem>0) = 0;
     Total_HFO = Total_HFO + abs(Diff_Tem);
end
I = 0 + K2*Total;

HFO = K2*Total_HFO;
end

function [I] = Maxima_Suppression(I)
[a,b] = size(I);
tem = zeros(a,b);
direction = zeros(a,b);

P1 = [-1 1];
P2 = [-1;1];
 
I1 = imfilter(I,P1);
I2 = imfilter(I,P2);


%Step1
for i = 1:a
    for j = 1:b
        tem(i,j) = atand(I2(i,j)/I1(i,j));
        if (tem(i,j)<67.5)&&(tem(i,j)>22.5)  
            direction(i,j) =  0;    
        elseif (tem(i,j)<22.5)&&(tem(i,j)>-22.5)
            direction(i,j) =  3;    
        elseif (tem(i,j)<-22.5)&&(tem(i,j)>-67.5)
            direction(i,j) =  2;    
        else
            direction(i,j) =  1;    
        end   
    end    
end

%Step2
%   2 1 0
%   3 X 3
%y  0 1 2

I5= zeros (a,b);
for i = 2:a-1
    for j = 2:b-1    
        if 0 == direction(i,j) %0
            if (I(i,j)>I(i-1,j+1) )&&( I(i,j)>I(i+1,j-1)  )
                I5(i,j) = I(i,j);
            else
                I5(i,j) = 0;
            end
        elseif 1 == direction(i,j) %1
            if ( I(i,j)>I(i-1,j) )&&( I(i,j)>I(i+1,j)  )
                I5(i,j) = I(i,j);
            else
                I5(i,j) = 0;
            end
        elseif 2 == direction(i,j) %2
            if ( I(i,j)>I(i-1,j-1) )&&( I(i,j)>I(i+1,j+1)  )
                I5(i,j) = I(i,j);
            else
                I5(i,j) = 0;
            end
        elseif 3 == direction(i,j) %3
            if ( I(i,j)>I(i,j+1) )&&( I(i,j)>I(i,j-1)  )
                I5(i,j) = I(i,j);
            else
                I5(i,j) = 0;
            end
        end        
    end
end
I = I5;
end

function [I5] = Hysteresis_Thresholding(I5)
I6 = im2double(I5);
% I6 = I5;

[a,b] = size(I5);
Upper_Threshold = 0.7896;
Lower_Threshold = 0.1762;

Fake_edge_numbers = 0;
True_edge_numbers = 0;

for i = 1:a
    for j = 1:b  
     %Strong edge
     if(I6(i,j)>Upper_Threshold)
         I6(i,j) = 100;
     end
     %Non edge
     if(I6(i,j)<Lower_Threshold)
         I6(i,j) = 0;
     end
    end
end

 %Find the weak edge
 %50 means labeled of weak edge
 range = 3;
 
 for i = 1:a
    for j = 1:b
     if(I6(i,j)~= 50 && I6(i,j) < Upper_Threshold && I6(i,j) > Lower_Threshold)
     %Deep First Search
     clear LIFO
     clear NODE
     LIFO(1,1) = i; %Initialization
     LIFO(1,2) = j;  
     LIFO_queue_number = 1;
     
     Node{1,1} = i;
     Node{1,2} = j; 
     Node_number = 1;
     
     loop = 1;
     
     while 1
     %Finish Check
     if(LIFO_queue_number == 0) %run out the queue
         disp('A fake edge was detected');
         Fake_edge_numbers = Fake_edge_numbers +1;
         
         %Set the weak edge to 0
         for p = 1:Node_number
             I6(Node{p,1},Node{p,2}) = 0;
         end
         break
     end
     
     if(loop == 0) %run out the queue
         disp('A weak edge is saved！');
         True_edge_numbers = True_edge_numbers + 1;
         %Set the weak edge to strong edge
         for p = 1:Node_number
             I6(Node{p,1},Node{p,2}) = 100;
         end
         break
     end
     
     %Find 8-Connectivity Successor
     I6(LIFO(LIFO_queue_number,1),LIFO(LIFO_queue_number,2)) = 50; %label
     
     Temporary_i = LIFO(LIFO_queue_number,1);
     Temporary_j = LIFO(LIFO_queue_number,2);
     
     LIFO(LIFO_queue_number,:) =[];  %remove the node in the LIFO queue
     LIFO_queue_number = LIFO_queue_number -1;
      for q = -range : range
          for w = -range : range
              
              
              if(Temporary_i+q<a && Temporary_j+w<b && Temporary_i+q>0 && Temporary_j+w>0)
                if(I6(Temporary_i+q,Temporary_j+w) == 100)
                  loop=0;
                end
   
              
              if(I6(Temporary_i+q,Temporary_j+w) ~= 0 && I6(Temporary_i+q,Temporary_j+w) ~= 50 && I6(Temporary_i+q,Temporary_j+w) ~= 100 && q~=0 && w~=0 ) %8-Connectivity
                  I6(Temporary_i+q,Temporary_j+w) = 50; %label
                  LIFO_queue_number = LIFO_queue_number +1;
                  LIFO(LIFO_queue_number,1) = Temporary_i+q;
                  LIFO(LIFO_queue_number,2) = Temporary_j+w;
                  
                  
                  Node_number = Node_number +1;
                  Node{Node_number,1} = Temporary_i+q;;
                  Node{Node_number,2} = Temporary_j+w;;
               end
             end
          end
      end

     end
     end
    end
 end
disp(['Fake edge: ' num2str(Fake_edge_numbers)]);
disp(['True edge: ' num2str(True_edge_numbers)]);
I5 = I6;
end
