%%%% 对图像做预处理 

close all; clear all;
addpath('./data')
% 读取拍摄到的图像
mark = imread('mark.png');
mark = im2double(mark);

ref = imread('ref_.png');
ref = im2double(ref);

% mark =sqrt(mark./ref);% Normalized hologram1
% MI = mark - ref;
MI = mark./ ref;

figure
imshow(MI, [0 2]);

% reshape the figure 确定剪裁参数
 stx_i = 1638;
 sty_i = 1349;
 edx_i = 3119;
 edy_i = 2832;
 MI_reshaped = MI(sty_i:edy_i,stx_i:edx_i);  % 振幅constrain


MI_reshaped = imresize(MI_reshaped, [1500 1500]);
%  
figure
imshow(MI_reshaped, []);

for num = 1:16 %设置最大迭代次数
    %当前应该用第几幅图
    nm = strcat('pp',num2str(num),'.png');%      
    % 读取拍摄到的图像
    OBJ = imread(nm);
    OBJ = im2double(OBJ);
    %归一化 并减去背景
%     OBJ=sqrt(OBJ./ref);% Normalized hologram1
% %     OBJ = (OBJ - min(min(OBJ)))./(max(max(OBJ))-min(min(OBJ)));
% %     ref = (ref - min(min(ref)))./(max(max(ref))-min(min(ref)));
%     OBJ_pro = OBJ - ref;
    OBJ_pro = OBJ ./ (ref+1);
    OBJ_reshaped = OBJ_pro(sty_i:edy_i,stx_i:edx_i);  % 振幅constrain
    OBJ_reshaped = imresize(OBJ_reshaped, [1500 1500]);
    OBJ_reshaped = OBJ_reshaped./(max(max(OBJ_reshaped)));% 归一化
    %%%%%%%%%%%%%%%%%%%%%%% 直接处理 %%%%%%%%%%%%%%%%%%%%%%%
    
    image_shift = 240; % 像平面 图像移动

    figure
    imshow(abs(OBJ_reshaped), [0 1]);   
    nm2 = strcat('reshaped_obj_',num2str(num),'.mat');% 
    save (nm2,'OBJ_reshaped');
    
end