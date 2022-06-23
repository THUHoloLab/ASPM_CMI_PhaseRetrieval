%%%% ��ͼ����Ԥ���� 

close all; clear all;
addpath('./data')
% ��ȡ���㵽��ͼ��
mark = imread('mark.png');
mark = im2double(mark);

ref = imread('ref_.png');
ref = im2double(ref);

% mark =sqrt(mark./ref);% Normalized hologram1
% MI = mark - ref;
MI = mark./ ref;

figure
imshow(MI, [0 2]);

% reshape the figure ȷ�����ò���
 stx_i = 1638;
 sty_i = 1349;
 edx_i = 3119;
 edy_i = 2832;
 MI_reshaped = MI(sty_i:edy_i,stx_i:edx_i);  % ���constrain


MI_reshaped = imresize(MI_reshaped, [1500 1500]);
%  
figure
imshow(MI_reshaped, []);

for num = 1:16 %��������������
    %��ǰӦ���õڼ���ͼ
    nm = strcat('pp',num2str(num),'.png');%      
    % ��ȡ���㵽��ͼ��
    OBJ = imread(nm);
    OBJ = im2double(OBJ);
    %��һ�� ����ȥ����
%     OBJ=sqrt(OBJ./ref);% Normalized hologram1
% %     OBJ = (OBJ - min(min(OBJ)))./(max(max(OBJ))-min(min(OBJ)));
% %     ref = (ref - min(min(ref)))./(max(max(ref))-min(min(ref)));
%     OBJ_pro = OBJ - ref;
    OBJ_pro = OBJ ./ (ref+1);
    OBJ_reshaped = OBJ_pro(sty_i:edy_i,stx_i:edx_i);  % ���constrain
    OBJ_reshaped = imresize(OBJ_reshaped, [1500 1500]);
    OBJ_reshaped = OBJ_reshaped./(max(max(OBJ_reshaped)));% ��һ��
    %%%%%%%%%%%%%%%%%%%%%%% ֱ�Ӵ��� %%%%%%%%%%%%%%%%%%%%%%%
    
    image_shift = 240; % ��ƽ�� ͼ���ƶ�

    figure
    imshow(abs(OBJ_reshaped), [0 1]);   
    nm2 = strcat('reshaped_obj_',num2str(num),'.mat');% 
    save (nm2,'OBJ_reshaped');
    
end