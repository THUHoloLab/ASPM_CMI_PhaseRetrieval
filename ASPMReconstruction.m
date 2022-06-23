clc;clear all; close all;
addpath('./data');
tic;
M =16; % CCD平面采样数
iter_max = 500; %最大迭代次数
wavelength = 532e-6;
factor=1;%缩放数据
pp = 3.8e-3/factor;
% phase_size = 500; % 测试区域大小
phase_size = 500*factor; % 测试区域大小

prop_d = 68.3; %衍射距离
prop_d = 66.3; %衍射距离

% 迭代时绘制RMSE曲线
RMSE_x = linspace(0,iter_max-1,iter_max);
RMSE_800iter = RMSE_x;

num = 0;

for prop_d= 64:1:64 %应该是64 64/63最好  2020 12 07
phi_v = 1.5; %之前的pattern 最大值是0.5 此处应为0.5的倍数 % 应为1.8 最好 2020 12 07
 disp(num2str(prop_d));

% 输入初始随机相位
U = randn(phase_size);
U = U-(min(min(U)));
U = U./(max(max(U)));
U = exp(1i*2*pi.*U); % 初始猜测的随机相位场
U = padarray(U, [phase_size phase_size]);
Ui = U;
for iter = 1:1:iter_max %设置最大迭代次数
    

 %% 连续读图
   num = mod(iter,M); %当前应该用第几幅图
    num = num+1;
    nm = strcat('modulation_',num2str(num),'.mat');% 调制相位
    nm2 = strcat('reshaped_obj_',num2str(num),'.mat');% 拍摄到的强度信息

   %% 
    load (nm,'modulation_grayscale'); 
%     modulation_grayscale = imresize(modulation_grayscale,factor,'nearest');
    modulation_grayscale = imresize(modulation_grayscale,factor);

    load (nm2,'OBJ_reshaped'); 
    OBJ_reshaped = imresize(OBJ_reshaped,factor);
%     OBJ_reshaped = imresize(OBJ_reshaped,factor,'nearest');


    RP_reshaped = modulation_grayscale.*phi_v*pi;
%     %% 7.5修改 SLM与CCD水平方向上应为镜像关系 由于ccd上下放反了 因而左右正常 上下翻转
%     if mod(iter+1,1)
%     RP_reshaped = imrotate(RP_reshaped,180);
%     end
% %     
    Int_CCD = OBJ_reshaped;
    
   %%%%%%%%%%%%%


    % reshape the modulation


    % 猜测随机相位场经过已知加载随机相位调制
    Ui_p = Ui; %上次迭代后的猜测相位场
    Ui = Ui.*exp(1i.*RP_reshaped);
    
    % 正传播至CCD平面
    Ui_ccd = ASM_diffraction(wavelength, Ui, prop_d,pp);
    
    % 利用强度约束 替换强度信息
    Ui_ccd_constrained = Int_CCD.^(0.5).*(exp(1i.*angle(Ui_ccd)));
    
    % 替换后反传播
    Ui_slm = ASM_diffraction(wavelength, Ui_ccd_constrained, -1*prop_d,pp );
    
    % 去掉调制相位 获得迭代后的猜测相位
    Ui = Ui_slm .* exp (-1*1i.*RP_reshaped);
%     Ui = Ui+((exp (-1*1i.*RP_reshaped)./(max(max(abs(exp (-1*1i.*RP_reshaped)))))).^2.*(Ui-Ui_slm)); %王冰洋文章里的更新公式 该公式有效
    
    
    
    % 限制区域
%     Constrain = zeros(1500);
%     Constrain(phase_size+1:2*phase_size,phase_size+1:2*phase_size)=Ui(phase_size+1:2*phase_size,phase_size+1:2*phase_size);
%     Ui =  Constrain;
    % 计算两次迭代之间的差
    ph_Ui = angle(Ui(phase_size+1:2*phase_size,phase_size+1:2*phase_size));% 信号区域相位
    uph_Ui = unwrap(ph_Ui);
%     ph_p_Ui = angle(Ui_p(phase_size+1:2*phase_size,phase_size+1:2*phase_size));% 前一次迭代的信号区域相位
    ph_p_Ui = angle(Ui_p(phase_size+1:2*phase_size,phase_size+1:2*phase_size));% 前一次迭代的信号区域相位

    uph_p_Ui = unwrap(ph_p_Ui);
    Err = abs(ph_Ui - ph_p_Ui);
    RMSE_err = (mean2(Err.^2)/pi).^0.5;
    
    RMSE_x(iter+1) = RMSE_err;
    disp(num2str(iter));

 
end

% % 剪裁 解包裹 从零开始
figure
subplot(1,2,1)
imshow (angle(Ui), [0 max(max(angle(Ui)))]);
title('全场相位distance')
subplot(1,2,2)
imshow(abs(Ui),[0 (max(max(abs(Ui))))^0.25]);
title('全场振幅')
% 
Ui_disp = Ui(phase_size+1:2*phase_size,phase_size+1:2*phase_size);
phase_Ui_disp = angle(Ui_disp);
% phase_Ui_disp = unwrap(angle(Ui_disp));
% phase_Ui_disp = phase_Ui_disp-min(min(phase_Ui_disp));



% figure
% subplot(1,2,1)
% imshow (angle(Ui_disp), [0 max(max(angle(Ui)))]);
% title('成像区域相位')
% 
% title(strcat('d=',num2str(prop_d),'phi_v',num2str(phi_v)));
% subplot(1,2,2)
% % imshow(unwrap(angle(Ui_disp)), [0 max(max(angle(Ui)))]);
% % title('Unwrapped GS reconstruction');
% 
% % figure 
% % mesh(phase_Ui_disp./pi);
% % title('')
figure
plot(RMSE_x);
title('RMSE between adjecnt iteration');
xlim([0 iter_max])
ylim([0.5 1.2])
end
% figure
% imshow(abs(Ui_disp),[]);
% title('成像区域振幅')
for recon_prop = -115:-1:-115

Ui_prop = ASM_diffraction(wavelength, Ui, recon_prop,pp);
Ui_prop_disp = Ui_prop(phase_size+1:2*phase_size,phase_size+1:2*phase_size);
phase_Ui_disp_prop = angle(Ui_prop_disp);
phase_Ui_disp = unwrap(angle(Ui_disp));
phase_Ui_disp_prop = phase_Ui_disp_prop-min(min(phase_Ui_disp_prop));

figure
subplot(1,2,1)
imshow (angle(Ui_prop_disp), [0 max(max(angle(Ui)))]);
subplot(1,2,2)
imshow (angle(Ui_prop), [0 max(max(angle(Ui)))])
title('after propogation');
title('成像区域相位')
end


% rec_filed_backp = Ui(501:1000,501:1000);
rec_filed_backp = Ui(phase_size+1:2*phase_size,phase_size+1:2*phase_size);

% rec_filed_backp = Ui(401:1100,401:1100);

b=padarray(rec_filed_backp,[phase_size phase_size]);
c=ASM_diffraction(wavelength,b,-115,pp);
d = angle(c(phase_size+1:2*phase_size,phase_size+1:2*phase_size));
figure
imshow(d, []);
colormap(parula);

rgb = ind2rgb(gray2ind(d,255),jet(255));
figure
imshow(rgb)
% figure
% imshow (real(Ui_prop_disp).^2, [0 max(max(real(Ui)))^0.1]);
% title('成像区域振幅')

toc;