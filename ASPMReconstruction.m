clc;clear all; close all;
addpath('./data');
tic;
M =16; % CCDƽ�������
iter_max = 500; %����������
wavelength = 532e-6;
factor=1;%��������
pp = 3.8e-3/factor;
% phase_size = 500; % ���������С
phase_size = 500*factor; % ���������С

prop_d = 68.3; %�������
prop_d = 66.3; %�������

% ����ʱ����RMSE����
RMSE_x = linspace(0,iter_max-1,iter_max);
RMSE_800iter = RMSE_x;

num = 0;

for prop_d= 64:1:64 %Ӧ����64 64/63���  2020 12 07
phi_v = 1.5; %֮ǰ��pattern ���ֵ��0.5 �˴�ӦΪ0.5�ı��� % ӦΪ1.8 ��� 2020 12 07
 disp(num2str(prop_d));

% �����ʼ�����λ
U = randn(phase_size);
U = U-(min(min(U)));
U = U./(max(max(U)));
U = exp(1i*2*pi.*U); % ��ʼ�²�������λ��
U = padarray(U, [phase_size phase_size]);
Ui = U;
for iter = 1:1:iter_max %��������������
    

 %% ������ͼ
   num = mod(iter,M); %��ǰӦ���õڼ���ͼ
    num = num+1;
    nm = strcat('modulation_',num2str(num),'.mat');% ������λ
    nm2 = strcat('reshaped_obj_',num2str(num),'.mat');% ���㵽��ǿ����Ϣ

   %% 
    load (nm,'modulation_grayscale'); 
%     modulation_grayscale = imresize(modulation_grayscale,factor,'nearest');
    modulation_grayscale = imresize(modulation_grayscale,factor);

    load (nm2,'OBJ_reshaped'); 
    OBJ_reshaped = imresize(OBJ_reshaped,factor);
%     OBJ_reshaped = imresize(OBJ_reshaped,factor,'nearest');


    RP_reshaped = modulation_grayscale.*phi_v*pi;
%     %% 7.5�޸� SLM��CCDˮƽ������ӦΪ�����ϵ ����ccd���·ŷ��� ����������� ���·�ת
%     if mod(iter+1,1)
%     RP_reshaped = imrotate(RP_reshaped,180);
%     end
% %     
    Int_CCD = OBJ_reshaped;
    
   %%%%%%%%%%%%%


    % reshape the modulation


    % �²������λ��������֪���������λ����
    Ui_p = Ui; %�ϴε�����Ĳ²���λ��
    Ui = Ui.*exp(1i.*RP_reshaped);
    
    % ��������CCDƽ��
    Ui_ccd = ASM_diffraction(wavelength, Ui, prop_d,pp);
    
    % ����ǿ��Լ�� �滻ǿ����Ϣ
    Ui_ccd_constrained = Int_CCD.^(0.5).*(exp(1i.*angle(Ui_ccd)));
    
    % �滻�󷴴���
    Ui_slm = ASM_diffraction(wavelength, Ui_ccd_constrained, -1*prop_d,pp );
    
    % ȥ��������λ ��õ�����Ĳ²���λ
    Ui = Ui_slm .* exp (-1*1i.*RP_reshaped);
%     Ui = Ui+((exp (-1*1i.*RP_reshaped)./(max(max(abs(exp (-1*1i.*RP_reshaped)))))).^2.*(Ui-Ui_slm)); %������������ĸ��¹�ʽ �ù�ʽ��Ч
    
    
    
    % ��������
%     Constrain = zeros(1500);
%     Constrain(phase_size+1:2*phase_size,phase_size+1:2*phase_size)=Ui(phase_size+1:2*phase_size,phase_size+1:2*phase_size);
%     Ui =  Constrain;
    % �������ε���֮��Ĳ�
    ph_Ui = angle(Ui(phase_size+1:2*phase_size,phase_size+1:2*phase_size));% �ź�������λ
    uph_Ui = unwrap(ph_Ui);
%     ph_p_Ui = angle(Ui_p(phase_size+1:2*phase_size,phase_size+1:2*phase_size));% ǰһ�ε������ź�������λ
    ph_p_Ui = angle(Ui_p(phase_size+1:2*phase_size,phase_size+1:2*phase_size));% ǰһ�ε������ź�������λ

    uph_p_Ui = unwrap(ph_p_Ui);
    Err = abs(ph_Ui - ph_p_Ui);
    RMSE_err = (mean2(Err.^2)/pi).^0.5;
    
    RMSE_x(iter+1) = RMSE_err;
    disp(num2str(iter));

 
end

% % ���� ����� ���㿪ʼ
figure
subplot(1,2,1)
imshow (angle(Ui), [0 max(max(angle(Ui)))]);
title('ȫ����λdistance')
subplot(1,2,2)
imshow(abs(Ui),[0 (max(max(abs(Ui))))^0.25]);
title('ȫ�����')
% 
Ui_disp = Ui(phase_size+1:2*phase_size,phase_size+1:2*phase_size);
phase_Ui_disp = angle(Ui_disp);
% phase_Ui_disp = unwrap(angle(Ui_disp));
% phase_Ui_disp = phase_Ui_disp-min(min(phase_Ui_disp));



% figure
% subplot(1,2,1)
% imshow (angle(Ui_disp), [0 max(max(angle(Ui)))]);
% title('����������λ')
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
% title('�����������')
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
title('����������λ')
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
% title('�����������')

toc;