clc; clear all; close all;
%天线组单信道测试代码

isCompensateAll  = 1;  % 1代表整个蓝牙包做频率补偿，0代表对每个天线进行频率补偿
M = 2;%天线组天线个数，如果测三天线只需要改为3
d = 6.25e-2;%天线阵列间距6.25cm
lambda = 3e8/2.4e9;%蓝牙传输频段为2.4G
theta = -60:0.5:60; %MUSIC谱峰搜索的范围
A = exp(-1i*2*pi*[0:M-1]'*d/lambda*sin(theta*pi/180));
gt=20;%测试角度
dir = 'G:\实验室事务\导航组\室内定位\蓝牙AOA\测试数据\190401美迪索科大厅的NOKIA改装天线测试\';
for kn=1:20
    iq=dlmread([dir,num2str(gt),'deg/',num2str(kn,'%2d'),'.dat'],'',1,0);
    if kn == 60
        a=1;
    end
    data = 1i*iq(1:2:end) + iq(2:2:end);
%     for k=1:10
%         x0 = [x0, [data((k-1)*N*3+(1:16)).'; data((k-1)*N*3+(17:32)).'; data((k-1)*N*3+(33:48)).']]; %第一行是天线1第二行是天线2第三行是天线3
%     end
    plot(abs(data))
%     data = data./abs(data) %归一化
    
    figure(111);
    plot(1:length(data),real(data),'r--*',1:length(data),imag(data),'b--*');grid;xlim([1,512]);%天线实部与虚部的波形图
%     subplot(212);plot(1:Ns,imag(x0(1,:)),'r--*',1:Ns,imag(x0(2,:)),'b--*',1:Ns,imag(x0(3,:)),'g--*');grid;xlim([1,Ns]);

    NFFT = 2048;      
    if isCompensateAll == 1
        [~,Fs]=max(abs(fft(data(1:end),NFFT)));
        Fs = (Fs -1)/NFFT*4e6;
        delta_f = 250e3-Fs;
        data = data.*exp(1i*2*pi*(delta_f)*(1:size(data,1))/4e6).';
        data2 =[];
    elseif isCompensateAll == 0
            if M == 2
                for i = 0:5
                    [~,Fs]=max(abs(fft(data(80*i+1:80*(i+1)),NFFT)));
                    Fs = (Fs -1)/NFFT*4e6;
                    delta_f = 250e3-Fs;
                    data2 =[data2 ;data(80*i+1:80*(i+1)).*exp(1i*2*pi*(delta_f)*(80*i+1:80*(i+1))/4e6).'];
                end
            elseif M == 3
                 %TODO
            end
     end
     
 %   data =data2;
    figure(113);
    plot((0:NFFT-1)/NFFT*4e6,abs(fft(data(1:end),NFFT)));% FFT的频谱图，横坐标表示的是采样频率，见采样定理（DSP，P116)
    figure;freqz(data(1:480),NFFT) %Matlab自有函数画的频谱图
    
    x0 = [];

    if M==3
    T = 160;
    offT = 20;
    N = 130;
    for k=1:1
        x0 =[x0,  [data((k-1)*T*3+(1:N)+offT).';
              data((k-1)*T*3+(1:N)+offT+T).';
              data((k-1)*T*3+(1:N)+offT+T).';]];  
    end
    elseif M==2
        T = 80;
    offT = 20;
    N = 50;
    
    for k=1:3
        x0 =[x0,  [data((k-1)*T*2+(1:N)+offT).';
              data((k-1)*T*2+(1:N)+offT+T).';]];  
    end
    
%     figure;   plot((0:NFFT-1)/NFFT*4e6,abs(fft(x0(1,1:end),NFFT)));
%     figure;   plot((0:NFFT-1)/NFFT*4e6,abs(fft(x0(2,1:end),NFFT)));
end
%     Ns = N*10;
%     figure(kn);
%     subplot(211);plot(1:Ns,real(x0(1,:)),'r--*',1:Ns,real(x0(2,:)),'b--*',1:Ns,real(x0(3,:)),'g--*');grid;xlim([1,Ns]);
%     subplot(212);plot(1:Ns,imag(x0(1,:)),'r--*',1:Ns,imag(x0(2,:)),'b--*',1:Ns,imag(x0(3,:)),'g--*');grid;xlim([1,Ns]);

%     figure(kn+20);
%     subplot(211);plot(1:M*N,real([x0(1,:),x0(2,:),x0(3,:)]),'r--*');grid;xlim([1,M*N]);
%     subplot(212);plot(1:M*N,imag([x0(1,:),x0(2,:),x0(3,:)]),'r--*');grid;xlim([1,M*N]);
%     
    Rs = x0*x0'/N;
%     plot(real(x0(1,:)),'r--*')
    %--- delay-and-sum Beamformer --------------------
   for kt=1:length(theta)
        dsb(kt) = abs(A(:,kt)'*Rs*A(:,kt));
    end
    figure(112);subplot(221);plot(theta,dsb);grid;title('CB');
    doas(kn,1) = find_doa_peaks(dsb,theta,1);

    
    %--- music --------------------
    [H,d]=eig(Rs);
    vH = A'*H(:,(1:M-1));
    music(:) = 1./sum(abs(vH).^2,2);
    figure(112);subplot(222);plot(theta,music);grid;title('MUSIC');
    doas(kn,2) = find_doa_peaks(music,theta,1);
 
%     ylim([thetas(1) thetas(end)]); zlim([0 1.1]); view(90,0);
   
    
    %--- mvdr --------------------
    Rsi = inv(Rs);
    for kt=1:length(theta)
        mvdr(kt) = 1/abs(A(:,kt)'*Rsi*A(:,kt));
    end
    figure(112);subplot(223);plot(theta,mvdr);grid;title('MVDR');
    doas(kn,3) = find_doa_peaks(mvdr,theta,1);
%     figure(112);subplot(223);plot(theta,mvdr(k1,:));grid;title('MVDR');
%     doas(kn,3) = find_doa_peaks(mvdr(k1,:),theta,1);
%     
 
%     [m1, im1] = max(dsb); 
%     [m2, im2] = max(m1);
%     DOAe = theta(im2)

    
%     [X,Y] = meshgrid(theta, Fs);
%     figure(2); subplot(221); mesh(X,Y,dsb/max(dsb(:))); title('DSB');
%     ylabel('Frequency'); xlabel('AoA (degree)');view([0,90]);
%     subplot(222); mesh(X,Y,music/max(music(:))); title('MUSIC');
%     ylabel('Frequency'); xlabel('AoA (degree)');view([0,90]);
%     subplot(223); mesh(X,Y,mvdr/max(mvdr(:))); title('MVDR');
%     ylabel('Frequency'); xlabel('AoA (degree)');view([0,90]);

%     return;
    
    %-----Sparse method, (does not need source number)------------
    sparse = samv(A, Rs, N);
    figure(112);subplot(224);plot(theta,sparse);grid;title('Sparse');
    doas(kn,4) = find_doa_peaks(sparse,theta,1);
   
end

for k=1:size(doas,2)
    RMSE(k) = sqrt(mean((doas(:,k) - mean(doas(:,k))).^2));
end

doas;
RMSE
MEAN = mean(doas,1)
Median = median(doas,1)

figure(1);
plot(1:kn,doas(:,1),'r--o',1:kn,doas(:,2),'b--o',1:kn,doas(:,3),'b--*',1:kn,doas(:,4),'g--^',1:kn,-ones(1,kn)*gt,'k-');grid;
legend('CB','MUSIC','MVDR','SAMV','Ground truth','Location','Best');xlim([1,kn+round(kn/3)]);
ylabel('AoA (deg)');
% title(sprintf('',));
