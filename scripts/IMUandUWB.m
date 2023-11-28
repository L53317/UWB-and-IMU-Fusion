%SINS与UWB组合导航
%SINS的数据采用IMU规范输出，UWB采用解算输出；
%简化版的SINS与UWB组合导航，采当地水平坐标系为导航坐标系，INS更新过程已简化：
%根据地球模型的高精度SINS组合系统，按照室内或者小范围环境的实际情况，简化状态转移矩阵；
%UWB位置解算采用LS算法；
%Kalman组合采用标准Kalman，采用松组合方式：状态向量为15维，观测量为6维。
%20190501.liu

clear;clc;
% %设置本文件所在路径为当前工作空间路径
[pathstr,namestr]=fileparts(mfilename('fullpath'));
cd(pathstr);%转到当前文件所在路径

glvs;  %建立SINS解算需要的全局变量
%生成轨迹上的IMU器件模拟信息
ts = 0.010; % imu采样时间
att0=[0;0;-90]*arcdeg/1; %初始姿态rad，与东北天的夹角，即载体坐标系（右前上）与导航坐标系（东北天）重合(北偏西为正，-π到π)
vn0=[0.2;0;0];         %初始速度m/s，东、北、天方向上的速度（需要对应初始姿态来设置,否则出错）
posimu0=[34.057*arcdeg;118.786*arcdeg;0];  %位置初值，采用纬经度和高度(m)
avpimu0=[att0;vn0;posimu0];  %组合成列向量，初始值信息
%%%%轨迹生成
%%%%    俯仰角速率  横滚就速率  方位角速率  纵向加速度  持续时间   %单位：°/s（若采用rad/s则后面不需再转换）,m/s,s
wat=[   0,          0,          0,          0,        7.5/0.2    %匀速或静止7.5/0.2=37.5s，运动采样数据共3750点
        0,          0,    180/(pi*2.5/0.2), 0,       pi*2.5/0.2  %半圆（半径2.5m,速度0.2m/s）运动，用时(pi*2.5/0.2)s
        0,          0,          0,          0,        7.5/0.2    %匀速或静止7.5/0.2=37.5s 
        ];                                  
wat(:,1:3)=wat(:,1:3)*arcdeg;    % 把角速度的单位由°/s 转换为 弧度/s
[attimu,vnimu,posimu] = trjprofile(att0,vn0,posimu0, wat,ts);%根据运动情况生成载体导航系下的真实的姿态（弧度）、速度和位置（弧度）
%[attimu,vnimu,posimu] = trjprofile(att0,vn0,posimu0, wat,ts,0);%不采用滤波器的版本
[wm,vm]=av2imu(attimu,vnimu,posimu,ts);%生成IMU的输出信息，wm是陀螺输出（真值），vm是加速度计输出（真值）

%组合导航解算
nn = 1;% 子样数
nts = nn*ts; % 进行一次解算需要的数据时间 
qnb0 = a2qua(att0);  %姿态初值，姿态角转四元数
pos0=[0;0;0];  %导航位置初值，采用东北天-xyz
avp0=[att0;vn0;pos0];
qnb = qnb0;  vn = vn0;  pos = pos0; % 姿态、速度和位置初始化
phi = [10; 20; 120]*arcmin;% 失准角(在陀螺好的情况下影响不大，主要是陀螺精度)
qnb = qaddphi(qnb, phi);  % 包含失准角的姿态四元数 

%IMU误差参数以及UWB误差参数
eb = [5.1;5.1;5.1]*dph; %陀螺常值零偏(°/h)，实验室陀螺实际值
%eb = [1;1;1]*dph;
web = [0.26;0.26;0.26]*dpsh; %角度随机游走系数(°/sqrt(h))，实验室陀螺实际值 
db = [80;90;100]*ug; % 加速度计常值偏值ug%wdb = [1;1;1]*ugpsHz;
wdb = [5;5;5]*ugpsHz; % 速度随机游走系数(1ug/sqrt(Hz)或=1*60*ug*[(m/s)/(sqrt(h)))
Qk = diag([web; wdb; zeros(9,1)])^2*nts;   % 系统误差协方差阵(半正定)，角度和速度随机游走
rk = [[0.15;0.15;0.2];[[0.025;0.025];0.25]]; % 观测值(IMU-UWB)误差:速度、位置误差（标准差）不需/Re
Rk = diag(rk)^2; %观测误差协方差阵(正定)    %注意预设误差大小将影响滤波效果，设大则重观测，设小重滤波值
P0 = diag([[1;1;10]*arcdeg; [0.1;0.1;0.1]; [[0.20;0.20];0.25]; [0.1;0.1;0.1]*dph; [100;100;100]*ug])^2;
%初始误差协方差阵，角度(°),速度(m/s),位置(m),陀螺仪(°/h),加速度计(10e-6*g=9.7803267714*10e-6m/s^2)
%X=[phi,dv,dP,epsilon,nabla];%位置误差向量顺序应该要看系统矩阵的定义，此处是NEU（纬度，经度和高度）顺序
%X=[phi_E,phi_N,phi_U,dv_E,dv_N,dv_U,dP_N,dP_E,dP_U,epsilon_x,epsilon_y,epsilon_z,nabla_x,nabla_y,nabla_z]
Hk = [zeros(6,3),eye(6),zeros(6)];   %观测矩阵，6*15(只用位置匹配时，应该是3*15)
%滤波器初始化
[kf_m, kf_n] = size(Hk);%计算维度
kf = kfinit_simple(Qk, Rk, P0, kf_n, Hk); % kf滤波器初始化，KF的各种量
% kf.Qk = Qk; %系统误差方差阵(半正定)
% kf.Rk = Rk; %观测误差协方差阵(正定)
% kf.Pk = P0; %初始误差
% kf.Xk = zeros(kf.n,1);  %状态估计量
% kf.Phikk_1 =zeros(15);  %状态转移矩阵（初始设为0）
% kf.Hk = Hk;             %量测矩阵
% kf.Tauk = eye(kf.n);  %系统误差的系数矩阵

%uwb位置（需要转换到导航坐标系）
[uwb(:,4:5),uwb(:,1:2)]=LSM2D(5);%东和北向位置和速度，LSM解算的UWB值位置(由于解算的速度误差可能比较大，可以不解算速度)
uwb(:,4)=uwb(:,4)-uwb(1,4);%转为相对第一点的偏移，若已经做好坐标系转换则不需要
uwb(:,5)=uwb(:,5)-uwb(1,5);
uwb(:,6)=ones(length(uwb(:,5)),1)*pos0(3)+rk(6).*randn(length(uwb(:,5)),1);%高度位置固定不变，只加误差；可根据设计的轨迹修改
uwb(:,3)=ones(length(uwb(:,1)),1)*vn0(3)+rk(3).*randn(length(uwb(:,5)),1); %高度速度固定不变，只加误差

lenimu = fix(length(wm(:,1))); %时长由陀螺仪(IMU)数据长度而定 length(wm(:,1))=11426
avp = zeros(lenimu, 10);xkpk = zeros(lenimu, 2*kf.n+1); kk = 1;  t = 0; % 记录结果存储空间预分配
for k=1:nn:lenimu 
    t = t + nts;
    %根据imu参数加噪声
    %[wm1, vm1] = imuadderr(wm(k:k+nn-1,:), vm(k:k+nn-1,:), eb, web, db, wdb, ts);%imu加噪声
    size_wm1 = size(wm(k:k+nn-1,:),1); 
    wm1 = wm(k:k+nn-1,:) + [ ts*eb(1) + sqrt(ts)*web(1)*randn(size_wm1,1), ...
                ts*eb(2) + sqrt(ts)*web(2)*randn(size_wm1,1), ...
                ts*eb(3) + sqrt(ts)*web(3)*randn(size_wm1,1) ];
    vm1 = vm(k:k+nn-1,:) + [ ts*db(1) + sqrt(ts)*wdb(1)*randn(size_wm1,1), ...
                ts*db(2) + sqrt(ts)*wdb(2)*randn(size_wm1,1), ...
                ts*db(3) + sqrt(ts)*wdb(3)*randn(size_wm1,1) ];

    %[qnb, vn, pos, eth] = insupdate(qnb, vn, pos, wm1, vm1, ts);%惯导更新
    %[qnb, vn, pos, eth] = insupdate_simple(qnb, vn, pos, wm1, vm1, ts);%惯导更新，这步或可使用kalman滤波更新
    [qnb, vn, pos] = insupdate_simple1(qnb, vn, pos, wm1, vm1, ts);%最终采用的简化惯导更新
    %kf.Phikk_1 = eye(15) + kfft15(eth, q2mat(qnb), sum(vm1,1)'/nts)*nts;%系统离散化P102，P155,P80
    kf.Phikk_1 = eye(15) + kfft15_simple(q2mat(qnb), sum(vm1,1)'/nts)*nts;%简化版系统离散化P102，P155,P80
    kf = kfupdate(kf);%时间更新,kf阶数为15阶
    if mod(t,0.1)<nts %UWB有数据，则进行Kalman更新（需要匹配数据20190426）
        uwb(round((t)/0.1+1),4:6)=pos0'+uwb(round((t)/0.1+1),4:6);
        kf = kfupdate(kf, [vn;pos]-uwb(round((t)/0.1+1),:)', 'M');  % 状态量更新
        vn(3) = vn(3) - kf.Xk(6);  kf.Xk(6) = 0; % 天向速度反馈
        vn(1:2) = vn(1:2) - kf.Xk(4:5);  kf.Xk(4:5) = 0; % 水平速度反馈使轨迹平滑，滤波后的值较准
        pos = pos - kf.Xk(7:9);  kf.Xk(7:9) = 0; % 位置反馈使误差不发散，滤波后的值较准
    end
    %avp(kk,:) = [qq2phi(qnb,qnb0); vn; pos; t]';%姿态角，速度和位置解算值（直接解算并反馈后的）。共10维
    qerr = qmul(qnb, [qnb0(1); -qnb0(2:4)]);%失准角误差=四元数计算值-四元数真值
    %phi = q2rv(qerr);%误差四元数转姿态角（等效旋转矢量）
    if qerr(1)<0,  qerr = -qerr;  end
    nmhalf = acos(qerr(1));  % 等效旋转矢量模值的一半
    if nmhalf>1e-20,   qerr_b = 2*nmhalf/sin(nmhalf);
    else             qerr_b = 2;    end
    phi = qerr_b*qerr(2:4);%转动的姿态角   % q = [ cos(|rv|/2); sin(|rv|/2)/|rv|*rv ];
    avp(kk,:) = [phi; vn; pos; t]';
    xkpk(kk,:) = [kf.Xk; diag(kf.Pk); t];   %状态量估计值以及估计误差协方差，加上时间。共31维
              %状态量估计值(15维)，估计误差(15维)
    kk = kk+1;
    if mod(t,100)<nts,  disp(fix(t));  end  % 显示进度
end
avp(kk:end,:) = [];  xkpk(kk:end,:) = [];  tt = avp(:,end);
dpos = [(avp(:,7)-avp(1,7)),(avp(:,8)-avp(1,8)),(avp(:,9)-avp(1,9))];%计算位置变化值
% 状态真值与估计效果对比图
mysubplot(321, tt, [avp(:,1:2),xkpk(:,1:2)]/arcmin, '\phi_E,\phi_N / \prime');%水平姿态角以及估计值
mysubplot(322, tt, [avp(:,3),xkpk(:,3)]/arcmin, '\phi_U / \prime');           %天向姿态角以及估计值
mysubplot(323, tt, [avp(:,4:6),xkpk(:,4:6)], '\deltav ^n / m/s');             %速度解算值及速度估计值
mysubplot(324, tt, [dpos,xkpk(:,7:9)],'\DeltaP / m');        %位置及位置估计值
mysubplot(325, tt, xkpk(:,10:12)/dph, '\epsilon / \circ/h'); %陀螺误差估计值
mysubplot(326, tt, xkpk(:,13:15)/ug, '\nabla / ug');         %加速度计误差估计值
subtitle('Positioning error and performance');%定位误差与效果

% 方差收敛图
pk = sqrt(xkpk(:,16:end-1));
mysubplot(321, tt, pk(:,1:2)/arcmin, '\phi_E,\phi_N / \prime');%水平姿态角误差
mysubplot(322, tt, pk(:,3)/arcmin, '\phi_U / \prime');         %天向姿态角误差
mysubplot(323, tt, pk(:,4:6), '\deltav ^n / m/s');             %速度误差
mysubplot(324, tt, pk(:,7:9), '\DeltaP / m');             %位置误差
mysubplot(325, tt, pk(:,10:12)/dph, '\epsilon / \circ/h');%陀螺估计误差
mysubplot(326, tt, pk(:,13:15)/ug, '\nabla / ug');        %加速度计估计误差
subtitle('Variance');%方差，协方差covariance

%弧度转角度再转距离值（km），转成以起始点的经纬度为原点（0,0）再显示
%dpos = lolahe2dis(avp(:,7:9));%计算位移的相对变化量
figure;                       %采用较为准确的地球半径参数
plot(dpos(:,1), dpos(:,2),'linewidth', 2 );
xlabel('East / m'); 
ylabel('North / m');
title('Relative displacements');%相对位移
hold on, plot(dpos(1,1),dpos(1,2) , 'ro');%标记起点
hold on, plot(dpos(end,1), dpos(end,2), 'g+');%标记终点

figure;     %位置收敛图
subplot 211;
plot(tt, dpos,'linewidth', 2);%解算的位置，转换成m为单位，地球半径Re采用球形模型
title('Relative displacements(feedback)');%相对位移
set(legend('E','N','U'),'Orientation','horizon','location','NorthEast','box','off');
xlabel('t/s');                                
ylabel('\DeltaP / m');
hold on
subplot 212;
plot(tt, [xkpk(:,7),xkpk(:,8),pk(:,9)],'linewidth', 2);%Kalman解算的状态值
xlabel('t/s');                                                         
ylabel('\DeltaXP / m');                                                                         
set(legend('E','N','U'),'Orientation','horizon','location','NorthEast','box','off');

%定位误差
Xr=Route2D_100Hz();%轨迹真值[x,y]=[E,N]，100Hz,共11427点
Xr=[Xr(:,1)-Xr(1,1),Xr(:,2)-Xr(1,2)];%转为相对第一点的变化量
Xr=Xr(2:end,:); % 去掉初始值,以适应dpos
ErrordposN=dpos(:,1)-Xr(:,1);%x方向误差
ErrordposE=dpos(:,2)-Xr(:,2);%y方向误差
Errordpos=sqrt(ErrordposN.^2+ErrordposE.^2);%欧式距离误差
RMSE_dpos=sqrt(sum((Errordpos).^2)/length(Errordpos));%求均方根误差
MeanVdpos=mean(ErrordposE);%按列求误差均值
StdVdpos=std(Errordpos,0);%求标准差值
Maxerrordpos=max(abs(Errordpos));%最大误差
%画总误差图
figure('units','normalized','position',[0.1,0.1,0.55,0.6]); 
plot(tt,Errordpos,'-');
ylabel('Error/m');
xlabel('t/s');
text(10,0.035,['Error: Mean= ',num2str(MeanVdpos),' Std= ',num2str(StdVdpos),...
    ' Max=',num2str(Maxerrordpos),' RMSE=',num2str(RMSE_dpos),' (m)']);
%gtext(['均方根误差= ',num2str(RMSE_dpos)]);
title('Integrated positioning by Kalman Method');
