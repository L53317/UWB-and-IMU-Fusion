function [Xp,v,ErrorX] = LSM2D(BSN)
%function [X ErrorX] = myLSM2(BSN, BS, R)
%最小二乘法定位，2维
%输入参数：节点数BSN，BS，R
%       - BSN为基站数目，3-6
%       - BS为(2,BSN)矩阵，为各个BS的坐标 x和y
%       - R为(BSN) 向量，为论文中的 r_i，即第 1，2,3,...BSN 个基站到MS的距离
%输出参数X:MS坐标(x,y),ErrorX

if nargin < 1 %判断是否传入BSN，用于日常调试本函数
    % 基站数目
    BSN = 5;
%else
end

BS = BSPosition2D(); %基站位置获取
BS = BS(:,1:BSN);  %取其中BSN个基站，(BS(1,i)为x坐标，(BS(2,i)为对应y坐标

% MS的实际位置
%按照0.2m/s,10Hz(采样时间0.1s)生成5*10m^2的数据,移动总路程(2*7.5+π*2.5)m，
%共计22.854m(1143+7=375+393+375+7)点
% y=0;%x=0-7.5
% y=2.5+sqrt(6.25-(x-7.5)^2);%x=7.5-10,这里应该有正负
% y=5;%x=0-7.5
%%%求解轨迹点
Xr=Route2D(); %生成轨迹矩阵X=[x y];共两列第一列为x第二列为y

%最小二乘法求最优解
for j=1:length(Xr(:,1))
    MS = [Xr(j,1),Xr(j,2)]; %原本的值
    %R0i是各个BS与MS的实际距离，无噪声
    for i = 1: BSN  %注意一般MATLAB都默认为行向量
        R0(i) = sqrt((BS(1,i) - MS(1))^2 + (BS(2,i) - MS(2))^2); 
    end
    Noise = 0.025; %单位m，测距标准差是0.025m，最大偏差值为0.1m
    for i = 1: BSN
        R(i) = R0(i) + Noise * randn(1); %需要每个噪声都不一样
    end    
    A=[2*(BS(1,1)-BS(1,BSN)) 2*(BS(2,1)-BS(2,BSN)); %采用BSN个基站的数据
       2*(BS(1,2)-BS(1,BSN)) 2*(BS(2,2)-BS(2,BSN))];
    B=[(BS(1,1)^2-BS(1,BSN)^2+BS(2,1)^2-BS(2,BSN)^2+R(BSN)^2-R(1)^2);
       (BS(1,2)^2-BS(1,BSN)^2+BS(2,2)^2-BS(2,BSN)^2+R(BSN)^2-R(2)^2)];
   if  BSN>3
   for i=(4-1):(BSN-1)
       A=[A; %采用BSN个基站的数据
          2*(BS(1,i)-BS(1,BSN)) 2*(BS(2,i)-BS(2,BSN))];
       B=[B;
          (BS(1,i)^2-BS(1,BSN)^2+BS(2,i)^2-BS(2,BSN)^2+R(BSN)^2-R(i)^2)];
   end
   end
    XL(j,:)=inv(A'*A)*A'*B; %最小二乘法求最优解
end

%最小二乘法误差
ErrorXLx=XL(:,1)-Xr(:,1);%x方向误差
ErrorXLy=XL(:,2)-Xr(:,2);%y方向误差
ErrorXL=sqrt(ErrorXLx.^2+ErrorXLy.^2);%欧式距离误差
RMSE_XL=sqrt(sum((ErrorXL).^2)/length(ErrorXL));%求均方根误差(2范数)
MeanVL=mean(ErrorXL);%按列求误差均值
StdVL=std(ErrorXL,0);%求标准差值
MaxerrorL=max(abs(ErrorXL));%最大误差
figure  %画误差图
plot(ErrorXL,'-');%画最小二乘法误差图
title([num2str(BSN),' anchors LS method positioning error']);
xlabel('Point index');%子图x轴，作为整张图的x轴
ylabel('y/m'); %子图y轴定义,作为全图y轴
text(10,0.005,['Error: Mean= ',num2str(MeanVL),' Std= ',num2str(StdVL),...
    ' Max=',num2str(MaxerrorL),' RMSE= ',num2str(RMSE_XL),' (m)']);
% %画轨迹图
% figure
% plot(BS(1,:),BS(2,:),'rd');%画基站图
% title('定位结果对比');
% xlabel('x方向/m');
% ylabel('y方向/m');
% hold on
% plot(Xr(:,1),Xr(:,2),'r-');%画实际轨迹
% hold on
% plot(XL(:,1),XL(:,2),'k.');    %画最小二乘法解算轨迹
% legend('基站','真实轨迹','LS');
% hold off;

%UWB速度计算(可以采用三次样条插值spline，会有光滑的速度且二阶导(加速度)连续,但会多段)
%或者求差分diff，这样就不会有连续的速度（本来就是离散点）
v0=[0.2,0]; %初始速度设定为0.2m/s
vt=diff(XL)/0.1;%取前一个时刻与该时刻的位置差分与时间之比作为当前时刻的速度
                         %diff(X)=[X(2:n,:) - X(1:n-1,:)];
v=[v0;vt];      %如此解算的速度及其不准确，不应采用（达到了0.2m/s的误差及以上）
% %计算速度误差，真值是0.2m/s（欧式距离）
% % stdv=std(v);%=[0.3,0.3]以上，主要是由于UWB解算有极大误差导致的
% % meanv=mean(abs(v(end-375:end,:)));%直线段的速度每个方向都在0.2-0.3m/s之间:速度不对
% vr=diff(Xr)/0.1;%计算速度真值
% vr=[v0;vr];  
% Errorvx=v(:,1)-vr(:,1);%x方向速度误差
% Errorvy=v(:,2)-vr(:,2);%y方向速度误差
% Errorv=sqrt(Errorvx.^2+Errorvy.^2);%欧式距离意义的速度误差
% RMSE_v=sqrt(sum((Errorv).^2)/length(Errorv));%求均方根误差(2范数)=0.435
% StdVLv=std(Errorv,0);%求误差误差标准差=0.209
% StdVLvx=std(Errorvx,0);%求速度误差标准差x方向=0.258
% StdVLvy=std(Errorvy,0);%求速度误差标准差y方向=0.350

%输出UWB解算位置
%XL=[XL(:,1)-XL(1,1),XL(:,2)-XL(1,2)];%输出时以第一点为起始点，东向和北向位置增量
                                     %不应在此转换，在外面转换
if 1==nargout %只输出位置
    Xp=XL;  %保存输出值,采用最小二乘法的值
elseif 2==nargout %输出速度和位置
    Xp=XL;  
    v=v;      %不推荐输出速度
else
    Xp=XL;  %保存输出值,采用最小二乘法的值
    v=v;    %保存速度
    ErrorX=ErrorXL;  %保存输出值误差,采用最小二乘法的值
end
