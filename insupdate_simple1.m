function [qnb, vn, pos] = insupdate_simple1(qnb, vn, pos, wm, vm, ts)  % 捷联惯导数值更新算法
%原型function [qnb, vn, pos, eth] = insupdate_simple(qnb, vn, pos, wm, vm, ts)
    %最终版的惯导更新程序
    %此程序中保留了qmul(),rv2q(),askew()等基本函数，去除地球参数
    nn = size(wm,1);  nts = nn*ts;
    %[phim, dvbm] = cnscl(wm, vm);  % 圆锥误差/划船误差补偿
    %               %（不补偿，10ms，1°圆锥半角->1.03°/h;划桨效应引起的速度误差量级未知）
    phim = sum(wm,1)';  dvbm = sum(vm,1)';%对各列求和(得到关于wm，vm的行向量再转置为列，得到角和速度增量)
    %eth = earth(pos, vn);  % 地球相关参数计算，不需要了
    
    %四元数法计算速度增量dvbm在导航系下的投影dvnm 
    dvnm_qo = qmul(qmul(qnb,[0;dvbm]),[qnb(1); -qnb(2:4)]); %
    dvnm=dvnm_qo(2:4,1);%三维矢量的四元数转换（实现速度增量的坐标系转换到n系）    
    %vn1 = vn + dvnm + eth.gcc*nts;  % 速度更新（由于导航系变化极小，可替换上式）
    vn1 = vn + dvnm+[0;0;-9.7803267]*nts;%最终可采用该式代替速度更新,g0=9.7803267m/s^2
    vn = (vn+vn1)/2;%采用该段时间的平均速度作为最终速度
    
    pos = pos + [vn(1);vn(2);vn(3)]*nts;  vn = vn1;  % 位置更新（改回东北天-xyz）
    % 姿态更新
    %qnb = qmul(rv2q(-eth.wnin*nts), qmul(qnb, rv2q(phim)));  % 姿态更新
    pos0=[34.057*pi/180;118.786*pi/180;0];%当地导航系所在的地理位置
    wie=7.292e-05;  %地球自转角速度
    wnie = wie*[0; cos(pos0(1)); sin(pos0(1))];%计算地球自转对姿态的影响（补偿后会好）
    qnb = qmul(rv2q(-wnie*nts), qmul(qnb, rv2q(phim)));  % 姿态更新
    %qnb =  qmul(qnb, rv2q(phim));%可以采用此代替上式，但是kf位置误差rk不能取太小，否则影响效果
    
    %四元数归一化
    %qnb = qnormlz(qnb);
    if qnb'*qnb<1e-6,  qnb = [1; 0; 0; 0];%判断是否正确
    else
        qnb = qnb/sqrt(qnb'*qnb); % 归一化
    end