function [wm, vm] = av2imu(att, vn, pos, ts)  % 由姿态、速度和位置生成惯性传感器信息
    wm0 = zeros(3,1); vm0 = wm0; I33 = eye(3);
    wm = att(2:end,:); vm = wm; %分配存储空间
    for k=2:length(att)
        eth = earth((pos(k-1,:)+pos(k,:))'/2, (vn(k-1,:)+vn(k,:))'/2);
        qbb = qmul(qmul(qconj(a2qua(att(k-1,:))),rv2q(eth.wnin*ts)),a2qua(att(k,:)));%两个数据间的四元数
        phim = q2rv(qbb);                 %计算两个数据间的旋转矢量
        wm1 = (I33+askew(1/12*wm0))\phim; %求陀螺仪仿真数据
        dvnsf = vn(k,:)'-vn(k-1,:)'-eth.gcc*ts;
        Cnb0 = a2mat(att(k-1,:)');
        vm1 = (I33+1/2*askew(1/6*wm0+wm1))\...
              (Cnb0'*(I33+askew(eth.wnin*ts/2))*dvnsf-1/12*cross(vm0,wm1));
        wm(k-1,:) = wm1';  vm(k-1,:) = vm1;  %存储陀螺仪和加速度计仿真数据
        wm0 = wm1; vm0 = vm1;
    end