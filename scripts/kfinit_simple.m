function kf = kfinit_simple(Qk, Rk, P0, Xn, Hk)
% Kalman filter initializes for structure array 'kf', this precedure 
% usually includs the setting of structure fields: Qk, Rk, Pxk, Hk.
%
% See also kfupdate.

kf.Qk = Qk; %系统误差方差阵(半正定)
kf.Rk = Rk; %观测误差协方差阵(正定)
kf.Pk = P0; %初始误差
kf.n = Xn;  %可以根据HK求出
[kf.m, ~] = size(Hk);
kf.Xk = zeros(kf.n,1);  %状态估计量kf.n
kf.Phikk_1 =zeros(kf.n);  %状态转移矩阵（初始设为0）,
kf.Hk = Hk;             %量测矩阵
kf.Tauk = eye(kf.n);  %系统误差的系数矩阵
end
