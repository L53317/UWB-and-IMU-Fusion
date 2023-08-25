function kf = kfupdate(kf, Zk, TimeMeasBoth)
%
    if nargin==1,      TimeMeasBoth = 'T';%只有状态值，只进行时间更新
    elseif nargin==2,  TimeMeasBoth = 'B';%有状态值,量测且无标识,完整KF    
    end
    if TimeMeasBoth=='T' || TimeMeasBoth=='B'     % 时间更新
        kf.Xkk_1 = kf.Phikk_1*kf.Xk; %状态量的一步预测
        kf.Pkk_1 = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Tauk*kf.Qk*kf.Tauk';%一步预测误差协方差阵
    else % TimeMeasBoth=='M'   %无量测更新
        kf.Xkk_1 = kf.Xk;      %应该是kf.Xk = kf.Xkk_1;才对吧，不过不会执行到此句
        kf.Pkk_1 = kf.Pk;
    end
    if TimeMeasBoth=='M' || TimeMeasBoth=='B'     % 量测更新
        kf.PXZkk_1 = kf.Pkk_1*kf.Hk';                %中间量，状态一步预测与量测一步预测的协方差矩阵
        kf.PZkk_1 = kf.Hk*kf.PXZkk_1 + kf.Rk;        %中间量，量测一步预测的均方误差矩阵
        kf.Kk = kf.PXZkk_1/kf.PZkk_1;                %增益矩阵
        kf.Xk = kf.Xkk_1 + kf.Kk*(Zk-kf.Hk*kf.Xkk_1);%状态量更新
        kf.Pk = kf.Pkk_1 - kf.Kk*kf.PZkk_1*kf.Kk';   % 估计误差协方差矩阵更新
                                    %kf.Pk = kf.Pkk_1 - kf.Kk*kf.Hk*kf.Pkk_1;%不应该是这条？
    else % TimeMeasBoth=='T'    % 无时间更新
        kf.Xk = kf.Xkk_1;
        kf.Pk = kf.Pkk_1;
    end
    kf.Pk = (kf.Pk+kf.Pk')/2;   % P阵对称化
