function [att, vn, pos] = trjprofile(att0, vn0, pos0, wat, ts,num)  % 生成轨迹姿态、速度和位置参数
    if nargin <6 %原来没有这个语句，增加该句是为了imu数据输出时延低，在短距离（几米）导航需考虑
        num=1;%说明需要进行滤波，num置为1。按原来的程序执行20190427
    end       %不进行滤波传入的参数应num=0
    len = fix(sum(wat(:,5))/ts);                       %计算数据总长度
    att = zeros(len, 3); vn = att; pos = att;  kk=1;   %分配存储空间
    att(1,:) = att0'; vn(1,:) = vn0'; pos(1,:) = pos0';%存储初值，转为行向量
    vb = a2mat(att0)'*vn0; vby = vb(2);   % 求纵向速度
    b = fir1(20, 0.01, 'low');%生成低通滤波器的参数（单位脉冲响应的序列）
    %N=20为滤波器阶数，Wn=0.01Hz为截止频率，输出b是N+1维行向量的滤波器系数。会造成时延。
    if 0==num %原来没有这个语句，增加该句是为了imu数据输出时延低，在短距离（几米）导航需考虑
        b =1;%此句替代前一句则不进行低通滤波(20190427.liu)
    end      %不进行滤波传入的参数应num=0
    b = b/sum(b); x = repmat([att0;vby]',length(b),1); % 低通滤波器，数据扩展
    %N=20为滤波器阶数，Wn=0.01Hz为截止频率，输出b是N+1维行向量的滤波器系数。会造成时延。N=1时不进行滤波
    for m=1:size(wat,1);
        watk = wat(m,:);
        for tk=ts:ts:(watk(5)+ts/10) %对相邻时间段的相应的数据进行低通滤波（姿态和速度低通滤波，使轨迹光滑）
            att0 = att0 + watk(1:3)'*ts;%计算姿态角，此处直接相加（只变一个角度，定轴转动所以可以如此）   
            vby = vby + watk(4)*ts;
            x = [x(2:end,:); [att0;vby]']; y = b*x;  % 低通滤波。会使轨迹稍微平滑
            att(kk+1,:) = y(1:3);                    % 存储低通滤波值
            vn(kk+1,:) = (a2mat(att(kk+1,:)')*[0;y(4);0])';  vn01 = (vn(kk,:)+vn(kk+1,:))/2;% 计算速度
            eth = earth(pos(kk,:)',vn01'); % 计算当地地球半径参数
            pos(kk+1,:) = pos(kk,:) + [vn01(2)/eth.RMh;vn01(1)/eth.clRNh;vn01(3)]'*ts;%速度积分，得到位置
            kk = kk+1;
        end
    end
    att(kk:end,:) = []; vn(kk:end,:) = []; pos(kk:end,:) = [];