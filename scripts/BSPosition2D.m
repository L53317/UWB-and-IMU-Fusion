function [BS]=BSPosition2D()
%基站的位置,最大值为6
%输入：无
%输出：基站位置坐标 BS：第一行为x，第二行为y；
% --- liu.20181212 ---

% 各个基站的位置，单位m  %为了减少Chan算法的奇值，使前四个基站不对称，
% BS = [0, 15, 15, 0,  7.5, 7.5,; %对称四个基站时，出现奇值，
%       0, 0,  10, 10, 5,   9];  %是否出现奇值主要看几何精度因子（GDOP）图
BS = [0, 15, 15, 7.5, 0,  7.5;   %都在轨迹外精度稍高（四基站时能看出）
      0, 0,  10, 9,   10, 5]; %四个基站时不能对称放置，否则易出现奇值  
                      %设置第五个基站位置未知，最后一个一般用不到
% figure  %画基站位置
% plot(BS(1,:),BS(2,:),'rd');%画基站图
% title('基站位置');
% xlabel('x方向/m');
% ylabel('y方向/m');
% legend('基站');