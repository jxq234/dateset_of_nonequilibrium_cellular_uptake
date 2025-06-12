clc
clear

xc = 20;
xmesh = [0:1:19, xc, xc, 21:1:40];
yinit = [0; 0; 0; 10; 10; 0];
sol = bvpinit(xmesh,yinit);

load('sol_initial.mat');



% 膜物理性质设置
kappa = 19.9;% kBT 弯曲刚度
a = 13.1;% sigma 球半径
w = 0.37;% kBT/sigma^2 黏附能
gamma = 2;% kBT/sigma 线张力
sigma = 0.26; % kBT/sigma^2 表面张力
F = 0; % kBT/sigma 活性力

% 对膜物理性质归一化
gamma = gamma*a/kappa;
sigma = sigma*a^2/kappa;
w = 2*w*10^2/kappa;
F = F*a/(pi*kappa);

alpha =  pi*9/16;
z_wrap = 1-cos(alpha);

p1 = [a F];
p2 = [alpha a gamma sigma F];

% 对自由膜部分求解
%sol = bvp5c(@(x,y,r) f(x,y,r,p1), @(YL,YR) bc(YL,YR,p2), combined_sol);
sol = bvp5c(@(x,y,r) f(x,y,r,p1), @(YL,YR) bc(YL,YR,p2), sol);

% 保存各步形状迭代求解
%save('sol_initial.mat', 'sol');
% 保存最终形状用于绘图
save('outcome.mat', 'sol');

% 包覆部分归一化能量和面积
A_wrap = 2*pi*a^2*z_wrap;
E_wrap = (-w+4+2*sigma)*z_wrap-F*z_wrap/a;
fprintf('A_wrap: %4.3f E_wrap: %4.3f \n',A_wrap,E_wrap);

% 输出无归一化表面张力
tension_raft = sol.y(6,1)*kappa/a^2;
tension_out = sol.y(6,length(sol.x))*kappa/a^2;
formatSpec = 'tension_raft: %4.4f tension_out: %4.4f \n';
fprintf(formatSpec,tension_raft, tension_out);


% 输出脂阀面积
A_raft = 0;
A_out = 0;
for si=1:length(sol.x)-1
    if sol.x(si) <= xc
        ds = sol.x(si+1) - sol.x(si);
        r = (sol.y(4,si) + sol.y(4,si+1))/2;
        A_raft = A_raft + 2*pi*r*ds;
    end
    if sol.x(si) > xc
        ds = sol.x(si+1) - sol.x(si);
        r = (sol.y(4,si) + sol.y(4,si+1))/2;
        A_out = A_out + 2*pi*r*ds;
    end
end
fprintf('A_raft: %4.3f A_out: %4.3f \n',A_raft+A_wrap,A_out);

% 输出各项能量和交界处半径
E_bend = 0;
E_ten = 0;
E_line_ten = 0;
W = 0;
s_raft = 0;
for si=1:length(sol.x)-1
    ds = sol.x(si+1) - sol.x(si);
    
    psi = (sol.y(1,si) + sol.y(1,si+1))/2;
    dpsi = (sol.y(2,si) + sol.y(2,si+1))/2;
    r = (sol.y(4,si) + sol.y(4,si+1))/2;
    sigma = (sol.y(6,si) + sol.y(6,si+1))/2;
    
    E_bend = E_bend + r*(dpsi+sin(psi)/r)^2*ds;
    E_ten = E_ten + 2*r*sigma/a^2*ds;
    W = W+F*sin(psi)/a*ds;
    
    if sol.x(si) == xc
        E_line_ten = 2*gamma/a*r; 
        boundary_psi = sol.y(1,si)/pi*180;
        r_sc = sol.y(4,si);
        s_raft = si;
    end
end
fprintf('boundary psi: %4.3f \n',boundary_psi);
fprintf('R(sc): %4.3f \n',r_sc);


E_total = E_bend + E_ten + E_line_ten-W;
formatSpec = 'E_bend: %4.3f E_ten: %4.3f E_line_ten: %4.3f E_free_total: %4.3f\n';
fprintf(formatSpec,E_bend, E_ten, E_line_ten, E_total);
fprintf('W: %4.3f \n',W);

% 验证4个边界哈密顿量的值
H_1L = sol.y(4,1)*sol.y(2,1)^2 - sin(sol.y(1,1))^2/sol.y(4,1) - 2/a^2*sol.y(4,1)*sol.y(6,1) + sol.y(3,1)*cos(sol.y(1,1));
s = length(sol.x);
H_2R = sol.y(4,s)*sol.y(2,s)^2 - sin(sol.y(1,s))^2/sol.y(4,s) - 2/a^2*sol.y(4,s)*sol.y(6,s) + sol.y(3,s)*cos(sol.y(1,s));
fprintf('H_1L: %4.3f H_2R: %4.3f \n',H_1L,H_2R);

A_all = A_wrap + A_raft;
E_all = E_wrap + E_total;
fprintf('E_all: %4.3f \n',E_all);
formatSpec = '%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.4f\t%4.4f\t%4.4f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t';
%fprintf(formatSpec,E_bend, E_ten, E_line_ten, E_total,boundary_psi,tension_raft,tension_out,A_raft,A_out,H_2R,r_sc,A_wrap,E_wrap,A_all,E_all);


% 对细胞膜形状作图
if 1

%plot(sol.y(4,:),sol.y(5,:),'--o')
plot(sol.y(4,1:s_raft),sol.y(5,1:s_raft),'linewidth',2,'color','r')
hold on
plot(sol.y(4,s_raft:length(sol.x)),sol.y(5,s_raft:length(sol.x)),'linewidth',2,'color','black')
hold on
plot(-sol.y(4,1:s_raft),sol.y(5,1:s_raft),'linewidth',2,'color','r')
hold on
plot(-sol.y(4,s_raft:length(sol.x)),sol.y(5,s_raft:length(sol.x)),'linewidth',2,'color','black')

    
axis([-35,35,-10,60])
end

function dydx = f(x,y,region,p) % equations being solved
a = p(1); % sigma 球半径
F = p(2); % 归一化活性力

dydx = zeros(6,1);

% y(1)---psi
% y(2)---dpsi/ds
% y(3)---lamdar
% y(4)---r
% y(5)---h
% y(6)---sigma

dydx(1) = y(2); 
%dydx(2) = sin(y(1))*cos(y(1))/y(4)^2  - y(2)*cos(y(1))/y(4) + y(3)*sin(y(1))/(2*y(4));
dydx(2) = sin(y(1))*cos(y(1))/y(4)^2  - y(2)*cos(y(1))/y(4) + y(3)*sin(y(1))/(2*y(4))-F*cos(y(1))/(2*a*y(4));
dydx(3) = y(2)^2 - sin(y(1))^2/y(4)^2 + 2*y(6)/a^2;
dydx(4) = cos(y(1));
dydx(5) = sin(y(1));
dydx(6) = 0;

end

function res = bc(YL,YR,p)

alpha = p(1); % 初始角
a = p(2);% sigma 球半径
gamma = p(3) ;% 归一化线张力
sigma = p(4) ;% 归一化表面张力
F = p(5);% 归一化活性力

res = [YL(1,1) - pi/2                % psi(s1) = alpha
       YL(4,1) - 6.79         % r(s1) = a*sin(alpha)
       YL(5,1)         % h(s1) = -a*cos(alpha)
       %YR(3,1)
       YR(6,1)-sigma

       
       YR(1,1) - YL(1,2)              % psi(sc-epsilon) = psi(s1+epsilon)
       YR(2,1) - YL(2,2)              % dpsi/ds(sc-epsilon) = dpsi/ds(s1+epsilon)
       YL(3,2) - YR(3,1) - 2*gamma/a  % dpsi/ds(sc-epsilon) = dpsi/ds(s1+epsilon)
       YR(4,1) - YL(4,2)              % r(sc-epsilon) = r(s1+epsilon)
       %YR(4,1) - 20
       YR(5,1) - YL(5,2)              % h(sc-epsilon) = h(s1+epsilon)
       (YL(6,2) - YR(6,1))*YR(4,1)-gamma*a*cos(YR(1,1)) %哈密顿量边界处连续
       
       YR(1,end) - pi/2               % psi(s2) = 0
       YR(3,end)
       %YR(6,end)-0.1
       
       %YR(2,end)
       %YR(6,end)-0.01                % sigma(s2) = sigma
       ];                    
end