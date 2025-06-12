clc
clear

N = 80;
ds = 0.5;
h = linspace(ds, ds, N);
s = linspace(ds, ds*N, N);

%% 给定初始形状 
C_s_s = linspace(0, 0, N);

phi = linspace(0, 0, N);
phi(end) = pi/2;
for i=1:N-1
    phi(end-i)=phi(end-i+1)-C_s_s(end-i+1)*ds;
end

r = linspace(0, 0, N);
r(1) = 10;
for i=1:N-1
    r(i+1)=r(i)+cos(phi(i))*ds;
end

z = linspace(0, 0, N);
for i=1:N-1
    z(i+1)=z(i)+sin(phi(i))*ds;
end

C_phi_phi = sin(phi)./r;
T_phi_s_phi = cos(phi)./r;

rou = linspace(0.99, 0.99, N);


Lb=20;
f = exp(-(s-(Lb)).^2/2);
df = zeros(1, N);
for i = 2:N-1
    df(i) = (f(i+1) - f(i-1)) / (2*h(i));
end
df(1) = (f(2) - f(1)) / (h(1));
df(N) = (f(N) - f(N-1)) / (h(N));


%% 主程序

kappa = 19.9;
miu = 100;
Ks = 25.9;
sigma = 1;
save_outcome = 1;

t_all = 100000;
dt = 0.005;

v_n = linspace(0, 0, N);
dv_n = linspace(0, 0, N);
d2v_n = linspace(0, 0, N);
v_s = linspace(0, 0, N);
dv_s = linspace(0, 0, N);

r_t = {r};
z_t = {z};
phi_t = {phi};
C_s_s_t = {C_s_s};
C_phi_phi_t = {C_phi_phi};
T_phi_s_phi_t = {T_phi_s_phi};
h_t = {h};
rou_t = {rou};


v_n_t = {v_n};
dv_n_t = {dv_n};
d2v_n_t = {d2v_n};
v_s_t = {v_s};
dv_s_t = {dv_s};

E_bend_t = linspace(0, 0, t_all);
E_ten_t = linspace(0, 0, t_all);
E_line_ten_t=linspace(0, 0, t_all);
W_t = linspace(0, 0, t_all);



for t = 1:t_all
    % 计算膜表面流
    [v_s,dv_s,v_n,dv_n,d2v_n] = v_one_step(kappa,miu,Ks,sigma,r,z,phi,h,rou,f,df,C_s_s,C_phi_phi,T_phi_s_phi,t);
    
    % 调用函数motion
    [r, phi, z, rou,h,C_s_s,C_phi_phi,T_phi_s_phi] = motion(r,z,phi, rou,h,C_s_s,C_phi_phi,T_phi_s_phi, v_n,dv_n,d2v_n,v_s,dv_s, dt);
    
    r_t{t} = r; 
    z_t{t} = z; 
    phi_t{t} = phi;
    h_t{t} = h;
    rou_t{t} = rou;
    C_s_s_t{t} = C_s_s;
    C_phi_phi_t{t} = C_phi_phi;
    T_phi_s_phi_t{t} = T_phi_s_phi;
    
    v_n_t{t} = v_n; 
    dv_n_t{t} = dv_n; 
    v_s_t{t} = v_s;
    dv_s_t{t} = dv_s; 
    
    if mod(t, 100) == 0
        fprintf('time: %4.3f \n',t);
    end
    
    
    if mod(t, 5000) == 0
        plot(r,z)
        %plot(rou)
        %plot(v_s)
        %plot(h)
        %plot(v_n)
        %plot(C_s_s)
        hold on
        
        fprintf('r_Lb: %4.3f \n',r(2*Lb+1));
        fprintf('r1: %4.3f \n',r(1));
        
        fprintf('rou_1: %4.3f \n',rou(1));
        fprintf('rou_2: %4.3f \n',rou(N));
        
        if save_outcome
        folderPath = './result';
        fileName = sprintf('r_t%d.mat', t);
        fullFilePath = fullfile(folderPath, fileName);
        save(fullFilePath, 'r');
        
        fileName = sprintf('z_t%d.mat', t);
        fullFilePath = fullfile(folderPath, fileName);
        save(fullFilePath, 'z');
        
        fileName = sprintf('rou_t%d.mat', t);
        fullFilePath = fullfile(folderPath, fileName);
        save(fullFilePath, 'rou');
        
        fileName = sprintf('vs_t%d.mat', t);
        fullFilePath = fullfile(folderPath, fileName);
        save(fullFilePath, 'v_s');
        end
        
    end
    
end
%% v_one_step
function [v_s,dv_s,v_n,dv_n,d2v_n] = v_one_step(kappa,miu,Ks,sigma,r,z,phi,h,rou,f,df,C_s_s,C_phi_phi,T_phi_s_phi,t)

    N = length(phi);
    
    drou = zeros(1, N);
    for i = 2:N-1
        drou(i) = (rou(i+1) - rou(i-1)) / (2*h(i));
    end
    drou(1) = (rou(2) - rou(1)) / (h(1));
    drou(N) = (rou(N) - rou(N-1)) / (h(N));

    dC_s_s = zeros(1, N);
    for i = 3:N-2
        dC_s_s(i) = (C_s_s(i-2) - 8*C_s_s(i-1) + 8*C_s_s(i+1) - C_s_s(i+2)) / (12*h(i));
    end
    dC_s_s(2) = (C_s_s(3) - C_s_s(1)) / (2*h(2));
    dC_s_s(1) = dC_s_s(2);
    dC_s_s(N-1) = (C_s_s(N) - C_s_s(N-2)) / (2*h(N-1));
    dC_s_s(N) = dC_s_s(N-1);
    
    
    d2C_s_s = zeros(1, N);
    for i = 3:N-2
        d2C_s_s(i) = (dC_s_s(i-2) - 8*dC_s_s(i-1) + 8*dC_s_s(i+1) - dC_s_s(i+2)) / (12*h(i));
    end
    d2C_s_s(2) = (dC_s_s(3) - dC_s_s(1)) / (2*h(2));
    d2C_s_s(1) = d2C_s_s(2);
    d2C_s_s(N-1) = (dC_s_s(N) - dC_s_s(N-2)) / (2*h(N-1));
    d2C_s_s(N) = d2C_s_s(N-1);
    

    d3C_s_s = zeros(1, N);
    for i = 3:N-2
        d3C_s_s(i) = (d2C_s_s(i-2) - 8*d2C_s_s(i-1) + 8*d2C_s_s(i+1) - d2C_s_s(i+2)) / (12*h(i));
    end
    d3C_s_s(2) = (d2C_s_s(3) - d2C_s_s(1)) / (2*h(2));
    d3C_s_s(1)=d3C_s_s(2);
    d3C_s_s(N-1) = (d2C_s_s(N) - d2C_s_s(N-2)) / (2*h(N-1));
    d3C_s_s(N) = d3C_s_s(N-1);
    
    
    dC_phi_phi = zeros(1, N);
    for i = 3:N-2
        dC_phi_phi(i) = (C_phi_phi(i-2) - 8*C_phi_phi(i-1) + 8*C_phi_phi(i+1) - C_phi_phi(i+2)) / (12*h(i));
    end
    dC_phi_phi(2) = (C_phi_phi(3) - C_phi_phi(1)) / (2*h(2));
    dC_phi_phi(1)=dC_phi_phi(2);
    dC_phi_phi(N-1) = (C_phi_phi(N) - C_phi_phi(N-2)) / (2*h(N-1));
    dC_phi_phi(N)=dC_phi_phi(N-1);
    
    
    d2C_phi_phi = zeros(1, N);
    for i = 3:N-2
        d2C_phi_phi(i) = (dC_phi_phi(i-2) - 8*dC_phi_phi(i-1) + 8*dC_phi_phi(i+1) - dC_phi_phi(i+2)) / (12*h(i));
    end
    d2C_phi_phi(2) = (dC_phi_phi(3) - dC_phi_phi(1)) / (2*h(2));
    d2C_phi_phi(1)=d2C_phi_phi(2);
    d2C_phi_phi(N-1) = (dC_phi_phi(N) - dC_phi_phi(N-2)) / (2*h(N-1));
    d2C_phi_phi(N) = d2C_phi_phi(N-1);

    
    d3C_phi_phi = zeros(1, N);
    for i = 3:N-2
        d3C_phi_phi(i) = (d2C_phi_phi(i-2) - 8*d2C_phi_phi(i-1) + 8*d2C_phi_phi(i+1) - d2C_phi_phi(i+2)) / (12*h(i));
    end
    d3C_phi_phi(2) = (d2C_phi_phi(3) - d2C_phi_phi(1)) / (2*h(2));
    d3C_phi_phi(1)= d3C_phi_phi(2);
    d3C_phi_phi(N-1) = (d2C_phi_phi(N) - d2C_phi_phi(N-2)) / (2*h(N-1));
    d3C_phi_phi(N) = d3C_phi_phi(N-1);

    
    C_k_k = C_s_s + C_phi_phi;
    dC_k_k = dC_s_s + dC_phi_phi;
    d2C_k_k = d2C_s_s + d2C_phi_phi;
    d3C_k_k = d3C_s_s + d3C_phi_phi;

    % 构建 main_matrix
    A1 = C_phi_phi.^2./(C_s_s.^2+C_phi_phi.^2);

    A2 = -T_phi_s_phi+(C_s_s.*dC_s_s+T_phi_s_phi.*C_s_s.*C_phi_phi+C_s_s.*dC_k_k)./(C_s_s.^2+C_phi_phi.^2)-2*C_s_s.^2.*(C_s_s.*dC_s_s+C_phi_phi.*dC_phi_phi)./(C_s_s.^2+C_phi_phi.^2).^2;

    A31 = (C_s_s.*T_phi_s_phi.*dC_phi_phi+C_s_s.*C_phi_phi.*(-C_s_s.*C_phi_phi-T_phi_s_phi.^2)+T_phi_s_phi.*C_phi_phi.*dC_k_k)./(C_s_s.^2+C_phi_phi.^2);
    A32 = -2*(C_s_s.*dC_s_s+C_phi_phi.*dC_phi_phi).*T_phi_s_phi.*C_s_s.*C_phi_phi./(C_s_s.^2+C_phi_phi.^2).^2;
    A3 = A31+A32+T_phi_s_phi.*T_phi_s_phi;

    A41 = Ks/(4*miu)*(1-rou.^2).*(C_s_s.*dC_k_k+C_k_k.*dC_k_k)./(C_s_s.^2+C_phi_phi.^2);
    A42 = -Ks/(4*miu)*(1-rou.^2)*2.*(C_s_s.*dC_s_s+C_phi_phi.*dC_phi_phi).*C_k_k.*C_s_s./(C_s_s.^2+C_phi_phi.^2).^2;
    A43 = Ks/(2*miu).*rou.*drou.*(1-C_s_s.*C_k_k./(C_s_s.^2+C_phi_phi.^2));
    
    A44 = sigma/(2*miu)*df.*C_s_s.*C_k_k./(C_s_s.^2+C_phi_phi.^2);
    A45 = sigma/(2*miu)*f.*(C_s_s.*dC_k_k+C_k_k.*dC_k_k)./(C_s_s.^2+C_phi_phi.^2);
    A46 = -sigma/(2*miu)*f*2.*(C_s_s.*dC_s_s+C_phi_phi.*dC_phi_phi).*C_k_k.*C_s_s./(C_s_s.^2+C_phi_phi.^2).^2;
    A4 = A41+A42+A43+A44+A45+A46;

    A51 = -kappa/(4*miu)*(C_s_s.*(C_phi_phi-C_s_s).^2.*dC_k_k+2*C_k_k.*C_s_s.*(C_phi_phi-C_s_s).*(dC_phi_phi-dC_s_s))./(C_s_s.^2+C_phi_phi.^2);
    A52 = kappa/(2*miu)*C_k_k.*C_s_s.*(C_phi_phi-C_s_s).^2.*(C_s_s.*dC_s_s+C_phi_phi.*dC_phi_phi)./(C_s_s.^2+C_phi_phi.^2).^2;

    A53 = -kappa/(2*miu)*C_s_s.*(d3C_k_k+dC_k_k.*(-C_phi_phi.*C_s_s-T_phi_s_phi.^2)+T_phi_s_phi.*d2C_k_k)./(C_s_s.^2+C_phi_phi.^2);
    A54 = kappa/(miu)*C_s_s.*(C_s_s.*dC_s_s+C_phi_phi.*dC_phi_phi).*(d2C_k_k+T_phi_s_phi.*dC_k_k)./(C_s_s.^2+C_phi_phi.^2).^2;
    
    A55 = -kappa/(4*miu)*C_k_k.*dC_k_k.*(C_phi_phi-C_s_s).^2./(C_s_s.^2+C_phi_phi.^2);
    A56 = -kappa/(2*miu)*dC_k_k.*(d2C_k_k+T_phi_s_phi.*dC_k_k)./(C_s_s.^2+C_phi_phi.^2);
    A5 = A51+A52+A53+A54+A55+A56;

    % Initialize diagonals for the sparse matrix
    main_diag = -2 + A3./A1 .* h.^2;
    off_diag_upper = 1 - (A2./A1 .* h) / 2;
    off_diag_lower = 1 + (A2./A1 .* h) / 2;

    % Initialize right-hand side
    rhs = -(A4+A5)./A1 .* h.^2;

    % Boundary conditions
    main_diag(1) = 1; 
    main_diag(end) = 1;

    off_diag_lower(end)=0;
    off_diag_lower(end-1)=0;
    off_diag_upper(1)=0;
    off_diag_upper(2)=0;

    rhs(1) = 0; 
    rhs(end) = 0;

    % Create the sparse matrix
    A_matrix = spdiags([off_diag_lower', main_diag', off_diag_upper'], [-1 0 1], N, N);

    % Solve the system
    v_s = A_matrix \ rhs';
    v_s = v_s';
    
    dv_s = zeros(1, N);
    for i = 2:N-1
        dv_s(i) = (v_s(i+1) - v_s(i-1)) / (2*h(i));
    end
    %dv_s(1) = dv_s(2);
    %dv_s(N) = dv_s(N-1);
    dv_s(1) = (v_s(2) - v_s(1)) / (h(1));
    dv_s(N) = (v_s(N) - v_s(N-1)) / (h(N));
    
    B1 = -C_s_s./(C_s_s.^2+C_phi_phi.^2);
    B2 = -T_phi_s_phi.*C_phi_phi./(C_s_s.^2+C_phi_phi.^2);
    B3 = -Ks/(4*miu).*C_k_k./(C_s_s.^2+C_phi_phi.^2).*(1-rou.^2)-sigma/(2*miu).*f.*C_k_k./(C_s_s.^2+C_phi_phi.^2)+kappa/(4*miu).*C_k_k./(C_s_s.^2+C_phi_phi.^2).*(C_phi_phi-C_s_s).^2+kappa/(2*miu).*(d2C_k_k+T_phi_s_phi.*dC_k_k)./(C_s_s.^2+C_phi_phi.^2);
    

    
    v_n = B1.*dv_s+B2.*v_s+B3;

    dv_n = zeros(1, N);
    for i = 2:N-1
        dv_n(i) = (v_n(i+1) - v_n(i-1)) / (2*h(i));
    end
    dv_n(1) = (v_n(2) - v_n(1)) / (h(1));
    dv_n(N) = (v_n(N) - v_n(N-1)) / (h(N));
    
    d2v_n = zeros(1, N);
    for i = 2:N-1
        d2v_n(i) = (dv_n(i+1) - dv_n(i-1)) / (2*h(i));
    end
    d2v_n(1) = (dv_n(2) - dv_n(1)) / (h(1));
    d2v_n(N) = (dv_n(N) - dv_n(N-1)) / (h(N));

end
%% motion
function [r_new,phi_new,z_new,rou_new,h_new,C_s_s_new,C_phi_phi_new,T_phi_s_phi_new] = motion(r,z,phi, rou,h,C_s_s,C_phi_phi,T_phi_s_phi, v_n,dv_n,d2v_n,v_s,dv_s, dt)

    N = length(r);
    
    drou = zeros(1, N);
    for i = 2:N-1
        drou(i) = (rou(i+1) - rou(i-1)) / (2*h(i));
    end
    drou(1) = (rou(2) - rou(1)) / (h(1));
    drou(N) = (rou(N) - rou(N-1)) / (h(N));
    
    C_k_k = C_s_s + C_phi_phi;
    
    k = 0.25; %0.25
    rou1 = 0.96;
    rou2 = 0.982;
    rou_0 = linspace(rou1, rou2, N);
    
    
    h_new = h + dt * C_s_s.* v_n;
    C_s_s_new = C_s_s - dt * (C_s_s.^2.*v_n+d2v_n);
    C_phi_phi_new = C_phi_phi - dt * (C_phi_phi.^2.*v_n+T_phi_s_phi.*dv_n);
    T_phi_s_phi_new = T_phi_s_phi  + dt * (-T_phi_s_phi.*C_phi_phi.*v_n + C_phi_phi.*dv_n);
    rou_new = rou + dt*(-v_s.*drou-rou.*dv_s-rou.*T_phi_s_phi.*v_s-rou.*C_k_k.*v_n-k*(rou-rou_0));

    
    phi_new = linspace(0, 0, N);
    phi_new(1) = phi(1) - dt * dv_n(1);
    for i = 1:(N-1)
        phi_new(i+1) = phi_new(i) + (C_s_s_new(i) * h_new(i)+C_s_s_new(i+1) * h_new(i+1))/2;
    end

    r_new = linspace(0, 0, N);
    z_new = linspace(0, 0, N);
    r_new(1) = r(1) + dt * v_n(1)*sin(phi(1));
    z_new(1) = z(1);
    for i = 1:(N-1)
        z_new(i+1) = z_new(i) + (sin(phi_new(i)) * h_new(i)+sin(phi_new(i+1)) * h_new(i+1))/2;
        r_new(i+1) = r_new(i) + (cos(phi_new(i)) * h_new(i)+cos(phi_new(i+1)) * h_new(i+1))/2;
    end
    
end

