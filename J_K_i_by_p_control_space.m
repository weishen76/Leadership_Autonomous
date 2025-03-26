
% 该代码可以通过调整pJ和pK来寻找自主型耦合的swarmalator集体行为状态 
% This code can find collective behavior states of self-coupling swarmalators by adjusting pJ and pK

clc,clear,close all
% 定义全局变量 | Define global variables
global N J K d

% 振子总数,固有频率，自推进速度 | Number of oscillators
N=500;

% 设置双δ分布的初始状态 | Set initial state for double delta distribution
pJ = 1;
pK = 0.7;
JA = 0.5;
JB = -1;
KA = 1;
KB = -1;

% 生成J和K的双δ分布 | Generate double delta distribution for J and K
J = floor(rand(N,1)-pJ)*(-1)*(JA-JB)+JB;
K = floor(rand(N,1)-pK)*(-1)*(KA-KB)+KB;

% 时间范围 | Time span
tspan = 0:0.1:500;

% 初始空间位置和相位 | Initial spatial positions and phases
initial_space = 2*rand(1, 2*N)-1; % 随机生成-1~1之间2N个随机数 | Randomly generate 2N numbers between -1 and 1
initial_theta = 2*pi*rand(1, N)-pi; % 生成-pi~pi之间N个随机数 | Generate N random numbers between -pi and pi
initial_space_theta = horzcat(initial_space, initial_theta); % 将两个向量拼接为一个向量 | Concatenate two vectors

% 定义距离函数，两点之间距离公式 | Define distance function between two points
d = @(x1,x2,y1,y2) sqrt((x1-x2).^2+(y1-y2).^2);

% 数值解常微分方程 | Numerical solution of ordinary differential equations
[t,space_theta_matrix] = ode45(@(t,space_theta) swarmalator_ode(t,space_theta),tspan,initial_space_theta);

figure;
% 为每个时刻，N个振子模拟系统行为 | Simulate system behavior for N oscillators at each moment
for i = 1:100:length(tspan)

    % 获得时刻i的所有振子空间位置和相位 | Get spatial positions and phases of all oscillators at moment i
    velocity = space_theta_matrix(i,:);
    space_x_velocity = velocity(1:2:2*N-1);
    space_y_velocity = velocity(2:2:2*N);  
    theta_velocity = velocity(2*N+1:1:3*N);
    
    % 限制θ在-pi到pi中，即对进行周期性循环，超出的自动加减2pi，使其一直在范围内 | Restrict θ to -pi to pi, periodic cycling
    theta_velocity = mod(theta_velocity + pi,2*pi) - pi; 
    
    % 绘制N个振子在当前时刻的位置 | Plot positions of N oscillators at current moment
    scatter(space_x_velocity,space_y_velocity,30,theta_velocity,'filled');  
    colormap(jet); % 设置颜色映射 | Set color mapping
    caxis([-pi pi]);   % 手动设置颜色条范围 | Manually set color bar range
    axis('auto'); % 自动设置坐标轴范围 | Auto set axis range
    axis equal;

    xlabel('x');
    ylabel('y');
    %grid on; % 显示网格线 | Show grid lines
    pause(0.1); % 暂停0.1秒 | Pause for 0.1 seconds
end

% 获得最后时刻的所有振子空间角 | Get spatial angles of all oscillators at last moment
space_angle = atan(space_y_velocity./space_x_velocity) + ((sign(space_x_velocity)-1)/2) * pi; 
space_angle = mod(space_angle + pi,2*pi) - pi; 

% 计算序参量S+和S-的值 | Calculate order parameters S+ and S-
Sup = abs((sum(exp(1i*(space_angle + theta_velocity))))/N);
Sdown  = abs((sum(exp(1i*(space_angle - theta_velocity))))/N);
S2 = (1/N)*max(abs(sum(exp(1i*(space_angle + 2.*theta_velocity)))),abs((sum(exp(1i*(space_angle - 2.*theta_velocity))))));
R  = abs((sum(exp(1i*theta_velocity)))/N);
R2  = abs((sum(exp(2*1i*theta_velocity)))/N);

% 显示计算结果 | Display calculation results
disp(['S+ = ',num2str(Sup)])
disp(['S- = ',num2str(Sdown)])
disp(['S2 = ',num2str(S2)])
disp(['R = ',num2str(R)])
disp(['R2 = ',num2str(R2)])

% swarmalator模型的ODEs | ODEs for swarmalator model
function dtheta_and_dspace_dt = swarmalator_ode(t,space_theta)
    global N J K d
    % 初始化返回向量 | Initialize return vector
    dtheta_and_dspace_dt = zeros(3*N,1);
  
    % 分别得到横坐标，纵坐标，和相位的向量 | Get vectors of x-coordinates, y-coordinates, and phases
    space_x = space_theta(1:2:2*N-1);
    space_y = space_theta(2:2:2*N);
    theta = space_theta(2*N+1:1:3*N);

    for i=1:N
        % 所有点与第i个点的距离组成的向量 | Vector of distances between all points and the i-th point
        distance = d(space_x,space_x(i),space_y,space_y(i));

        % 将方程级数部分看成的向量 | Vector of equation series part
        vq_space_x = ((space_x-space_x(i))./distance).*(1+J(i).*cos(theta-theta(i)))-(space_x-space_x(i))./(distance.^2);
        vq_space_y = ((space_y-space_y(i))./distance).*(1+J(i).*cos(theta-theta(i)))-(space_y-space_y(i))./(distance.^2);
        vq_theta = K(i).*(sin(theta-theta(i)))./distance;

        % 将向量中所有为NaN的项改为0 | Replace NaN items in vector with 0
        vq_space_x(isnan(vq_space_x)) = 0;
        vq_space_y(isnan(vq_space_y)) = 0;
        vq_theta(isnan(vq_theta)) = 0;

        % 计算方程右式 | Calculate right side of equation
        dtheta_and_dspace_dt(2*i-1) = (1/N)*sum(vq_space_x);
        dtheta_and_dspace_dt(2*i) = (1/N)*sum(vq_space_y);
        dtheta_and_dspace_dt(i+2*N) = (1/N)*sum(vq_theta);

    end
    % 显示当前时间 | Display current time
    disp(['t = ' ,num2str(t)]);
end
