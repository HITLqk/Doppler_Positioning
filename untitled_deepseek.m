clc; clear; close all;

%% 参数定义
c = 299792458; % 光速 (m/s)
%f0 = 1575.42e6; % GPS L1载波频率 (Hz)
f0 = 26.5e9;
lambda = c/f0; % 波长 (m)
w_earth = 7.292115e-5; % 地球自转角速度 (rad/s)
a = 6378137.0; % WGS-84椭球长半轴 (m)
e = 0.0818191908426; % WGS-84椭球偏心率

%% 真实状态定义 (ECEF坐标系)
x_true = struct(...
    'position', [ -2314263.0; 4561865.0;  3814862.0],... % ECEF位置 (m)
    'velocity', [      -2.3;     -14.2;      29.7],... % ECEF速度 (m/s)
    'clock_bias', 0.0003,... % 接收机时钟偏差 (s)
    'clock_bias_rate', 1.2e-9); % 时钟偏差率 (s/s)

%% 初始猜测 (带较大误差)
x0 = struct(...
    'position', x_true.position + [150e3; -80e3; 50e3],... % 初始位置偏差约150km
    'velocity', [0; 0; 0],...
    'clock_bias', 0,... 
    'clock_bias_rate', 0);

%% 生成8颗卫星的轨道参数 (LLA及变化率)
% 格式: [Lat(deg), Lon(deg), Alt(km), LatRate(deg/s), LonRate(deg/s), AltRate(km/s)]
sat_lla = [...
    35.1,  -98.2, 1200,  0.003,  0.015, 0.002;...
    42.3, -112.5, 1150, -0.002,  0.012, 0.001;...
    28.7,  -85.6, 1250,  0.004,  0.018, 0.003;...
    51.2,  -72.3, 1100, -0.001,  0.010, 0.002;...
   -33.5,  150.8, 1300,  0.002, -0.014, 0.001;...
    15.9,   45.6, 1150,  0.003,  0.016, 0.002;...
   -12.7,  -68.9, 1200, -0.002,  0.013, 0.003;...
    60.8,  -45.1, 1050,  0.001,  0.011, 0.001];

%% 将卫星参数转换为ECEF坐标系下的位置和速度
num_sats = size(sat_lla,1);
sat_ecef = struct('position', zeros(3,num_sats), 'velocity', zeros(3,num_sats));

for i = 1:num_sats
    % 提取LLA参数
    lat = deg2rad(sat_lla(i,1)); 
    lon = deg2rad(sat_lla(i,2));
    alt = sat_lla(i,3)*1e3; % km转m
    
    % LLA转ECEF位置
    N = a / sqrt(1 - e^2*sin(lat)^2);
    X = (N + alt) * cos(lat) * cos(lon);
    Y = (N + alt) * cos(lat) * sin(lon);
    Z = (N*(1 - e^2) + alt) * sin(lat);
    sat_ecef.position(:,i) = [X; Y; Z];
    
    % LLA变化率转ECEF速度
    dlat = deg2rad(sat_lla(i,4)); % lat rate (rad/s)
    dlon = deg2rad(sat_lla(i,5)); % lon rate (rad/s)
    dalt = sat_lla(i,6)*1e3;      % alt rate (m/s)
    
    % 微分转换矩阵 (LLA变化率 -> ECEF速度)
    dx_dlat = -sin(lat)*cos(lon)*(N + alt) - cos(lat)*cos(lon)*(N*e^2*sin(lat)*cos(lat))/sqrt(1 - e^2*sin(lat)^2);
    dx_dlon = -Y;
    dx_dh = cos(lat)*cos(lon);
    
    dy_dlat = -sin(lat)*sin(lon)*(N + alt) - cos(lat)*sin(lon)*(N*e^2*sin(lat)*cos(lat))/sqrt(1 - e^2*sin(lat)^2);
    dy_dlon = X;
    dy_dh = cos(lat)*sin(lon);
    
    dz_dlat = cos(lat)*(N*(1 - e^2) + alt) - sin(lat)*(N*(1 - e^2)*e^2*sin(lat)/sqrt(1 - e^2*sin(lat)^2));
    dz_dlon = 0;
    dz_dh = sin(lat);
    
    Vx = dx_dlat*dlat + dx_dlon*dlon + dx_dh*dalt;
    Vy = dy_dlat*dlat + dy_dlon*dlon + dy_dh*dalt;
    Vz = dz_dlat*dlat + dz_dlon*dlon + dz_dh*dalt;
    
    sat_ecef.velocity(:,i) = [Vx; Vy; Vz] + cross([0;0;w_earth], sat_ecef.position(:,i));
end

%% 计算真实多普勒频移（带噪声）
sigma_doppler = 0.01; % 多普勒测量噪声标准差 (m/s)
doppler_meas = zeros(num_sats,1);

for i = 1:num_sats
    % 卫星到接收机的相对位置矢量
    r_rel = x_true.position - sat_ecef.position(:,i);
    range = norm(r_rel);
    los = r_rel / range; % 视线单位向量
    
    % 相对速度（接收机速度 - 卫星速度 + 地球自转补偿）
    v_rel = x_true.velocity - sat_ecef.velocity(:,i) + cross([0;0;w_earth], r_rel);
    
    % 多普勒频移模型（式13简化版）
    doppler_true = (-los'*v_rel + c*x_true.clock_bias_rate)/lambda;
    
    % 添加高斯噪声
    doppler_meas(i) = doppler_true + sigma_doppler/lambda*randn;
end

%% 结果输出
disp('=== 真实状态 ===');
disp(['位置 (m): ', num2str(x_true.position')]);
disp(['速度 (m/s): ', num2str(x_true.velocity')]);
disp(['时钟偏差 (s): ', num2str(x_true.clock_bias)]);
disp(['时钟偏差率 (s/s): ', num2str(x_true.clock_bias_rate)]);

disp('=== 卫星ECEF位置 (m) ===');
disp(sat_ecef.position);

disp('=== 生成的多普勒测量值 (Hz) ===');
disp(doppler_meas);

%% 函数：LLA到ECEF转换（备选）
function [X,Y,Z] = lla2ecef(lat, lon, alt)
    a = 6378137.0; e = 0.0818191908426;
    N = a / sqrt(1 - e^2*sin(lat)^2);
    X = (N + alt) * cos(lat) * cos(lon);
    Y = (N + alt) * cos(lat) * sin(lon);
    Z = (N*(1 - e^2) + alt) * sin(lat);
end

