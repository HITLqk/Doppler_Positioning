%% 清除环境
clear; clc;

%% 1. 参数设定与真实状态构造
% WGS-84参数
a = 6378137;                     % 长半轴 [m]
f = 1/298.257223563;             % 扁率
e2 = 2*f - f^2;                  % 第一偏心率平方

% 常量
c = 3e8;                         % 光速 [m/s]
lambda = 0.19029367;              % 载波波长 (例如GPS L1) [m]

% 构造真实状态 x_true (7维： r(3); v(3); d_deltaR/dt )
r_true = [6370e3; 0; 0];           % 接收机位置 [m]
v_true = [10; 5; 0];               % 接收机速度 [m/s]
dDeltaRdt_true = 1e-8;             % 接收机钟偏率 [s/s]
x_true = [r_true; v_true; dDeltaRdt_true];...

% 初始猜测 x0 (在真实状态上加扰动)
x0 = x_true + [100; -50; 80; 2; -3; 1; 5e-9];

%% 2. 卫星数据输入与转换
% 输入8颗卫星的地理参数: Lat(deg), Lon(deg), Alt(km), LatRate(deg/s), LonRate(deg/s), AltRate(km/s)
satData = [...
    30,  120, 0.8,  0.001,  -0.001,  0.0001;
    31,  121, 1.0, -0.001,   0.0005,  0.0002;
    29,  119, 0.9,  0.0005, -0.0005,  0.0001;
    32,  122, 1.2, -0.0008,  0.0008,  0.0003;
    28,  118, 0.7,  0.0012, -0.0007,  0.0001;
    33,  123, 1.1, -0.0011,  0.0006,  0.0002;
    27,  117, 0.6,  0.0009, -0.0009,  0.0001;
    34,  124, 1.3, -0.001,   0.0010,  0.0003];

lat_deg   = satData(:,1);    % [deg]
lon_deg   = satData(:,2);    % [deg]
alt_km    = satData(:,3);    % [km]
latRate_degPerSec = satData(:,4);  % [deg/s]
lonRate_degPerSec = satData(:,5);  % [deg/s]
altRate_kmPerSec  = satData(:,6);  % [km/s]

% 转换单位
lat_rad = deg2rad(lat_deg);
lon_rad = deg2rad(lon_deg);
alt = alt_km * 1000;  % [m]
latRate_radPerSec = deg2rad(latRate_degPerSec);
lonRate_radPerSec = deg2rad(lonRate_degPerSec);
altRate = altRate_kmPerSec * 1000; % [m/s]

nSat = length(lat_rad);
r_sat = zeros(3, nSat);
v_sat = zeros(3, nSat);
for j = 1:nSat
    r_sat(:,j) = geodetic2ecef(lat_rad(j), lon_rad(j), alt(j), a, e2);
    v_sat(:,j) = geodeticRates2ecef(lat_rad(j), lon_rad(j), alt(j), ...
                                    latRate_radPerSec(j), lonRate_radPerSec(j), altRate(j), a, e2);
end

% 卫星钟偏率 (假设已知，这里取非零示例)
dot_delta_sat = 1e-10 * ones(nSat,1);

%% 3. 构造测量方程组
% 对于每颗卫星 j，方程为：
% f_j(x) = ( (r - r_j)^T/(||r - r_j||) )*(v - v_j) + c*x(7) - c*dot_delta_sat(j) + lambda*D_j = 0.
% 其中 D_j 为多普勒频移测量（单位：Hz），这里我们利用真实状态生成观测值

doppler_meas = zeros(nSat,1);
for j = 1:nSat
    d_vec = x_true(1:3) - r_sat(:,j);
    hat_rho = d_vec / norm(d_vec);
    relative_v = x_true(4:6) - v_sat(:,j);
    term = hat_rho' * relative_v + c * x_true(7) - c * dot_delta_sat(j);
    doppler_meas(j) = -term / lambda;
end

%% 4. 牛顿迭代求解
maxIter = 50;
tol = 1e-8;
x_est = x0;
err_history = zeros(maxIter,1);

for k = 1:maxIter
    % 计算残差向量 f(x)
    f_vec = dopplerResidual7(x_est, r_sat, v_sat, doppler_meas, dot_delta_sat, c, lambda);
    
    % 判断收敛
    if norm(f_vec) < tol
        fprintf('迭代收敛，共迭代 %d 次。\n', k);
        break;
    end
    
    % 计算雅可比矩阵
    J = dopplerJacobian7(x_est, r_sat, v_sat, dot_delta_sat, c, lambda);
    
    % 牛顿更新：J * Delta_x = -f_vec
    Delta_x = - (J \ f_vec);
    x_est = x_est + Delta_x;
    err_history(k) = norm(Delta_x);
end

%% 输出结果
fprintf('真实状态 x_true:\n');
disp(x_true);
fprintf('初始猜测 x0:\n');
disp(x0);
fprintf('最终估计状态 x_est:\n');
disp(x_est);
fprintf('最终残差 norm: %e\n', norm(dopplerResidual7(x_est, r_sat, v_sat, doppler_meas, dot_delta_sat, c, lambda)));

figure;
plot(err_history(1:k),'-o');
xlabel('迭代次数');
ylabel('状态增量范数');
title('牛顿迭代收敛情况');

%% --- 函数定义 ---
function f = dopplerResidual7(x, r_sat, v_sat, D, dot_delta_sat, c, lambda)
    % 输入 x: 7x1状态向量 [r(3); v(3); d_deltaR/dt]
    % r_sat: 3xnSat, v_sat: 3xnSat, D: nSatx1, dot_delta_sat: nSatx1
    nSat = size(r_sat,2);
    f = zeros(nSat,1);
    r = x(1:3);
    v = x(4:6);
    dDeltaRdt = x(7);
    
    for j = 1:nSat
        d_vec = r - r_sat(:,j);
        norm_d = norm(d_vec);
        hat_rho = d_vec / norm_d;
        % 测量模型：hat_rho^T (v - v_sat) + c*dDeltaRdt - c*dot_delta_sat(j) + lambda*D(j) = 0
        f(j) = hat_rho'*(v - v_sat(:,j)) + c*dDeltaRdt - c*dot_delta_sat(j) + lambda*D(j);
    end
end

function J = dopplerJacobian7(x, r_sat, v_sat, dot_delta_sat, c, lambda)
    % 计算状态向量 x = [r(3); v(3); d_deltaR/dt] 对残差 f 的雅可比矩阵 J (nSat x 7)
    nSat = size(r_sat,2);
    J = zeros(nSat,7);
    r = x(1:3);
    v = x(4:6);
    dDeltaRdt = x(7);
    
    for j = 1:nSat
        d_vec = r - r_sat(:,j);
        norm_d = norm(d_vec);
        hat_rho = d_vec / norm_d;
        % 对 r 的偏导：
        % ∂hat_rho/∂r = I/norm_d - (d_vec*d_vec')/norm_d^3
        dHat_dr = (eye(3)/norm_d) - (d_vec*d_vec')/(norm_d^3);
        % 对 r： = dHat_dr*(v - v_sat)
        J(j,1:3) = (dHat_dr*(v - v_sat(:,j)))';
        % 对 v 的偏导： hat_rho'
        J(j,4:6) = hat_rho';
        % 对 dDeltaRdt 的偏导： c
        J(j,7) = c;
        % 注意：本模型中 x 中不涉及钟偏 delta_R
    end
end

%% --- 地理坐标转换函数 ---
function r_ecef = geodetic2ecef(lat, lon, alt, a, e2)
    % 将地理坐标转换为ECEF坐标 (WGS-84)
    N = a ./ sqrt(1 - e2 * sin(lat).^2);
    x = (N + alt) .* cos(lat) .* cos(lon);
    y = (N + alt) .* cos(lat) .* sin(lon);
    z = ((1 - e2) * N + alt) .* sin(lat);
    r_ecef = [x; y; z];
end

%% --- 地理速率转换为ECEF速度 ---
function v_ecef = geodeticRates2ecef(lat, lon, alt, latRate, lonRate, altRate, a, e2)
    % 近似将地理速率转换为ECEF速度
    N = a/sqrt(1 - e2*sin(lat)^2);
    dN_dlat = a*e2*sin(lat)*cos(lat)/( (1 - e2*sin(lat)^2)^(3/2) );
    % 对 x = (N+alt)*cos(lat)*cos(lon)
    dx_dlat = dN_dlat*cos(lat)*cos(lon) - (N+alt)*sin(lat)*cos(lon);
    dx_dlon = -(N+alt)*cos(lat)*sin(lon);
    dx_dalt = cos(lat)*cos(lon);
    % 对 y = (N+alt)*cos(lat)*sin(lon)
    dy_dlat = dN_dlat*cos(lat)*sin(lon) - (N+alt)*sin(lat)*sin(lon);
    dy_dlon = (N+alt)*cos(lat)*cos(lon);
    dy_dalt = cos(lat)*sin(lon);
    % 对 z = ((1-e2)*N+alt)*sin(lat)
    dz_dlat = (1-e2)*dN_dlat*sin(lat) + ((1-e2)*N+alt)*cos(lat);
    dz_dlon = 0;
    dz_dalt = sin(lat);
    
    vx = dx_dlat*latRate + dx_dlon*lonRate + dx_dalt*altRate;
    vy = dy_dlat*latRate + dy_dlon*lonRate + dy_dalt*altRate;
    vz = dz_dlat*latRate + dz_dlon*lonRate + dz_dalt*altRate;
    v_ecef = [vx; vy; vz];
end
