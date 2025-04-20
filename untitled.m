%% 1. 未知数设定与真实状态构造
clear; clc;

% WGS-84参数
a = 6378137;                     % 长半轴 [m]
f = 1/298.257223563;             % 扁率
e2 = 2*f - f^2;                  % 第一偏心率平方

% 常量
c = 3e8;                         % 光速 [m/s]
lambda = 0.19029367;              % 载波波长（例如GPS L1）[m]

% 构造真实状态 x_true = [r; delta_R; v; d_deltaR/dt]
% r: ECEF位置 (m), delta_R: 钟偏 (s), v: ECEF速度 (m/s), d_deltaR/dt: 钟偏率 (s/s)
% 此处示例：接收机位于地球表面某点（近似）
r_true = [6370e3; 0; 0];           % [m]
deltaR_true = 0.001;               % [s]
v_true = [10; 5; 0];               % [m/s]
dDeltaRdt_true = 1e-8;             % [s/s]
x_true = [r_true; deltaR_true; v_true; dDeltaRdt_true];

% 初始猜测 x0（用于迭代求解，扰动真实状态）
x0 = x_true + [100; -50; 80; 0.0005; 2; -3; 1; 5e-9];

%% 2. 卫星数据输入与转换
% 输入8颗卫星的地理参数：Lat(deg), Lon(deg), Alt(km), Lat Rate(deg/s), Lon Rate(deg/s), Alt Rate(km/s)
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

% 转换单位：纬度、经度转换为弧度；海拔从km转换为m；速率转换为适当单位
lat_rad = deg2rad(lat_deg);
lon_rad = deg2rad(lon_deg);
alt = alt_km * 1000;  % [m]
latRate_radPerSec = deg2rad(latRate_degPerSec);
lonRate_radPerSec = deg2rad(lonRate_degPerSec);
altRate = altRate_kmPerSec * 1000; % [m/s]

% 对8颗卫星，转换为ECEF位置和速度
nSat = length(lat_rad);
r_sat = zeros(3, nSat);
v_sat = zeros(3, nSat);
% 这里采用自定义函数进行地理->ECEF转换
for j = 1:nSat
    r_sat(:,j) = geodetic2ecef(lat_rad(j), lon_rad(j), alt(j), a, e2);
    % 卫星速度转换（简单近似）：假设速度由地理速率线性转换得到
    v_sat(:,j) = geodeticRates2ecef(lat_rad(j), lon_rad(j), alt(j), ...
                                    latRate_radPerSec(j), lonRate_radPerSec(j), altRate(j), a, e2);
end

%% 3. 根据公式(13)生成多普勒频移测量
% 公式(13)（简化形式）：
% -lambda*D_j = hat_rho_j^T (v - v_j) + c*d_deltaR_R - c*d_delta_j
% 假设卫星钟偏率 d_delta_j/dT_j 已知，这里为了简化，取为0
% 则简化为：
% -lambda*D_j = hat_rho_j^T (v - v_j) + c*d_deltaR_R
% 其中：v为接收机速度 (x_true(5:7))，d_deltaR_R为接收机钟偏率 (x_true(8))
doppler_meas = zeros(nSat,1);
for j = 1:nSat
    % 计算接收机真实位置与卫星位置之差
    d_vec = x_true(1:3) - r_sat(:,j);
    % 单位视线向量
    hat_rho = d_vec / norm(d_vec);
    % 计算相对速度分量沿视线方向： (v - v_sat)
    relative_v = x_true(5:7) - v_sat(:,j);
    % 计算右边项
    term = hat_rho' * relative_v + c * x_true(8);  % 此处忽略卫星钟偏率项
    % 根据公式，得到多普勒测量（Hz）
    doppler_meas(j) = -term / lambda;
end

%% 输出结果
fprintf('真实状态 x_true:\n');
disp(x_true);
fprintf('初始猜测 x0:\n');
disp(x0);
fprintf('卫星ECEF位置 (m):\n');
disp(r_sat);
fprintf('卫星ECEF速度 (m/s):\n');
disp(v_sat);
fprintf('生成的8颗卫星多普勒频移 (Hz):\n');
disp(doppler_meas);


% 假设已经定义了 x_true, r_sat, v_sat, doppler_meas, c, lambda
residuals = dopplerResidual(x_true, r_sat, v_sat, doppler_meas, c, lambda);
disp('各卫星对应的残差（应接近0）：');
disp(residuals);




%% --- 地理坐标转换函数 ---
function r_ecef = geodetic2ecef(lat, lon, alt, a, e2)
    % 将地理坐标转换为ECEF坐标 (WGS-84)
    % lat, lon: 弧度制
    % alt: 高度 [m]
    % a: 长半轴 [m]
    % e2: 第一偏心率平方
    N = a ./ sqrt(1 - e2 * sin(lat).^2);
    x = (N + alt) .* cos(lat) .* cos(lon);
    y = (N + alt) .* cos(lat) .* sin(lon);
    z = ((1 - e2) * N + alt) .* sin(lat);
    r_ecef = [x; y; z];
end

%% --- 地理速率转换为ECEF速度 ---
function v_ecef = geodeticRates2ecef(lat, lon, alt, latRate, lonRate, altRate, a, e2)
    % 近似将地理速率转换为ECEF速度
    % lat, lon: 弧度制
    % alt: 高度 [m]
    % latRate, lonRate: [rad/s]
    % altRate: [m/s]
    % a, e2: WGS-84参数
    N = a/sqrt(1 - e2*sin(lat)^2);
    dN_dlat = a*e2*sin(lat)*cos(lat)/( (1 - e2*sin(lat)^2)^(3/2) );
    % 对x = (N+alt)*cos(lat)*cos(lon)
    dx_dlat = (dN_dlat)*cos(lat)*cos(lon) - (N+alt)*sin(lat)*cos(lon);
    dx_dlon = -(N+alt)*cos(lat)*sin(lon);
    dx_dalt = cos(lat)*cos(lon);
    % 对y = (N+alt)*cos(lat)*sin(lon)
    dy_dlat = (dN_dlat)*cos(lat)*sin(lon) - (N+alt)*sin(lat)*sin(lon);
    dy_dlon = (N+alt)*cos(lat)*cos(lon);
    dy_dalt = cos(lat)*sin(lon);
    % 对z = ((1-e2)*N+alt)*sin(lat)
    dz_dlat = ((1-e2)*dN_dlat)*sin(lat) + ((1-e2)*N+alt)*cos(lat);
    dz_dlon = 0;
    dz_dalt = sin(lat);
    
    vx = dx_dlat*latRate + dx_dlon*lonRate + dx_dalt*altRate;
    vy = dy_dlat*latRate + dy_dlon*lonRate + dy_dalt*altRate;
    vz = dz_dlat*latRate + dz_dlon*lonRate + dz_dalt*altRate;
    v_ecef = [vx; vy; vz];
end


function f = dopplerResidual(x, r_sat, v_sat, D, c, lambda)
    % dopplerResidual: 构造基于多普勒测量的残差向量
    % 输入：
    %   x      : 8x1状态向量 [r(1:3); delta_R; v(1:3); d_deltaR/dt]
    %            其中 r 为接收机位置, v 为接收机速度, d_deltaR/dt 为钟偏率
    %   r_sat  : 3xnSat 卫星位置矩阵 (ECEF, 单位: m)
    %   v_sat  : 3xnSat 卫星速度矩阵 (ECEF, 单位: m/s)
    %   D      : nSatx1 多普勒频移测量 (单位: Hz)
    %   c      : 光速 (m/s)
    %   lambda : 载波波长 (m)
    %
    % 输出：
    %   f      : nSatx1 残差向量，每个分量对应一个卫星的方程 f_j(x)=0
    
    nSat = size(r_sat, 2);
    f = zeros(nSat, 1);
    
    % x 的结构:
    % x(1:3) -> 接收机位置 r
    % x(4)   -> 接收机钟偏 delta_R (本方程中不直接出现)
    % x(5:7) -> 接收机速度 v
    % x(8)   -> 接收机钟偏率 (d_deltaR/dt)
    
    for j = 1:nSat
        d_vec = x(1:3) - r_sat(:, j);     % 接收机与第j颗卫星之间的向量
        norm_d = norm(d_vec);
        unit_d = d_vec / norm_d;           % 单位视线向量 hat_rho_j
        % 根据公式构造残差：
        % hat_rho_j^T (v - v_j) + c * (d_deltaR/dt) + lambda * D_j = 0
        f(j) = unit_d' * (x(5:7) - v_sat(:, j)) + c * x(8) + lambda * D(j);
    end
end
