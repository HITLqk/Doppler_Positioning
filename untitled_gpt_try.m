%function demoDopplerGroundReceiverNewSatData_GeodeticOutput()
    %% ========== 1. 参数与常量定义 ==========
    clear; clc;
    
    % WGS-84参数
    a = 6378137;                     % 长半轴 [m]
    f = 1/298.257223563;             % 扁率
    e2 = 2*f - f^2;                  % 第一偏心率平方

    % 常量
    c = 3e8;                         % 光速 [m/s]
    lambda = 0.19029367;              % 载波波长 (例如GPS L1) [m]

    %% ========== 2. 接收机状态设定 ==========
    % 地面接收机地理坐标: lat = 45.75 deg, lon = 126.65 deg, alt = 0 km
    lat_rec_deg = 45.75;
    lon_rec_deg = 126.65;
    alt_rec_km  = 0;
    
    lat_rec_deg_guess = 44.75;
    lon_rec_deg_guess = 125.65;
    alt_rec_km_guess  = 0;

    % 转换为弧度和米
    lat_rec = deg2rad(lat_rec_deg);
    lon_rec = deg2rad(lon_rec_deg);
    alt_rec = alt_rec_km * 1000;  % [m]
    
    lat_rec_guess = deg2rad(lat_rec_deg_guess);
    lon_rec_guess = deg2rad(lon_rec_deg_guess);
    alt_rec_guess = alt_rec_km_guess * 1000;  % [m]

    % 计算接收机ECEF位置
    r_true = geodetic2ecef(lat_rec, lon_rec, alt_rec, a, e2);
    % 接收机速度：设为零（静止）
    v_true = [0; 0; 0];
    % 假设接收机钟偏率 (dδ_R/dt)
    dDeltaRdt_true = 1e-8;  % [s/s]
    % 状态向量 x_true 为 7维： [r (3); v (3); d_deltaR/dt (1)]
    x_true = [r_true; v_true; dDeltaRdt_true];
    
    % 初始猜测 x0 (在真实状态上加扰动)
    %x0 = x_true + [100; -50; 80; 2; -3; 1; 5e-9];
    x0 = [geodetic2ecef(lat_rec_guess, lon_rec_guess, alt_rec_guess, a, e2); [0; 0; 0]; dDeltaRdt_true];
    
    %% ========== 3. 卫星数据输入与转换 ==========
    % 卫星数据: 每行: [Lat(deg), Lon(deg), Alt(km), LatRate(deg/s), LonRate(deg/s), AltRate(km/s)]
    satData = [...
        51.979,   93.150,  552.679853,   0.014623,   0.094727,   0.002798;  % STARLINK-1008_44714
        43.716,  105.748,  551.083427,   0.035137,   0.067749,   0.007329;  % STARLINK-1039_44744
        37.589,  115.247,  549.894428,   0.041150,   0.055722,   0.007486;  % STARLINK-1193_45100
        23.995,  120.646,  547.522829,  -0.047565,   0.040975,  -0.004785;  % STARLINK-1582_46043
        31.153,  133.213,  548.216734,  -0.044950,   0.047234,  -0.006859;  % STARLINK-1292_45394
        37.890,  144.494,  549.371205,  -0.040929,   0.056219,  -0.007745;  % STARLINK-1300_45374
        48.505,  135.964,  346.164370,   0.028053,   0.085285,   0.006675;  % STARLINK-1170_45073
        52.683,  139.925,  430.119274,  -0.009879,   0.100688,  -0.002644]; % STARLINK-1474_45738
    
    lat_deg = satData(:,1);
    lon_deg = satData(:,2);
    alt_km  = satData(:,3);
    latRate_degPerSec = satData(:,4);
    lonRate_degPerSec = satData(:,5);
    altRate_kmPerSec  = satData(:,6);
    
    % 单位转换: 纬度和经度转换为弧度；alt从 km 转为 m
    lat_rad = deg2rad(lat_deg);
    lon_rad = deg2rad(lon_deg);
    alt_m = alt_km * 1000;
    latRate_radPerSec = deg2rad(latRate_degPerSec);
    lonRate_radPerSec = deg2rad(lonRate_degPerSec);
    altRate_mPerSec   = altRate_kmPerSec * 1000;
    
    nSat = size(satData,1);
    r_sat = zeros(3, nSat);
    v_sat = zeros(3, nSat);
    for j = 1:nSat
        r_sat(:,j) = geodetic2ecef(lat_rad(j), lon_rad(j), alt_m(j), a, e2);
        v_sat(:,j) = geodeticRates2ecef(lat_rad(j), lon_rad(j), alt_m(j), ...
                                         latRate_radPerSec(j), lonRate_radPerSec(j), altRate_mPerSec(j), a, e2);
    end
    
    % 卫星钟偏率：取0 (简化)
    dot_delta_sat = zeros(nSat,1);
    
    %% ========== 4. 根据简化多普勒公式生成观测值 ==========
    % 对于第 j 颗卫星：
    % f_j(x_true) = hat_rho_j'*(v_true - v_sat(:,j)) + c*dDeltaRdt_true - c*dot_delta_sat(j) + lambda*D_j = 0
    % 故：
    % D_j = -[hat_rho_j'*(v_true - v_sat(:,j)) + c*dDeltaRdt_true - c*dot_delta_sat(j)]/lambda
    doppler_meas = zeros(nSat,1);
    for j = 1:nSat
        d_vec = x_true(1:3) - r_sat(:,j);
        norm_d = norm(d_vec);
        hat_rho = d_vec / norm_d;
        relative_v = x_true(4:6) - v_sat(:,j);
        term = hat_rho' * relative_v + c*x_true(7) - c*dot_delta_sat(j);
        doppler_meas(j) = - term / lambda;
    end
    
    fprintf('生成的多普勒测量值 (Hz):\n');
    disp(doppler_meas);












    %% ========== 5. 构造测量方程组 ==========
    % 对于每颗卫星 j，测量方程为：
    % f_j(x) = ( (r - r_sat(:,j))^T / ||r - r_sat(:,j)|| )*(v - v_sat(:,j)) + c*dDeltaRdt - c*dot_delta_sat(j) + lambda*D_j = 0
    % 其中 x为7x1状态向量 [r(3); v(3); d_deltaR/dt]
    res = dopplerResidual7(x_true, r_sat, v_sat, doppler_meas, dot_delta_sat, c, lambda);
    fprintf('真实状态下各卫星的残差 (应接近0):\n');
    disp(res);
    
    %% ========== 6. 牛顿迭代求解 ==========
    % doIteration = true;
    % if doIteration
        maxIter = 50;
        tol = 1e-8;
        x_est = x0;
        err_history = zeros(maxIter,1);
        
        for k = 1:maxIter
            f_vec = dopplerResidual7(x_est, r_sat, v_sat, doppler_meas, dot_delta_sat, c, lambda);
            if norm(f_vec) < tol
                fprintf('迭代收敛，共迭代 %d 次。\n', k);
                break;
            end
            J = dopplerJacobian7(x_est, r_sat, v_sat, dot_delta_sat, c, lambda);
            Delta_x = - (J \ f_vec);
            x_est = x_est + Delta_x;
            err_history(k) = norm(Delta_x);
        end
        
        %% 将初始猜测和最终估计的接收机位置转换为地理坐标
        [lat0, lon0, alt0] = ecef2geodetic(x0(1:3), a, f);
        [lat_est, lon_est, alt_est] = ecef2geodetic(x_est(1:3), a, f);
        [lat_true, lon_true, alt_true] = ecef2geodetic(x_true(1:3), a, f);
        
        fprintf('真实状态 (geodetic: lat,lon,alt[km]):\n');
        fprintf('  Lat = %.6f deg, Lon = %.6f deg, Alt = %.3f km\n', lat_true, lon_true, alt_true/1000);
        fprintf('初始猜测 (geodetic):\n');
        fprintf('  Lat = %.6f deg, Lon = %.6f deg, Alt = %.3f km\n', lat0, lon0, alt0/1000);
        fprintf('最终估计 (geodetic):\n');
        fprintf('  Lat = %.6f deg, Lon = %.6f deg, Alt = %.3f km\n', lat_est, lon_est, alt_est/1000);
        fprintf('最终残差 norm: %e\n', norm(dopplerResidual7(x_est, r_sat, v_sat, doppler_meas, dot_delta_sat, c, lambda)));
        
        figure;
        plot(err_history(1:k),'-o');
        xlabel('迭代次数');
        ylabel('状态增量范数');
        title('牛顿迭代收敛情况');
    % end
%end











%% ----------------- 辅助函数 -----------------
function r_ecef = geodetic2ecef(lat, lon, alt, a, e2)
    % 将地理坐标 (lat, lon, alt) 转换为ECEF (WGS-84)
    % lat, lon: 弧度制, alt: m
    N = a ./ sqrt(1 - e2*sin(lat).^2);
    x = (N + alt).*cos(lat).*cos(lon);
    y = (N + alt).*cos(lat).*sin(lon);
    z = ((1 - e2)*N + alt).*sin(lat);
    r_ecef = [x; y; z];
end

function [lat, lon, alt] = ecef2geodetic(r, a, f)
    % ecef2geodetic: 将ECEF坐标转换为WGS-84地理坐标 (lat, lon, alt)
    % 输入：
    %   r: 3x1向量 [x; y; z] (m)
    %   a: 长半轴 (m)
    %   f: 扁率
    % 输出：
    %   lat, lon: 单位为度, alt: 单位为 m
    e2 = 2*f - f^2;
    x = r(1); y = r(2); z = r(3);
    lon = atan2(y, x);
    p = sqrt(x.^2 + y.^2);
    % 初始纬度估计
    lat = atan2(z, p*(1 - e2));
    lat0 = 0;
    % 迭代求解
    while abs(lat - lat0) > 1e-12
        lat0 = lat;
        N = a / sqrt(1 - e2*sin(lat)^2);
        alt = p / cos(lat) - N;
        lat = atan2(z, p*(1 - e2*(N/(N+alt))));
    end
    % 转为度
    lat = rad2deg(lat);
    lon = rad2deg(lon);
end

function v_ecef = geodeticRates2ecef(lat, lon, alt, latRate, lonRate, altRate, a, e2)
    % 将地理速率转换为ECEF速度 (m/s)
    N = a/sqrt(1 - e2*sin(lat)^2);
    dN_dlat = a*e2*sin(lat)*cos(lat)/( (1 - e2*sin(lat)^2)^(3/2) );
    
    dx_dlat = dN_dlat*cos(lat)*cos(lon) - (N+alt)*sin(lat)*cos(lon);
    dx_dlon = -(N+alt)*cos(lat)*sin(lon);
    dx_dalt = cos(lat)*cos(lon);
    
    dy_dlat = dN_dlat*cos(lat)*sin(lon) - (N+alt)*sin(lat)*sin(lon);
    dy_dlon = (N+alt)*cos(lat)*cos(lon);
    dy_dalt = cos(lat)*sin(lon);
    
    dz_dlat = (1-e2)*dN_dlat*sin(lat) + ((1-e2)*N+alt)*cos(lat);
    dz_dlon = 0;
    dz_dalt = sin(lat);
    
    vx = dx_dlat*latRate + dx_dlon*lonRate + dx_dalt*altRate;
    vy = dy_dlat*latRate + dy_dlon*lonRate + dy_dalt*altRate;
    vz = dz_dlat*latRate + dz_dlon*lonRate + dz_dalt*altRate;
    v_ecef = [vx; vy; vz];
end








function f = dopplerResidual7(x, r_sat, v_sat, D, dot_delta_sat, c, lambda)
    % 计算7维状态 x = [r(3); v(3); d_deltaR/dt] 对每颗卫星测量方程的残差向量 (nSat x 1)
    nSat = size(r_sat,2);
    f = zeros(nSat,1);
    r = x(1:3);
    v = x(4:6);
    dDeltaRdt = x(7);
    
    for j = 1:nSat
        d_vec = r - r_sat(:,j);
        norm_d = norm(d_vec);
        hat_rho = d_vec / norm_d;
        relative_v = v - v_sat(:,j);
        f(j) = hat_rho' * relative_v + c*dDeltaRdt - c*dot_delta_sat(j) + lambda*D(j);
    end
end

function J = dopplerJacobian7(x, r_sat, v_sat, dot_delta_sat, c, lambda)
    % 计算状态 x = [r(3); v(3); d_deltaR/dt] 对残差 f 的雅可比 (nSat x 7)
    nSat = size(r_sat,2);
    J = zeros(nSat,7);
    r = x(1:3);
    v = x(4:6);
    dDeltaRdt = x(7);
    
    for j = 1:nSat
        d_vec = r - r_sat(:,j);
        norm_d = norm(d_vec);
        hat_rho = d_vec / norm_d;
        dHat_dr = (eye(3)/norm_d) - (d_vec*d_vec')/(norm_d^3);
        J(j,1:3) = (dHat_dr*(v - v_sat(:,j)))';
        J(j,4:6) = hat_rho';
        J(j,7) = c;
    end
end
