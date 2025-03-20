use crate::eph::eph_base::{CS_AGX, h2g, j1_j2};
use crate::eph::eph_base::{
    CS_AU, CS_BA, CS_BA2, CS_GS, CS_R_EAR, CS_XX_HH, PI2, llr_conv, llr2xyz, mqc, p_gst2, parallax,
    rad2mrad, rad2rrad, rad2str,
};
use crate::eph::nutation::Nutation;
use crate::eph::prece::Prece;
use crate::eph::xl::XL;
pub use std::f64::consts::PI;

/// 计算行星距角
/// # Arguments
/// * `xt` - 行星序号
/// * `t` - 儒略世纪数
/// * `jing` - 精度控制参数
///   - 0: 低精度
///   - 1: 中精度
///   - >=2: 高精度(补光行时)
pub fn xing_jj(xt: usize, t: f64, jing: i32) -> f64 {
    // 初始化地球和行星的坐标
    let (mut a, mut z) = ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);

    // 低精度计算
    if jing == 0 {
        a = XL::p_coord(0, t, 10, 10, 10); // 地球
        z = XL::p_coord(xt, t, 10, 10, 10); // 行星
        z = h2g(&z, &a); // 转到地心坐标
    }
    // 中精度计算
    else if jing == 1 {
        a = XL::p_coord(0, t, 60, 60, 60); // 地球
        z = XL::p_coord(xt, t, 60, 60, 60); // 行星
        z = h2g(&z, &a); // 转到地心坐标
    }
    // 高精度计算(补光行时)
    else {
        a = XL::p_coord(0, t - a[2] * CS_AGX, -1, -1, -1); // 地球
        z = XL::p_coord(xt, t - z[2] * CS_AGX, -1, -1, -1); // 行星
        z = h2g(&z, &a); // 转到地心坐标
    }

    // 转换太阳坐标
    a[0] += PI;
    a[1] = -a[1];

    // 返回两个天体的夹角
    j1_j2(z[0], z[1], a[0], a[1])
}

/// 计算内行星大距
/// 大距计算超底速算法
/// # Arguments
/// * `xt` - 行星编号(1=水星,2=金星)
/// * `t` - 儒略世纪数TD
/// * `dx` - true表示东大距,false表示西大距
/// # Returns
/// * (t, r) - (大距时间,角距离)
pub fn da_ju(xt: usize, mut t: f64, dx: bool) -> (f64, f64) {
    // 行星周期和修正参数
    let (a, c) = match xt {
        1 => (
            // 水星
            115.8774777586 / 36525.0,
            [2.0, 0.2, 0.01, 46.0, 87.0],
        ),
        2 => (
            // 金星
            583.9213708245 / 36525.0,
            [4.0, 0.2, 0.01, 382.0, 521.0],
        ),
        _ => panic!("Invalid planet number"),
    };

    // 确定大距时间基准值
    let b = if dx {
        c[3] / 36525.0 // 东大距
    } else {
        c[4] / 36525.0 // 西大距
    };

    // 计算大距平时间
    t = b + a * ((t - b) / a + 0.5).floor();

    // 三次迭代提高精度
    let mut r1 = 0.0;
    let mut r2;
    let mut r3 = 0.0;

    for i in 0..3 {
        let dt = c[i] / 36525.0;
        r1 = xing_jj(xt, t - dt, i as i32);
        r2 = xing_jj(xt, t, i as i32);
        r3 = xing_jj(xt, t + dt, i as i32);
        t += (r1 - r3) / (r1 + r3 - 2.0 * r2) * dt / 2.0;
    }

    // 最终修正
    r2 = xing_jj(xt, t, 2); // 使用最后一次的精度
    r2 += (r1 - r3) / (r1 + r3 - 2.0 * r2) * (r3 - r1) / 8.0;

    (t, r2)
}

/// 计算行星的视坐标
/// # Arguments
/// * `xt` - 行星序号
/// * `t` - 儒略世纪数
/// * `n` - 计算精度
/// * `gxs` - 光行时修正
/// # Returns
/// 行星视赤道坐标'\['经度,纬度,距离'\]'
pub fn xing_liu0(xt: usize, t: f64, n: i32, gxs: f64) -> [f64; 3] {
    // 计算黄赤交角
    let mut e = Prece::hcjj(t);

    // 计算地球坐标
    let a = XL::p_coord(0, t - gxs, n, n, n);

    // 计算行星坐标
    let mut z = XL::p_coord(xt, t - gxs, n, n, n);

    // 转到地心坐标
    z = h2g(&z, &a);

    // 如果计算了光行时,需要补充章动
    if gxs != 0.0 {
        let zd = Nutation::nutation2(t); // 章动计算
        z[0] += zd[0]; // 补充黄经章动
        e += zd[1]; // 补充交角章动
    }

    // 转换到赤道坐标系
    llr_conv(&z, e)
}

/// 计算行星留(顺留逆留)
/// # Arguments
/// * `xt` - 行星序号
/// * `t` - 儒略世纪数
/// * `sn` - true表示顺留,false表示逆留
/// # Returns
/// 留的时刻(儒略世纪数)
pub fn xing_liu(xt: usize, mut t: f64, sn: bool) -> f64 {
    // 先求冲(下合)
    let hh = CS_XX_HH[xt - 1] as f64 / 36525.0; // 会合周期

    // 计算行星相对地球的黄经平速度
    let mut v = PI2 / hh;
    if xt > 2 {
        v = -v;
    }

    // 水星的平速度与真速度相差较多,多算几次
    for _ in 0..6 {
        t -= rad2rrad(XL::xl0_calc(xt, 0, t, 8) - XL::xl0_calc(0, 0, t, 8)) / v;
    }

    // 时间步长参数
    let tt = [5.0 / 36525.0, 1.0 / 36525.0, 0.5 / 36525.0, 2e-6, 2e-6];

    // 行星特征周期
    let tc = [17.4, 28.0, 52.0, 82.0, 86.0, 88.0, 89.0, 90.0];
    let tc = tc[xt - 1] / 36525.0;

    // 根据顺留逆留调整时间
    if sn {
        if xt > 2 {
            t -= tc;
        } else {
            t += tc;
        }
    } else if xt > 2 {
            t += tc;
    } else {
            t -= tc;
        
    }

    // 四次迭代提高精度
    let mut y2 = [0.0, 0.0, 0.0];
    for i in 0..4 {
        let dt = tt[i];
        let mut n = 8;
        let mut g = 0.0;

        if i >= 3 {
            g = y2[2] * CS_AGX;
            n = -1;
        }

        let y1 = xing_liu0(xt, t - dt, n, g);
        y2 = xing_liu0(xt, t, n, g);
        let y3 = xing_liu0(xt, t + dt, n, g);

        t += (y1[0] - y3[0]) / (y1[0] + y3[0] - 2.0 * y2[0]) * dt / 2.0;
    }

    t
}

/// 计算月亮行星视赤经差
/// # Arguments
/// * `xt` - 行星序号
/// * `t` - 儒略世纪数TD
/// * `n` - 项数
/// * `e` - 交角
/// * `g` - '\['光行时1,光行时2,章动经度,章动交角'\]'
/// # Returns
/// '\['赤经差,纬度差,光行时比率1,光行时比率2'\]'
pub fn xing_mp(xt: usize, t: f64, n: i32, e: f64, g: &[f64]) -> Vec<f64> {
    // 计算地球坐标
    let a = XL::p_coord(0, t - g[1], n, n, n);

    // 计算行星坐标
    let mut p = XL::p_coord(xt, t - g[1], n, n, n);

    // 计算月亮坐标
    let mut m = XL::m_coord(t - g[0], n, n, n);

    // 转到地心坐标
    p = h2g(&p, &a);

    // 补充章动
    m[0] += g[2];
    p[0] += g[2];

    // 转到赤道坐标
    m = llr_conv(&m, e + g[3]);
    p = llr_conv(&p, e + g[3]);

    // 计算赤经差及光行时
    vec![
        rad2rrad(m[0] - p[0]),
        m[1] - p[1],
        m[2] / CS_GS / 86400.0 / 36525.0,
        p[2] / CS_GS / 86400.0 / 36525.0 * CS_AU,
    ]
}

/// 计算行星合月(视赤经)
/// # Arguments
/// * `xt` - 行星序号
/// * `t` - 儒略世纪TD
/// # Returns
/// (时间,纬度差)
pub fn xing_hy(xt: usize, mut t: f64) -> (f64, f64) {
    // 迭代求精确合月时间
    let mut g = vec![0.0; 4];
    let mut d = Vec::new();

    for _ in 0..3 {
        d = xing_mp(xt, t, 8, 0.4091, &g);
        t -= d[0] / 8192.0;
    }

    // 计算黄赤交角
    let e = Prece::hcjj(t);

    // 计算章动
    let zd = Nutation::nutation2(t);

    // 更新光行时和章动参数
    g = vec![d[2], d[3], zd[0], zd[1]];

    // 计算速度
    let d = xing_mp(xt, t, 8, e, &g);
    let d2 = xing_mp(xt, t + 1e-6, 8, e, &g);
    let v = (d2[0] - d[0]) / 1e-6;

    // 再次迭代求精
    let d = xing_mp(xt, t, 30, e, &g);
    t -= d[0] / v;
    let d = xing_mp(xt, t, -1, e, &g);
    t -= d[0] / v;

    (t, d[1])
}

/// 计算行星太阳视黄经差与w0的差
/// # Arguments
/// * `xt` - 行星序号
/// * `t` - 儒略世纪数TD
/// * `n` - 精度项数
/// * `w0` - 黄经差基准值
/// * `ts` - 太阳光行时
/// * `tp` - 行星光行时
/// # Returns
/// '\['赤经差,纬度差,光行时s,光行时p'\]'
pub fn xing_sp(xt: usize, t: f64, n: i32, w0: f64, ts: f64, tp: f64) -> Vec<f64> {
    // 计算地球、行星、太阳坐标
    let a = XL::p_coord(0, t - tp, n, n, n); // 地球
    let mut p = XL::p_coord(xt, t - tp, n, n, n); // 行星
    let mut s = XL::p_coord(0, t - ts, n, n, n); // 太阳

    // 太阳坐标修正
    s[0] += PI;
    s[1] = -s[1];

    // 转到地心坐标
    p = h2g(&p, &a);

    // 返回赤经差及光行时
    vec![
        rad2rrad(p[0] - s[0] - w0),
        p[1] - s[1],
        s[2] * CS_AGX,
        p[2] * CS_AGX,
    ]
}

/// 计算行星合日或冲日
/// # Arguments
/// * `xt` - 行星序号
/// * `t` - 儒略世纪数TD
/// * `f` - true表示求冲(或下合),false表示求合
/// # Returns
/// (时间,纬度差)
pub fn xing_hr(xt: usize, mut t: f64, f: bool) -> (f64, f64) {
    let dt = 2e-5;

    // 设置合(冲)参数
    let (w0, w1) = if f {
        // 求冲或下合
        let w1 = if xt > 2 { PI } else { 0.0 }; // 地心黄经差(冲)
        (0.0, w1) // 日心黄经差
    } else {
        // 求合或上合
        (PI, 0.0) // 合时日心黄经差180度,地心黄经差0
    };

    // 计算行星相对地球的黄经平速度
    let mut v = PI2 / CS_XX_HH[xt - 1] as f64 * 36525.0;
    if xt > 2 {
        v = -v;
    }

    // 水星的平速度与真速度相差较多,多算几次
    for _ in 0..6 {
        t -= rad2rrad(XL::xl0_calc(xt, 0, t, 8) - XL::xl0_calc(0, 0, t, 8) - w0) / v;
    }

    // 严格计算
    let mut a = xing_sp(xt, t, 8, w1, 0.0, 0.0);
    let b = xing_sp(xt, t + dt, 8, w1, 0.0, 0.0);

    // 计算速度
    let v = (b[0] - a[0]) / dt;

    // 迭代求精
    a = xing_sp(xt, t, 40, w1, a[2], a[3]);
    t -= a[0] / v;

    a = xing_sp(xt, t, -1, w1, a[2], a[3]);
    t -= a[0] / v;

    (t, a[1])
}

/// 计算行星位置 todo: 改成非字符串返回？
/// # Arguments
/// * `xt` - 天体编号(1-9行星,10为月亮)
/// * `jd` - 儒略日(力学时)
/// * `l` - 地理经度
/// * `fa` - 地理纬度
/// # Returns
/// 位置信息字符串
pub fn xing_x(xt: usize, jd: f64, l: f64, fa: f64) -> String {
    // 基本参数计算
    let t = jd / 36525.0;
    let zd = Nutation::nutation2(t);
    let (dl, de) = (zd[0], zd[1]); // 章动 
    let e = Prece::hcjj(t) + de; // 真黄赤交角
    let gst_ping = p_gst2(jd); // 平恒星时
    let gst = gst_ping + dl * e.cos(); // 真恒星时

    let mut s = String::new();
    let mut z;
    let a;
    let z2;
    let mut a2;
    let (ra, rb, rc);
    let mut rfn = 8;

    if xt == 10 {
        // 月亮计算
        rfn = 2;
        // 求光行时并精确求出地月距
        a = XL::e_coord(t, 15, 15, 15); // 地球
        z = XL::m_coord(t, 1, 1, -1); // 月亮
        ra = z[2];

        let t2 = t - ra * CS_AGX / CS_AU; // 光行时计算

        // 求视坐标
        a2 = XL::e_coord(t2, 15, 15, 15); // 地球
        z = XL::m_coord(t2, -1, -1, -1); // 月亮
        rc = z[2];

        // 求光行距
        a2 = h2g(&a, &a2);
        a2[2] *= CS_AU;
        z2 = h2g(&z, &a2);
        rb = z2[2];

        // 地心黄道及赤道坐标
        z[0] = rad2mrad(z[0] + dl);
        s.push_str(&format!(
            "视黄经 {} 视黄纬 {} 地心距 {}\n",
            rad2str(z[0], 0),
            rad2str(z[1], 0),
            format!("{:.1$}", ra, rfn)
        ));

        z = llr_conv(&z, e); // 转到赤道坐标
        s.push_str(&format!(
            "视赤经 {} 视赤纬 {} 光行距 {}\n",
            rad2str(z[0], 1),
            rad2str(z[1], 0),
            format!("{:.1$}", rb, rfn)
        ));
    } else if xt < 10 {
        // 行星和太阳计算
        a = XL::p_coord(0, t, -1, -1, -1); // 地球
        z = XL::p_coord(xt, t, -1, -1, -1); // 行星
        z[0] = rad2mrad(z[0]);
        s.push_str(&format!(
            "黄经一 {} 黄纬一 {} 向径一 {}\n",
            rad2str(z[0], 0),
            rad2str(z[1], 0),
            format!("{:.1$}", z[2], rfn)
        ));

        // 地心黄道
        z = h2g(&z, &a);
        ra = z[2]; // 地心距
        let t2 = t - ra * CS_AGX; // 光行时

        // 重算坐标
        a2 = XL::p_coord(0, t2, -1, -1, -1); // 地球
        z2 = XL::p_coord(xt, t2, -1, -1, -1); // 行星
        z = h2g(&z2, &a);
        rb = z[2]; // 光行距
        z = h2g(&z2, &a2);
        rc = z[2]; // 视距
        z[0] = rad2mrad(z[0] + dl); // 补章动

        s.push_str(&format!(
            "视黄经 {} 视黄纬 {} 地心距 {}\n",
            rad2str(z[0], 0),
            rad2str(z[1], 0),
            format!("{:.1$}", ra, rfn)
        ));
        z = llr_conv(&z, e); // 转到赤道坐标
        s.push_str(&format!(
            "视赤经 {} 视赤纬 {} 光行距 {}\n",
            rad2str(z[0], 1),
            rad2str(z[1], 0),
            format!("{:.1$}", rb, rfn)
        ));
    } else {
        panic!("xt number error!")
    }

    // 计算时角和视差
    let sj = rad2rrad(gst + l - z[0]); // 得到天体时角
    parallax(&mut z, sj, fa, 0.0); // 视差修正
    s.push_str(&format!(
        "站赤经 {} 站赤纬 {} 视距离 {}\n",
        rad2str(z[0], 1),
        rad2str(z[1], 0),
        format!("{:.1$}", rc, rfn)
    ));

    // 转换到地平坐标
    z[0] += PI / 2.0 - gst - l;
    z = llr_conv(&z, PI / 2.0 - fa);
    z[0] = rad2mrad(PI / 2.0 - z[0]);

    // 大气折射修正
    if z[1] > 0.0 {
        z[1] += mqc(z[1]);
    }

    s.push_str(&format!(
        "方位角 {} 高度角 {}\n",
        rad2str(z[0], 0),
        rad2str(z[1], 0)
    ));
    s.push_str(&format!(
        "恒星时 {}(平) {}(真)\n",
        rad2str(rad2mrad(gst_ping), 1),
        rad2str(rad2mrad(gst), 1)
    ));
    s
}

//========日月食计算使用的一些函数=============

/// 交点信息结构体
#[derive(Default)]
pub struct EllPoint {
    /// 判别式,<0表示无解
    pub d: f64,
    /// 交点x坐标
    pub x: f64,
    /// 交点y坐标
    pub y: f64,
    /// 交点z坐标
    pub z: f64,
    /// x1到交点的距离
    pub r1: f64,
    /// x2到交点的距离
    pub r2: f64,
}

/// 求空间两点连线与地球的交点(靠近点x1的交点)
/// # Arguments
/// * `x1,y1,z1` - 第一个点的空间直角坐标
/// * `x2,y2,z2` - 第二个点的空间直角坐标
/// * `e` - 地球偏心率
/// * `r` - 地球赤道半径
/// # Returns
/// 交点信息结构体
pub fn line_ell(x1: f64, y1: f64, z1: f64, x2: f64, y2: f64, z2: f64, e: f64, r: f64) -> EllPoint {
    let mut p = EllPoint::default();

    // 计算方向向量
    let dx = x2 - x1;
    let dy = y2 - y1;
    let dz = z2 - z1;
    let e2 = e * e;

    // 求解二次方程
    let a = dx * dx + dy * dy + dz * dz / e2;
    let b = x1 * dx + y1 * dy + z1 * dz / e2;
    let c = x1 * x1 + y1 * y1 + z1 * z1 / e2 - r * r;

    p.d = b * b - a * c;
    if p.d < 0.0 {
        return p; // 判别式小于0无解
    }

    // 只求靠近x1的交点
    let d = p.d.sqrt();
    let d = if b < 0.0 { -d } else { d };
    let t = (-b + d) / a;

    // 计算交点坐标
    p.x = x1 + dx * t;
    p.y = y1 + dy * t;
    p.z = z1 + dz * t;

    // 计算到两端点的距离
    let r = (dx * dx + dy * dy + dz * dz).sqrt();
    p.r1 = r * t.abs();
    p.r2 = r * (t - 1.0).abs();

    p
}

/// 交点的经纬度等信息
#[derive(Default)]
pub struct BesselPoint {
    /// 判别式,<0表示无解
    pub d: f64,
    /// 交点x坐标
    pub x: f64,
    /// 交点y坐标
    pub y: f64,
    /// 交点z坐标
    pub z: f64,
    /// x1到交点的距离
    pub r1: f64,
    /// x2到交点的距离
    pub r2: f64,

    /// 地理经度
    pub j: f64,
    /// 地理纬度
    pub w: f64,
}

/// 贝塞尔坐标系下的交点计算
/// # Arguments
/// * `x1,y1,z1` - 第一个点坐标
/// * `x2,y2,z2` - 第二个点坐标
/// * `e` - 地球偏心率
/// * `r` - 地球赤道半径
/// * `i` - 贝塞尔坐标参数
/// # Returns
/// 交点的经纬度等信息
pub fn line_ear2(
    x1: f64,
    y1: f64,
    z1: f64,
    x2: f64,
    y2: f64,
    z2: f64,
    e: f64,
    r: f64,
    i: &[f64],
) -> BesselPoint {
    let mut p = BesselPoint::default();

    let p1 = i[1].cos();
    let q = i[1].sin();

    // 坐标变换
    // let x1 = x1;
    let y1 = p1 * y1 - q * z1;
    let z1 = q * y1 + p1 * z1;

    // let x2 = x2;
    let y2 = p1 * y2 - q * z2;
    let z2 = q * y2 + p1 * z2;

    // 求与椭球的交点
    let ep = line_ell(x1, y1, z1, x2, y2, z2, e, r);
    p.d = ep.d;
    p.x = ep.x;
    p.y = ep.y;
    p.z = ep.z;
    p.r1 = ep.r1;
    p.r2 = ep.r2;

    // 默认值
    p.j = 100.0;
    p.w = 100.0;

    if p.d < 0.0 {
        return p;
    }

    // 计算经纬度
    p.j = rad2rrad(ep.y.atan2(ep.x) + i[0] - i[2]);
    p.w = (ep.z / e / e / (ep.x * ep.x + ep.y * ep.y).sqrt()).atan();

    p
}

/// 在分点坐标中求空间两点连线与地球的交点
pub fn line_ear(p: &[f64], q: &[f64], gst: f64) -> BesselPoint {
    let p_xyz = llr2xyz(p);
    let q_xyz = llr2xyz(q);

    // 求与椭球的交点
    let r = line_ell(
        p_xyz[0], p_xyz[1], p_xyz[2], q_xyz[0], q_xyz[1], q_xyz[2], CS_BA, CS_R_EAR,
    );

    let mut bp =BesselPoint { d: r.d, ..Default::default() };
    // bp.d = r.d;

    if r.d < 0.0 {
        bp.j = 100.0;
        bp.w = 100.0;
        return bp;
    }

    // 计算地理坐标
    bp.w = (r.z / CS_BA2 / (r.x * r.x + r.y * r.y).sqrt()).atan();
    bp.j = rad2rrad(r.y.atan2(r.x) - gst);

    bp
}

/// 椭圆与圆的交点计算结果
#[derive(Default)]
pub struct OvalCircleIntersect {
    /// 交点个数
    pub n: i32,
    /// 第一个交点坐标
    pub a: [f64; 2],
    /// 第二个交点坐标
    pub b: [f64; 2],
}

/// 计算椭圆与圆的交点
/// # Arguments
/// * `r` - 椭圆长半径
/// * `ba` - 椭圆短半轴与长半轴的比值
/// * `r2` - 圆半径
/// * `x0` - 圆心x坐标
/// * `y0` - 圆心y坐标
pub fn cir_ovl(r: f64, ba: f64, r2: f64, x0: f64, y0: f64) -> OvalCircleIntersect {
    let mut re = OvalCircleIntersect::default();

    // 计算圆心到原点距离
    let d = (x0 * x0 + y0 * y0).sqrt();
    let sin_b = y0 / d;
    let cos_b = x0 / d;

    // 计算余弦值
    let cos_a = (r * r + d * d - r2 * r2) / (2.0 * d * r);
    if cos_a.abs() > 1.0 {
        // 无解
        re.n = 0;
        return re;
    }

    let sin_a = (1.0 - cos_a * cos_a).sqrt();
    let ba2 = ba * ba;

    // 计算两个交点
    for k in [-1.0, 1.0].iter() {
        let s = cos_a * sin_b + sin_a * cos_b * k;
        let g = r - s * s * (1.0 / ba2 - 1.0) / 2.0;

        let cos_a = (g * g + d * d - r2 * r2) / (2.0 * d * g);
        if cos_a.abs() > 1.0 {
            re.n = 0;
            return re;
        }

        let sin_a = (1.0 - cos_a * cos_a).sqrt();
        let c = cos_a * cos_b - sin_a * sin_b * k;
        let s = cos_a * sin_b + sin_a * cos_b * k;

        if *k == 1.0 {
            re.a = [g * c, g * s];
        } else {
            re.b = [g * c, g * s];
        }
    }

    re.n = 2;
    re
}

/// 直线与椭圆相交的计算结果
#[derive(Default)]
pub struct LineOvalIntersect {
    /// 交点个数
    pub n: i32,
    /// 第一个交点
    pub a: [f64; 2],
    /// 第二个交点
    pub b: [f64; 2],
    /// 到第一个交点的距离
    pub r1: f64,
    /// 到第二个交点的距离
    pub r2: f64,
}

/// 计算直线与椭圆的交点
/// # Arguments
/// * `x1,y1` - 直线上一点的坐标
/// * `dx,dy` - 直线的方向向量
/// * `r` - 椭圆长半径
/// * `ba` - 椭圆短半轴与长半轴的比值
pub fn line_ovl(x1: f64, y1: f64, dx: f64, dy: f64, r: f64, ba: f64) -> LineOvalIntersect {
    let mut p = LineOvalIntersect::default();

    // 计算二次方程系数
    let f = ba * ba;
    let a = dx * dx + dy * dy / f;
    let b = x1 * dx + y1 * dy / f;
    let c = x1 * x1 + y1 * y1 / f - r * r;

    // 求判别式
    let d = b * b - a * c;
    if d < 0.0 {
        p.n = 0;
        return p; // 判别式小于0无解
    }

    // 求解参数方程
    p.n = if d == 0.0 { 1 } else { 2 };
    let d = d.sqrt();
    let t1 = (-b + d) / a;
    let t2 = (-b - d) / a;

    // 计算交点坐标
    p.a = [x1 + dx * t1, y1 + dy * t1];
    p.b = [x1 + dx * t2, y1 + dy * t2];

    // 计算到交点的距离
    let l = (dx * dx + dy * dy).sqrt();
    p.r1 = l * t1.abs();
    p.r2 = l * t2.abs();

    p
}
