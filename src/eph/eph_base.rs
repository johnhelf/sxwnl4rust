use crate::eph::delta_t::dt_t;
use crate::eph::prece::Prece;
use crate::eph::xl::XL;
pub use std::f64::consts::PI;

// 天文常数
pub const CS_R_EAR: f64 = 6378.1366; // 地球赤道半径(千米)
pub const CS_R_EAR_A: f64 = 0.99834 * CS_R_EAR; // 平均半径
pub const CS_BA: f64 = 0.99664719; // 地球极赤半径比
pub const CS_BA2: f64 = CS_BA * CS_BA; // 地球极赤半径比的平方
pub const CS_AU: f64 = 1.49597870691e8; // 天文单位长度(千米)
pub const CS_SIN_P: f64 = CS_R_EAR / CS_AU; // sin(太阳视差)
// pub const CS_PI: f64 = CS_SIN_P.asin(); // 太阳视差
pub const CS_GS: f64 = 299792.458; // 光速(千米/秒)
pub const CS_AGX: f64 = CS_AU / CS_GS / 86400.0 / 36525.0; // 每天文单位的光行时间(儒略世纪)

// 行星会合周期
pub const CS_XX_HH: [i32; 8] = [116, 584, 780, 399, 378, 370, 367, 367];

// 天体名称
pub const XX_NAME: [&str; 9] = [
    "地球",
    "水星",
    "金星",
    "火星",
    "木星",
    "土星",
    "天王星",
    "海王星",
    "冥王星",
];

// 角度转换常数
pub const RAD: f64 = 180.0 * 3600.0 / PI; // 每弧度的角秒数
pub const RADD: f64 = 180.0 / PI; // 每弧度的度数
pub const PI2: f64 = 2.0 * PI;
pub const PI_2: f64 = PI / 2.0;
pub const J2000: f64 = 2451545.0;

// 月食参数
pub const CS_K: f64 = 0.2725076; // 月亮与地球的半径比(用于半影计算)
pub const CS_K2: f64 = 0.2722810; // 月亮与地球的半径比(用于本影计算)
pub const CS_K0: f64 = 109.1222; // 太阳与地球的半径比
pub const CS_S_MOON: f64 = CS_K * CS_R_EAR * 1.0000036 * RAD; // 用于月亮视半径计算
pub const CS_S_MOON2: f64 = CS_K2 * CS_R_EAR * 1.0000036 * RAD; // 用于月亮视半径计算
pub const CS_S_SUN: f64 = 959.64; // 用于太阳视半径计算

// 太阳视差
pub fn cs_pi() -> f64 {
    CS_SIN_P.asin()
}

// 数学工具函数
// 注:大部分数学函数可以直接使用Rust标准库,不需要转换
pub fn int2(v: f64) -> f64 {
    v.floor()
}

// pub fn mod2(v: f64, n: f64) -> f64 {
//     ((v % n) + n) % n
// }

//=================================角度格式化=======================================

/// 将弧度转为字符串,ext为小数保留位数
/// tim=0 输出格式示例: -23°59' 48.23"
/// tim=1 输出格式示例: 18h 29m 44.52s
pub fn rad2str_e(mut d: f64, tim: i32, ext: i32) -> String {
    let mut s = String::from(" ");
    let (mut w1, mut w2, mut w3) = ("°", "'", "\"");

    if d < 0.0 {
        d = -d;
        s = String::from("-");
    }

    if tim == 1 {
        d *= 12.0 / PI;
        w1 = "h ";
        w2 = "m";
        w3 = "s";
    } else {
        d *= 180.0 / PI;
    }

    let mut a = d.floor();
    d = (d - a) * 60.0;
    let mut b = d.floor();
    d = (d - b) * 60.0;
    let mut c = d.floor();

    let q = 10.0_f64.powi(ext);
    d = ((d - c) * q + 0.5).floor();

    if d >= q {
        d -= q;
        c += 1.0;
    }
    if c >= 60.0 {
        c -= 60.0;
        b += 1.0;
    }
    if b >= 60.0 {
        b -= 60.0;
        a += 1.0;
    }

    let a_str = format!("   {}", a as i32);
    let b_str = format!("0{}", b as i32);
    let c_str = format!("0{}", c as i32);
    let d_str = format!("00000{}", d as i32);

    s.push_str(&a_str[a_str.len() - 3..]);
    s.push_str(w1);
    s.push_str(&b_str[b_str.len() - 2..]);
    s.push_str(w2);
    s.push_str(&c_str[c_str.len() - 2..]);

    if ext > 0 {
        s.push('.');
        s.push_str(&d_str[d_str.len() - ext as usize..]);
        s.push_str(w3);
    }

    s
}

/// 将弧度转为字符串,保留2位小数
pub fn rad2str(d: f64, tim: i32) -> String {
    rad2str_e(d, tim, 2)
}

/// 将弧度转为字符串,精确到分
/// 输出格式示例: -23°59'
pub fn rad2str2(mut d: f64) -> String {
    let mut s = String::from("+");
    if d < 0.0 {
        d = -d;
        s = String::from("-");
    }

    d *= 180.0 / PI;
    let a = d.floor();
    let b = ((d - a) * 60.0 + 0.5).floor();
    let mut a = a as i32;
    let mut b = b as i32;

    if b >= 60 {
        b -= 60;
        a += 1;
    }

    let a_str = format!("   {}", a);
    let b_str = format!("0{}", b);

    s.push_str(&a_str[a_str.len() - 3..]);
    s.push('°');
    s.push_str(&b_str[b_str.len() - 2..]);
    s.push('\'');

    s
}

/// 秒转为分秒
/// fx为小数点位数
/// fs为1转为"分秒"格式,为2转为"m s"格式,否则转为角度分秒格式
pub fn m2fm(mut v: f64, fx: usize, fs: i32) -> String {
    let mut gn = String::new();
    if v < 0.0 {
        v = -v;
        gn = String::from("-");
    }

    let f = (v / 60.0).floor();
    let m = v - f * 60.0;

    match fs {
        0 => format!("{}{}'{}\"", gn, f as i32, format!("{:.1$}", m, fx)),
        1 => format!("{}{}分{}秒", gn, f as i32, format!("{:.1$}", m, fx)),
        2 => format!("{}{}m{}s", gn, f as i32, format!("{:.1$}", m, fx)),
        _ => format!("{}{}'{}\"", gn, f as i32, format!("{:.1$}", m, fx)),
    }
}

/// 字符串转弧度
/// f=1表示输入的s为时分秒
pub fn str2rad(s: &str, f: bool) -> f64 {
    let mut fh = 1.0;
    let f = if f { 15.0 } else { 1.0 };

    if s.contains('-') {
        fh = -1.0;
    }

    let s = s
        .replace(&['h', 'm', 's', '-', '°', '\'', '"'][..], " ")
        .split_whitespace()
        .take(3)
        .map(|x| x.parse::<f64>().unwrap_or(0.0))
        .collect::<Vec<f64>>();

    fh * (s[0] * 3600.0 + s[1] * 60.0 + s[2]) / RAD * f
}

//=================================数学工具=========================================
/// 对超过0-2PI的角度转为0-2PI
pub fn rad2mrad(mut v: f64) -> f64 {
    v %= 2.0 * PI;
    if v < 0.0 {
        return v + 2.0 * PI;
    }
    v
}

/// 对超过-PI到PI的角度转为-PI到PI
pub fn rad2rrad(mut v: f64) -> f64 {
    v %= 2.0 * PI;
    if v <= -PI {
        return v + 2.0 * PI;
    }
    if v > PI {
        return v - 2.0 * PI;
    }
    v
}

/// 临界余数(a与最近的整倍数b相差的距离)
pub fn mod2(a: f64, b: f64) -> f64 {
    let c = (a + b) % b;
    if c > b / 2.0 { c - b } else { c }
}

/// 球面转直角坐标
pub fn llr2xyz(jw: &[f64]) -> [f64; 3] {
    let (j, w, r) = (jw[0], jw[1], jw[2]);
    [r * w.cos() * j.cos(), r * w.cos() * j.sin(), r * w.sin()]
}

/// 直角坐标转球面
pub fn xyz2llr(xyz: &[f64]) -> [f64; 3] {
    let (x, y, z) = (xyz[0], xyz[1], xyz[2]);
    let r = (x * x + y * y + z * z).sqrt();
    [rad2mrad(y.atan2(x)), (z / r).asin(), r]
}

/// 球面坐标旋转
/// 黄道赤道坐标变换,赤到黄E取负
pub fn llr_conv(jw: &[f64], e: f64) -> [f64; 3] {
    let (j, w) = (jw[0], jw[1]);
    let mut r = [
        (j.sin() * e.cos() - w.tan() * e.sin()).atan2(j.cos()),
        (e.cos() * w.sin() + e.sin() * w.cos() * j.sin()).asin(),
        jw[2],
    ];
    r[0] = rad2mrad(r[0]);
    r
}

/// 赤道坐标转为地平坐标
pub fn cd2dp(z: &[f64], l: f64, fa: f64, gst: f64) -> [f64; 3] {
    let mut a = [z[0] + PI / 2.0 - gst - l, z[1], z[2]];
    a = llr_conv(&a, PI / 2.0 - fa);
    a[0] = rad2mrad(PI / 2.0 - a[0]);
    a
}

/// 求角度差
pub fn j1_j2(j1: f64, w1: f64, j2: f64, w2: f64) -> f64 {
    let dj = rad2rrad(j1 - j2);
    let dw = w1 - w2;

    if dj.abs() < 1.0 / 1000.0 && dw.abs() < 1.0 / 1000.0 {
        let dj = dj * ((w1 + w2) / 2.0).cos();
        return (dj * dj + dw * dw).sqrt();
    }

    (w1.sin() * w2.sin() + w1.cos() * w2.cos() * dj.cos()).acos()
}

/// 日心球面转地心球面
/// z: 星体球面坐标
/// a: 地球球面坐标
/// 本函数是通用的球面坐标中心平移函数,行星计算中将反复使用
pub fn h2g(z: &[f64], a: &[f64]) -> [f64; 3] {
    let a = llr2xyz(a); // 地球
    let mut z = llr2xyz(z); // 星体
    z[0] -= a[0];
    z[1] -= a[1];
    z[2] -= a[2];
    xyz2llr(&z)
}

/// 视差角(不是视差)
pub fn shi_cha_j(gst: f64, l: f64, fa: f64, j: f64, w: f64) -> f64 {
    let h = gst + l - j; // 天体的时角
    rad2mrad(h.sin().atan2(fa.tan() * w.cos() - w.sin() * h.cos()))
}

//=================================蒙气改正=========================================
/// 大气折射修正 (h为真高度)
pub fn mqc(h: f64) -> f64 {
    0.0002967 / ((h + 0.003138 / (h + 0.08919)).tan())
}

/// 大气折射修正 (ho为视高度)
pub fn mqc2(ho: f64) -> f64 {
    -0.0002909 / ((ho + 0.002227 / (ho + 0.07679)).tan())
}

//=================================视差改正=========================================
/// 视差修正
/// z: 赤道坐标
/// h: 时角
/// fa: 地理纬度
/// high: 海拔(千米)
pub fn parallax(z: &mut [f64], h: f64, fa: f64, high: f64) -> &mut [f64] {
    // 视差修正
    let mut dw = 1.0;
    if z[2] < 500.0 {
        dw = CS_AU;
    }
    z[2] *= dw;

    let f = CS_BA;
    let u = (f * fa.tan()).atan();
    let g = z[0] + h;

    // 站点与地地心向径的赤道投影长度
    let r0 = CS_R_EAR * u.cos() + high * fa.cos();
    // 站点与地地心向径的轴向投影长度
    let z0 = CS_R_EAR * u.sin() * f + high * fa.sin();

    let x0 = r0 * g.cos();
    let y0 = r0 * g.sin();

    // 转换到直角坐标
    let mut s = llr2xyz(z);

    // 视差修正
    s[0] -= x0;
    s[1] -= y0;
    s[2] -= z0;

    // 转回球面坐标
    let s = xyz2llr(&s);

    z[0] = s[0];
    z[1] = s[1];
    z[2] = s[2] / dw;

    z
}

//=============================一些天文基本问题=====================================
/// 返回朔日的编号
/// jd应在朔日附近,允许误差数天
pub fn suo_n(jd: f64) -> i32 {
    ((jd + 8.0) / 29.5306).floor() as i32
}

/// 太阳光行差
/// t是世纪数
pub fn gxc_sun_lon(t: f64) -> f64 {
    // 平近点角
    let v = -0.043126 + 628.301955 * t - 0.000002732 * t * t;
    let e = 0.016708634 - 0.000042037 * t - 0.0000001267 * t * t;

    // 黄经光行差
    -20.49552 * (1.0 + e * v.cos()) / RAD
}

/// 黄纬光行差
pub fn gxc_sun_lat(_t: f64) -> f64 {
    0.0
}

/// 月球经度光行差
/// 误差0.07"
pub fn gxc_moon_lon(_t: f64) -> f64 {
    -3.4e-6
}

/// 月球纬度光行差
/// 误差0.006"
pub fn gxc_moon_lat(t: f64) -> f64 {
    0.063 * (0.057 + 8433.4662 * t + 0.000064 * t * t).sin() / RAD
}

/// 计算格林尼治平恒星时
/// # Arguments
/// * `t` - 2000年首起算的日数(UT)
/// * `dt` - deltatT(日),精度要求不高时dt可取值为0
/// # Returns
/// * 格林尼治平恒星时(不含赤经章动及非多项式部分)
///   即格林尼治子午圈的平春风点起算的赤经
pub fn p_gst(t: f64, dt: f64) -> f64 {
    let t = (t + dt) / 36525.0;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;

    PI2 * (0.7790572732640 + 1.00273781191135448 * t)
        + (0.014506 + 4612.15739966 * t + 1.39667721 * t2 - 0.00009344 * t3 + 0.00001882 * t4) / RAD
}

/// 计算平恒星时
/// # Arguments
/// * `jd` - 力学时J2000起算日数
/// # Returns
/// * 平恒星时
pub fn p_gst2(jd: f64) -> f64 {
    let dt = dt_t(jd);
    p_gst(jd - dt, dt)
}

/// 太阳升降计算
/// # Arguments
/// * `jd` - 儒略日(须接近L当地平午UT)
/// * `l` - 地理经度
/// * `fa` - 地理纬度
/// * `sj` - -1表示日出,1表示日落
/// # Returns
/// * 返回格林尼治UT
pub fn sun_sheng_j(mut jd: f64, l: f64, fa: f64, sj: i32) -> f64 {
    // 初始化
    jd = (jd + 0.5).floor() - l / PI2;

    for _ in 0..2 {
        // 计算黄赤交角
        let t = jd / 36525.0;
        let e = (84381.4060 - 46.836769 * t) / RAD;

        // 儒略世纪年数,力学时
        let t = t + (32.0 * (t + 1.8) * (t + 1.8) - 20.0) / 86400.0 / 36525.0;

        // 计算太阳黄经
        let j = (48950621.66 + 6283319653.318 * t + 53.0 * t * t - 994.0
            + 334166.0 * (4.669257 + 628.307585 * t).cos()
            + 3489.0 * (4.6261 + 1256.61517 * t).cos()
            + 2060.6 * (2.67823 + 628.307585 * t).cos() * t)
            / 10000000.0;

        // 太阳黄经正余弦值
        let sin_j = j.sin();
        let cos_j = j.cos();

        // 恒星时(子午圈位置)
        let gst = (0.7790572732640 + 1.00273781191135448 * jd) * PI2
            + (0.014506 + 4612.15739966 * t + 1.39667721 * t * t) / RAD;

        // 太阳赤经
        let a = sin_j * e.cos();
        let a = a.atan2(cos_j);

        // 太阳赤纬
        let d = (e.sin() * sin_j).asin();

        // 太阳在地平线上的cos(时角)计算
        let cos_h0 = ((-50.0 * 60.0 / RAD).sin() - fa.sin() * d.sin()) / (fa.cos() * d.cos());

        if cos_h0.abs() >= 1.0 {
            return 0.0;
        }

        // (升降时角-太阳时角)/太阳速度
        jd += rad2rrad(sj as f64 * cos_h0.acos() - (gst + l - a)) / 6.28;
    }

    jd
}

/// 时差计算(高精度)
/// # Arguments
/// * `t` - 力学时儒略世纪数
/// # Returns
/// * 时差(单位:周/天)
pub fn pty_zty(t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    // 计算黄经 L
    let mut l = (1753470142.0 + 628331965331.8 * t + 5296.74 * t2 + 0.432 * t3
        - 0.1124 * t4
        - 0.00009 * t5)
        / 1000000000.0
        + PI
        - 20.5 / RAD;

    // 计算章动
    let dl = -17.2 * (2.1824 - 33.75705 * t).sin() / RAD; // 黄经章
    let de = 9.2 * (2.1824 - 33.75705 * t).cos() / RAD; // 交角章

    // 计算真黄赤交角
    let e = Prece::hcjj(t) + de;

    // 计算地球坐标
    let mut z = [
        XL::xl0_calc(0, 0, t, 50) + PI + gxc_sun_lon(t) + dl,
        -(2796.0 * (3.1987 + 8433.46616 * t).cos()
            + 1016.0 * (5.4225 + 550.75532 * t).cos()
            + 804.0 * (3.88 + 522.3694 * t).cos())
            / 1000000000.0,
        0.0,
    ];

    // 转换到赤道坐标
    z = llr_conv(&z, e);
    z[0] -= dl * e.cos();

    // 计算时差
    l = rad2rrad(l - z[0]);

    // 返回结果(单位:周/天)
    l / PI2
}

/// 时差计算(低精度)
/// # Arguments
/// * `t` - 力学时儒略世纪数
/// # Returns
/// * 时差(单位:周/天)
///
/// 误差约在1秒以内
pub fn pty_zty2(t: f64) -> f64 {
    // 计算黄经
    let l = (1753470142.0 + 628331965331.8 * t + 5296.74 * t * t) / 1000000000.0 + PI;

    // 计算黄赤交角
    let e = (84381.4088 - 46.836051 * t) / RAD;

    // 地球坐标
    let mut z = [XL::xl0_calc(0, 0, t, 5) + PI, 0.0, 0.0];

    // 转换到赤道坐标
    z = llr_conv(&z, e);

    // 计算时差并转换单位
    let l = rad2rrad(l - z[0]);
    l / PI2 // 单位是周(天)
}
