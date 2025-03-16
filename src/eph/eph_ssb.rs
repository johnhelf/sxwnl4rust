use crate::eph::eph_base::{
    CS_AU, CS_GS, PI, RAD, llr_conv, mqc, p_gst2, rad2mrad, rad2str_e, str2rad,J2000, //,rad2rrad
};

use crate::eph::nutation::Nutation;
use crate::eph::prece::Prece;
use crate::eph::xl::XL;
use crate::eph::jd::JD;

/****************************
算法取自《天文算法》
****************************/

/// 地球的SSB速度计算表
const EV_TAB: &[i32] = &[
    // 金 地 火 木 土   v1x   v1y   v1z   v2x   v2y   v2z
    1, 1, 0, 0, 0, 0, -1719914, -25, 25, 1578089, 10, 684185, //01
    2, 1, 0, 0, 0, 0, 6434, 28007, 25697, -5904, 11141, -2559, 1, 3, 0, 0, 0, 0, 715, 0, 6, -657,
    -15, -282, 1, 7, 0, 0, 0, 0, 715, 0, 0, -656, 0, -285, 3, 1, 0, 0, 0, 0, 486, -236, -216, -446,
    -94, -193, 1, 4, 0, 0, 0, 0, 159, 0, 2, -147, -6, -61, //06
    1, 10, 0, 0, 0, 0, 0, 0, 0, 26, 0, -59, 1, 7, 1, 9, 0, 0, 39, 0, 0, -36, 0, -16, 2, 3, 0, 0, 0,
    0, 33, -10, -9, -30, -5, -13, 2, 1, -1, 3, 0, 0, 31, 1, 1, -28, 0, -12, 3, 1, -8, 2, 3, 3, 8,
    -28, 25, 8, 11, 3, //11
    5, 1, -8, 2, 3, 3, 8, -28, -25, -8, -11, -3, 2, 0, -1, 1, 0, 0, 21, 0, 0, -19, 0, -8, 1, 0, 0,
    0, 0, 0, -19, 0, 0, 17, 0, 8, 1, 5, 0, 0, 0, 0, 17, 0, 0, -16, 0, -7, 1, 1, -2, 3, 0, 0, 16, 0,
    0, 15, 1, 7, //16
    1, 6, 0, 0, 0, 0, 16, 0, 1, -15, -3, -6, 1, 1, 1, 3, 0, 0, 11, -1, -1, -10, -1, -5, 2, 0, -2,
    1, 0, 0, 0, -11, -10, 0, -4, 0, 1, 1, -1, 3, 0, 0, -11, -2, -2, 9, -1, 4, 4, 1, 0, 0, 0, 0, -7,
    -8, -8, 6, -3, 3, //21
    3, 1, -2, 3, 0, 0, -10, 0, 0, 9, 0, 4, 1, 0, -2, 1, 0, 0, -9, 0, 0, -9, 0, -4, 2, 0, -3, 1, 0,
    0, -9, 0, 0, -8, 0, -4, 2, 4, 0, 0, 0, 0, 0, -9, -8, 0, -3, 0, 2, 0, -4, 1, 0, 0, 0, -9, 8, 0,
    3, 0, //26
    3, 1, -2, 2, 0, 0, 8, 0, 0, -8, 0, -3, 1, 7, 2, 8, -1, 9, 8, 0, 0, -7, 0, -3, 8, 0, -12, 1, 0,
    0, -4, -7, -6, 4, -3, 2, 8, 0, -14, 1, 0, 0, -4, -7, 6, -4, 3, -2, 2, 2, 0, 0, 0, 0, -6, -5,
    -4, 5, -2, 2, //31
    3, 0, -4, 1, 0, 0, -1, -1, -2, -7, 1, -4, 2, 1, -2, 3, 0, 0, 4, -6, -5, -4, -2, -2, 3, 0, -3,
    1, 0, 0, 0, -7, -6, 0, -3, 0, 2, 1, -2, 2, 0, 0, 5, -5, -4, -5, -2, -2, 1, 7, -2, 8, 0, 0, 5,
    0, 0, -5, 0, -2, //36
    1, 1, 0, 0, 0, 0, -2, 0, -13, 156, 32, -358, // 泊松项开始
    2, 1, 0, 0, 0, 0, 141, -107, -95, -130, -48, -55, 3, 1, 0, 0, 0, 0, -5, -4, -4, 5, 0, 0,
];

/// 计算地球的SSB速度
/// # Arguments
/// * `t` - 儒略世纪数
/// # Returns
/// * `[f64; 3]` - 速度分量(AU/世纪),平均误差4*10^-8AU/日
pub fn ev_ssb(t: f64) -> [f64; 3] {
    // 天体平均轨道要素
    let j = [
        3.1761467 + 1021.3285546 * t, // 金星
        1.7534703 + 628.3075849 * t,  // 地球
        6.2034809 + 334.0612431 * t,  // 火星
        0.5995465 + 52.9690965 * t,   // 木星
        0.8740168 + 21.3299095 * t,   // 土星
        5.4812939 + 7.4781599 * t,    // 天王星
        5.3118863 + 3.8133036 * t,    // 海王星
        3.8103444 + 8399.6847337 * t, // 月球 L'
        5.1984667 + 7771.3771486 * t, // 月球 D
        2.3555559 + 8328.6914289 * t, // 月球 M'
        1.6279052 + 8433.4661601 * t, // 月球 F
    ];

    let mut v = [0.0, 0.0, 0.0]; // 速度分量
    let mut tn = 1.0; // 时间因子

    // 周期项计算
    for i in 0..39 {
        // 泊松项使用t
        if i >= 36 {
            tn = t;
        }

        // 计算角度
        let k = i * 12;
        let c = EV_TAB[k] as f64 * j[EV_TAB[k + 1] as usize]
            + EV_TAB[k + 2] as f64 * j[EV_TAB[k + 3] as usize]
            + EV_TAB[k + 4] as f64 * j[EV_TAB[k + 5] as usize];

        // 计算正弦和余弦值
        let (s, c) = c.sin_cos();

        // 累加速度分量
        v[0] += (EV_TAB[k + 6] as f64 * s + EV_TAB[k + 7] as f64 * c) * tn;
        v[1] += (EV_TAB[k + 8] as f64 * s + EV_TAB[k + 9] as f64 * c) * tn;
        v[2] += (EV_TAB[k + 10] as f64 * s + EV_TAB[k + 11] as f64 * c) * tn;
    }

    // 单位转换为AU/世纪
    const SCALE: f64 = 0.00036525;
    v[0] *= SCALE;
    v[1] *= SCALE;
    v[2] *= SCALE;

    v
}

/// 地心SSB位置计算表
const EP_TAB: &[(i32, f64, f64)] = &[
    // 振幅, 相位(弧度), 频率(弧度/世纪)
    (999829, 1.753486, 6283.07585),
    (999892, 0.182659, 6283.07585),
    (8353, 1.710345, 12566.1517),
    (24427, 3.141593, 0.0),
    (5611, 0.0, 0.0),
    (8353, 0.139529, 12566.1517),
    (105, 1.667226, 18849.22755),
    (105, 0.096417, 18849.22755),
    // ...其他周期项...
    // 木星引力项
    (-5196635 / 1048, 0.599451, 529.6909651),
    (-5195200 / 1048, 5.312032, 529.6909651),
    // 土星引力项
    (-9516383 / 3504, 0.874414, 213.2990954),
    (-9529869 / 3504, 5.586006, 213.2990954),
    // 海王星引力项
    (-30058900 / 19434, 5.312113, 38.1330356),
    (-30060560 / 19434, 3.740863, 38.1330356),
    // 天王星引力项
    (-19173710 / 22927, 5.481334, 74.7815986),
    (-19165180 / 22927, 3.910457, 74.7815986),
];

/// 计算地心SSB坐标
/// # Arguments
/// * `t` - 儒略世纪数
/// # Returns
/// * `[f64; 3]` - 坐标分量(单位:AU),精度4位有效数字
pub fn ep_ssb(t: f64) -> [f64; 3] {
    let t = t / 10.0;
    let mut x = 0.0;
    let mut y = 0.0;
    let mut z = 0.0;

    // 计算周期项
    for &(amp, phase, freq) in EP_TAB {
        let angle = phase + freq * t;
        let cos_angle = angle.cos();

        let term = amp as f64 * cos_angle;
        match EP_TAB.iter().position(|&p| p == (amp, phase, freq)) {
            Some(i) if i % 2 == 0 => x += term,
            Some(_) => y += term,
            None => (),
        }
    }

    // 添加时间项修正
    x += t
        * (1234.0
            + 515.0 * (6.002663 + 12566.1517 * t).cos()
            + 13.0 * (5.959431 + 18849.22755 * t).cos()
            + 11.0 * (2.015542 + 6283.07585 * t).cos());

    y += t
        * (930.0
            + 515.0 * (4.431805 + 12566.1517 * t).cos()
            + 13.0 * (4.388605 + 18849.22755 * t).cos());

    z += t
        * (54.0
            + 2278.0 * (3.413725 + 6283.07585 * t).cos()
            + 19.0 * (3.370613 + 12566.15170 * t).cos());

    // 单位转换
    x /= 1000000.0;
    y /= 1000000.0;
    z /= 1000000.0;

    // 黄赤转换
    let e = -84381.448 / RAD; // 黄赤交角
    [x, z * e.sin() + y * e.cos(), z * e.cos() - y * e.sin()]
}

/// 引力偏转计算
///
/// # Arguments
/// * `z` - 天体赤道坐标 [赤经, 赤纬, 距离]
/// * `a` - 太阳赤道坐标 [赤经, 赤纬, 距离]
///
/// # Returns
/// * 修正后的天体赤道坐标 [赤经, 赤纬, 距离]
pub fn ylpz(z: &[f64; 3], a: &[f64; 3]) -> [f64; 3] {
    let mut r = [z[0], z[1], z[2]];

    // 计算角度差
    let d = z[0] - a[0];

    // 角度差的余弦
    let d_cos = z[1].sin() * a[1].sin() + z[1].cos() * a[1].cos() * d.cos();

    // 引力偏转系数
    // let d_factor = 0.00407 * (1.0 / (1.0 - d_cos) + d_cos / 2.0) / RAD;
    // 错误：系数计算有误，应该检查分母是否为零
    let d_factor = if d_cos.abs() >= 0.9999999 {
        0.0 // 当天体与太阳重合或相对时，不计算引力偏转
    } else {
        0.00407 * (1.0 / (1.0 - d_cos) + d_cos / 2.0) / RAD
    };

    // 修正赤经
    // r[0] += d_factor * (a[1].cos() * d.sin() / z[1].cos());
    // 错误：当天体接近天极时，cos(赤纬)接近零，需要避免除以零
    if z[1].cos().abs() > 1e-9 {
        r[0] += d_factor * (a[1].cos() * d.sin() / z[1].cos());
    }

    // 修正赤纬
    r[1] += d_factor * (z[1].sin() * a[1].cos() * d.cos() - a[1].sin() * z[1].cos());

    // 标准化赤经到 [0, 2π)
    r[0] = rad2mrad(r[0]);

    r
}

//=================================低精度周年光行差=================================
/// 光行差计算常数
#[derive(Debug)]
pub struct GxcConst {
    /// 光行差常数
    pub k: f64,
    /// 太阳的几何平黄经
    pub l: f64,
    /// 近点
    pub p: f64,
    /// 离心率
    pub e: f64,
}
// impl SSBPosition {}

/// 取恒星光行差计算相关的常数
/// # Arguments
/// * `t` - 儒略世纪数
/// # Returns
/// * `GxcConst` - 光行差计算常数
pub fn get_gxc_const(t: f64) -> GxcConst {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    GxcConst {
        // 光行差常数
        k: 20.49552 / RAD,

        // 太阳的几何平黄经
        l: (280.4664567 + 36000.76982779 * t + 0.0003032028 * t2 + t3 / 49931000.0
            - t4 / 153000000.0
            - 5e-12 * t5)
            * PI
            / 180.0,

        // 近点
        p: (102.93735 + 1.71946 * t + 0.00046 * t2) * PI / 180.0,

        // 离心率
        e: 0.016708634 - 0.000042037 * t - 0.0000001267 * t2,
    }
}

/// 光行差黄道修正(恒星周年光行差)
/// # Arguments
/// * `t` - 儒略世纪数
/// * `a` - 天体黄道坐标[经度,纬度,距离]
/// # Returns
/// * `[f64; 3]` - 修正后的黄道坐标
pub fn hd_zn_gxc(t: f64, a: &[f64; 3]) -> [f64; 3] {
    let r = get_gxc_const(t);
    let dl = r.l - a[0];
    let dp = r.p - a[0];

    [
        // 黄经修正
        a[0] - r.k * (dl.cos() - r.e * dp.cos()) / a[1].cos(),
        // 黄纬修正
        a[1] - r.k * a[1].sin() * (dl.sin() - r.e * dp.sin()),
        // 距离不变
        a[2],
    ]
}

/// 计算周年光行差对赤道坐标的影响值
/// # Arguments
/// * `t` - 儒略世纪数(J2000起算)
/// * `a` - 天体赤道坐标[赤经,赤纬,距离]
/// # Returns
/// * `[f64; 3]` - 修正后的赤道坐标
pub fn cd_zn_gxc(t: f64, a: &[f64; 3]) -> [f64; 3] {
    let r = get_gxc_const(t);

    // 黄赤交角,取0.409也可以,对结果影响不很大
    let e = Prece::hcjj(t);

    let (sin_l, cos_l) = r.l.sin_cos();
    let cos_e = e.cos();
    let tan_e = e.tan();

    let (sin_r, cos_r) = a[0].sin_cos();
    let (sin_d, cos_d) = a[1].sin_cos();

    let tcss = tan_e * cos_d - sin_r * sin_d;

    // 修正值
    let mut b = [
        // 赤经周年光行差修正
        a[0] - r.k * (cos_r * cos_e * cos_l + sin_r * sin_l) / cos_d,
        // 赤纬周年光行差修正
        a[1] - r.k * (cos_l * cos_e * tcss + cos_r * sin_d * sin_l),
        // 距离不变
        a[2],
    ];

    // e项修正(考虑非正圆运动的误差)
    let (sin_p, cos_p) = r.p.sin_cos();
    b[0] += r.e * r.k * (cos_r * cos_e * cos_p + sin_r * sin_p) / cos_d;
    b[1] += r.e * r.k * (cos_p * cos_e * tcss + cos_r * sin_d * sin_p);

    b
}

//=================================周年视差或光行差=================================
/// 严格的恒星视差或光行差改正
/// # Arguments
/// * `z` - 某时刻天体赤道坐标(球面坐标),含自行但不含章动和光行差 [α,δ,r]
/// * `v` - 同时刻地球赤道坐标(直角坐标) [x,y,z]
/// * `f` - 修正类型:
///   - false: 进行光行差修正,v须为SSB速度,返回z+v向量(z向径为光速)
///   - true: 做周年视差修正,v须为SSB位置,返回z-v向量(z向径为距离)
/// # Returns
/// * `[f64; 3]` - 修正后的赤道坐标
/// # Notes
/// z和v应统一使用J2000赤道坐标系
pub fn sc_gxc(z: &[f64; 3], v: &[f64; 3], f: bool) -> [f64; 3] {
    let mut r = *z;

    // 光速(AU每儒略世纪)或视差因子
    let c = if f == false {
        CS_GS / CS_AU * 86400.0 * 36525.0 // 光速
    } else {
        -z[2] // 视差因子
    };

    // 计算三角函数值
    let (sin_j, cos_j) = z[0].sin_cos();
    let (sin_w, cos_w) = z[1].sin_cos();

    // 赤经修正
    r[0] += rad2mrad((v[1] * cos_j - v[0] * sin_j) / cos_w / c);

    // 赤纬修正
    r[1] += (v[2] * cos_w - (v[0] * cos_j + v[1] * sin_j) * sin_w) / c;

    r
}

/// 一次计算周年视差和光行差
/// (AI生成的)
/// # Arguments
/// * `z` - 天体赤道坐标 [α,δ,r]
/// * `t` - 儒略世纪数
/// # Returns
/// * `[f64; 3]` - 修正后的赤道坐标
pub fn annual_aberration(z: &[f64; 3], t: f64) -> [f64; 3] {
    // 获取地球SSB位置和速度
    let p = ep_ssb(t); // SSB位置
    let v = ev_ssb(t); // SSB速度

    // 先做视差修正
    let mut r = sc_gxc(z, &p, true);

    // 再做光行差修正
    r = sc_gxc(&r, &v, false);

    r
}

/// 计算太阳J2000球面坐标
/// # Arguments
/// * `t` - 儒略世纪数  
/// * `prec` - 计算精度,取20表示精确到0.1角秒
/// # Returns
/// * `[f64; 3]` - 太阳J2000坐标[经度,纬度,距离]
pub fn sun_2000(t: f64, prec: i32) -> [f64; 3] {
    // 获取地球Date坐标(相对于瞬时黄道)
    let mut a = XL::e_coord(t, prec, prec, prec);

    // 计算太阳坐标
    a[0] += PI; // 太阳黄经=地球黄经+180度
    a[1] = -a[1]; // 太阳黄纬=地球黄纬取反

    // 转换到J2000坐标系
    let a = Prece::hdllr_d2j(t, &a, "P03");

    [a[0], a[1], a[2]]
}

//=================================恒星星历计算=====================================

/// 星座信息
#[derive(Debug)]
pub struct Constellation {
    /// 中文名
    pub name_cn: &'static str,
    /// 缩写
    pub abbr: &'static str,
    /// 面积(平方度)
    pub area: f64,
    /// 赤经(时分)
    pub ra: &'static str,
    /// 赤纬(度分)
    pub dec: &'static str,
    /// 象限角
    pub quadrant: &'static str,
    /// 星座族
    pub family: &'static str,
    /// 英文名
    pub name_en: &'static str,
}

/// 88星座数据
pub const CONSTELLATIONS: &[Constellation] = &[
    Constellation {
        name_cn: "仙女座",
        abbr: "And",
        area: 722.278,
        ra: "0 48.46",
        dec: "37 25.91",
        quadrant: "NQ1",
        family: "英仙",
        name_en: "Andromeda",
    },
    Constellation {
        name_cn: "唧筒座",
        abbr: "Ant",
        area: 238.901,
        ra: "10 16.43",
        dec: "-32 29.01",
        quadrant: "SQ2",
        family: "拉卡伊",
        name_en: "Antlia",
    },
    Constellation {
        name_cn: "唧筒座",
        abbr: "Ant",
        area: 238.901,
        ra: "10 16.43",
        dec: "-32 29.01",
        quadrant: "SQ2",
        family: "拉卡伊",
        name_en: "Antlia",
    },
    Constellation {
        name_cn: "天燕座",
        abbr: "APS",
        area: 206.327,
        ra: "16 08.65",
        dec: "-75 18.00",
        quadrant: "SQ3",
        family: "拜耳",
        name_en: "Apus",
    },
    Constellation {
        name_cn: "宝瓶座",
        abbr: "Aqr",
        area: 979.854,
        ra: "22 17.38",
        dec: "-10 47.35",
        quadrant: "SQ4",
        family: "黄道",
        name_en: "Aquarius",
    },
    Constellation {
        name_cn: "天鹰座",
        abbr: "Aql",
        area: 652.473,
        ra: "19 40.02",
        dec: "3 24.65",
        quadrant: "NQ4",
        family: "武仙",
        name_en: "Aquila",
    },
    Constellation {
        name_cn: "天坛座",
        abbr: "Ara",
        area: 237.057,
        ra: "17 22.49",
        dec: "-56 35.30",
        quadrant: "SQ3",
        family: "武仙",
        name_en: "Ara",
    },
    Constellation {
        name_cn: "白羊座",
        abbr: "Ari",
        area: 441.395,
        ra: "2 38.16",
        dec: "20 47.54",
        quadrant: "NQ1",
        family: "黄道",
        name_en: "Aries",
    },
    Constellation {
        name_cn: "御夫座",
        abbr: "Aur",
        area: 657.438,
        ra: "6 04.42",
        dec: "42 01.68",
        quadrant: "NQ2",
        family: "英仙",
        name_en: "Auriga",
    },
    Constellation {
        name_cn: "牧夫座",
        abbr: "Boo",
        area: 906.831,
        ra: "14 42.64",
        dec: "31 12.16",
        quadrant: "NQ3",
        family: "大熊",
        name_en: "Bootes",
    },
    Constellation {
        name_cn: "雕具座",
        abbr: "Cae",
        area: 124.865,
        ra: "4 42.27",
        dec: "-37 52.90",
        quadrant: "SQ1",
        family: "拉卡伊",
        name_en: "Caelum",
    },
    Constellation {
        name_cn: "鹿豹座",
        abbr: "Cam",
        area: 756.828,
        ra: "8 51.37",
        dec: "69 22.89",
        quadrant: "NQ2",
        family: "大熊",
        name_en: "Camelopardalis",
    },
    Constellation {
        name_cn: "巨蟹座",
        abbr: "Cnc",
        area: 505.872,
        ra: "8 38.96",
        dec: "19 48.35",
        quadrant: "NQ2",
        family: "黄道",
        name_en: "Cancer",
    },
    Constellation {
        name_cn: "猎犬座",
        abbr: "CVn",
        area: 465.194,
        ra: "13 06.96",
        dec: "40 06.11",
        quadrant: "NQ3",
        family: "大熊",
        name_en: "Canes Venatici",
    },
    Constellation {
        name_cn: "大犬座",
        abbr: "CMa",
        area: 380.118,
        ra: "6 49.74",
        dec: "-22 08.42",
        quadrant: "SQ2",
        family: "猎户",
        name_en: "Canis Major",
    },
    Constellation {
        name_cn: "小犬座",
        abbr: "CMi",
        area: 183.367,
        ra: "7 39.17",
        dec: "6 25.63",
        quadrant: "NQ2",
        family: "猎户",
        name_en: "Canis Minor",
    },
    Constellation {
        name_cn: "摩羯座",
        abbr: "CAP",
        area: 413.947,
        ra: "21 02.93",
        dec: "-18 01.39",
        quadrant: "SQ4",
        family: "黄道",
        name_en: "Capricornus",
    },
    Constellation {
        name_cn: "船底座",
        abbr: "Car",
        area: 494.184,
        ra: "8 41.70",
        dec: "-63 13.16",
        quadrant: "SQ2",
        family: "幻之水",
        name_en: "Carina",
    },
    Constellation {
        name_cn: "仙后座",
        abbr: "Cas",
        area: 598.407,
        ra: "1 19.16",
        dec: "62 11.04",
        quadrant: "NQ1",
        family: "英仙",
        name_en: "Cassiopeia",
    },
    Constellation {
        name_cn: "半人马座",
        abbr: "Cen",
        area: 1060.422,
        ra: "13 04.27",
        dec: "-47 20.72",
        quadrant: "SQ3",
        family: "武仙",
        name_en: "Centaurus",
    },
    Constellation {
        name_cn: "仙王座",
        abbr: "Cep",
        area: 587.787,
        ra: "22 00.00",
        dec: "71 00.51",
        quadrant: "NQ4",
        family: "英仙",
        name_en: "Cepheus",
    },
    Constellation {
        name_cn: "鲸鱼座",
        abbr: "Cet",
        area: 1231.411,
        ra: "1 40.10",
        dec: "-7 10.76",
        quadrant: "SQ1",
        family: "英仙",
        name_en: "Cetus",
    },
    Constellation {
        name_cn: "堰蜒座",
        abbr: "Cha",
        area: 131.592,
        ra: "10 41.53",
        dec: "-79 12.30",
        quadrant: "SQ2",
        family: "拜耳",
        name_en: "Chamaeleon",
    },
    Constellation {
        name_cn: "圆规座",
        abbr: "Cir",
        area: 93.353,
        ra: "14 34.54",
        dec: "-63 01.82",
        quadrant: "SQ3",
        family: "拉卡伊",
        name_en: "Circinus",
    },
    Constellation {
        name_cn: "天鸽座",
        abbr: "Col",
        area: 270.184,
        ra: "5 51.76",
        dec: "-35 05.67",
        quadrant: "SQ1",
        family: "幻之水",
        name_en: "Columba",
    },
    Constellation {
        name_cn: "后发座",
        abbr: "Com",
        area: 386.475,
        ra: "12 47.27",
        dec: "23 18.34",
        quadrant: "NQ3",
        family: "大熊",
        name_en: "Coma Berenices",
    },
    Constellation {
        name_cn: "南冕座",
        abbr: "CrA",
        area: 127.696,
        ra: "18 38.79",
        dec: "-41 08.85",
        quadrant: "SQ4",
        family: "武仙",
        name_en: "Corona Australis",
    },
    Constellation {
        name_cn: "北冕座",
        abbr: "CrB",
        area: 178.710,
        ra: "15 50.59",
        dec: "32 37.49",
        quadrant: "NQ3",
        family: "大熊",
        name_en: "Corona Borealis",
    },
    Constellation {
        name_cn: "乌鸦座",
        abbr: "Crv",
        area: 183.801,
        ra: "12 26.52",
        dec: "-18 26.20",
        quadrant: "SQ3",
        family: "武仙",
        name_en: "Corvus",
    },
    Constellation {
        name_cn: "巨爵座",
        abbr: "Crt",
        area: 282.398,
        ra: "11 23.75",
        dec: "-15 55.74",
        quadrant: "SQ2",
        family: "武仙",
        name_en: "Crater",
    },
    Constellation {
        name_cn: "南十字座",
        abbr: "Cru",
        area: 68.447,
        ra: "12 26.99",
        dec: "-60 11.19",
        quadrant: "SQ3",
        family: "武仙",
        name_en: "Crux",
    },
    Constellation {
        name_cn: "天鹅座",
        abbr: "Cyg",
        area: 803.983,
        ra: "20 35.28",
        dec: "44 32.70",
        quadrant: "NQ4",
        family: "武仙",
        name_en: "Cygnus",
    },
    Constellation {
        name_cn: "海豚座",
        abbr: "Del",
        area: 188.549,
        ra: "20 41.61",
        dec: "11 40.26",
        quadrant: "NQ4",
        family: "幻之水",
        name_en: "Delphinus",
    },
    Constellation {
        name_cn: "剑鱼座",
        abbr: "Dor",
        area: 179.173,
        ra: "5 14.51",
        dec: "-59 23.22",
        quadrant: "SQ1",
        family: "拜耳",
        name_en: "Dorado",
    },
    Constellation {
        name_cn: "天龙座",
        abbr: "Dra",
        area: 1082.952,
        ra: "15 08.64",
        dec: "67 00.40",
        quadrant: "NQ3",
        family: "大熊",
        name_en: "Draco",
    },
    Constellation {
        name_cn: "小马座",
        abbr: "Equ",
        area: 71.641,
        ra: "21 11.26",
        dec: "7 45.49",
        quadrant: "NQ4",
        family: "幻之水",
        name_en: "Equuleus",
    },
    Constellation {
        name_cn: "波江座",
        abbr: "Eri",
        area: 1137.919,
        ra: "3 18.02",
        dec: "-28 45.37",
        quadrant: "SQ1",
        family: "幻之水",
        name_en: "Eridanus",
    },
    Constellation {
        name_cn: "天炉座",
        abbr: "For",
        area: 397.502,
        ra: "2 47.88",
        dec: "-31 38.07",
        quadrant: "SQ1",
        family: "拉卡伊",
        name_en: "Fornax",
    },
    Constellation {
        name_cn: "双子座",
        abbr: "Gem",
        area: 513.761,
        ra: "7 04.24",
        dec: "22 36.01",
        quadrant: "NQ2",
        family: "黄道",
        name_en: "Gemini",
    },
    Constellation {
        name_cn: "天鹤座",
        abbr: "Gru",
        area: 365.513,
        ra: "22 27.39",
        dec: "-46 21.11",
        quadrant: "SQ4",
        family: "拜耳",
        name_en: "Grus",
    },
    Constellation {
        name_cn: "武仙座",
        abbr: "Her",
        area: 1225.148,
        ra: "17 23.16",
        dec: "27 29.93",
        quadrant: "NQ3",
        family: "武仙",
        name_en: "Hercules",
    },
    Constellation {
        name_cn: "时钟座",
        abbr: "Hor",
        area: 248.885,
        ra: "3 16.56",
        dec: "-53 20.18",
        quadrant: "SQ1",
        family: "拉卡伊",
        name_en: "Horologium",
    },
    Constellation {
        name_cn: "长蛇座",
        abbr: "Hya",
        area: 1302.844,
        ra: "11 36.73",
        dec: "-14 31.91",
        quadrant: "SQ2",
        family: "武仙",
        name_en: "Hydra",
    },
    Constellation {
        name_cn: "水蛇座",
        abbr: "Hyi",
        area: 243.035,
        ra: "2 20.65",
        dec: "-69 57.39",
        quadrant: "SQ1",
        family: "拜耳",
        name_en: "Hydrus",
    },
    Constellation {
        name_cn: "印第安座",
        abbr: "Ind",
        area: 294.006,
        ra: "21 58.33",
        dec: "-59 42.40",
        quadrant: "SQ4",
        family: "拜耳",
        name_en: "Indus",
    },
    Constellation {
        name_cn: "蝎虎座",
        abbr: "Lac",
        area: 200.688,
        ra: "22 27.68",
        dec: "46 02.51",
        quadrant: "NQ4",
        family: "英仙",
        name_en: "Lacerta",
    },
    Constellation {
        name_cn: "狮子座",
        abbr: "Leo",
        area: 946.964,
        ra: "10 40.03",
        dec: "13 08.32",
        quadrant: "NQ2",
        family: "黄道",
        name_en: "Leo",
    },
    Constellation {
        name_cn: "小狮座",
        abbr: "LMi",
        area: 231.956,
        ra: "10 14.72",
        dec: "32 08.08",
        quadrant: "NQ2",
        family: "大熊",
        name_en: "Leo Minor",
    },
    Constellation {
        name_cn: "天兔座",
        abbr: "Lep",
        area: 290.291,
        ra: "5 33.95",
        dec: "-19 02.78",
        quadrant: "SQ1",
        family: "猎户",
        name_en: "Lepus",
    },
    Constellation {
        name_cn: "天秤座",
        abbr: "Lib",
        area: 538.052,
        ra: "15 11.96",
        dec: "-15 14.08",
        quadrant: "SQ3",
        family: "黄道",
        name_en: "Libra",
    },
    Constellation {
        name_cn: "豺狼座",
        abbr: "Lup",
        area: 333.683,
        ra: "15 13.21",
        dec: "-42 42.53",
        quadrant: "SQ3",
        family: "武仙",
        name_en: "Lupus",
    },
    Constellation {
        name_cn: "天猫座",
        abbr: "Lyn",
        area: 545.386,
        ra: "7 59.53",
        dec: "47 28.00",
        quadrant: "NQ2",
        family: "大熊",
        name_en: "Lynx",
    },
    Constellation {
        name_cn: "天琴座",
        abbr: "Lyr",
        area: 286.476,
        ra: "18 51.17",
        dec: "36 41.36",
        quadrant: "NQ4",
        family: "武仙",
        name_en: "Lyra",
    },
    Constellation {
        name_cn: "山案座",
        abbr: "Men",
        area: 153.484,
        ra: "5 24.90",
        dec: "-77 30.24",
        quadrant: "SQ1",
        family: "拉卡伊",
        name_en: "Mensa",
    },
    Constellation {
        name_cn: "显微镜座",
        abbr: "Mic",
        area: 209.513,
        ra: "20 57.88",
        dec: "-36 16.49",
        quadrant: "SQ4",
        family: "拉卡伊",
        name_en: "Microscopium",
    },
    Constellation {
        name_cn: "麒麟座",
        abbr: "Mon",
        area: 481.569,
        ra: "7 03.63",
        dec: "0 16.93",
        quadrant: "NQ2",
        family: "猎户",
        name_en: "Monoceros",
    },
    Constellation {
        name_cn: "苍蝇座",
        abbr: "Mus",
        area: 138.355,
        ra: "12 35.28",
        dec: "-70 09.66",
        quadrant: "SQ3",
        family: "拜耳",
        name_en: "Musca",
    },
    Constellation {
        name_cn: "矩尺座",
        abbr: "Nor",
        area: 165.290,
        ra: "15 54.18",
        dec: "-51 21.09",
        quadrant: "SQ3",
        family: "拉卡伊",
        name_en: "Norma",
    },
    Constellation {
        name_cn: "南极座",
        abbr: "Oct",
        area: 291.045,
        ra: "23 00.00",
        dec: "-82 09.12",
        quadrant: "SQ4",
        family: "拉卡伊",
        name_en: "Octans",
    },
    Constellation {
        name_cn: "蛇夫座",
        abbr: "Oph",
        area: 948.340,
        ra: "17 23.69",
        dec: "-7 54.74",
        quadrant: "SQ3",
        family: "武仙",
        name_en: "Ophiuchus",
    },
    Constellation {
        name_cn: "猎户座",
        abbr: "Ori",
        area: 594.120,
        ra: "5 34.59",
        dec: "5 56.94",
        quadrant: "NQ1",
        family: "猎户",
        name_en: "Orion",
    },
    Constellation {
        name_cn: "孔雀座",
        abbr: "Pav",
        area: 377.666,
        ra: "19 36.71",
        dec: "-65 46.89",
        quadrant: "SQ4",
        family: "拜耳",
        name_en: "Pavo",
    },
    Constellation {
        name_cn: "飞马座",
        abbr: "Peg",
        area: 1120.794,
        ra: "22 41.84",
        dec: "19 27.98",
        quadrant: "NQ4",
        family: "英仙",
        name_en: "Pegasus",
    },
    Constellation {
        name_cn: "英仙座",
        abbr: "Per",
        area: 614.997,
        ra: "3 10.50",
        dec: "45 00.79",
        quadrant: "NQ1",
        family: "英仙",
        name_en: "Perseus",
    },
    Constellation {
        name_cn: "凤凰座",
        abbr: "Phe",
        area: 469.319,
        ra: "0 55.91",
        dec: "-48 34.84",
        quadrant: "SQ1",
        family: "拜耳",
        name_en: "Phoenix",
    },
    Constellation {
        name_cn: "绘架座",
        abbr: "Pic",
        area: 246.739,
        ra: "5 42.46",
        dec: "-53 28.45",
        quadrant: "SQ1",
        family: "拉卡伊",
        name_en: "Pictor",
    },
    Constellation {
        name_cn: "双鱼座",
        abbr: "Psc",
        area: 889.417,
        ra: "0 28.97",
        dec: "13 41.23",
        quadrant: "NQ1",
        family: "黄道",
        name_en: "Pisces",
    },
    Constellation {
        name_cn: "南鱼座",
        abbr: "PsA",
        area: 245.375,
        ra: "22 17.07",
        dec: "-30 38.53",
        quadrant: "SQ4",
        family: "幻之水",
        name_en: "Piscis Austrinus",
    },
    Constellation {
        name_cn: "船尾座",
        abbr: "Pup",
        area: 673.434,
        ra: "7 15.48",
        dec: "-31 10.64",
        quadrant: "SQ2",
        family: "幻之水",
        name_en: "Puppis",
    },
    Constellation {
        name_cn: "罗盘座",
        abbr: "Pyx",
        area: 220.833,
        ra: "8 57.16",
        dec: "-27 21.10",
        quadrant: "SQ2",
        family: "幻之水",
        name_en: "Pyxis",
    },
    Constellation {
        name_cn: "网罟座",
        abbr: "Ret",
        area: 113.936,
        ra: "3 55.27",
        dec: "-59 59.85",
        quadrant: "SQ1",
        family: "拉卡伊",
        name_en: "Reticulum",
    },
    Constellation {
        name_cn: "天箭座",
        abbr: "Sge",
        area: 79.932,
        ra: "19 39.05",
        dec: "18 51.68",
        quadrant: "NQ4",
        family: "武仙",
        name_en: "Sagitta",
    },
    Constellation {
        name_cn: "人马座",
        abbr: "Sgr",
        area: 867.432,
        ra: "19 05.94",
        dec: "-28 28.61",
        quadrant: "SQ4",
        family: "黄道",
        name_en: "Sagittarius",
    },
    Constellation {
        name_cn: "天蝎座",
        abbr: "Sco",
        area: 496.783,
        ra: "16 53.24",
        dec: "-27 01.89",
        quadrant: "SQ3",
        family: "黄道",
        name_en: "Scorpius",
    },
    Constellation {
        name_cn: "玉夫座",
        abbr: "Scl",
        area: 474.764,
        ra: "0 26.28",
        dec: "-32 05.30",
        quadrant: "SQ1",
        family: "拉卡伊",
        name_en: "Sculptor",
    },
    Constellation {
        name_cn: "盾牌座",
        abbr: "Sct",
        area: 109.114,
        ra: "18 40.39",
        dec: "-9 53.32",
        quadrant: "SQ4",
        family: "武仙",
        name_en: "Scutum",
    },
    Constellation {
        name_cn: "巨蛇座",
        abbr: "Ser",
        area: 636.928,
        ra: "16 57.04",
        dec: "6 07.32",
        quadrant: "NQ3",
        family: "武仙",
        name_en: "Serpens",
    },
    Constellation {
        name_cn: "六分仪座",
        abbr: "Sex",
        area: 313.515,
        ra: "10 16.29",
        dec: "-2 36.88",
        quadrant: "SQ2",
        family: "武仙",
        name_en: "Sextans",
    },
    Constellation {
        name_cn: "金牛座",
        abbr: "Tau",
        area: 797.249,
        ra: "4 42.13",
        dec: "14 52.63",
        quadrant: "NQ1",
        family: "黄道",
        name_en: "Taurus",
    },
    Constellation {
        name_cn: "望远镜座",
        abbr: "Tel",
        area: 251.512,
        ra: "19 19.54",
        dec: "-51 02.21",
        quadrant: "SQ4",
        family: "拉卡伊",
        name_en: "Telescopium",
    },
    Constellation {
        name_cn: "三角座",
        abbr: "Tri",
        area: 131.847,
        ra: "2 11.07",
        dec: "31 28.56",
        quadrant: "NQ1",
        family: "英仙",
        name_en: "Triangulum",
    },
    Constellation {
        name_cn: "南三角座",
        abbr: "TrA",
        area: 109.978,
        ra: "16 04.95",
        dec: "-65 23.28",
        quadrant: "SQ3",
        family: "武仙",
        name_en: "Triangulum Australe",
    },
    Constellation {
        name_cn: "杜鹃座",
        abbr: "Tuc",
        area: 294.557,
        ra: "23 46.64",
        dec: "-65 49.80",
        quadrant: "SQ4",
        family: "拜耳",
        name_en: "Tucana",
    },
    Constellation {
        name_cn: "大熊座",
        abbr: "UMa",
        area: 1279.660,
        ra: "11 18.76",
        dec: "50 43.27",
        quadrant: "NQ2",
        family: "大熊",
        name_en: "Ursa Major",
    },
    Constellation {
        name_cn: "小熊座",
        abbr: "UMi",
        area: 255.864,
        ra: "15 00.00",
        dec: "77 41.99",
        quadrant: "NQ3",
        family: "大熊",
        name_en: "Ursa Minor",
    },
    Constellation {
        name_cn: "船帆座",
        abbr: "Vel",
        area: 499.649,
        ra: "9 34.64",
        dec: "-47 10.03",
        quadrant: "SQ2",
        family: "幻之水",
        name_en: "Vela",
    },
    Constellation {
        name_cn: "室女座",
        abbr: "Vir",
        area: 1294.428,
        ra: "13 24.39",
        dec: "-4 09.51",
        quadrant: "SQ3",
        family: "黄道",
        name_en: "Virgo",
    },
    Constellation {
        name_cn: "飞鱼座",
        abbr: "Vol",
        area: 141.354,
        ra: "7 47.73",
        dec: "-69 48.07",
        quadrant: "SQ2",
        family: "拜耳",
        name_en: "Volans",
    },
    Constellation {
        name_cn: "狐狸座",
        abbr: "Vul",
        area: 268.165,
        ra: "20 13.88",
        dec: "24 26.56",
        quadrant: "NQ4",
        family: "武仙",
        name_en: "Vulpecula",
    },
];

/// 恒星库数据
pub const HXK: &[&str] = &[
    // 库0
    "库0#* 0 01 57.620,- 6 00 50.68, 0.0031, -0.041, 0.008, 4.37 ,星1630 ,Psc 30 M3#* 0 03 44.391,-17 20 09.58, 0.0020, -0.007, 0.014, 4.55 ,星905  ,Cet 2  B9#* 0 05 20.142,- 5 42 27.45,-0.0009,  0.089, 0.025, 4.61 ,星1002 ,Psc 33 K1#* 0 08 23.260, 29 05 25.54, 0.0104, -0.163, 0.034, 2.07 ,星1    ,And α B9#* 0 09 10.686, 59 08 59.19, 0.0681, -0.180, 0.060, 2.28 ,星2    ,Cas β F2#* 0 10 19.247, 46 04 20.17, 0.0005,  0.001, 0.003, 5.01 ,星4    ,And 22 F2#* 0 11 34.421,-27 47 59.06, 0.0003,  0.016, 0.006, 5.41 ,星5    ,Scl κ2 K2#* 0 11 44.010,-35 07 59.24, 0.0138,  0.115, 0.046, 5.24 ,星6    ,Scl θ F3#* 0 13 14.154, 15 11 00.93, 0.0003, -0.008, 0.010, 2.83 ,星7    ,Peg γ B2#* 0 14 36.165, 20 12 24.12, 0.0064, -0.001, 0.010, 4.79 ,星1004 ,Peg χ M2#* 0 17 05.500, 38 40 53.87,-0.0046, -0.013, 0.013, 4.61 ,N30    ,And θ A2#* 0 18 19.658, 36 47 06.79,-0.0055, -0.042, 0.023, 4.51 ,星1005 ,And σ A2#* 0 18 38.258, 31 31 02.01, 0.0044, -0.004, 0.006, 5.88 ,星1006 ,Pi 0h38 A0#* 0 19 25.676,- 8 49 26.14,-0.0010, -0.037, 0.011, 3.56 ,星9    ,Cet ι K2#* 0 20 35.863,  8 11 24.96,-0.0003,  0.010, 0.008, 5.38 ,星1008 ,Psc 41 K3#* 0 21 07.270, 37 58 06.95, 0.0049, -0.039, 0.020, 5.16 ,星1009 ,And ρ F5#* 0 24 47.506, 61 49 51.80, 0.0018, -0.002, 0.004, 5.38 ,GC     ,Cas 12 B9#* 0 25 24.210,  1 56 22.87,-0.0010, -0.013, 0.006, 5.77 ,星1010 ,Psc 44 G5#* 0 25 45.092,-77 15 15.30, 0.6689,  0.323, 0.134, 2.82 ,星11   ,Hyi β G2#* 0 26 17.052,-42 18 21.55, 0.0210, -0.354, 0.042, 2.40 ,星12   ,Phe α K0#",
    // 库1
    "库1#* 0 48 22.978,  5 16 50.19,  0.0507,-1.141,0.134, 5.74,星1019,G.Psc 96 K2#* 0 26 17.052,-42 18 21.55,  0.0210,-0.354,0.042, 2.40,星12  ,Phe α   K0#* 2 36 00.049,- 7 49 53.77, -0.0022,-0.060,0.006, 5.53,星1074,Cet 80   M0#* 2 35 52.472,  5 35 35.67, -0.0019,-0.024,0.009, 4.87,星1072,Cet υ   G8#*18 36 27.834,  9 07 20.98, -0.0001,-0.132,0.026, 5.38,星1484,Oph 9  F5#*18 36 56.338, 38 47 01.29,  0.0172, 0.288,0.129, 0.03,星699, Lyr α A0#*18 37 54.426,-21 23 51.81, -0.0001,-0.066,0.010, 5.93,星1485,Sgr 83 A5#*18 42 16.428,- 9 03 09.14,  0.0005,-0.002,0.017, 4.70,星1486,Sct δ F2#*18 43 31.254,- 8 16 30.76,  0.0013, 0.008,0.006, 4.88,星702, Sct ε G8#",
];

/// 星座库检索
///
/// # Arguments
/// * `key` - 星座缩写，如 "And", "Psc" 等
///
/// # Returns
/// * 包含该星座信息的字符串
pub fn sch_hxk(key: &str) -> String {
    let mut result = String::new();

    // 遍历所有子库
    for &catalog in HXK {
        let n0 = catalog.find('#').unwrap_or(0);
        let mut n1 = n0;

        loop {
            n1 = match catalog[n1 + 1..].find(key) {
                Some(pos) => n1 + 1 + pos,
                None => break,
            };

            // 找到该记录的结束位置
            let n2 = match catalog[n1..].find('#') {
                Some(pos) => n1 + pos,
                None => catalog.len(),
            };

            // 找到该记录的开始位置
            let n3 = match catalog[..n1].rfind('#') {
                Some(pos) => pos,
                None => n0,
            };

            result.push_str(&catalog[n3..n2]);
        }
    }

    // 提取星座中心位置
    for constellation in CONSTELLATIONS {
        if constellation.abbr == key {
            let mut a = constellation.ra.to_string();
            let mut b = constellation.dec.to_string();

            // 处理赤经
            let ra_parts: Vec<&str> = a.split_whitespace().collect();
            if ra_parts.len() >= 2 {
                let minutes = ra_parts[1].parse::<f64>().unwrap_or(0.0) * 0.6;
                a = format!("{} {:.1}", ra_parts[0], minutes);
            }

            // 处理赤纬
            let dec_parts: Vec<&str> = b.split_whitespace().collect();
            if dec_parts.len() >= 2 {
                let minutes = dec_parts[1].parse::<f64>().unwrap_or(0.0) * 0.6;
                b = format!("{} {:.1}", dec_parts[0], minutes);
            }

            let center_info = format!(
                "#*{},{},0,0,0,0.0,中心{}方,{}",
                a, b, constellation.name_cn, constellation.name_en
            );

            result = center_info + &result;
            break;
        }
    }

    result
}

/// 恒星数据结构
#[derive(Debug, Clone)]
pub struct StarData {
    pub ra: f64,          // 赤经(弧度)
    pub dec: f64,         // 赤纬(弧度)
    pub pm_ra: f64,       // 赤经世纪自行(弧度/世纪)
    pub pm_dec: f64,      // 赤纬世纪自行(弧度/世纪)
    pub parallax: f64,    // 视差(弧度)
    pub magnitude: f64,   // 星等
    pub info: String,     // 星座等信息
    pub spectrum: String, // 光谱类型
}

/// 提取并格式化恒星库
///
/// # Arguments
/// * `s` - 恒星库字符串
/// * `all` - 是否提取全部数据，true表示全部取出
///
/// # Returns
/// * 恒星数据向量
pub fn get_hxk(s: &str, all: bool) -> Vec<f64> {
    let mut result = Vec::new();

    // 预处理字符串
    let s = s.replace('\r', "");
    let s = s.replace('\n', "#");

    // 去除第一行
    let s = match s.find('#') {
        Some(pos) => &s[pos + 1..],
        None => return result,
    };

    // 替换分隔符
    let s = s.replace("#+", ",");
    let s = s.replace(", +", ",");

    // 分割数据
    let parts: Vec<&str> = s.split(',').collect();

    // let mut k = 0;
    let mut i = 0;

    while i + 7 < parts.len() {
        if parts[i].is_empty() || parts[i].len() < 5 {
            i += 8;
            continue;
        }

        // 如果不是星号开头且不是全部提取，则跳过
        if !parts[i].starts_with('*') && !all {
            i += 8;
            continue;
        }

        // 去除星号
        let star_data = parts[i].trim_start_matches('*');

        // 转换数据
        result.push(str2rad(star_data, true)); // 赤经
        result.push(str2rad(parts[i + 1], false)); // 赤纬
        result.push(parts[i + 2].parse::<f64>().unwrap_or(0.0) / (PI / 180.0) * 15.0); // 赤经世纪自行
        result.push(parts[i + 3].parse::<f64>().unwrap_or(0.0) / (PI / 180.0)); // 赤纬世纪自行
        result.push(parts[i + 4].parse::<f64>().unwrap_or(0.0) / (PI / 180.0)); // 视差
        result.push(parts[i + 5].parse::<f64>().unwrap_or(0.0)); // 星等
        result.push(0.0); // 占位，原JS中存储字符串，这里用0.0占位
        result.push(0.0); // 占位，原JS中存储字符串，这里用0.0占位

        // k += 8;
        i += 8;
    }

    result
}

/// 提取并格式化恒星库(返回结构体版本)
pub fn get_hxk_struct(s: &str, all: bool) -> Vec<StarData> {
    let mut result = Vec::new();

    // 预处理字符串
    let s = s.replace('\r', "");
    let s = s.replace('\n', "#");

    // 去除第一行
    let s = match s.find('#') {
        Some(pos) => &s[pos + 1..],
        None => return result,
    };

    // 替换分隔符
    let s = s.replace("#+", ",");
    let s = s.replace(", +", ",");

    // 分割数据
    let parts: Vec<&str> = s.split(',').collect();

    let mut i = 0;

    while i + 7 < parts.len() {
        if parts[i].is_empty() || parts[i].len() < 5 {
            i += 8;
            continue;
        }

        // 如果不是星号开头且不是全部提取，则跳过
        if !parts[i].starts_with('*') && !all {
            i += 8;
            continue;
        }

        // 去除星号
        let star_data = parts[i].trim_start_matches('*');

        // 转换数据
        let star = StarData {
            ra: str2rad(star_data, true),
            dec: str2rad(parts[i + 1], false),
            pm_ra: parts[i + 2].parse::<f64>().unwrap_or(0.0) / (PI / 180.0) * 15.0,
            pm_dec: parts[i + 3].parse::<f64>().unwrap_or(0.0) / (PI / 180.0),
            parallax: parts[i + 4].parse::<f64>().unwrap_or(0.0) / (PI / 180.0),
            magnitude: parts[i + 5].parse::<f64>().unwrap_or(0.0),
            info: parts[i + 6].to_string(),
            spectrum: parts[i + 7].trim_end_matches('#').to_string(),
        };

        result.push(star);
        i += 8;
    }

    result
}

/// 恒星计算函数
///
/// # Arguments
/// * `t` - 儒略世纪TD
/// * `stars` - 恒星数据数组
/// * `q` - 章动周期天数(为0不限制)
/// * `lx` - 计算类型(0:视赤经赤纬, 1:站心坐标, 2:平赤经赤纬)
/// * `l` - 地理经度(弧度，东经为正)
/// * `fa` - 地理纬度(弧度，北纬为正)
///
/// # Returns
/// * 计算结果字符串
pub fn hx_calc(t: f64, stars: &[f64], q: i32, lx: usize, l: f64, fa: f64) -> String {
    let mut result = String::new();
    // 设置标题
    let s0 = match lx {
        0 => "视赤经 视赤纬".to_string(),
        1 => "站心坐标".to_string(),
        2 => "平赤经 平赤纬".to_string(),
        _ => "未知坐标".to_string(),
    };

    // 计算所需的天文参数
    let mut d = [0.0, 0.0];
    let mut e = 0.0;
    let mut v = [0.0, 0.0, 0.0];
    let mut p = [0.0, 0.0, 0.0];
    let mut a = [0.0, 0.0, 0.0];
    let gst_p;
    let mut gst = 0.0;

    if lx == 0 || lx == 1 {
        d = Nutation::nutation(t, q); // 章动
        e = Prece::hcjj(t); // 黄赤交角
        v = ev_ssb(t); // 地球SSB直角速度(光行差使用的)
        p = ep_ssb(t); // 地球SSB直角位置(视差修正使用的)
        a = sun_2000(t, 20); // 太阳J2000球面坐标(引力偏转用的)
        a = llr_conv(&a, 84381.406 / RAD); // 太阳赤道坐标
        gst_p = p_gst2(t * 36525.0); // 平恒星时
        gst = gst_p + d[0] * e.cos(); // 真恒星时
    }

    // 遍历所有恒星
    for i in (0..stars.len()).step_by(8) {
        // 获取恒星信息
        let star_info = format!(
            "{} {} {} ",
            stars[i + 6],
            stars[i + 7],
            // get_star_info(stars, i + 6),
            // get_star_info(stars, i + 7),
            stars[i + 5]
        );

        // 计算J2000赤经赤纬(含自行)
        let mut z = [
            stars[i] + stars[i + 2] * t * 100.0,     // J2000赤经含自行
            stars[i + 1] + stars[i + 3] * t * 100.0, // J2000赤纬含自行
            if stars[i + 4] != 0.0 {
                1.0 / stars[i + 4]
            } else {
                1e11
            }, // 距离
        ];

        // 标准化赤经
        z[0] = rad2mrad(z[0]);

        if lx == 0 || lx == 1 {
            // 引力偏转修正
            z = ylpz(&z, &a);

            // 周年视差修正
            z = sc_gxc(&z, &p, true);

            // 光行差修正
            z = sc_gxc(&z, &v, false);

            // 转到当日赤道(岁差修正)
            z = Prece::cdllr_j2d(t, &z, "P03");

            // 章动修正
            z = Nutation::cd_nutation(&z, e, d[0], d[1]);

            if lx == 1 {
                // 站心坐标
                // let sj = rad2rrad(gst + l - z[0]); // 得到此刻天体时角 注：没用到
                z[0] += PI / 2.0 - gst - l; // 转到相对于地平赤道分点的赤道坐标
                let zz = llr_conv(&z, PI / 2.0 - fa); // 恒星地平坐标
                z = [zz[0],zz[1],zz[2]];
                z[0] = rad2mrad(PI / 2.0 - z[0]); // 方位角,高度角

                // 大气折射修正
                if z[1] > 0.0 {
                    z[1] += mqc(z[1]);
                }
            }
        } else if lx == 2 {
            // 转到当日赤道(岁差修正)
            z = Prece::cdllr_j2d(t, &z, "P03");
        }

        // 格式化输出
        let pos_str = match lx {
            0 | 2 => format!(
                "{} {}\n",
                rad2str_e(z[0], 1, 3),
                rad2str_e(z[1], 0, 2)
            ),
            _ => format!(
                "{} {}\n",
                rad2str_e(z[0], 0, 2),
                rad2str_e(z[1], 0, 2)
            ),
        };

        result.push_str(&(star_info + &pos_str));
    }

    // 添加时间信息和标题
    let header = format!("{} TD {}\n",JD::new().jd2str(t * 36525.0 + J2000), s0);

    header + &result + "\n"
}
