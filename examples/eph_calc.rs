//! 天文历法计算示例
//!
//! 本示例展示如何使用sxwnl4rust库进行天文历法计算，
//! 包括太阳位置、月亮位置、行星位置等计算。

use std::f64::consts::PI;
use sxwnl4rust::eph::{
    delta_t::dt_t,
    eph::{xing_jj, xing_liu0},
    eph_base::{J2000, llr_conv, rad2str_e},
    jd::JD,
    nutation::Nutation,
    prece::Prece,
    xl::XL,
};

fn main() {
    println!("天文历法计算示例");

    // 创建日期计算实例
    let jd_obj = JD::new();

    // 计算儒略日（2023年5月1日12时）
    let d = 1.0 + 12.0 / 24.0; // 1日12时
    let jd = jd_obj.jd(2025, 5, d);
    println!("儒略日: {:.5}", jd);

    // 计算力学时
    let dt = dt_t(jd);
    let t = (jd + dt - J2000) / 36525.0; // 儒略世纪数
    println!("力学时差: {:.5}秒", dt * 86400.0);

    // 计算太阳位置
    let sun = XL::p_coord(0, t, -1, -1, -1);
    // 转换太阳坐标（地心黄道坐标）
    let sun_lon = sun[0] + PI; // 太阳黄经 = 地球黄经 + 180°
    let sun_lat = -sun[1]; // 太阳黄纬 = -地球黄纬
    let sun_r = sun[2]; // 日地距离

    // 计算黄赤交角
    let e = Prece::hcjj(t);
    // 计算章动
    let nut = Nutation::nutation2(t);
    let lon_nut = nut[0]; // 黄经章动
    let obl_nut = nut[1]; // 交角章动

    // 转换到赤道坐标
    let sun_eq = llr_conv(&[sun_lon + lon_nut, sun_lat, sun_r], e + obl_nut);

    println!("\n太阳位置:");
    println!("  黄经: {:.4}°", sun_lon * 180.0 / PI);
    println!("  黄纬: {:.4}°", sun_lat * 180.0 / PI);
    println!("  地心距离: {:.8} AU", sun_r);
    println!("  视赤经: {}", rad2str_e(sun_eq[0], 1, 2));
    println!("  视赤纬: {}", rad2str_e(sun_eq[1], 0, 2));

    // 计算月亮位置
    let moon = XL::m_coord(t, -1, -1, -1);

    // 转换到赤道坐标
    let moon_eq = llr_conv(&[moon[0] + lon_nut, moon[1], moon[2]], e + obl_nut);

    // 计算月相（月日黄经差）
    let phase_angle = (moon[0] - sun_lon).abs();
    let phase = (1.0 - phase_angle.cos()) / 2.0;

    println!("\n月亮位置:");
    println!("  黄经: {:.4}°", moon[0] * 180.0 / PI);
    println!("  黄纬: {:.4}°", moon[1] * 180.0 / PI);
    println!("  地心距离: {:.8} 千米", moon[2]);
    println!("  视赤经: {}", rad2str_e(moon_eq[0], 1, 2));
    println!("  视赤纬: {}", rad2str_e(moon_eq[1], 0, 2));
    println!("  月相: {:.1}%", phase * 100.0);

    // 计算行星位置
    println!("\n行星位置:");
    let planets = [
        "水星",
        "金星",
        "火星",
        "木星",
        "土星",
        "天王星",
        "海王星",
        "冥王星",
    ];

    for (i, name) in planets.iter().enumerate() {
        let planet_idx = i + 1;
        if planet_idx <= 8 {
            // 计算行星视坐标
            let planet = xing_liu0(planet_idx, t, -1, 0.0);

            println!("{}:", name);
            println!("  视赤经: {}", rad2str_e(planet[0], 1, 2));
            println!("  视赤纬: {}", rad2str_e(planet[1], 0, 2));
            println!("  地心距离: {:.8} AU", planet[2]);

            // 计算行星与太阳的角距离
            let angle = xing_jj(planet_idx, t, 2);
            println!("  与太阳角距: {:.2}°", angle * 180.0 / PI);
        }
    }
}
