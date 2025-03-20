// 声明 eph 模块
pub mod eph;

use crate::eph::{delta_t::dt_t, eph_base::pty_zty, xl::XL};

fn main() {
    println!("寿星天文历法 Rust 实现");

    // 示例:计算 2024 年春分点
    let jd = 2460034.5; // 2024年春分

    // 计算 ΔT
    let dt = dt_t(jd);
    println!("ΔT = {:.6} 日", dt);

    // 计算时差
    let t = (jd - 2451545.0) / 36525.0; // 转换为儒略世纪数
    let eq = pty_zty(t);
    println!("时差 = {:.6} 日", eq);

    // 计算太阳黄经
    let sun_lon = XL::xl0_calc(0, 0, t, 50);
    println!("太阳黄经 = {:.6}°", sun_lon.to_degrees());

    // 计算月亮位置
    let moon = XL::m_coord(t, 50, 50, 50);
    println!("月球位置:");
    println!("  黄经 = {:.6}°", moon[0].to_degrees());
    println!("  黄纬 = {:.6}°", moon[1].to_degrees());
    println!("  地月距离 = {:.6} AU", moon[2]);
}
