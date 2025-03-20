//! 寿星天文历法 Rust 实现 - 天文计算模块
//!
//! 本模块包含各种天文计算所需的子模块和函数

// 声明并导出子模块
pub mod delta_t; // ΔT 计算模块
pub mod eph;
pub mod eph_base; // 基础天文计算模块
pub mod jd;
pub mod msc;
pub mod szj;
pub mod xl; // 行星和月亮位置计算模块 // 太阳系质心计算模块

// 可选：添加其他子模块
pub mod eph_ssb;
pub mod nutation; // 章动计算模块
pub mod prece; // 岁差计算模块

// 可选：重新导出常用函数，使调用更方便
// pub use self::delta_t::dt_t;
// pub use self::eph_base::eph_base;
// pub use self::xl::XL;

/// 模块版本信息
pub const VERSION: &str = "0.1.0";
