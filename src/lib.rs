//! sxwnl4rust - 寿星天文历的Rust移植版
//!
//! 本库提供了中国传统历法和天文历法的计算功能，
//! 包括农历日期转换、节气计算、干支纪年等功能。

// 导出农历相关模块
pub mod lunar;

// 导出天文历相关模块
pub mod eph;

// 导出功能函数模块
pub mod funcs;