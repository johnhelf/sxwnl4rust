//! sxwnl4rust - 寿星天文历的Rust移植版
//!
//! 本库提供了中国传统历法和天文历法的计算功能，
//! 包括农历日期转换、节气计算、干支纪年等功能。

// 导出农历相关模块
// #[path] 指定大小写正确的路径（源码目录为 src/Lunar/），兼容 Linux 大小写敏感文件系统
#[path = "Lunar/mod.rs"]
pub mod lunar;

// 导出天文历相关模块
pub mod eph;

// 导出功能函数模块
pub mod funcs;