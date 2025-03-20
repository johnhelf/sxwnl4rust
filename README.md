# sxwnl4rust - 寿星天文历的Rust移植版

## 项目简介

sxwnl4rust 是在 Claude-3.7-sonnet 的辅助下，从许剑伟老师的寿星天文历（万年历）JavaScript 版本移植到 Rust 的项目。原项目是一款精准的跨度大的天文历法工具，可作为一般的实用日历工具，对史学工作者、历算工作者、天文爱好者均有较大的参考价值。

原项目地址：[https://github.com/sxwnl/sxwnl](https://github.com/sxwnl/sxwnl)

## 项目状态

⚠️ **注意事项**：
- 本项目目前仍在开发和调整中
- 不能保证所有计算结果的准确性
- 部分功能可能尚未完全实现或与原版有所差异

## 实现差异说明

由于原项目是用 JavaScript 编写的，在移植过程中存在一些实现上的差异：

1. 原 JS 版本中许多方法直接输出格式化的字符串，而在 Rust 版本中更倾向于返回结构化数据
2. 部分算法逻辑进行了调整以适应 Rust 的所有权和借用规则
3. 一些全局变量和函数被重构为结构体和方法
4. 数据结构和接口设计更符合 Rust 的惯用法

## 基本用法

### 农历日期转换

```rust

use chrono::{TimeZone, Utc};

fn main() {
    // 创建农历计算实例
    let lunar = sxwnl4rust::lunar::lunar::LunarCalendar::default();

    // 设置公历日期（年、月、日）
    let date = Utc.with_ymd_and_hms(2025, 5, 1, 12, 0, 0).single().unwrap();

    // 获取农历信息
    let lunar_day = lunar.get_day_info(&date, 116.4); // 使用北京经度

    // 输出农历日期信息
    println!("公历日期: 2025年5月1日");
    println!(
        "农历日期: {}年{}月{}",
        lunar_day.lunar_year, lunar_day.lunar_month_name, lunar_day.lunar_day_name
    );
    println!("干支纪年: {}", lunar_day.gan_zhi_year);
    println!("干支纪月: {}", lunar_day.gan_zhi_month);
    println!("干支纪日: {}", lunar_day.gan_zhi_day);
    println!("生肖: {}", lunar_day.sheng_xiao);
}
```

更多示例请参考 `examples` 目录下的示例代码。

## 目录结构

```
sxwnl4rust/
├── Cargo.toml              # 项目配置文件
├── README.md               # 项目说明文档
├── src/
│   ├── lib.rs              # 库入口文件
│   ├── lunar/              # 农历相关模块
│   │   ├── mod.rs          # 农历模块入口
│   │   ├── lunar.rs        # 农历核心算法
│   │   ├── solar_terms.rs  # 农历节气计算
│   │   └── niao_hao.rs     # 中国古代年号
│   ├── eph/                # 天文历相关模块
│   │   ├── mod.rs          # 天文历模块入口
│   │   ├── eph.rs          # 天文历核心算法
│   │   ├── msc.rs          # 日月计算
│   │   ├── eph_base.rs     # 基础天文算法
│   │   ├── jd.rs           # 儒略日转换
│   │   ├── delta_t.rs      # 黄道坐标转换
│   │   ├── eph_ssb.rs      # 地球的SSB速度计算
│   │   ├── nutation.rs     # 章动计算
│   │   ├── prece.rs        # 岁差表
│   │   ├── szj.rs          # 日月升降计算
│   │   └── xl.rs           # 星历计算
│   ├── funcs/              # 其它功能模块
│   │   ├── mod.rs          # 工具模块入口
│   └── └── festival.rs     # 节日计算模块
├── examples/               # 示例代码
│   ├── eph_calc.rs         # 天文历法计算示例
│   ├── lunar_date.rs       # 农历日期转换示例
│   └── fest_terms.rs       # 节日和节气信息示例

```

## 致谢

感谢许剑伟老师开发的原版寿星天文历，为中国传统历法的研究和应用做出了重要贡献。本项目旨在将这些宝贵的算法移植到 Rust 生态系统中，使更多的开发者能够在现代应用中使用这些历法算法。

## 版权说明

本项目是对原寿星天文历的移植，遵循原项目的开源精神。如在您的软件中使用了本程序的核心算法及数据，建议在您的软件中申明"算法来源于寿星天文历的Rust移植版"。