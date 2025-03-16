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
use sxwnl4rust::lunar::Lunar;

fn main() {
    // 创建农历计算实例
    let mut lunar = Lunar::new();
    
    // 设置公历日期（年、月、日）
    lunar.set_day(2023, 5, 1);
    
    // 获取农历信息
    let lunar_info = lunar.get_lunar_info();
    println!("农历日期: {}年{}月{}", lunar_info.lunar_year_name, lunar_info.lunar_month_name, lunar_info.lunar_day_name);
    println!("干支纪年: {}", lunar_info.gz_year);
    println!("干支纪月: {}", lunar_info.gz_month);
    println!("干支纪日: {}", lunar_info.gz_day);
    
    // 获取节气信息
    let jq_info = lunar.get_jieqi_info();
    if jq_info.is_jieqi {
        println!("今天是节气: {}", jq_info.name);
    }
    
    // 获取下一个节气
    let next_jq = lunar.get_next_jieqi();
    println!("下一个节气: {}, 时间: {}", next_jq.name, next_jq.time);
}
```

### 天文历计算

```rust
use sxwnl4rust::eph::Eph;
use sxwnl4rust::util::jd::{ date_to_jd, jd_to_date };

fn main() {
    // 创建天文历实例
    let mut eph = Eph::new();
    
    // 计算指定日期的儒略日
    let jd = date_to_jd(2023, 5, 1, 12.0);
    
    // 计算太阳位置
    let sun_pos = eph.get_sun_position(jd);
    println!("太阳位置: 赤经 {:.4}°, 赤纬 {:.4}°", sun_pos.ra, sun_pos.dec);
    
    // 计算月亮位置
    let moon_pos = eph.get_moon_position(jd);
    println!("月亮位置: 赤经 {:.4}°, 赤纬 {:.4}°", moon_pos.ra, moon_pos.dec);
    println!("月相: {:.2}%", moon_pos.phase * 100.0);
    
    // 计算行星位置
    let planets = eph.get_planet_positions(jd);
    for (name, pos) in planets {
        println!("{}: 赤经 {:.4}°, 赤纬 {:.4}°", name, pos.ra, pos.dec);
    }
}
```

### 日食和月食计算

```rust
use sxwnl4rust::eph::Eph;
use chrono::{NaiveDate, Datelike};

fn main() {
    let mut eph = Eph::new();
    
    // 计算指定年份的日食
    let year = 2023;
    let solar_eclipses = eph.get_solar_eclipses(year);
    
    println!("{}年的日食:", year);
    for eclipse in solar_eclipses {
        let date = NaiveDate::from_ymd_opt(
            eclipse.date.year, 
            eclipse.date.month as u32, 
            eclipse.date.day as u32
        ).unwrap();
        
        println!("日期: {}, 类型: {}, 最大食分: {:.2}", 
                 date.format("%Y-%m-%d"), 
                 eclipse.eclipse_type, 
                 eclipse.magnitude);
    }
    
    // 计算指定年份的月食
    let lunar_eclipses = eph.get_lunar_eclipses(year);
    
    println!("\n{}年的月食:", year);
    for eclipse in lunar_eclipses {
        let date = NaiveDate::from_ymd_opt(
            eclipse.date.year, 
            eclipse.date.month as u32, 
            eclipse.date.day as u32
        ).unwrap();
        
        println!("日期: {}, 类型: {}, 食分: {:.2}", 
                 date.format("%Y-%m-%d"), 
                 eclipse.eclipse_type, 
                 eclipse.magnitude);
    }
}
```

### 二十四节气计算

```rust
use sxwnl4rust::lunar::Lunar;
use chrono::NaiveDateTime;

fn main() {
    let mut lunar = Lunar::new();
    
    // 计算2023年的二十四节气
    let year = 2023;
    let jieqi_list = lunar.get_year_jieqi(year);
    
    println!("{}年二十四节气:", year);
    for jq in jieqi_list {
        // 将儒略日转换为日期时间
        let dt = NaiveDateTime::from_timestamp_opt(
            (jq.time - 2440587.5) * 86400.0 as i64, 
            0
        ).unwrap();
        
        println!("{}: {}", jq.name, dt.format("%Y-%m-%d %H:%M:%S"));
    }
}
```

更多示例请参考 `examples` 目录下的示例代码（还没上传，后续添加）。

## 目录结构

```
sxwnl4rust/
├── Cargo.toml           # 项目配置文件
├── README.md            # 项目说明文档
├── src/
│   ├── lib.rs           # 库入口文件
│   ├── lunar/           # 农历相关模块
│   │   ├── mod.rs       # 农历模块入口
│   │   ├── lunar.rs     # 农历核心算法
│   │   ├── jd.rs        # 儒略日转换
│   │   └── ...
│   ├── eph/             # 天文历相关模块
│   │   ├── mod.rs       # 天文历模块入口
│   │   ├── eph.rs       # 天文历核心算法
│   │   ├── msc.rs       # 日月食计算
│   │   ├── eph0.rs      # 基础天文算法
│   │   └── ...
│   ├── util/            # 工具函数
│   │   ├── mod.rs       # 工具模块入口
│   │   ├── constants.rs # 常量定义
│   │   └── ...
│   └── ...
├── examples/            # 示例代码
│   ├── basic_usage.rs   # 基本用法示例
│   └── ...
└── tests/               # 测试代码
    ├── lunar_tests.rs   # 农历测试
    ├── eph_tests.rs     # 天文历测试
    └── ...
```

## 致谢

感谢许剑伟老师开发的原版寿星天文历，为中国传统历法的研究和应用做出了重要贡献。本项目旨在将这些宝贵的算法移植到 Rust 生态系统中，使更多的开发者能够在现代应用中使用这些历法算法。

## 版权说明

本项目是对原寿星天文历的移植，遵循原项目的开源精神。如在您的软件中使用了本程序的核心算法及数据，建议在您的软件中申明"算法来源于寿星天文历的Rust移植版"。