//! 农历日期转换示例
//! 
//! 本示例展示如何使用sxwnl4rust库进行公历日期到农历日期的转换，
//! 以及获取节气、干支纪年等信息。

use chrono::{TimeZone, Utc};

fn main() {
    println!("农历日期转换示例");
    
    // 创建农历计算实例
    let lunar = sxwnl4rust::lunar::lunar::LunarCalendar::default();
    
    // 设置公历日期（年、月、日）
    let date = Utc.with_ymd_and_hms(2025, 5, 1, 12, 0, 0).single().unwrap();
    
    // 获取农历信息
    let lunar_day = lunar.get_day_info(&date, 116.4); // 使用北京经度
    
    // 输出农历日期信息
    println!("公历日期: 2023年5月1日");
    println!("农历日期: {}年{}月{}", 
             lunar_day.lunar_year, 
             lunar_day.lunar_month_name, 
             lunar_day.lunar_day_name);
    println!("干支纪年: {}", lunar_day.gan_zhi_year);
    println!("干支纪月: {}", lunar_day.gan_zhi_month);
    println!("干支纪日: {}", lunar_day.gan_zhi_day);
    println!("生肖: {}", lunar_day.sheng_xiao);
    
    // 获取节气信息
    if !lunar_day.jie_qi.is_empty() {
        println!("今天是节气: {}", lunar_day.jie_qi);
        println!("节气时刻: {}", lunar_day.jie_qi_str);
    }
    
    // 获取月相信息
    println!("月相: {}", lunar_day.moon_phase);
    println!("月相时刻: {}", lunar_day.moon_phase_str);
    
    // 计算不同年份的农历信息
    println!("\n其他年份的农历信息:");
    for year in [2000, 2008, 2020, 2024] {
        let date = Utc.with_ymd_and_hms(year, 1, 1, 12, 0, 0).single().unwrap();
        let day = lunar.get_day_info(&date, 116.4);
        println!("{0}年1月1日 - 农历: {1}年{2}月{3} {4}年 生肖: {5}", 
                 year, 
                 day.lunar_year, 
                 day.lunar_month_name, 
                 day.lunar_day_name, 
                 day.gan_zhi_year, 
                 day.sheng_xiao);
    }
}