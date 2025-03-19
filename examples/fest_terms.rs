
use sxwnl4rust::lunar::lunar::{JIE_QI_NAMES, LunarCalendar};
use sxwnl4rust::eph::jd::JD;
use sxwnl4rust::lunar::solar_terms::SolarTerms;
// use sxwnl4rust::funcs::day::Day;
// use sxwnl4rust::eph::xl::XL;
// use std::f64::consts::PI;
use sxwnl4rust::eph::eph_base::J2000;

fn main() {
    println!("节日和节气信息示例");
    
    let mut lunar = LunarCalendar::default();
    
    // 修正JD初始化方式
    println!("\n2025年节气信息:");
    let year = 2025;
    let mut jd = JD::new();
    jd.year = year;
        jd.month = 1;
        jd.day = 1;
        jd.hour = 12.0;
        jd.minute = 0.0;
        jd.second = 0.1;
    
   
    //月份是0-11

    // if month == 12 {
    //     jd.year = year + 1;
    //     jd.month = 1;
    // } else {
        jd.month = jd.month + 1; 
    // }

     // 计算月首的儒略日
     let bd0 = jd.to_jd() as i32 - J2000 as i32; //公历月首,中午

    // 计算全年节气
    let mut ssq = SolarTerms::new();
   
    ssq.calc_year(bd0 as f64);
    
    for (i, name) in JIE_QI_NAMES.iter().enumerate() {
        let jq_time = ssq.solar_terms[i];
        let mut jd = JD::new();
        jd.set_from_jd(jq_time);  // 使用set_from_jd方法
        println!("  {}: {}年{}月{}日 {:02}:{:02}:{:02}", 
            name,
            jd.year,
            jd.month,
            jd.day,
            jd.hour as i32,
            jd.minute as i32,
            jd.second as i32
        );
    }
    
    // 计算2025年1月的月历（使用库方法）
    println!("\n2025年1月月历:");
    let month_calendar = lunar.yue_li_calc(year, 1);
    
    // 显示节日信息（通过FestivalInfo结构）
    println!("\n2025年1月节日和节气:");
    for day in &month_calendar.days_info {
        let festival = &day.festival;
        println!("  {}/{}/{} {}月{}", year, 1, day.day.d, day.lunar_month_name, day.lunar_day_name);
        
        if !festival.major.is_empty() {
            println!("    重要节日: {}", festival.major.join(", "));
        }
        if !festival.minor.is_empty() {
            println!("    一般节日: {}", festival.minor.join(", "));
        }
        if !festival.other.is_empty() {
            println!("    其他纪念日: {}", festival.other.join(", "));
        }
        if !day.jie_qi.is_empty() {
            println!("    节气: {}", day.jie_qi);
        }
    }
    
    // 计算特定节日的正确方式
    println!("\n特定日期信息:");
    
    // 1. 获取2025年春节（通过月历数据）
    let spring_festival = month_calendar.days_info.iter()
        .find(|d| d.lunar_month_name == "正" && d.lunar_day_name == "初一")
        .unwrap();
    println!("  2025年春节: {}/{}/{}", year,spring_festival.day.m+1, spring_festival.day.d);
    
    // 2. 获取清明节（通过节气数据）
    // 修正清明节日期获取方式
    let qingming = ssq.solar_terms[7];
    let mut jd_qm = JD::new();
    jd_qm.set_from_jd(qingming);
    println!("  2025年清明节: {}/{}", jd_qm.month, jd_qm.day);
    
    // 3. 通过月历获取元旦信息
    let new_year_day = &month_calendar.days_info[0];
    println!("  2025年元旦: {:?}", new_year_day.festival.major);
}