use chrono::{Datelike, Timelike}; //DateTime, TimeZone, Utc, 
use sxwnl4rust::lunar::lunar::LunarCalendar;

#[test]
fn test_from_lunar_date() {
    let lunar_calendar = LunarCalendar::default();
    
    // 测试用例1: 2023年正月初一 (非闰月) - 农历新年
    let result = lunar_calendar.from_lunar_date(2023, 1, false, 1, 12, 0, 0.0, 120.0, 8.0);
    //打印结果
    println!("{:?}", result);
    assert!(result.is_some());
    let date_time = result.unwrap();
    assert_eq!(date_time.year(), 2023);
    assert_eq!(date_time.month(), 1); // 2023年农历正月初一是2023年1月22日
    assert_eq!(date_time.day(), 22);
    
    // 测试用例2: 2023年闰二月初一 - 测试闰月
    let result = lunar_calendar.from_lunar_date(2023, 2, true, 1, 12, 0, 0.0, 120.0, 8.0);
    assert!(result.is_some());
    let date_time = result.unwrap();
    assert_eq!(date_time.year(), 2023);
    assert_eq!(date_time.month(), 3); // 2023年闰二月初一是2023年3月22日
    assert_eq!(date_time.day(), 22);
    
    // 测试用例3: 2022年腊月三十 (除夕)
    let result = lunar_calendar.from_lunar_date(2022, 12, false, 30, 12, 0, 0.0, 120.0, 8.0);
    assert!(result.is_some());
    let date_time = result.unwrap();
    assert_eq!(date_time.year(), 2023);
    assert_eq!(date_time.month(), 1);
    assert_eq!(date_time.day(), 21); // 2022年腊月三十是2023年1月21日
    
    // 测试用例4: 2020年闰四月十五 - 另一个闰月测试
    let result = lunar_calendar.from_lunar_date(2020, 4, true, 15, 12, 0, 0.0, 120.0, 8.0);
    assert!(result.is_some());
    let date_time = result.unwrap();
    assert_eq!(date_time.year(), 2020);
    assert_eq!(date_time.month(), 6);
    assert_eq!(date_time.day(), 6); // 2020年闰四月十五是2020年6月6日
    
    // 测试用例5: 测试时分秒和经度参数
    let result = lunar_calendar.from_lunar_date(2023, 5, false, 5, 18, 30, 45.0, 116.4, 8.0);
    assert!(result.is_some());
    let date_time = result.unwrap();
    assert_eq!(date_time.hour(), 18);
    assert_eq!(date_time.minute(), 45); //北京东经116.4度，比北京时间（120.0）晚14分24秒
    assert_eq!(date_time.second(), 9);
    
    // 测试用例6: 边界测试 - 无效年份
    let result = lunar_calendar.from_lunar_date(3001, 1, false, 1, 12, 0, 0.0, 120.0, 8.0);
    assert!(result.is_none());
    
    // 测试用例7: 边界测试 - 无效月份
    let result = lunar_calendar.from_lunar_date(2023, 13, false, 1, 12, 0, 0.0, 120.0, 8.0);
    assert!(result.is_none());
    
    // 测试用例8: 边界测试 - 无效日期
    let result = lunar_calendar.from_lunar_date(2023, 1, false, 31, 12, 0, 0.0, 120.0, 8.0);
    assert!(result.is_none());
    
    // 测试用例9: 边界测试 - 无效小时
    let result = lunar_calendar.from_lunar_date(2023, 1, false, 1, 24, 0, 0.0, 120.0, 8.0);
    assert!(result.is_none());
    
    // 测试用例10: 边界测试 - 无效分钟
    let result = lunar_calendar.from_lunar_date(2023, 1, false, 1, 12, 60, 0.0, 120.0, 8.0);
    assert!(result.is_none());
    
    // 测试用例11: 边界测试 - 无效秒数
    let result = lunar_calendar.from_lunar_date(2023, 1, false, 1, 12, 0, 61.0, 120.0, 8.0);
    assert!(result.is_none());
    
    // 测试用例12: 测试不存在的闰月
    let result = lunar_calendar.from_lunar_date(2023, 3, true, 1, 12, 0, 0.0, 120.0, 8.0); // 2023年没有闰三月
    assert!(result.is_none());
    
    // 测试用例13: 测试不同时区
    let result = lunar_calendar.from_lunar_date(2023, 1, false, 1, 12, 0, 0.0, 120.0, 0.0); // UTC时区
    assert!(result.is_some());
    let date_time = result.unwrap();
    // 验证时区转换是否正确
    assert_eq!(date_time.hour(), 4); // 北京时间12点对应UTC时间4点
}