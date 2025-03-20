use chrono::{DateTime, NaiveDateTime, Utc}; //Datelike,Timelike,
// use sxwnl4rust::eph::jd::JD;
use sxwnl4rust::lunar::lunar::LunarCalendar;
// use sxwnl4rust::eph::eph_base::J2000;

fn main() {
    // 解析命令行参数
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("用法: {} <日期时间(YYYY-MM-DD HH:MM:SS)> <经度>", args[0]);
        std::process::exit(1);
    }

    // 解析日期时间
    let dt_str = &args[1];
    let dt = NaiveDateTime::parse_from_str(dt_str, "%Y-%m-%d %H:%M:%S")
        .expect("日期时间格式应为YYYY-MM-DD HH:MM:SS");
    // let dt_utc = DateTime::<Utc>::from_utc(dt, Utc);
    let dt_utc = DateTime::from_naive_utc_and_offset(dt, Utc);

    // 解析经度
    let longi = args[2].parse::<f64>().expect("经度应为浮点数");


    // let mut jd = JD::new();
    // jd.year = dt.year();
    // jd.month = dt.month() as i32;
    // jd.day = dt.day() as i32;
    // jd.hour = dt.hour() as f64;
    // jd.minute = dt.minute() as f64;
    // jd.second = dt.second() as f64;
    
    let calendar = LunarCalendar::default();

    let ba_zi_d = calendar.ming_li_ba_zi_detail(&dt_utc, longi,-8); //北京时间= -8

    // 格式化输出结果
    println!("输入时间：{}\n", dt);
    // println!("真太阳时：{:.2}\n", jd);
    println!("八字结果：");
    println!("年柱：{}", ba_zi_d.bazi.year);
    println!("月柱：{}", ba_zi_d.bazi.month);
    println!("日柱：{}", ba_zi_d.bazi.day);
    println!("时柱：{}", ba_zi_d.bazi.hour);
}