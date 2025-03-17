use crate::eph::eph_base::{J2000, PI,PI2,pty_zty2};
use chrono::{DateTime, Datelike, TimeZone, Utc,Timelike};
// use lazy_static::lazy_static;
use crate::eph::delta_t::dt_t;
use crate::eph::jd::JD;
use crate::lunar::nian_hao::{get_nian_hao};
use crate::eph::xl::XL;
use crate::funcs::festival::FestivalInfo;
use crate::lunar::solar_terms::SolarTerms;
use std::array;

/// 农历年份信息
#[derive(Default)]
pub struct LunarYear {
    /// 闰月位置
    pub leap: i32,
    /// 各月名称
    pub month_names: Vec<String>,
    /// 中气表,其中立秋的儒略日用于计算三伏
    pub solar_terms: Vec<f64>,
    /// 合朔表
    pub new_moons: Vec<f64>,
    /// 各月大小
    pub month_days: Vec<i32>,
    /// 年计数
    pub year_count: Vec<i32>,
}

#[derive(Debug, Default, Clone)]
pub struct Day {
    pub d0: i32,    //儒略日,北京时12:00
    pub di: i32,    //公历月内日序数
    pub y: i32,     //公历年
    pub m: i32,     //公历月
    pub dn: i32,    //公历月天数
    pub week0: i32, //月首的星期
    pub week: i32,  //当前日的星期
    pub weeki: i32, //本日所在的周序号
    pub weekn: i32, //本月的总周数
    pub d: i32,     //公历日
}

/// 农历日期信息
#[derive(Debug, Default, Clone)]
pub struct LunarDay {
    pub day: Day,

    pub ldi: i32, //距农历月首的编移量,0对应初一

    /// 距冬至天数
    pub cur_dz: i32,
    /// 距夏至天数  
    pub cur_xz: i32,
    /// 距立秋天数
    pub cur_lq: i32,
    /// 距芒种天数
    pub cur_mz: i32,
    /// 距小暑天数
    pub cur_xs: i32,

    /// 农历日名称(如:初一)
    pub lunar_day_name: String,
    /// 农历月名称
    pub lunar_month_name: String,
    /// 是否闰月
    pub is_leap_month: bool,
    /// 月大小(30/29天)
    pub lunar_month_days: i32,
    /// 下个月名称
    pub next_month_name: String,

    /// 节气名称
    pub jie_qi: String,

    pub jie_qi_time: f64,

    pub jie_qi_str: String,

    /// 月相名称
    pub moon_phase: String,

    pub moon_phase_time: f64,

    pub moon_phase_str: String,

    /// 年号纪年
    pub nian_hao: String,

    pub lunar_year: i32, //农历纪年(10进制,1984年起算)

    pub lunar_year0: i32, //农历纪年(10进制,1984年起算)

    //黄帝纪年
    pub lunar_year1: String,

    /// 干支纪年(立春)
    pub gan_zhi_year: String,
    //干支纪年(正月)
    pub gan_zhi_year1: String,
    //农历纪月
    pub lunar_month: i32,
    //干支纪月
    pub gan_zhi_month: String,

    pub gan_zhi_day: String,

    /// 生肖
    pub sheng_xiao: String,

    /// 星座
    pub constellation: String,

    /// 添加回历字段
    pub hui_li: Option<HuiLiDay>,

    /// 节日
    pub festival: FestivalInfo,
}

#[derive(Debug, Default)]
pub struct LunarSimpleDay {
    pub day: Day,

    pub ldi: i32, //距农历月首的编移量,0对应初一

    /// 农历日名称(如:初一)
    pub lunar_day_name: String,
    /// 农历月名称
    pub lunar_month_name: String,
    /// 是否闰月
    pub is_leap_month: bool,
    /// 月大小(30/29天)
    pub lunar_month_days: i32,
    /// 下个月名称
    pub next_month_name: String,

    /// 干支纪年(立春)
    pub gan_zhi_year: String,
    //干支纪年(正月)
    pub gan_zhi_year1: String,
    //农历纪月
    pub lunar_month: i32,
    //干支纪月
    pub gan_zhi_month: String,

    pub gan_zhi_day: String,

    /// 生肖
    pub sheng_xiao: String,
}

// /// 年号信息
// #[derive(Debug)]
// pub struct YearInfo {
//     /// 起始年(公元)
//     start_year: i32,
//     /// 使用年数
//     years: i32,
//     /// 已用年数
//     used_years: i32,
//     /// 朝代
//     dynasty: String,
//     /// 朝号
//     reign_title: String,
//     /// 皇帝
//     emperor: String,
//     /// 年号
//     year_name: String,
// }

/// 回历日期信息
#[derive(Debug, Default, Clone)]
pub struct HuiLiDay {
    /// 回历年
    pub year: i32,
    /// 回历月
    pub month: i32,
    /// 回历日
    pub day: i32,
}

/// 农历月历
#[derive(Debug)]
pub struct LunarMonth {
    /// 本月第一天的星期
    pub first_day_week: i32,
    /// 本月天数
    pub days: i32,
    /// 本月各天信息
    pub days_info: Vec<LunarDay>,
    /// 闰月位置(如果有)
    pub leap_month: Option<i32>,
}

pub struct JulianDate {
    pub year: i32,
    pub month: i32,
    pub day: i32,
    pub hour: i32,
    pub minute: i32,
    pub second: f64,
}

/// 农历计算器
#[derive(Default,Debug)]
pub struct LunarCalendar {
    pub gan: String,
    pub zhi: String,
    pub sheng_xiao: String,
    pub nian_hao: String,
}

pub const GAN: [&'static str; 10] = ["甲", "乙", "丙", "丁", "戊", "己", "庚", "辛", "壬", "癸"];
pub const ZHI: [&'static str; 12] = [
    "子", "丑", "寅", "卯", "辰", "巳", "午", "未", "申", "酉", "戌", "亥",
];
pub const SHENG_XIAO: [&'static str; 12] = [
    "鼠", "牛", "虎", "兔", "龙", "蛇", "马", "羊", "猴", "鸡", "狗", "猪",
];
pub const CONSTELLATION: [&'static str; 12] = [
    "摩羯", "水瓶", "双鱼", "白羊", "金牛", "双子", "巨蟹", "狮子", "处女", "天秤", "天蝎", "射手",
];
pub const MOON_PHASE_NAMES: [&'static str; 4] = ["朔", "上弦", "望", "下弦"];
pub const JIE_QI_NAMES: [&'static str; 24] = [
    "冬至", "小寒", "大寒", "立春", "雨水", "惊蛰", "春分", "清明", "谷雨", "立夏", "小满", "芒种",
    "夏至", "小暑", "大暑", "立秋", "处暑", "白露", "秋分", "寒露", "霜降", "立冬", "小雪", "大雪",
];
pub const YUE_MING: [&str; 12] = [
    "十一", "十二", "正", "二", "三", "四", "五", "六", "七", "八", "九", "十",
];

// 添加农历月名
pub const YUE_MC: [&'static str; 12] = [
    "正", "二", "三", "四", "五", "六", "七", "八", "九", "十", "冬", "腊",
];

// 添加农历日名
pub const RI_MC: [&'static str; 31] = [
    "初一", "初二", "初三", "初四", "初五", "初六", "初七", "初八", "初九", "初十", "十一", "十二",
    "十三", "十四", "十五", "十六", "十七", "十八", "十九", "二十", "廿一", "廿二", "廿三", "廿四",
    "廿五", "廿六", "廿七", "廿八", "廿九", "三十", "卅一",
];

/// 八字信息
#[derive(Default, Debug)]
pub struct BaZiInfo {
    /// 年干支
    pub year: String,
    /// 月干支  
    pub month: String,
    /// 日干支
    pub day: String,
    /// 时辰干支
    pub hour: String,
    /// 当前时辰
    pub cur_hour: String,
    /// 全天时辰表
    pub time_table: [String; 13],
}

/// 八字信息
#[derive(Default, Debug)]
pub struct BaZiDetail {

    pub date : DateTime<Utc>,

    pub year: i32,
    pub month: i32,
    pub day: i32,
    pub hour: i32,
    pub minute: i32,
    pub second: f64,

    pub day_info : LunarSimpleDay,
    pub bazi : BaZiInfo,
}

impl LunarCalendar {
    /// 创建农历计算器实例
    // pub fn new() -> Self {
    //     Self {
    //     //     gan: ["甲","乙","丙","丁","戊","己","庚","辛","壬","癸"],
    //     //     zhi: ["子","丑","寅","卯","辰","巳","午","未","申","酉","戌","亥"],
    //     //     sheng_xiao: ["鼠","牛","虎","兔","龙","蛇","马","羊","猴","鸡","狗","猪"],
    //     //     constellation: ["摩羯","水瓶","双鱼","白羊","金牛","双子",
    //     //                   "巨蟹","狮子","处女","天秤","天蝎","射手"],
    //     //     moon_phase_names: ["朔","上弦","望","下弦"],
    //     //     jie_qi_names: ["冬至","小寒","大寒","立春","雨水","惊蛰",
    //     //                   "春分","清明","谷雨","立夏","小满","芒种",
    //     //                   "夏至","小暑","大暑","立秋","处暑","白露",
    //     //                   "秋分","寒露","霜降","立冬","小雪","大雪"],
    //     }
    // }

    // 精确节气计算
    pub fn qi_accurate(&self, w: f64) -> f64 { //w 参数是节气的黄经值（太阳黄经
        // 调用天文算法计算节气时刻
        let t = XL::s_alon_t(w) * 36525.0;
        // 转换为UT1时间
        t - dt_t(t) + 8.0 / 24.0 // 转为北京时间(+8h)
    }

    // 精确朔望计算
    pub fn so_accurate(&self, w: f64) -> f64 {
        // 调用天文算法计算月相时刻
        let t = XL::m_lon_t(w) * 36525.0;
        // 转换为UT1时间
        t - dt_t(t) + 8.0 / 24.0 // 转为北京时间(+8h)
    }

    // 寻找最近的节气
    pub fn qi_accurate2(&self, jd: f64) -> f64 {
        let d = PI / 12.0; // 每个节气间隔15度

        // 估算节气时刻
        let w = ((jd + 293.0) / 365.2422 * 24.0).floor() * d;

        // 精确计算
        let a = self.qi_accurate(w);

        // 寻找最近的节气
        if (a - jd) > 5.0 {
            self.qi_accurate(w - d)
        } else if (a - jd) < -5.0 {
            self.qi_accurate(w + d)
        } else {
            a
        }
    }

    // 寻找最近的朔望
    pub fn so_accurate2(&self, jd: f64) -> f64 {
        // 估算朔望时刻
        self.so_accurate(((jd + 8.0) / 29.5306).floor() * 2.0 * PI)
    }

    /// 获取年号
    pub fn get_nian_hao(&self, year: i32) -> String {
        get_nian_hao(year)
    }

    /// 计算回历日期
    pub fn calc_hui_li(&self, jd: f64) -> HuiLiDay {
        // 以下算法使用Excel测试得到,主要关心年临界与月临界
        let d = jd + 503105.0; // 调整至回历纪元

        // 周期计算
        let z = (d / 10631.0).floor(); // 10631为一周期(30年)
        let mut d = d - z * 10631.0;

        // 年份计算,加0.5作用是保证闰年正确(一周中的闰年是第2,5,7,10,13,16,18,21,24,26,29年)
        let y = ((d + 0.5) / 354.366).floor();
        d -= (y * 354.366 + 0.5).floor();

        // 月份计算,分子加0.11,分母加0.01的作用是第354或355天的月份保持为12月
        let m = ((d + 0.11) / 29.51).floor();
        d -= (m * 29.5 + 0.5).floor();

        HuiLiDay {
            year: (z * 30.0 + y + 1.0) as i32,
            month: (m + 1.0) as i32,
            day: (d + 1.0) as i32,
        }
    }

    /// 计算指定公历年月的三合历
    /// # Arguments
    /// * `year` - 公历年份
    /// * `month` - 公历月份
    /// # Returns
    /// * `LunarMonth` - 月历信息
    pub fn yue_li_calc(&mut self, year: i32, month: i32) -> LunarMonth {
        // 基本参数计算
        let mut jd = JD::new();
        // jd.jd(year, month, 1, 12, 0, 0.1);
        jd.year = year;
        jd.month = month;
        jd.day = 1;
        jd.hour = 12.0;
        jd.minute = 0.0;
        jd.second = 0.1;
        // jd.h = 12;
        // jd.m = 0;
        // jd.s = 0.1;
        // jd.y = year;
        // jd.m = month;
        // jd.d = 1;

        // 计算月首的儒略日
        let bd0 = jd.to_jd() as i32 - J2000 as i32; //公历月首,中午

        // 计算本月天数
        if month == 12 {
            jd.year = year + 1;
            jd.month = 1;
        } else {
            jd.month = month + 1;
        }
        let bdn = jd.to_jd() as i32 - J2000 as i32 - bd0; //本月天数(公历)

        // 计算月首星期
        let w0 = (bd0 + J2000 as i32 + 1 + 7000000) % 7;

        // 生成LunarMonth结构
        let mut lunar_month = LunarMonth {
            first_day_week: w0,
            days: bdn,
            days_info: Vec::with_capacity(bdn as usize),
            leap_month: None,
        };

        // 计算该年的农历信息
        let cycle_year = year - 1984 + 12000;
        self.gan = GAN[(cycle_year % 10) as usize].to_string();
        self.zhi = ZHI[(cycle_year % 12) as usize].to_string();
        // let gan_zhi = format!("{}{}", self.gan, self.zhi);
        self.sheng_xiao = SHENG_XIAO[(cycle_year % 12) as usize].to_string();
        self.nian_hao = get_nian_hao(year);

        // 初始化农历历法计算器
        let mut ssq = SolarTerms::new();

        // 计算每日信息
        for i in 0..bdn {
            let mut day = LunarDay::default();

            // 1. 公历信息
            day.day.d0 = bd0 + i; //儒略日,北京时12:00
            day.day.di = i; // 月内序号
            day.day.y = year;
            day.day.m = month;
            day.day.dn = bdn;
            day.day.week0 = w0;
            day.day.week = (w0 + i) % 7;
            day.day.weeki = (w0 + i) / 7;
            day.day.weekn = (w0 + bdn - 1) / 7 + 1;
            jd.set_from_jd(day.day.d0 as f64 + J2000);
            day.day.d = jd.day;

            // 2. 农历信息
            // 如果超出计算农历范围则重新计算
            if ssq.solar_terms.is_empty()
                || (day.day.d0 as f64) < ssq.solar_terms[0]
                || (day.day.d0 as f64) >= ssq.solar_terms[24]
            {
                ssq.calc_year(day.day.d0 as f64);
            }

            // 找农历月
            let mut mk = ((day.day.d0 as f64 - ssq.new_moons[0]) / 30.0) as usize;
            if mk < 13 && ssq.new_moons[mk + 1] <= day.day.d0  as f64 {
                mk += 1;
            }

            // 农历月内日期
            day.ldi = (day.day.d0 as f64 - ssq.new_moons[mk]) as i32;
            day.lunar_day_name = RI_MC[day.ldi as usize].to_string();

            // 距节气天数
            day.cur_dz = ((day.day.d0 as f64) - ssq.solar_terms[0]) as i32; // 距冬至天数
            day.cur_xz = ((day.day.d0 as f64) - ssq.solar_terms[12]) as i32; // 距夏至天数
            day.cur_lq = ((day.day.d0 as f64) - ssq.solar_terms[15]) as i32; // 距立秋天数  
            day.cur_mz = ((day.day.d0 as f64) - ssq.solar_terms[11]) as i32; // 距芒种天数
            day.cur_xs = ((day.day.d0 as f64) - ssq.solar_terms[13]) as i32; // 距小暑天数

            // 月的信息(月名称、大小及闰状况)
            if day.day.d0 == ssq.new_moons[mk] as i32 || i == 0 {
                day.lunar_month_name = ssq.month_names[mk].clone();
                day.lunar_month_days = ssq.month_days[mk];
                day.is_leap_month = ssq.leap == mk as i32;
                day.next_month_name = if mk < 13 {
                    ssq.month_names[mk + 1].clone()
                } else {
                    "未知".to_string()
                };
            } else {
                let prev = &lunar_month.days_info[(i - 1) as usize];
                day.lunar_month_name = prev.lunar_month_name.clone();
                day.lunar_month_days = prev.lunar_month_days;
                day.is_leap_month = prev.is_leap_month;
                day.next_month_name = prev.next_month_name.clone();
            }

            // 计算节气
            let mut qk = ((day.day.d0 as f64 - ssq.solar_terms[0] - 7.0) / 15.2184) as i32; // 注意:改为i32
            if qk < 23 && day.day.d0 as f64 >= ssq.solar_terms[qk as usize + 1] {
                qk += 1; // 增加这一行,与JS版本保持一致
            }

            // 判断是否恰好是节气日
            if (day.day.d0 as f64 - ssq.solar_terms[qk as usize]).abs() < 0.5 {
                // 允许0.5天误差
                day.jie_qi = JIE_QI_NAMES[qk as usize].to_string();
            } else {
                day.jie_qi = String::new(); // 不是节气日则清空
            }

            // 干支纪年处理
            // 以立春为界定年首
            let tmp = if (day.day.d0 as f64) < ssq.solar_terms[3] {
                -365.0
            } else {
                0.0
            };
            let mut d = ssq.solar_terms[3] + tmp + 365.25 * 16.0 - 35.0;

            day.lunar_year = (d / 365.2422 + 0.5).floor() as i32; // 农历纪年(1984年起算)

            // 以正月初一定年首
            d = ssq.new_moons[2]; // 一般第3个月为春节
            // 找春节
            for j in 0..14 {
                if ssq.month_names[j] != "正" || (ssq.leap == j as i32 && j > 0) {
                    continue;
                }
                d = ssq.new_moons[j];
                if (day.day.d0 as f64) < d {
                    d -= 365.0;
                    break;
                }
            }

            // 计算该年春节与1984年平均春节(立春附近)相差天数
            d = d + 5810.0;
            day.lunar_year0 = (d / 365.2422 + 0.5).floor() as i32;

            // 干支纪年
            let d1 = day.lunar_year + 12000;
            day.gan_zhi_year = format!("{}{}", GAN[d1 as usize % 10], ZHI[d1 as usize % 12]);

            let d2 = day.lunar_year0 + 12000;
            day.gan_zhi_year1 = format!("{}{}", GAN[d2 as usize % 10], ZHI[d2 as usize % 12]);

            // 黄帝纪年
            day.lunar_year1 = (day.lunar_year0 + 1984 + 2698).to_string();

            // 纪月处理
            // 1998年12月7(大雪)开始连续进行节气计数,0为甲子
            let mut mk = ((day.day.d0 as f64 - ssq.solar_terms[0]) / 30.43685) as usize;
            if mk < 12 && day.day.d0 as f64 >= ssq.solar_terms[2 * mk + 1] {
                mk += 1;
            }

            // 相对于1998年12月7(大雪)的月数
            d = mk as f64 + ((ssq.solar_terms[12] + 390.0) / 365.2422).floor() * 12.0 + 900000.0;

            day.lunar_month = (d as usize % 12) as i32;
            day.gan_zhi_month = format!("{}{}", GAN[d as usize % 10], ZHI[d as usize % 12]);

            // 纪日,2000年1月7日起算
            d = day.day.d0 as f64 - 6.0 + 9000000.0;
            day.gan_zhi_day = format!("{}{}", GAN[d as usize % 10], ZHI[d as usize % 12]);

            // 星座
            mk = ((day.day.d0 as f64 - ssq.solar_terms[0] - 15.0) / 30.43685) as usize;
            if mk < 11 && day.day.d0 as f64 >= ssq.solar_terms[2 * mk + 2] {
                mk += 1;
            }
            day.constellation = format!("{}座", CONSTELLATION[(mk + 12) % 12]);

            // 回历
            day.hui_li = Some(self.calc_hui_li(day.day.d0 as f64));

            let mut festival_info = FestivalInfo::default();
            // 节日
            day.festival = festival_info.calc_festivals(&day);

            lunar_month.days_info.push(day);
        } //for

        //以下是月相与节气的处理

        let jd2 = bd0 as f64 + dt_t(bd0 as f64) - 8.0 / 24.0; // 转换为力学时
        // 月相查找
        let mut w = XL::ms_alon(jd2 / 36525.0, 10, 3);
        w = ((w - 0.78) / PI * 2.0).floor() * PI / 2.0;

        loop {
            // 计算月相时刻
            let d = self.so_accurate(w);
            let d_int = (d + 0.5).floor() as i32;

            // 计算月相序号(朔望凸亏)
            let xn = ((w / PI2 * 4.0 + 4000000.01).floor() as i32 % 4) as usize;

            w += PI2 / 4.0;

            // 超出本月范围则结束
            if d_int >= bd0 + bdn {
                break;
            }

            // 未到本月则继续
            if d_int < bd0 {
                continue;
            }

            // 更新日对象的月相信息
            let mut day = lunar_month.days_info[(d_int - bd0) as usize].clone() as LunarDay;
            day.moon_phase = MOON_PHASE_NAMES[xn].to_string();
            day.moon_phase_time = d;
            day.moon_phase_str = jd.time_str(d);

            // 相邻日期超出本月则结束
            if d_int + 5 >= bd0 + bdn {
                break;
            }
        }

        // 节气查找
        w = XL::s_alon(jd2 / 36525.0, 3);
        w = ((w - 0.13) / PI2 * 24.0).floor() * PI2 / 24.0;

        loop {
            // 计算节气时刻
            let d = self.qi_accurate(w);
            let d_int = (d + 0.5).floor() as i32;

            // 计算节气序号
            let xn = ((w / PI2 * 24.0 + 24000006.01).floor() as i32 % 24) as usize;

            w += PI2 / 24.0;

            // 超出本月范围则结束
            if d_int >= bd0 + bdn {
                break;
            }

            // 未到本月则继续
            if d_int < bd0 {
                continue;
            }

            // 更新日对象的节气信息
            let mut day = lunar_month.days_info[(d_int - bd0) as usize].clone() as LunarDay;
            day.jie_qi = JIE_QI_NAMES[xn].to_string();
            day.jie_qi_time = d;
            day.jie_qi_str = jd.time_str(d);

            // 相邻日期超出本月则结束
            if d_int + 12 >= bd0 + bdn {
                break;
            }
        }

        lunar_month
    }

    /// 从公历日期获取当日详细信息,考虑时间和经度因素
    /// # Arguments
    /// * `date` - 公历日期时间(UTC时间)
    /// * `longitude` - 经度(度,东经为正)
    /// # Returns
    /// * `LunarDay` - 农历日期信息
    pub fn get_day_info(&self, date: &DateTime<Utc>, longitude: f64) -> LunarDay {
        // 1. 考虑经度确定真太阳时
        let time_zone = 8.0; // 北京时间为东8区
        let local_time_offset = (longitude - time_zone * 15.0) / 15.0; // 经度修正(小时)

        // 2. 生成儒略日,考虑经度时差
        let mut jd = JD::new();
        jd.year = date.year();
        jd.month = date.month() as i32;
        jd.day = date.day() as i32;
        jd.hour = date.hour() as f64;
        jd.minute = date.minute() as f64;
        jd.second = date.second() as f64;
        // jd.jd(
        //     date.year(),
        //     date.month(),
        //     date.day(),
        //     date.hour() as i32,
        //     date.minute() as i32,
        //     date.second() as f64,
        // );

        // 转换为力学时
        let jde = jd.to_jd() + dt_t(jd.to_jd()) / 24.0;
        // 转换为北京时
        let bj_time = jde - 8.0 / 24.0;
        // 经度修正
        let local_jd = bj_time + local_time_offset / 24.0;

        // 3. 获取儒略日整数部分作为日期
        let d0 = (local_jd.floor() as f64 - J2000).floor() ;

        // 4. 构造基础公历信息
        let mut day = LunarDay::default();
        day.day.d0 = d0 as i32;
        day.day.y = date.year();
        day.day.m = date.month() as i32;
        day.day.d = date.day() as i32;

        // 计算本月天数和周信息
        let first_day = date.with_day(1).unwrap();
        let mut next_month = first_day + chrono::Duration::days(32);
        next_month = next_month.with_day(1).unwrap();

        day.day.dn = (next_month - first_day).num_days() as i32;
        day.day.week0 = (first_day.weekday().num_days_from_sunday() as i32 + 7) % 7;
        day.day.week = (date.weekday().num_days_from_sunday() as i32 + 7) % 7;
        day.day.weeki = ((date.day() as i32 + day.day.week0 - 1) / 7) as i32;
        day.day.weekn = ((day.day.dn + day.day.week0 - 1) / 7 + 1) as i32;

        // 3. 初始化农历历法计算器
        let mut ssq = SolarTerms::new();
        ssq.calc_year(d0 as f64);

        // 4. 计算农历月
        let mut mk = ((d0 as f64 - ssq.new_moons[0]) / 30.0) as usize;
        if mk < 13 && ssq.new_moons[mk + 1] <= d0 {
            mk += 1;
        }

        // 5. 设置农历信息
        day.ldi = (d0 - ssq.new_moons[mk]) as i32;
        day.lunar_day_name = RI_MC[day.ldi as usize].to_string();
        day.lunar_month_name = ssq.month_names[mk].clone();
        day.lunar_month_days = ssq.month_days[mk];
        day.is_leap_month = ssq.leap == (mk as i32);
        day.next_month_name = if mk < 13 {
            ssq.month_names[mk + 1].clone()
        } else {
            "未知".to_string()
        };

        // 6. 计算距节气天数
        day.cur_dz = (d0 - ssq.solar_terms[0]) as i32; // 距冬至天数
        day.cur_xz = (d0 - ssq.solar_terms[12]) as i32; // 距夏至天数
        day.cur_lq = (d0 - ssq.solar_terms[15]) as i32; // 距立秋天数
        day.cur_mz = (d0 - ssq.solar_terms[11]) as i32; // 距芒种天数
        day.cur_xs = (d0 - ssq.solar_terms[13]) as i32; // 距小暑天数

        // 7. 计算干支纪年
        let tmp = if (d0 as f64) < ssq.solar_terms[3] { -365.0 } else { 0.0 };
        let mut d = ssq.solar_terms[3] + tmp + 365.25 * 16.0 - 35.0;

        day.lunar_year = (d / 365.2422 + 0.5).floor() as i32;

        d = ssq.new_moons[2];
        for j in 0..14 {
            if ssq.month_names[j] != "正" || (ssq.leap == j as i32 && j > 0) {
                continue;
            }
            d = ssq.new_moons[j];
            if (d0 as f64) < d {
                d -= 365.0;
                break;
            }
        }

        d = d + 5810.0;
        day.lunar_year0 = (d / 365.2422 + 0.5).floor() as i32;

        // 8. 设置干支纪年
        let d1 = day.lunar_year + 12000;
        day.gan_zhi_year = format!("{}{}", GAN[d1 as usize % 10], ZHI[d1 as usize % 12]);

        let d2 = day.lunar_year0 + 12000;
        day.gan_zhi_year1 = format!("{}{}", GAN[d2 as usize % 10], ZHI[d2 as usize % 12]);

        // 9. 设置其他属性
        day.lunar_year1 = format!("{}", day.lunar_year0 + 1984 + 2698);
        day.sheng_xiao = SHENG_XIAO[d2 as usize % 12].to_string();
        day.nian_hao = self.get_nian_hao(date.year());

        // 10. 计算干支纪月
        let mut mk = ((d0 as f64 - ssq.solar_terms[0]) / 30.43685) as usize;
        if mk < 12 && d0 as f64 >= ssq.solar_terms[2 * mk + 1] {
            mk += 1;
        }

        d = mk as f64 + ((ssq.solar_terms[12] + 390.0) / 365.2422).floor() * 12.0 + 900000.0;
        day.lunar_month = (d as usize % 12) as i32;
        day.gan_zhi_month = format!("{}{}", GAN[d as usize % 10], ZHI[d as usize % 12]);

        // 11. 计算干支纪日
        d = d0 as f64 - 6.0 + 9000000.0;
        day.gan_zhi_day = format!("{}{}", GAN[d as usize % 10], ZHI[d as usize % 12]);

        // 12. 计算星座
        mk = ((d0 as f64 - ssq.solar_terms[0] - 15.0) / 30.43685) as usize;
        if mk < 11 && d0 as f64 >= ssq.solar_terms[2 * mk + 2] {
            mk += 1;
        }
        day.constellation = format!("{}座", CONSTELLATION[(mk + 12) % 12]);

        // 13. 计算回历
        day.hui_li = Some(self.calc_hui_li(d0 as f64));

        let mut festival_info = FestivalInfo::default();
        // 14. 计算节日
        day.festival = festival_info.calc_festivals(&day);

        // 在计算节气和月相时使用精确时刻
        let mk = ((d0 as f64 - ssq.solar_terms[0] - 7.0) / 15.2184) as i32;
        if mk < 23 && local_jd >= ssq.solar_terms[mk as usize + 1] {
            let diff = local_jd - ssq.solar_terms[mk as usize + 1];
            // 如果时间差小于1天,则认为是节气
            if diff.abs() < 1.0 {
                day.jie_qi = JIE_QI_NAMES[mk as usize].to_string();
                day.jie_qi_time = ssq.solar_terms[mk as usize + 1];
                day.jie_qi_str = jd.time_str(day.jie_qi_time);
            }
        }

        // // 月相同理,使用精确时刻比较
        // let moon_phase = ((local_jd - ssq.new_moons[0]) / 29.5306).floor();
        // let phase_time = ssq.new_moons[moon_phase as usize];
        // if (local_jd - phase_time).abs() < 1.0 {
        //     day.moon_phase = MOON_PHASE_NAMES[(moon_phase % 4.0) as usize].to_string();
        //     day.moon_phase_time = phase_time;
        //     day.moon_phase_str = jd.time_str(day.moon_phase_time);
        // }

        // 月相查找 - 使用与 yue_li_calc 相同的算法
        let jd2 = local_jd + dt_t(local_jd); // 转换为力学时
        let t = jd2 / 36525.0; // 转换为儒略世纪数
        let w = XL::ms_alon(t, 10, 3); // 计算月日视黄经差
        let xn = ((w / PI2 * 4.0 + 4000000.01).floor() as i32 % 4) as usize; // 计算月相序号

        // 设置月相信息
        day.moon_phase = MOON_PHASE_NAMES[xn].to_string();

        // 计算精确的月相时刻 (可选)
        // 这里可以使用 so_accurate 方法计算精确的月相时刻
        // 类似于 yue_li_calc 中的实现
        let w_phase = ((w - 0.78) / PI * 2.0).floor() * PI / 2.0;
        let phase_time = self.so_accurate(w_phase);
        if (local_jd - phase_time).abs() < 1.0 {
            day.moon_phase_time = phase_time;
            day.moon_phase_str = jd.time_str(day.moon_phase_time);
        }

        day
    }


    /// 从公历日期获取农历简单信息
    /// # Arguments 
    /// * `date` - 公历日期时间(UTC时间)
    /// * `longitude` - 经度(度,东经为正)
    /// # Returns
    /// * `LunarSimpleDay` - 简单农历日期信息
    pub fn get_simple_day(&self, date: &DateTime<Utc>, longitude: f64) -> LunarSimpleDay {
        // 1. 考虑经度确定真太阳时
        let time_zone = 8.0; // 北京时间为东8区
        let local_time_offset = (longitude - time_zone * 15.0) / 15.0;

        // 2. 生成儒略日
        let mut jd = JD::new();
        jd.year = date.year();
        jd.month = date.month() as i32;
        jd.day = date.day() as i32;
        jd.hour = date.hour() as f64;
        jd.minute = date.minute() as f64;
        jd.second = date.second() as f64;
        // jd.jd(
        //     date.year(),
        //     date.month(),
        //     date.day(),
        //     date.hour() as i32, 
        //     date.minute() as i32,
        //     date.second() as f64
        // );

        // 转换为力学时并进行经度修正
        let jde = jd.to_jd() + dt_t(jd.to_jd()) / 24.0;
        let bj_time = jde - 8.0 / 24.0;
        let local_jd = bj_time + local_time_offset / 24.0;
        let d0 = local_jd.floor() - J2000;

        // 3. 构造基础信息
        let mut day = LunarSimpleDay::default();
        day.day.d0 = d0 as i32;
        day.day.y = date.year();
        day.day.m = date.month() as i32;
        day.day.d = date.day() as i32;

        // 计算月历信息
        let first_day = date.with_day(1).unwrap();
        let mut next_month = first_day + chrono::Duration::days(32);
        next_month = next_month.with_day(1).unwrap();

        day.day.dn = (next_month - first_day).num_days() as i32;
        day.day.week0 = (first_day.weekday().num_days_from_sunday() as i32 + 7) % 7;
        day.day.week = (date.weekday().num_days_from_sunday() as i32 + 7) % 7;
        day.day.weeki = ((date.day() as i32 + day.day.week0 - 1) / 7) as i32;
        day.day.weekn = ((day.day.dn + day.day.week0 - 1) / 7 + 1) as i32;

        // 4. 计算农历信息
        let mut ssq = SolarTerms::new();
        ssq.calc_year(d0 as f64);

        // 5. 找农历月
        let mut mk = ((d0 as f64 - ssq.new_moons[0]) / 30.0) as usize;
        if mk < 13 && ssq.new_moons[mk + 1] <= d0 {
            mk += 1;
        }

        // 6. 设置月历信息
        day.ldi = (d0 - ssq.new_moons[mk]) as i32;
        day.lunar_day_name = RI_MC[day.ldi as usize].to_string();
        day.lunar_month_name = ssq.month_names[mk].clone();
        day.lunar_month_days = ssq.month_days[mk];
        day.is_leap_month = ssq.leap == (mk as i32);
        day.next_month_name = if mk < 13 {
            ssq.month_names[mk + 1].clone()
        } else {
            "未知".to_string()
        };

        // 7. 计算农历年
        let tmp = if (d0 as f64) < ssq.solar_terms[3] { -365.0 } else { 0.0 };
        let mut d = ssq.solar_terms[3] + tmp + 365.25 * 16.0 - 35.0;
        let lunar_year = (d / 365.2422 + 0.5).floor() as i32;

        // 8. 计算节气年
        d = ssq.new_moons[2];
        for j in 0..14 {
            if ssq.month_names[j] != "正" || (ssq.leap == j as i32 && j > 0) {
                continue;
            }
            d = ssq.new_moons[j];
            if (d0 as f64) < d {
                d -= 365.0;
                break;
            }
        }
        d = d + 5810.0;
        let lunar_year0 = (d / 365.2422 + 0.5).floor() as i32;

        // 9. 设置干支纪年
        let d1 = lunar_year + 12000;
        day.gan_zhi_year = format!("{}{}", GAN[d1 as usize % 10], ZHI[d1 as usize % 12]);

        let d2 = lunar_year0 + 12000;
        day.gan_zhi_year1 = format!("{}{}", GAN[d2 as usize % 10], ZHI[d2 as usize % 12]);

        // 10. 计算农历月及其干支
        let mut mk = ((d0 as f64 - ssq.solar_terms[0]) / 30.43685) as usize;
        if mk < 12 && d0 as f64 >= ssq.solar_terms[2 * mk + 1] {
            mk += 1;
        }
        d = mk as f64 + ((ssq.solar_terms[12] + 390.0) / 365.2422).floor() * 12.0 + 900000.0;
        day.lunar_month = (d as usize % 12) as i32;
        day.gan_zhi_month = format!("{}{}", GAN[d as usize % 10], ZHI[d as usize % 12]);

        // 11. 计算日干支
        d = d0 as f64 - 6.0 + 9000000.0;
        day.gan_zhi_day = format!("{}{}", GAN[d as usize % 10], ZHI[d as usize % 12]);

        // 12. 设置生肖
        day.sheng_xiao = SHENG_XIAO[d2 as usize % 12].to_string();

        day
    }

    /// 命理八字计算
    /// # Arguments
    /// * `jd` - 格林尼治UT(J2000起算)
    /// * `longitude` - 本地经度
    /// # Returns
    /// * `BaZiInfo` - 八字信息
    pub fn ming_li_ba_zi(&self, jd: f64, longitude: f64) -> BaZiInfo {
        // 转换为力学时
        let jd2 = jd + dt_t(jd);

        // 计算太阳视黄经
        let w = XL::s_alon(jd2 / 36525.0, -1);

        // 1984年立春起算的节气数(不含中气)
        let k = ((w / PI2 * 360.0 + 45.0 + 15.0 * 360.0) / 30.0).floor() as i32;

        // 计算本地真太阳时
        let mut jd = jd + pty_zty2(jd2 / 36525.0) + longitude / PI / 2.0;

        // 转为前一日23点起算(原jd为本日中午12点起算)
        jd += 13.0 / 24.0;

        // 计算日数与时辰
        let d = jd.floor();
        let sc = ((jd - d) * 12.0).floor() as usize;

        // 计算年月日时的干支
        let mut bazi = BaZiInfo::default();

        // 年干支
        let v = (k / 12 + 6000000) as usize;
        bazi.year = format!("{}{}", GAN[v % 10], ZHI[v % 12]);

        // 月干支
        let v = (k + 2 + 60000000) as usize;
        bazi.month = format!("{}{}", GAN[v % 10], ZHI[v % 12]);

        // 日干支
        let v = ((d - 6.0) as i32 + 9000000) as usize;
        bazi.day = format!("{}{}", GAN[v % 10], ZHI[v % 12]);

        // 时干支
        let mut v = ((d - 1.0) as i32 * 12 + 90000000 + sc as i32) as usize;
        bazi.hour = format!("{}{}", GAN[v % 10], ZHI[v % 12]);

         // 生成全天时辰表
         v -= sc;
         let mut time_table = array::from_fn(|_| String::default());
         for i in 0..13 {
             let gz = format!("{}{}", GAN[(v + i) % 10], ZHI[(v + i) % 12]);
             
             // 存储时辰
             time_table[i] = gz.to_string();
             
             // 记录当前时辰
             if i == sc {
                 bazi.cur_hour = gz;
             }
         }
         bazi.time_table = time_table;

        bazi
    }


    /// 计算指定时间的详细八字信息
    /// # Arguments
    /// * `date` - 公历日期时间(UTC时间)
    /// * `longitude` - 经度(度,东经为正)
    /// # Returns
    /// * `BaZiDetail` - 详细八字信息
    pub fn ming_li_ba_zi_detail(&self, date: &DateTime<Utc>, longitude: f64) -> BaZiDetail {
        let mut detail = BaZiDetail::default();
        
        // 1. 设置基本时间信息
        detail.date = *date;
        detail.year = date.year();
        detail.month = date.month() as i32;
        detail.day = date.day() as i32;
        detail.hour = date.hour() as i32;
        detail.minute = date.minute() as i32;
        detail.second = date.second() as f64;

        // 2. 获取农历日期信息
        detail.day_info = self.get_simple_day(date, longitude);

        // 3. 计算八字信息
        // 转换为儒略日(J2000起算)
        let mut jd = JD::new();
        jd.year = detail.year;
        jd.month = detail.month as i32;
        jd.day = detail.day as i32;
        jd.hour = detail.hour as f64;
        jd.minute = detail.minute as f64;
        jd.second = detail.second as f64;
        // jd.jd(
        //     detail.year,
        //     detail.month as i32,
        //     detail.day as i32,
        //     detail.hour,
        //     detail.minute,
        //     detail.second
        // );
        let j2000_jd = jd.to_jd() - J2000 as f64;
        
        // 计算八字
        detail.bazi = self.ming_li_ba_zi(j2000_jd, longitude);

        detail
    }


    /// 从农历日期获取公历日期
    /// # Arguments
    /// * `year` - 农历年(1900-2100)
    /// * `month` - 农历月(1-12)
    /// * `is_leap` - 是否闰月
    /// * `day` - 农历日(1-30)
    /// * `hour` - 小时(0-23)
    /// * `minute` - 分钟(0-59)  
    /// * `second` - 秒(0-59.999)
    /// * `longitude` - 经度(度,东经为正)
    /// # Returns
    /// * `Option<DateTime<Utc>>` - 公历日期时间,如果参数无效则返回None
    pub fn from_lunar_date(
        &self,
        year: i32,
        month: i32,
        is_leap: bool,
        day: i32,
        hour: i32,
        minute: i32,
        second: f64,
        longitude: f64
    ) -> Option<DateTime<Utc>> {
        // 1. 参数检查
        if !(1900..=2100).contains(&year) 
            || !(1..=12).contains(&month)
            || !(1..=30).contains(&day)
            || !(0..=23).contains(&hour)
            || !(0..=59).contains(&minute)
            || !(0.0..=60.0).contains(&second) {
            return None;
        }

        // 2. 定位农历年
        let mut ssq = SolarTerms::new();
        let year_start = ((year - 2000) as f64 * 365.2422).floor() as i32;
        
        // 计算该年的信息
        ssq.calc_year(year_start as f64);

        // 3. 找到对应的农历月
        let mut month_index = 0;
        let mut found = false;
        
        for i in 0..14 {
            if ssq.month_names[i] == YUE_MC[month as usize - 1] {
                // 检查是否是闰月
                if is_leap {
                    if ssq.leap == (i as i32) {
                        month_index = i;
                        found = true;
                        break;
                    }
                } else {
                    if ssq.leap != (i as i32) {
                        month_index = i;
                        found = true;
                        break;
                    }
                }
            }
        }

        if !found {
            return None;
        }

        // 4. 计算日期对应的儒略日
        let jd = ssq.new_moons[month_index] + day as f64 - 1.0;  // 减1是因为农历初一对应HS值

        // 5. 加上时分秒
        let time_offset = hour as f64 / 24.0 
                       + minute as f64 / 1440.0 
                       + second / 86400.0;

        // 6. 考虑经度修正
        let time_zone = 8.0; // 北京时间为东8区
        let local_time_offset = (longitude - time_zone * 15.0) / 15.0 / 24.0;

        let final_jd = jd + J2000 + time_offset - local_time_offset;

        // 7. 转换为UTC日期时间
        let mut jd_obj = JD::new();
        if jd_obj.set_from_jd(final_jd) {
            Utc.with_ymd_and_hms(
                jd_obj.year,
                jd_obj.month as u32,
                jd_obj.day as u32,
                jd_obj.hour as u32,
                jd_obj.minute as u32,
                jd_obj.second as u32
            )
            // .and_hms_milli_opt(
            //     ((jd_obj.second - jd_obj.second.floor()) * 1000.0) as u32
            // )
            .single()
        } else {
            None
        }
    }

    //下面都是AI生成的辅助代码，暂时不考虑

    // /// 计算二十四节气时刻表
    // fn calc_terms_table(&self, t: f64) -> Vec<f64> {
    //     let mut terms = Vec::with_capacity(25);
    //     for i in 0..25 {
    //         let w = 2.0 * PI * i as f64 / 24.0;
    //         terms.push(eph::solar::term_time(t, w));
    //     }
    //     terms
    // }

    // /// 获取农历日名称(如初一、初二等)
    // fn get_day_name(&self, jd: f64) -> String {
    //     // 农历日名
    //     const RMC: [&str; 31] = [
    //         "初一","初二","初三","初四","初五","初六","初七","初八","初九","初十",
    //         "十一","十二","十三","十四","十五","十六","十七","十八","十九","二十",
    //         "廿一","廿二","廿三","廿四","廿五","廿六","廿七","廿八","廿九","三十","卅一"
    //     ];

    //     let t = (jd - J2000) / 36525.0;
    //     let moon_phase = eph::moon::phase(t);
    //     let days_since_new_moon = ((moon_phase / (2.0 * PI)) * 29.53).round() as usize;

    //     RMC[days_since_new_moon % 30].to_string()
    // }

    // /// 计算农历月份信息
    // fn calc_lunar_month(&self, jd: f64, terms: &[f64]) -> (String, bool, i32, String) {
    //     // 月名
    //     const YMC: [&str; 12] = [
    //         "正","二","三","四","五","六",
    //         "七","八","九","十","冬","腊"
    //     ];

    //     // 计算月首
    //     let mut month_start = self.find_month_start(jd);
    //     let is_leap = self.check_leap_month(month_start, terms);

    //     // 计算月大小
    //     let next_month_start = self.find_month_start(month_start + 35.0); // 约35天后必是下月
    //     let month_days = ((next_month_start - month_start).round()) as i32;

    //     // 确定月名
    //     let month_index = self.calc_month_index(month_start, terms);
    //     let month_name = YMC[month_index].to_string();

    //     // 下月名
    //     let next_month_index = (month_index + 1) % 12;
    //     let next_month = YMC[next_month_index].to_string();

    //     (month_name, is_leap, month_days, next_month)
    // }

    // /// 查找月首(定朔)
    // fn find_month_start(&self, jd: f64) -> f64 {
    //     let t = (jd - J2000) / 36525.0;
    //     let moon_lon = eph::moon::mean_longitude(t);
    //     let sun_lon = eph::solar::mean_longitude(t);

    //     // 月日视黄经差约0时为朔
    //     let delta = moon_lon - sun_lon;
    //     let correction = delta / (2.0 * PI) * 29.53; // 转为日数

    //     jd - correction
    // }

    // /// 检查是否闰月
    // fn check_leap_month(&self, month_start: f64, terms: &[f64]) -> bool {
    //     // 计算本月中气
    //     let mid_month = month_start + 15.0;
    //     let no_zhong_qi = !terms.iter().any(|&t| {
    //         (t - mid_month).abs() < 0.5 // 允许0.5天误差
    //     });

    //     // 无中气且在闰月年份范围内
    //     no_zhong_qi && self.is_leap_year(month_start)
    // }

    // /// 计算月序(子正寅月)
    // fn calc_month_index(&self, month_start: f64, terms: &[f64]) -> usize {
    //     // 以节气推算月序
    //     let t = (month_start - J2000) / 36525.0;
    //     let sun_lon = eph::solar::apparent_longitude(t);

    //     ((sun_lon / (2.0 * PI) * 12.0).round() as usize + 2) % 12
    // }

    // /// 判断是否闰年
    // fn is_leap_year(&self, jd: f64) -> bool {
    //     let year = ((jd - J2000) / 365.2422).floor() + 2000.0;

    //     // 农历闰年规则
    //     let cycle_year = ((year - 1984.0) % 60.0).floor();
    //     [2,5,7,10,13,16,18,21,24,26,29].contains(&(cycle_year as i32))
    // }

    // /// 获取生肖
    // fn get_sheng_xiao(&self, gz_year: &str) -> String {
    //     let zhi = gz_year.chars().nth(1).unwrap();
    //     let idx = zhi.iter().position(|&x| x.chars().next().unwrap() == zhi).unwrap();
    //     sheng_xiao[idx].to_string()
    // }

    // /// 计算星座
    // fn calc_constellation(&self, jd: f64) -> String {
    //     let t = (jd - J2000) / 36525.0;
    //     let sun_lon = eph::solar::apparent_longitude(t);
    //     let idx = ((sun_lon / (2.0 * PI) * 12.0).floor() as usize + 9) % 12;
    //     constellation[idx].to_string()
    // }

    // /// 格式化输出回历日期
    // pub fn format_hui_li(&self, hui_li: &HuiLiDay) -> String {
    //     format!("回历{}年{}月{}日", hui_li.year, hui_li.month, hui_li.day)
    // }

    // /// 格式化输出
    // pub fn format_day(&self, lunar_day: &LunarDay, solar_date: &DateTime<Utc>) -> String {
    //     let festivals = self.calc_festivals(lunar_day, solar_date);

    //     let mut result = String::new();

    //     // 1. 公历日期
    //     result.push_str(&format!("{}年{}月{}日 星期{}\n",
    //         solar_date.year(),
    //         solar_date.month(),
    //         solar_date.day(),
    //         "日一二三四五六"[solar_date.weekday().num_days_from_sunday() as usize..].chars().next().unwrap()
    //     ));

    //     // 2. 农历信息
    //     result.push_str(&format!("农历{}年{}{}月{}\n",
    //         lunar_day.gan_zhi_year,
    //         if lunar_day.is_leap_month { "闰" } else { "" },
    //         lunar_day.lunar_month_name,
    //         lunar_day.lunar_day_name
    //     ));

    //     // 3. 节日信息
    //     if !festivals.major.is_empty() {
    //         result.push_str("节日: ");
    //         result.push_str(&festivals.major.join("、"));
    //         result.push('\n');
    //     }

    //     // 4. 其他信息
    //     if !lunar_day.jie_qi.is_empty() {
    //         result.push_str(&format!("节气: {}\n", lunar_day.jie_qi));
    //     }

    //     // 添加回历输出
    //     if let Some(hui_li) = &lunar_day.hui_li {
    //         result.push_str(&format!("回历: {}\n",
    //             self.format_hui_li(hui_li)));
    //     }

    //     result
    // }

    // /// 完善干支纪年计算
    // fn calc_gz_year(&self, jd: f64, spring_jd: f64) -> String {
    //     // 以立春为界定年首
    //     let year = if jd < spring_jd {
    //         ((jd - J2000) / 365.2422).floor() + 2000.0 - 1.0
    //     } else {
    //         ((jd - J2000) / 365.2422).floor() + 2000.0
    //     };

    //     let cycle_year = ((year - 1864.0) % 60.0) as usize;
    //     let gan_idx = cycle_year % 10;
    //     let zhi_idx = cycle_year % 12;

    //     format!("{}{}", gan[gan_idx], zhi[zhi_idx])
    // }

    // /// 添加农历纪月计算
    // fn calc_month_gan_zhi(&self, jd: f64) -> String {
    //     // 以节气推算月序,1998年12月7(大雪)开始连续进行节气计数
    //     let start_jd = J2000 + 365.0 * 28.0; // 1998-12-07
    //     let months = ((jd - start_jd) / 29.5306).floor() as i32;

    //     let gan_idx = (months % 10) as usize;
    //     let zhi_idx = (months % 12) as usize;

    //     format!("{}{}", gan[gan_idx], zhi[zhi_idx])
    // }

    // /// 添加农历纪日计算
    // fn calc_day_gan_zhi(&self, jd: f64) -> String {
    //     // 以2000年1月7日为甲子日
    //     let days = (jd - J2000 + 6.0).floor() as i32;

    //     let gan_idx = ((days % 10) + 10) as usize % 10;
    //     let zhi_idx = ((days % 12) + 12) as usize % 12;

    //     format!("{}{}", gan[gan_idx], zhi[zhi_idx])
    // }

    // /// 从公历日期获取农历日期
    // /// # Arguments
    // /// * `date` - 公历日期时间
    // /// # Returns
    // /// * `LunarDay` - 农历日期信息
    // pub fn from_solar_date(&self, date: &DateTime<Utc>) -> LunarDay {
    //     // 转换为儒略日
    //     let jd = JD::new();
    //     jd.jd(date.year(), date.month(), date.day());
    //     // 加上时分秒的小数部分
    //     jd += date.hour()/24.0 + date.minute()/1440.0 + date.second()/86400.0;

    //     // 调用已有方法计算农历信息
    //     self.calc_day(jd)
    // }
}
