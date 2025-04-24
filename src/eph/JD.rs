use crate::eph::eph_base::J2000;

/// 日期结构体
#[derive(Debug, Clone, Default)]
pub struct JD {
    pub year: i32,   // 年
    pub month: i32,  // 月
    pub day: i32,    // 日
    pub hour: f64,   // 时
    pub minute: f64, // 分
    pub second: f64, // 秒
}

const WEEKS: [&str; 8] = ["日", "一", "二", "三", "四", "五", "六", "七"];

impl JD {
    /// 创建新的日期实例
    pub fn new(y:i32,m:i32,d:i32,h:f64,i:f64,s:f64) -> Self {
        JD {
            year: y,
            month: m,
            day: d,
            hour: h,
            minute: i,
            second: s,
        }
    }
    
    pub fn default() -> Self {
        JD {
            year: 2000,
            month: 1,
            day: 1,
            hour: 12.0,
            minute: 0.0,
            second: 0.0,
        }
    }

    /// 公历转儒略日
    pub fn jd(&self, mut y: i32, mut m: i32, d: f64) -> f64 {
        let mut n = 0.0;
        let mut g = 0;

        // 判断是否为格里高利历日1582*372+10*31+15(1582年10月15日)
        if (y * 372 + m * 31 + d.floor() as i32) >= 588829 {
            g = 1;
        }

        // 1-2月当作上一年的13-14月
        if m <= 2 {
            m += 12;
            y -= 1;
        }

        // 加百年闰
        if g == 1 {
            let n_tmp = (y as f64 / 100.0).floor();
            n = 2.0 - n_tmp + (n_tmp / 4.0).floor();
        }

        // 计算儒略日
        (365.25 * (y + 4716) as f64).floor() + (30.6001 * (m + 1) as f64).floor() + d + n - 1524.5
    }

    /// 儒略日数转日期
    pub fn dd(&self, jd: f64) -> Vec<f64> {
        let jd = jd + 0.5;
        let z = jd.floor();
        let f = jd - z;

        let z = z as i32;
        let mut a = z;

        if z >= 2299161 {
            let aa = ((z as f64 - 1867216.25) / 36524.25).floor() as i32;
            a = z + 1 + aa - aa / 4;
        }

        let b = a + 1524;
        let c = ((b as f64 - 122.1) / 365.25).floor() as i32; //年数
        let d = (365.25 * c as f64).floor() as i32;  
        let e = ((b - d) as f64 / 30.6001).floor() as i32; //月数

        let dd = (b - d - (30.6001 * e as f64).floor() as i32) as f64; //日数
        let mm;
        let yy;
        if e > 13 {
            mm = e - 13;
            yy = c - 4715;
        }else{
            mm = e -1;
            yy = c - 4716;
        }
        // let mm = if e < 14 { e - 1 } else { e - 13 };
        // let yy = if mm > 2 { c - 4716 } else { c - 4715 };
        // let yy = c;

        //日的小数转为时分秒
        let mut hhh = f * 24.0;
        let mut mmm = (hhh - hhh.floor()) * 60.0;
        let mut sss = ((mmm - mmm.floor()) * 60.0).ceil();
        hhh = hhh.floor();
        mmm = mmm.floor();

        if sss == 60.0 {
            mmm += 1.0;
            sss = 0.0;
            if mmm == 60.0 {
                hhh += 1.0;
                mmm = 0.0;
            }
        }
        //  F*=24; r.h=int2(F); F-=r.h;
        //  F*=60; r.m=int2(F); F-=r.m;
        //  F*=60; r.s=F;

        vec![yy as f64, mm as f64, dd + f, hhh, mmm, sss]
    }

    /// 日期转字符串
    pub fn dd2str(&self, r: &[f64]) -> String {
        let y = r[0].floor() as i32;
        let m = r[1].floor() as i32;
        let d = r[2].floor() as i32;
        let s = r[2] - d as f64;

        let h = (s * 24.0).floor() as i32;
        let min = ((s * 24.0 - h as f64) * 60.0).floor() as i32;
        let sec = ((s * 24.0 - h as f64) * 60.0 - min as f64) * 60.0;

        format!(
            "{:04}-{:02}-{:02} {:02}:{:02}:{:02}",
            y,
            m,
            d,
            h,
            min,
            sec.floor() as i32
        )
    }

    /// 儒略日转字符串
    pub fn jd2str(&self, jd: f64) -> String {
        self.dd2str(&self.dd(jd))
    }

    /// 当前时间转儒略日
    pub fn to_jd(&self) -> f64 {
        let d = self.day as f64 + (self.hour + (self.minute + self.second / 60.0) / 60.0) / 24.0;
        self.jd(self.year, self.month, d)
    }

    /// 儒略日转当前时间
    pub fn set_from_jd(&mut self, jd: f64) -> bool {
        let r = self.dd(jd + J2000);
        if r[1] > 12.0 || r[1] <= 0.0 || r[2].floor() > 31.0 || r[2] <= 0.0 {
            return false;
        }
        self.year = r[0] as i32;
        self.month = r[1] as i32;
        self.day = r[2].floor() as i32;
        // let s = r[2] - self.day as f64;

        // self.hour = (s * 24.0).floor();
        // self.minute = ((s * 24.0 - self.hour) * 60.0).floor();
        // self.second = (((s * 24.0 - self.hour) * 60.0 - self.minute) * 60.0).floor();

        self.hour = r[3];
        self.minute = r[4];
        self.second = r[5];

        true
    }

    /// 获取时间字符串
    pub fn time_str(&self, jd: f64) -> String {
        let r = self.dd(jd);
        let d = r[2].floor();
        let s = r[2] - d;
        let h = (s * 24.0).floor();
        let m = ((s * 24.0 - h) * 60.0).floor();
        format!("{:02}:{:02}", h as i32, m as i32)
    }

    /// 获取星期
    pub fn get_week(&self, jd: f64) -> String {
        let w = (jd + 0.5).floor() as i64 % 7;
        WEEKS[w as usize].to_string()
    }

    /// 获取指定年月的第n个星期w的儒略日数
    pub fn nn_week(&self, y: i32, m: i32, n: i32, w: i32) -> f64 {
        let jd = self.jd(y, m, 1.0);
        let w0 = (jd + 0.5).floor() as i32 % 7;
        let mut r = w - w0;
        if r < 0 {
            r += 7;
        }
        r += 7 * (n - 1);
        jd + r as f64
    }
}
