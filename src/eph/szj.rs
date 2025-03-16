use std::f64::consts::PI;
use crate::eph::{ 
    eph_base::{CS_R_EAR, RAD, PI2, rad2rrad, llr_conv, p_gst,mod2},
    xl::XL,
    jd::JD,
    delta_t::dt_t,
    prece::Prece,
};

/// 日月升降计算结构体
#[derive(Default)]
pub struct SZJ {
    /// 站点地理经度,向东测量为正
    pub l: f64,
    /// 站点地理纬度
    pub fa: f64,
    /// TD-UT
    pub dt: f64, 
    /// 黄赤交角
    pub e: f64,
    /// 多天的升中降结果
    pub rts: Vec<RiseTransitSet>,
}

/// 多天升中降数据结构
#[derive(Default)]
pub struct RiseTransitSet {
    /// 月升时间
    pub ms: String,
    /// 月中天时间
    pub mz: String,
    /// 月落时间 
    pub mj: String,
    /// 日出时间
    pub s: String,
    /// 日中天时间
    pub z: String,
    /// 日落时间
    pub j: String,
    /// 民用晨光时间
    pub c: String,
    /// 民用昏光时间
    pub h: String,
    /// 光照时间 
    pub ch: String,
    /// 昼长时间
    pub sj: String,
}

impl SZJ {
    /// 创建新的SZJ实例
    pub fn new() -> Self {
        Self {
            l: 0.0,
            fa: 0.0,
            dt: 0.0,
            e: 0.409092614,
            rts: Vec::new(),
        }
    }

    /// 计算时角
    /// # Arguments
    /// * `h` - 地平纬度 
    /// * `w` - 赤纬
    /// # Returns
    /// 时角(弧度)
    pub fn get_h(&self, h: f64, w: f64) -> f64 {
        let c = (h.sin() - self.fa.sin() * w.sin()) / 
               (self.fa.cos() * w.cos());
        if c.abs() > 1.0 {
            PI 
        } else {
            c.acos()
        }
    }

    /// 月亮坐标计算
    /// 章动同时影响恒星时和天体坐标,所以不计算章动。返回时角及赤经纬
    /// # Arguments
    /// * `jd` - 儒略日
    /// * `h0` - 是否计算升起时角
    /// * `r` - 结果容器
    pub fn m_coord(&self, jd: f64, h0: bool, r: &mut Coords) {
        // 低精度月亮赤经纬
        let mut z = XL::m_coord((jd + self.dt)/36525.0, 40, 30, 8);
        
        // 转为赤道坐标
        z = llr_conv(&z, self.e);
        
        // 得到此刻天体时角
        r.h = rad2rrad(p_gst(jd, self.dt) + self.l - z[0]);

        // 升起对应的时角
        if h0 {
            r.h0 = self.get_h(
                0.7275 * CS_R_EAR / z[2] - 34.0 * 60.0 / RAD,
                z[1]
            );
        }
    }

    /// 月亮到中升降时刻计算
    /// # Arguments
    /// * `jd` - 儒略日
    /// # Returns
    /// 月亮升中降时刻
    pub fn mt(&mut self, mut jd: f64) -> MoonTransit {
        self.dt = dt_t(jd);
        self.e = Prece::hcjj(jd/36525.0);
        
        // 查找最靠近当日中午的月上中天
        jd -= mod2(0.1726222 + 0.966136808032357 * jd 
            - 0.0366 * self.dt + self.l/PI2, 1.0);

        let sv = PI2 * 0.966;
        let mut r = MoonTransit::new(jd);
        let mut coords = Coords::default();
        
        // 月亮坐标
        self.m_coord(jd, true, &mut coords);
        r.h = coords.h; 
        r.h0 = coords.h0;
        
        // 计算升中降时刻
        r.s += (-r.h0 - r.h)/sv;  // 升
        r.j += (r.h0 - r.h)/sv;   // 降
        r.z += (0.0 - r.h)/sv;    // 中天
        r.x += (PI - r.h)/sv;     // 下中天

        // 多次迭代提高精度
        self.m_coord(r.s, true, &mut coords);
        r.s += rad2rrad(-coords.h0 - coords.h)/sv;
        
        self.m_coord(r.j, true, &mut coords);  
        r.j += rad2rrad(coords.h0 - coords.h)/sv;
        
        self.m_coord(r.z, false, &mut coords);
        r.z += rad2rrad(0.0 - coords.h)/sv;
        
        self.m_coord(r.x, false, &mut coords);
        r.x += rad2rrad(PI - coords.h)/sv;

        r
    }

    /// 计算太阳坐标和时角
    /// 章动同时影响恒星时和天体坐标,所以不计算章动
    /// 
    /// # Arguments
    /// * `jd` - 儒略日数
    /// * `xm` - 模式选择(1-4或10,用于不同深度的地平以下角度计算)
    /// * `coords` - 存储计算结果的结构体
    pub fn s_coord(&self, jd: f64, xm: i32, coords: &mut SolarCoords) {
        // 太阳坐标(修正了光行差)
        let mut z = [
            XL::e_lon((jd + self.dt)/36525.0, 5) + PI - 20.5/RAD,
            0.0,
            1.0
        ];
        
        // 转为赤道坐标
        z = llr_conv(&z, self.e);
        
        // 得到此刻天体时角
        coords.h = rad2rrad(p_gst(jd, self.dt) + self.l - z[0]);
        
        // 计算不同深度的地平角
        if xm == 10 || xm == 1 {
            coords.h1 = self.get_h(-50.0 * 60.0/RAD, z[1]); // 地平以下50分
        }
        if xm == 10 || xm == 2 {  
            coords.h2 = self.get_h(-6.0 * 3600.0/RAD, z[1]); // 地平以下6度
        }
        if xm == 10 || xm == 3 {
            coords.h3 = self.get_h(-12.0 * 3600.0/RAD, z[1]); // 地平以下12度
        }
        if xm == 10 || xm == 4 {
            coords.h4 = self.get_h(-18.0 * 3600.0/RAD, z[1]); // 地平以下18度
        }
    }

    /// 计算太阳到中升降时刻
    /// # Arguments
    /// * `jd` - 当地中午12点时间对应的2000年首起算的格林尼治时间UT
    /// # Returns
    /// 太阳升降时刻数据
    pub fn st(&mut self, mut jd: f64) -> SunTransit {
        self.dt = dt_t(jd);
        self.e = Prece::hcjj(jd/36525.0);
        
        // 查找最靠近当日中午的日上中天
        jd -= mod2(jd + self.l/PI2, 1.0);

        let sv = PI2;
        let mut r = SunTransit::new(jd);
        let mut coords = SolarCoords::default();

        // 太阳坐标
        self.s_coord(jd, 10, &mut coords);
        r.h = coords.h;
        r.h1 = coords.h1;
        r.h2 = coords.h2;
        r.h3 = coords.h3;
        r.h4 = coords.h4;

        // 计算各种时刻
        r.s += (-r.h1 - r.h)/sv;  // 升起
        r.j += (r.h1 - r.h)/sv;   // 降落
        r.c += (-r.h2 - r.h)/sv;  // 民用晨
        r.hh += (r.h2 - r.h)/sv;   // 民用昏
        r.c2 += (-r.h3 - r.h)/sv; // 航海晨  
        r.hh2 += (r.h3 - r.h)/sv;  // 航海昏
        r.c3 += (-r.h4 - r.h)/sv; // 天文晨
        r.hh3 += (r.h4 - r.h)/sv;  // 天文昏
        r.z += (0.0 - r.h)/sv;    // 中天
        r.x += (PI - r.h)/sv;     // 下中天

        // 精确修正
        self.s_coord(r.s, 1, &mut coords);
        r.s += rad2rrad(-coords.h1 - coords.h)/sv;
        if coords.h1 == PI {
            r.sm.push_str("无升起.");
        }

        self.s_coord(r.j, 1, &mut coords);
        r.j += rad2rrad(coords.h1 - coords.h)/sv;
        if coords.h1 == PI {
            r.sm.push_str("无降落.");
        }

        self.s_coord(r.c, 2, &mut coords);
        r.c += rad2rrad(-coords.h2 - coords.h)/sv;
        if coords.h2 == PI {
            r.sm.push_str("无民用晨."); 
        }

        self.s_coord(r.hh, 2, &mut coords);
        r.hh += rad2rrad(coords.h2 - coords.h)/sv;
        if coords.h2 == PI {
            r.sm.push_str("无民用昏.");
        }

        self.s_coord(r.c2, 3, &mut coords);
        r.c2 += rad2rrad(-coords.h3 - coords.h)/sv;
        if coords.h3 == PI {
            r.sm.push_str("无航海晨.");
        }

        self.s_coord(r.hh2, 3, &mut coords);
        r.hh2 += rad2rrad(coords.h3 - coords.h)/sv;
        if coords.h3 == PI {
            r.sm.push_str("无航海昏.");
        }

        self.s_coord(r.c3, 4, &mut coords);
        r.c3 += rad2rrad(-coords.h4 - coords.h)/sv; 
        if coords.h4 == PI {
            r.sm.push_str("无天文晨.");
        }

        self.s_coord(r.hh3, 4, &mut coords);
        r.hh3 += rad2rrad(coords.h4 - coords.h)/sv;
        if coords.h4 == PI {
            r.sm.push_str("无天文昏.");
        }

        self.s_coord(r.z, 0, &mut coords);
        r.z += (0.0 - coords.h)/sv;

        self.s_coord(r.x, 0, &mut coords);
        r.x += rad2rrad(PI - coords.h)/sv;

        r
    }

    /// 计算多天的升中降时刻
    /// # Arguments
    /// * `jd` - 当地起始略日(中午时刻)
    /// * `n` - 计算天数
    /// * `jdl` - 地理经度
    /// * `wdl` - 地理纬度  
    /// * `sq` - 时区
    pub fn calc_rts(&mut self, jd: f64, n: i32, jdl: f64, wdl: f64, mut sq: f64) {
        // 初始化rts数组
        if self.rts.is_empty() {
            for _ in 0..31 {
                self.rts.push(RiseTransitSet::default());
            }
        }

        // 设置站点参数
        self.l = jdl;
        self.fa = wdl;
        sq /= 24.0;

        // 初始化日期
        for i in 0..n {
            let r = &mut self.rts[i as usize];
            r.ms = "--:--:--".to_string();
            r.mz = "--:--:--".to_string();
            r.mj = "--:--:--".to_string();
        }

        let jd_struct = JD::new(); 

        // 计算多天升中降时刻
        for i in -1..=n {
            // 太阳
            if i >= 0 && i < n {
                let r = self.st(jd + i as f64 + sq);
                let idx = i as usize;
                self.rts[idx].s = jd_struct.time_str(r.s - sq);  // 升
                self.rts[idx].z = jd_struct.time_str(r.z - sq);  // 中
                self.rts[idx].j = jd_struct.time_str(r.j - sq);  // 降
                self.rts[idx].c = jd_struct.time_str(r.c - sq);  // 晨
                self.rts[idx].h = jd_struct.time_str(r.h - sq);  // 昏

                // 光照时间,time_str()内部+0.5,所以这里补上-0.5
                self.rts[idx].ch = jd_struct.time_str(r.h - r.c - 0.5);
                
                // 昼长 
                self.rts[idx].sj = jd_struct.time_str(r.j - r.s - 0.5);
            }

            // 月亮
            let r = self.mt(jd + i as f64 + sq);

            // 月升
            let c = ((r.s - sq + 0.5).floor() - jd) as i32;
            if c >= 0 && c < n {
                self.rts[c as usize].ms = jd_struct.time_str(r.s - sq);
            }

            // 月中天
            let c = ((r.z - sq + 0.5).floor() - jd) as i32;
            if c >= 0 && c < n {
                self.rts[c as usize].mz = jd_struct.time_str(r.z - sq);
            }

            // 月落
            let c = ((r.j - sq + 0.5).floor() - jd) as i32;
            if c >= 0 && c < n {
                self.rts[c as usize].mj = jd_struct.time_str(r.j - sq);
            }
        }
    }
}

/// 天体坐标结果
#[derive(Default)] 
pub struct Coords {
    /// 时角
    h: f64,
    /// 升起时角
    h0: f64,
}

/// 太阳坐标计算结果
#[derive(Default)]
pub struct SolarCoords {
    /// 时角
    pub h: f64,
    /// 地平以下50分时的时角
    pub h1: f64,
    /// 地平以下6度时的时角
    pub h2: f64,
    /// 地平以下12度时的时角 
    pub h3: f64,
    /// 地平以下18度时的时角
    pub h4: f64,
}

/// 月亮升中降时刻
#[derive(Default)]
pub struct MoonTransit {
    /// 月升时刻
    pub s: f64,
    /// 月落时刻  
    pub j: f64,
    /// 上中天时刻
    pub z: f64,
    /// 下中天时刻
    pub x: f64,
    /// 时角 
    pub h: f64,
    /// 升起时角
    pub h0: f64,
}

impl MoonTransit {
    /// 创建新实例
    pub fn new(jd: f64) -> Self {
        Self {
            s: jd,
            j: jd, 
            z: jd,
            x: jd,
            h: 0.0,
            h0: 0.0,
        }
    }
}


/// 太阳升降时刻数据
#[derive(Default)]
pub struct SunTransit {
    /// 升起时刻
    pub s: f64,
    /// 降落时刻
    pub j: f64, 
    /// 中天时刻
    pub z: f64,
    /// 下中天时刻
    pub x: f64,
    /// 民用晨时刻
    pub c: f64,
    /// 民用昏时刻 
    pub hh: f64,
    /// 航海晨时刻
    pub c2: f64,
    /// 航海昏时刻
    pub hh2: f64,
    /// 天文晨时刻
    pub c3: f64,
    /// 天文昏时刻
    pub hh3: f64,
    /// 时角
    pub h: f64,
    /// 地平高度角(50分)
    pub h1: f64,
    /// 地平高度角(6度)
    pub h2: f64,
    /// 地平高度角(12度)
    pub h3: f64,
    /// 地平高度角(18度)
    pub h4: f64,
    /// 特殊情况说明
    pub sm: String,
}

impl SunTransit {
    /// 创建新实例
    pub fn new(jd: f64) -> Self {
        Self {
            s: jd,
            j: jd,
            z: jd,
            x: jd,
            c: jd,
            hh: jd, 
            c2: jd,
            hh2: jd,
            c3: jd,
            hh3: jd,
            h: 0.0,
            h1: 0.0,
            h2: 0.0,
            h3: 0.0,
            h4: 0.0,
            sm: String::new(),
        }
    }
}