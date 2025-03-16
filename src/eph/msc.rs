
use crate::eph::delta_t::dt_t;
use crate::eph::nutation::Nutation;
use crate::eph::prece::Prece;
use crate::eph::eph_base::{CS_AU, PI,PI2,RAD,J2000,CS_S_MOON2, p_gst,gxc_moon_lon,gxc_moon_lat,llr_conv,rad2mrad,parallax,mqc,CS_S_MOON,CS_R_EAR,CS_R_EAR_A,suo_n,CS_BA,CS_BA2,CS_K,CS_K0,CS_K2,xyz2llr,rad2rrad,sun_sheng_j,shi_cha_j,p_gst2,gxc_sun_lon,gxc_sun_lat,llr2xyz};
use crate::eph::xl::XL;
use crate::eph::eph::{line_ear,line_ear2,line_ovl,cir_ovl};


    /// 已知t1时刻星体位置、速度，求x*x+y*y=r*r时,t的值
    /// # Arguments
    /// * `g` - 位置速度结构
    /// * `v` - x方向速度
    /// * `u` - y方向速度
    /// * `r` - 半径
    /// * `n` - 正负号标志
    pub fn line_t(g: &PosVel, v: f64, u: f64, r: f64, n: bool) -> f64 {
        let b = g.y * v - g.x * u;
        let a = u * u + v * v;
        let b2 = u * b;
        let c = b * b - r * r * v * v;
        let d = b2 * b2 - a * c;

        if d < 0.0 {
            return 0.0;
        }

        let d = d.sqrt();
        let d = if !n { -d } else { d };
        
        g.t + ((-b2 + d)/a - g.x)/v
    }

///太阳月亮计算类计算相关数据结构
pub struct MSC {
    // 基本参数
    pub t: f64,   // 力学时
    pub l: f64,   // 地理经度
    pub fa: f64,  // 地理纬度
    pub dt: f64,  // TD-UT
    pub jd: f64,  // UT
    pub dl: f64,  // 黄经章动
    pub de: f64,  // 交角章动
    pub e: f64,   // 真黄赤交角
    pub gst: f64, // 真恒星时

    // 月亮参数 
    pub m_hj: f64,  // 月球视黄经
    pub m_hw: f64,  // 月球视黄纬
    pub m_r: f64,   // 地月质心距
    pub m_cj: f64,  // 月球视赤经
    pub m_cw: f64,  // 月球赤纬
    pub m_shi_j: f64, // 月亮时角
    pub m_cj2: f64,  // 视差修正后赤经
    pub m_cw2: f64,  // 视差修正后赤纬 
    pub m_r2: f64,   // 视差修正后地月距
    pub m_dj: f64,   // 月亮方位角
    pub m_dw: f64,   // 月亮高度角
    pub m_pj: f64,   // 大气折射修正后方位角
    pub m_pw: f64,   // 大气折射修正后高度角

    // 太阳参数
    pub s_hj: f64,  // 太阳视黄经 
    pub s_hw: f64,  // 太阳视黄纬
    pub s_r: f64,   // 日地质心距
    pub s_cj: f64,  // 太阳视赤经
    pub s_cw: f64,  // 太阳赤纬
    pub s_shi_j: f64, // 太阳时角
    pub s_cj2: f64,  // 视差修正后赤经
    pub s_cw2: f64,  // 视差修正后赤纬
    pub s_r2: f64,   // 视差修正后日地距
    pub s_dj: f64,   // 太阳方位角
    pub s_dw: f64,   // 太阳高度角 
    pub s_pj: f64,   // 大气折射修正后方位角
    pub s_pw: f64,   // 大气折射修正后高度角

    // 其他参数
    pub sc: f64,    // 时差(日)
    pub pty: f64,   // 平太阳时
    pub zty: f64,   // 真太阳时
    pub m_rad: f64, // 月亮视半径(角秒)
    pub s_rad: f64, // 太阳视半径(角秒)
    pub e_m_rad: f64, // 月亮地心视半径(角秒)
    pub e_shadow: f64, // 地本影半径(角秒)
    pub e_shadow2: f64, // 地半影半径(角秒)
    pub m_ill: f64,  // 月相照明比
    pub zx_j: f64,   // 中心食经度  
    pub zx_w: f64,   // 中心食纬度
}


//太阳月亮计算类
impl MSC {
    /// 创建新的日月食计算实例
    pub fn new() -> Self {
        Self {
            t: 0.0, l: 0.0, fa: 0.0, dt: 0.0, jd: 0.0,
            dl: 0.0, de: 0.0, e: 0.0, gst: 0.0,
            m_hj: 0.0, m_hw: 0.0, m_r: 0.0, m_cj: 0.0, m_cw: 0.0,
            m_shi_j: 0.0, m_cj2: 0.0, m_cw2: 0.0, m_r2: 0.0,
            m_dj: 0.0, m_dw: 0.0, m_pj: 0.0, m_pw: 0.0,
            s_hj: 0.0, s_hw: 0.0, s_r: 0.0, s_cj: 0.0, s_cw: 0.0,
            s_shi_j: 0.0, s_cj2: 0.0, s_cw2: 0.0, s_r2: 0.0,
            s_dj: 0.0, s_dw: 0.0, s_pj: 0.0, s_pw: 0.0,
            sc: 0.0, pty: 0.0, zty: 0.0, m_rad: 0.0, s_rad: 0.0,
            e_m_rad: 0.0, e_shadow: 0.0, e_shadow2: 0.0,
            m_ill: 0.0, zx_j: 0.0, zx_w: 0.0,
        }
    }

    /// 太阳月亮计算
    pub fn calc(&mut self, t: f64, l: f64, fa: f64, high: f64) {
        // 保存基本参数
        self.t = t;
        self.l = l; 
        self.fa = fa;

        // 计算力学时和世界时
        self.dt = dt_t(t); 
        self.jd = t - self.dt;
        
        // 计算章动和交角
        let t = t/36525.0;
        let zd = Nutation::nutation2(t);
        self.dl = zd[0];  // 黄经章动
        self.de = zd[1];  // 交角章动
        self.e = Prece::hcjj(t) + self.de; // 真黄赤交角
        
        // 计算真恒星时
        self.gst = p_gst(self.jd, self.dt) + 
                   self.dl * self.e.cos(); 

        // 计算月亮坐标
        let mut z = XL::m_coord(t, -1, -1, -1);
        z[0] = rad2mrad(z[0] + gxc_moon_lon(t) + self.dl);
        z[1] += gxc_moon_lat(t);
        self.m_hj = z[0];
        self.m_hw = z[1]; 
        self.m_r = z[2];

        // 转换到赤道坐标
        z = llr_conv(&z, self.e);
        self.m_cj = z[0];
        self.m_cw = z[1];

        // 计算月亮时角
        self.m_shi_j = rad2mrad(self.gst + l - z[0]);
        if self.m_shi_j > PI {
            self.m_shi_j -= PI2;
        }

        // 视差修正
        parallax(&mut z, self.m_shi_j, fa, high);
        self.m_cj2 = z[0];
        self.m_cw2 = z[1];
        self.m_r2 = z[2];

        // 转换到地平坐标
        z[0] += PI/2.0 - self.gst - l;
        z = llr_conv(&z, PI/2.0 - fa);
        z[0] = rad2mrad(PI/2.0 - z[0]);
        self.m_dj = z[0];
        self.m_dw = z[1];
        
        // 大气折射修正
        if z[1] > 0.0 {
            z[1] += mqc(z[1]);
        }
        self.m_pj = z[0];
        self.m_pw = z[1];

        // 类似地计算太阳坐标...
        // TODO: 实现太阳坐标计算部分

        // 计算其他参数
        self.calc_other_params(t);
    }

    /// 计算其他参数
    fn calc_other_params(&mut self, t: f64) {
        // 计算时差
        let t = t/10.0;
        let t2 = t*t;
        let t3 = t2*t;
        let t4 = t3*t;
        let t5 = t4*t;

        let lon = (1753470142.0 + 6283319653318.0*t + 
                  529674.0*t2 + 432.0*t3 - 1124.0*t4 - 
                  9.0*t5)/1.0e9 + PI - 20.5/RAD;

        let mut lon = rad2mrad(lon - (self.s_cj - 
                     self.dl * self.e.cos()));
        if lon > PI {
            lon -= PI2;
        }
        self.sc = lon/PI2;

        // 计算真太阳时和平太阳时
        self.pty = self.jd + self.l/PI2;
        self.zty = self.jd + self.l/PI2 + self.sc;

        // 计算视半径等参数
        self.m_rad = CS_S_MOON/self.m_r2;
        self.s_rad = 959.63/self.s_r2;
        self.e_m_rad = CS_S_MOON/self.m_r;
        
        // 计算本影和半影
        self.e_shadow = (CS_R_EAR_A/self.m_r*RAD - 
                       (959.63-8.794)/self.s_r)*51.0/50.0;
        self.e_shadow2 = (CS_R_EAR_A/self.m_r*RAD + 
                        (959.63+8.794)/self.s_r)*51.0/50.0;

        // 计算月相照明
        self.m_ill = XL::moon_ill(t);

        // 计算中心食
        if (rad2rrad(self.m_cj - self.s_cj)).abs() < 
            50.0/180.0*PI {
            let pp = line_ear(
                &[self.m_cj, self.m_cw, self.m_r],
                &[self.s_cj, self.s_cw, self.s_r * CS_AU],
                self.gst
            );
            self.zx_j = pp.j;
            self.zx_w = pp.w;
        } else {
            self.zx_j = 100.0;
            self.zx_w = 100.0;
        }
    }
}




/// 月食快速计算器
#[derive(Default)]
pub struct YsPL {
    /// 月食时间序列(食甚,初亏,复圆,半影食始,半影食终,食既,生光)
    pub lt: Vec<f64>,
    /// 食分
    pub sf: f64,
    /// 食类型
    pub lx: String,
}

/// 位置速度结构 
#[derive(Default)]
pub struct PosVel {
    x: f64,
    y: f64,
    t: f64,
    mr: f64,

    er: f64,
    er2: f64,
    e_m_rad: f64, // 月亮地心视半径(角秒)
    e_shadow: f64,
    e_shadow2: f64,
}

impl PosVel {
    pub fn new(x: f64, y: f64, t: f64) ->Self {
        Self {
            x: x,
            y: y,
            t: t,
            mr: 0.0,
            er: 0.0,
            er2: 0.0,
            e_m_rad: 0.0,
            e_shadow: 0.0,
            e_shadow2: 0.0,
        }
    }
}

//月食批量快速计算
impl YsPL {
    /// 创建新的月食计算器实例
    pub fn new() -> Self {
        Self {
            lt: vec![0.0; 7],
            sf: 0.0, 
            lx: String::new()
        }
    }

    // /// 已知t1时刻星体位置、速度，求x*x+y*y=r*r时,t的值
    // fn line_t(&self, g: &PosVel, v: f64, u: f64, r: f64, n: bool) -> f64 {
    //     let b = g.y * v - g.x * u;
    //     let a = u * u + v * v;
    //     let b2 = u * b;
    //     let c = b * b - r * r * v * v;
    //     let d = b2 * b2 - a * c;
        
    //     if d < 0.0 {
    //         return 0.0;
    //     }
        
    //     let d = d.sqrt();
    //     let d = if !n { -d } else { d };
    //     g.t + ((-b2 + d)/a - g.x)/v
    // }

    /// 日月黄经纬差转为日面中心直角坐标(用于月食)
    fn lec_xy(&self, jd: f64) -> PosVel {
        let mut re = PosVel::default();
        let t = jd/36525.0;
        
        // 太阳月亮黄道坐标
        let mut zs = XL::e_coord(t, -1, -1, -1); // 地球坐标
        zs[0] = rad2mrad(zs[0] + PI + gxc_sun_lon(t));
        zs[1] = -zs[1] + gxc_sun_lat(t); // 补上太阳光行差
        
        let mut zm = XL::m_coord(t, -1, -1, -1); // 月球坐标 
        zm[0] = rad2mrad(zm[0] + gxc_moon_lon(t));
        zm[1] += gxc_moon_lat(t); // 补上月球光行差
        
        // 视半径
        re.e_m_rad = CS_S_MOON/zm[2]; // 月亮地心视半径(角秒)
        re.e_shadow = (CS_R_EAR_A/zm[2]*RAD-(959.63-8.794)/zs[2])*51.0/50.0; 
        re.e_shadow2 = (CS_R_EAR_A/zm[2]*RAD+(959.63+8.794)/zs[2])*51.0/50.0;
        
        re.x = rad2rrad(zm[0] + PI - zs[0]) * ((zm[1]-zs[1])/2.0).cos();
        re.y = zm[1] + zs[1];
        re.mr = re.e_m_rad/RAD;
        re.er = re.e_shadow/RAD;
        re.er2 = re.e_shadow2/RAD; 
        re.t = jd;
        
        re
    }

    /// 月食的食甚计算(jd为近朔的力学时,误差几天不要紧)
    pub fn lec_max(&mut self, mut jd: f64) {
        // 初始化
        self.lt.iter_mut().for_each(|x| *x = 0.0);
        self.sf = 0.0;
        self.lx.clear();

        // 低精度的朔(误差10分钟)
        jd = XL::ms_alon_t2(((jd-4.0)/29.5306).floor() * PI2 + PI) * 36525.0;

        // 求极值(平均误差数秒)
        let u = -18461.0 * (0.057109 + 0.23089571958*jd).sin() * 0.23090/RAD;
        let v = (XL::m_v(jd/36525.0) - XL::e_v(jd/36525.0))/36525.0;
        
        let mut g = self.lec_xy(jd);
        jd -= (g.y*u + g.x*v)/(u*u + v*v);

        // 精密求极值
        let dt = 60.0/86400.0;
        g = self.lec_xy(jd);
        let g2 = self.lec_xy(jd + dt);
        
        let u = (g2.y - g.y)/dt;
        let v = (g2.x - g.x)/dt;
        let dt = -(g.y*u + g.x*v)/(u*u + v*v);
        jd += dt;

        // 求直线到影子中心的最小值 
        let x = g.x + dt*v;
        let y = g.y + dt*u;
        let rmin = (x*x + y*y).sqrt();

         // 判断食相
         if rmin <= g.mr + g.er {
            // 偏食计算
            self.lt[1] = jd;
            self.lx = "偏".to_string(); 
            self.sf = (g.mr + g.er - rmin)/g.mr/2.0;

            // 计算初亏
            self.lt[0] = line_t(&g, v, u, g.mr + g.er, false);
            let mut g3 = self.lec_xy(self.lt[0]);
            self.lt[0] = line_t(&g3, v, u, g3.mr + g3.er, false);

            // 计算复圆
            self.lt[2] = line_t(&g, v, u, g.mr + g.er, true);
            g3 = self.lec_xy(self.lt[2]);
            self.lt[2] = line_t(&g3, v, u, g3.mr + g3.er, true);
        }

        // 半影食计算
        if rmin <= g.mr + g.er2 {
            // 计算半影食始
            self.lt[3] = line_t(&g, v, u, g.mr + g.er2, false);
            let mut g3 = self.lec_xy(self.lt[3]);
            self.lt[3] = line_t(&g3, v, u, g3.mr + g3.er2, false);

            // 计算半影食终
            self.lt[4] = line_t(&g, v, u, g.mr + g.er2, true);
            g3 = self.lec_xy(self.lt[4]);
            self.lt[4] = line_t(&g3, v, u, g3.mr + g3.er2, true);
        }

        // 全食计算
        if rmin <= g.er - g.mr {
            self.lx = "全".to_string();
            
            // 计算食既
            self.lt[5] = line_t(&g, v, u, g.er - g.mr, false);
            let mut g3 = self.lec_xy(self.lt[5]);
            self.lt[5] = line_t(&g3, v, u, g3.er - g3.mr, false);

            // 计算生光
            self.lt[6] = line_t(&g, v, u, g.er - g.mr, true);
            g3 = self.lec_xy(self.lt[6]);
            self.lt[6] = line_t(&g3, v, u, g3.er - g3.mr, true);
        }
    }
}

/// 日食计算结果
pub struct EclipseResult {
    /// 朔时刻（J2000起算的儒略日）
    pub jd: f64,
    /// 朔时刻（同jd）
    pub jd_suo: f64,
    /// 切迹情况
    pub ac: i32,
    /// 日食类型(N无食,P偏食,T全食,A环食)
    pub lx: String,
}

/// 快速日食计算器
pub fn ec_fast(jd: f64) -> EclipseResult {
    let mut re = EclipseResult {
        jd: 0.0,
        jd_suo: 0.0,
        ac: 1,
        lx: "N".to_string()
    };

    // 合朔时的日月黄经差
    let w = ((jd + 8.0)/29.5306).floor() * PI2;

    // 合朔时间计算(2000前后4000年误差1小时内)
    let mut t = (w + 1.08472)/7771.37714500204; // 平朔时间
    re.jd =  t * 36525.0;
    re.jd_suo = re.jd;

    let mut t2 = t*t;
    let mut t3 = t2*t;
    let t4 = t3*t;
    let l = (93.2720993 + 483202.0175273*t - 0.0034029*t2 
             - t3/3526000.0 + t4/863310000.0)/180.0*PI;

    // 一般大于21度已不可能有食
    if l.sin().abs() > 0.4 {
        return re;
    }

    // 精确合朔时间计算
    t -= (-0.0000331*t*t + 0.10976 * (0.785 + 8328.6914*t).cos())/7771.0;
    // let t2 = t*t;
    
    // 计算合朔时黄经差L
    let mut l = -1.084719 + 7771.377145013*t - 0.0000331*t2;
    l += (22640.0 * (0.785 + 8328.6914*t + 0.000152*t2).cos()
        + 4586.0 * (0.19 + 7214.063*t - 0.000218*t2).cos()
        + 2370.0 * (2.54 + 15542.754*t - 0.000070*t2).cos()
        + 769.0 * (3.1 + 16657.383*t).cos()
        + 666.0 * (1.5 + 628.302*t).cos()
        + 412.0 * (4.8 + 16866.93*t).cos()
        + 212.0 * (4.1 - 1114.63*t).cos()
        + 205.0 * (0.2 + 6585.76*t).cos()
        + 192.0 * (4.9 + 23871.45*t).cos()
        + 165.0 * (2.6 + 14914.45*t).cos()
        + 147.0 * (5.5 - 7700.39*t).cos()
        + 125.0 * (0.5 + 7771.38*t).cos()
        + 109.0 * (3.9 + 8956.99*t).cos()
        + 55.0 * (5.6 - 1324.18*t).cos()
        + 45.0 * (0.9 + 25195.62*t).cos()
        + 40.0 * (3.8 - 8538.24*t).cos()
        + 38.0 * (4.3 + 22756.82*t).cos()
        + 36.0 * (5.5 + 24986.07*t).cos()
        - 6893.0 * (4.669257 + 628.3076*t).cos()
        - 72.0 * (4.6261 + 1256.62*t).cos()
        - 43.0 * (2.67823 + 628.31*t).cos() * t
        + 21.0) / RAD;

    // 迭代计算精确合朔时间    
    t += (w - l) / (7771.38 
        - 914.0 * (0.7848 + 8328.691425*t + 0.0001523*t2).sin()
        - 179.0 * (2.543 + 15542.7543*t).sin()
        - 160.0 * (0.1874 + 7214.0629*t).sin());

    re.jd =  t * 36525.0;
    re.jd_suo = re.jd;

    // 计算月球黄纬、距离和速度

    //纬 52,15 (角秒)
    t2=t*t/10000.0;
    t3=t2*t/10000.0;
    let mut mb = 
        18461.0 * (0.0571 + 8433.46616*t - 0.640*t2 - 1.0*t3).cos()
        + 1010.0 * (2.413 + 16762.1576*t + 0.88*t2 + 25.0*t3).cos()
        + 1000.0 * (5.440 - 104.7747*t + 2.16*t2 + 26.0*t3).cos()
        + 624.0 * (0.915 + 7109.2881*t + 0.0*t2 + 7.0*t3).cos()
        + 199.0 * (1.82 + 15647.529*t - 2.8*t2 - 19.0*t3).cos()
        + 167.0 * (4.84 - 1219.403*t - 1.5*t2 - 18.0*t3).cos()
        + 117.0 * (4.17 + 23976.220*t - 1.3*t2 + 6.0*t3).cos()
        + 62.0 * (4.8 + 25090.849*t + 2.0*t2 + 50.0*t3).cos()
        + 33.0 * (3.3 + 15437.980*t + 2.0*t2 + 32.0*t3).cos()
        + 32.0 * (1.5 + 8223.917*t + 4.0*t2 + 51.0*t3).cos()
        + 30.0 * (1.0 + 6480.986*t + 0.0*t2 + 7.0*t3).cos()
        + 16.0 * (2.5 - 9548.095*t - 3.0*t2 - 43.0*t3).cos()
        + 15.0 * (0.2 + 32304.912*t + 0.0*t2 + 31.0*t3).cos()
        + 12.0 * (4.0 + 7737.590*t).cos()
        + 9.0 * (1.9 + 15019.227*t).cos()
        + 8.0 * (5.4 + 8399.709*t).cos() 
        + 8.0 * (4.2 + 23347.918*t).cos()
        + 7.0 * (4.9 - 1847.705*t).cos()
        + 7.0 * (3.8 - 16133.856*t).cos()
        + 7.0 * (2.7 + 14323.351*t).cos();

    mb /= RAD;

    
     // 距离(106, 23千米)
     let mut mr = 385001.0 
        + 20905.0 * (5.4971 + 8328.691425*t + 1.52*t2 + 25.0*t3).cos()
        + 3699.0 * (4.900 + 7214.06287*t - 2.18*t2 - 19.0*t3).cos()
        + 2956.0 * (0.972 + 15542.75429*t - 0.66*t2 + 6.0*t3).cos()
        + 570.0 * (1.57 + 16657.3828*t + 3.0*t2 + 50.0*t3).cos()
        + 246.0 * (5.69 - 1114.6286*t - 3.7*t2 - 44.0*t3).cos()
        + 205.0 * (1.02 + 14914.4523*t - 1.0*t2 + 6.0*t3).cos()
        + 171.0 * (3.33 + 23871.4457*t + 1.0*t2 + 31.0*t3).cos()
        + 152.0 * (4.94 + 6585.761*t - 2.0*t2 - 19.0*t3).cos()
        + 130.0 * (0.74 - 7700.389*t - 2.0*t2 - 25.0*t3).cos()
        + 109.0 * (5.20 + 7771.377*t).cos()
        + 105.0 * (2.31 + 8956.993*t + 1.0*t2 + 25.0*t3).cos()
        + 80.0 * (5.38 - 8538.241*t + 2.8*t2 + 26.0*t3).cos()
        + 49.0 * (6.24 + 628.302*t).cos()
        + 35.0 * (2.7 + 22756.817*t - 3.0*t2 - 13.0*t3).cos()
        + 31.0 * (4.1 + 16171.056*t - 1.0*t2 + 6.0*t3).cos()
        + 24.0 * (1.7 + 7842.365*t - 2.0*t2 - 19.0*t3).cos()
        + 23.0 * (3.9 + 24986.074*t + 5.0*t2 + 75.0*t3).cos()
        + 22.0 * (0.4 + 14428.126*t - 4.0*t2 - 38.0*t3).cos()
        + 17.0 * (2.0 + 8399.679*t).cos();
    
    mr /= 6378.1366; // 转换单位

    // 转换时间尺度
    let t = jd/365250.0;
    let t2 = t*t;
    let t3 = t2*t;

    let mut sr = 10001399.0 // 日地距离基值
        + 167070.0 * (3.098464 + 6283.07585*t).cos()
        + 1396.0 * (3.0552 + 12566.1517*t).cos()
        + 10302.0 * (1.10749 + 6283.07585*t).cos() * t
        + 172.0 * (1.064 + 12566.152*t).cos() * t
        + 436.0 * (5.785 + 6283.076*t).cos() * t2
        + 14.0 * (4.27 + 6283.08*t).cos() * t3;
    
    // 转换到地球半径单位
    sr *= 1.49597870691/6378.1366*10.0;

    // t = jd/36525.0;
    // 月日黄经差速度
    let vl = (7771.0 
        - 914.0 * (0.785 + 8328.6914*t).sin()
        - 179.0 * (2.543 + 15542.7543*t).sin()
        - 160.0 * (0.187 + 7214.0629*t).sin()
    ) / 36525.0;

    // 月亮黄纬速度
    let vb = (-755.0 * (0.057 + 8433.4662*t).sin()
        - 82.0 * (2.413 + 16762.1576*t).sin()
    ) / 36525.0;

    // 月亮径向速度
    let vr = (-27299.0 * (5.497 + 8328.691425*t).sin()
        - 4184.0 * (4.900 + 7214.06287*t).sin()
        - 7204.0 * (0.972 + 15542.75429*t).sin()
    ) / 36525.0;
    

    // 几何参数计算
    //gm伽马值,smR日月距
    let gm = mr * mb.sin() * vl / (vb*vb + vl*vl).sqrt();
    let smr = sr - mr;
    let mk = 0.2725076;
    let sk = 109.1222;

    // 计算本影半影锥角和半径
    let f1 = (sk + mk)/smr; //tanf1半影锥角
    let r1 = mk + f1*mr; // 半影半径
    let f2 = (sk - mk)/smr; //tanf2本影锥角
    let r2 = mk - f2*mr; // 本影半径

    let b = 0.9972;
    let agm = gm.abs();
    let ar2 = r2.abs();
    let fh2 = mr - mk/f2;
    let h = if agm < 1.0 { (1.0 - gm*gm).sqrt() } else { 0.0 };

    // 基本类型判断
    if fh2 < h {
        re.lx = "T".to_string();
    } else {
        re.lx = "A".to_string();
    }

    // 计算特征量
    let ls1 = agm - (b + r1);  // 无食分界
    let ls2 = agm - (b + ar2); // 偏食分界 
    let ls3 = agm - b;         // 无中心食分界
    let ls4 = agm - (b - ar2); // 有中心食分界

    // 判断切迹情况 
    if ls1.abs() < 0.016 { re.ac = 0; }
    if ls2.abs() < 0.016 { re.ac = 0; }
    if ls3.abs() < 0.016 { re.ac = 0; }
    if ls4.abs() < 0.016 { re.ac = 0; }

    // 判断食相类型
    if ls1 > 0.0 {
        re.lx = "N".to_string();          // 无日食
    } else if ls2 > 0.0 {
        re.lx = "P".to_string();          // 偏食
    } else if ls3 > 0.0 { 
        re.lx.push('0');                  // 无中心
    } else if ls4 > 0.0 {
        re.lx.push('1');                  // 有中心(本影未全部进入)
    } else {
        // 本影全进入情况
        if (fh2 - h).abs() < 0.019 {
            re.ac = 0;
        }
        if fh2.abs() < h {
            let dr = vr * h / vl / mr;
            let h1 = mr - dr - mk/f2;     // 入点影锥z坐标
            let h2 = mr + dr - mk/f2;     // 出点影锥z坐标
            
            // 根据入点出点位置判断环食类型
            if h1 > 0.0 && h2 > 0.0 {
                re.lx = "H".to_string();   // 环全环
            } else if h1 > 0.0 {
                re.lx = "H3".to_string();  // 环全全
            } else if h2 > 0.0 {
                re.lx = "H2".to_string();  // 全全环
            }

            // 判断临界情况
            if h1.abs() < 0.019 { re.ac = 0; }
            if h2.abs() < 0.019 { re.ac = 0; }
        }
    }

    re 
}



/// 太阳系计算结构
pub struct RsGS {
    /// 日月赤道坐标插值表
    zs: Vec<f64>,
    /// 插值点之间的时间间距
    zdt: f64,
    /// 插值表中心时间
    zjd: f64,
    /// TD-UT时差
    dt: f64,
    /// 半影锥角
    tanf1: f64,
    /// 本影锥角
    tanf2: f64,
    /// 太阳视半径
    srad: f64,
    /// 贝圆极赤比 
    bba: f64,
    /// 黄交线与赤交线的夹角
    bhc: f64,
    /// 地月距
    dyj: f64
}


/// 速度参数结构
#[derive(Default)]
pub struct VelocityParams {
    /// 地方x向速度
    pub vx: f64,
    /// 地方y向速度  
    pub vy: f64,
    /// 相对x向速度
    pub vx_rel: f64,
    /// 相对y向速度
    pub vy_rel: f64,
    /// 合成速度
    pub v: f64,
}

/// 影锥参数结构
#[derive(Default, Clone)] 
pub struct ShadowParams {
    /// 半影半径
    pub r1: f64,
    /// 本影半径
    pub r2: f64,
    /// 本影半径绝对值
    pub ar2: f64,
    /// 食分
    pub sf: f64,
}

/// 日食特征数据结构
#[derive(Default)]
pub struct EclipseFeature {
    pub jd_suo: f64,     // 朔时刻
    pub jd: f64,         // 中点时刻
    pub d_t: f64,        // deltaT
    pub ds: f64,         // 黄交线与赤交线夹角
    pub vx: f64,         // 影速x分量
    pub vy: f64,         // 影速y分量  
    pub ax: f64,         // 加速度x分量
    pub ay: f64,         // 加速度y分量
    pub v: f64,          // 速度
    pub k: f64,          // 斜率
    pub xc: f64,         // 中点x坐标
    pub yc: f64,         // 中点y坐标 
    pub zc: f64,         // 中点z坐标
    pub d: f64,          // 到圆心距离
    pub i: Vec<f64>,     // 贝塞尔坐标参数
    pub gk1: [f64;3],   // 中心始参数
    pub gk2: [f64;3],   // 中心终参数
    pub gk3: [f64;3],   // 偏食始参数 
    pub gk4: [f64;3],   // 偏食终参数
    pub gk5: Vec<f64>,   // 视午参数
    pub zx_j: f64,       // 中心食经度
    pub zx_w: f64,       // 中心食纬度
    pub sf: f64,         // 食分
    pub sdp: Vec<f64>,   // 太阳地平坐标
    pub lx: String,      // 食型
    pub dw: f64,         // 食带宽度
    pub tt: f64,         // 持续时间
}

/// 界线数据结构
#[derive(Default)]
pub struct BoundaryData {
    /// 存储界线点的数组
    pub points: Vec<f64>,
    /// 当前是否有解
    pub f: i32,
    /// 上次是否有解
    pub f2: Option<i32>,
}

/// 日食界线数据结构
#[derive(Default)]
pub struct EclipseLines {
    /// 初亏复圆界线
    pub p1: Vec<f64>,
    pub p2: Vec<f64>,
    pub p3: Vec<f64>, 
    pub p4: Vec<f64>,

    /// 日出日没食甚线
    pub q1: Vec<f64>,
    pub q2: Vec<f64>,
    pub q3: Vec<f64>,
    pub q4: Vec<f64>,

    /// 南北界线
    pub l0: Vec<f64>, // 中心线
    pub l1: BoundaryData, // 半影北界
    pub l2: BoundaryData, // 半影南界
    pub l3: BoundaryData, // 本影北界
    pub l4: BoundaryData, // 本影南界
    pub l5: BoundaryData, // 0.5食分北界 
    pub l6: BoundaryData  // 0.5食分南界
}

/// 日食界线数据结构(第二种)
#[derive(Default)]
pub struct EclipseLines2 {
    /// 本影界线
    pub p1: Vec<f64>,
    /// 半影界线
    pub p2: Vec<f64>,
    /// 晨昏圈
    pub p3: Vec<f64>,
}

/// 日食界线数据结构(第三种)
#[derive(Default, Clone)]
pub struct EclipseLineData {
    /// 时间(儒略日)
    pub time: f64,
    /// 五条界线的点坐标[经度,纬度]
    pub points: Vec<[f64; 2]>,
}

impl RsGS {
    /// 创建新实例
    pub fn new() -> Self {
        Self {
            zs: Vec::new(),
            zdt: 0.04,
            zjd: 0.0,
            dt: 0.0,
            tanf1: 0.0046,
            tanf2: 0.0045,
            srad: 0.0046,
            bba: 1.0,
            bhc: 0.0,
            dyj: 23500.0
        }
    }

    /// 创建插值表(根数表)
    /// # Arguments
    /// * `jd` - 儒略日
    /// * `n` - 插值点数量(2/3/7分别对应低/中/高精度)
    pub fn init(&mut self, mut jd: f64, n: usize) {
        // 复用已有插值表
        if suo_n(jd) == suo_n(self.zjd) && self.zs.len() == n * 9 {
            return;
        }

        // 清空插值表
        self.zs.clear();
        
        // 计算低精度的朔(误差10分钟)
        jd = XL::ms_alon_t2(((jd - 4.0)/29.5306).floor() * PI2 + PI) * 36525.0;
        self.zjd = jd;
        
        // 计算TD-UT时差
        self.dt = dt_t(jd);

        // 计算章动和交角
        let t = jd/36525.0;
        let zd = Nutation::nutation2(t);
        let e = Prece::hcjj(t) + zd[1]; 

        // 创建插值点
        for i in 0..n {
            // 计算插值点时间
            let t = (self.zjd + (i as f64 - n as f64/2.0 + 0.5) * self.zdt) / 36525.0;
            
            // 计算地球和月球坐标
            let (mut s, mut m) = match n {
                7 => (XL::e_coord(t,-1,-1,-1), XL::m_coord(t,-1,-1,-1)), // 高精度
                3 => (XL::e_coord(t,65,65,65), XL::m_coord(t,-1,150,150)), // 中精度
                _ => (XL::e_coord(t,20,20,20), XL::m_coord(t,30,30,30))  // 低精度
            };

            // 补充光行差和章动
            s[0] += zd[0] + gxc_sun_lon(t) + PI;
            s[1] = -s[1] + gxc_sun_lat(t);
            m[0] += zd[0] + gxc_moon_lon(t);
            m[1] += gxc_moon_lat(t);

            // 转换到赤道坐标
            s = llr_conv(&s, e);
            m = llr_conv(&m, e);
            s[2] *= CS_AU;

            // 确保数据连续性
            let k = i * 9;
            if i > 0 {
                if s[0] < self.zs[k-9] {
                    s[0] += PI2;
                }
                if m[0] < self.zs[k-6] {
                    m[0] += PI2;  
                }
            }

            // 存入插值表
            self.zs.extend_from_slice(&[
                s[0], s[1], s[2],
                m[0], m[1], m[2]
            ]);

            // 计算贝塞尔坐标z轴
            let s_xyz = llr2xyz(&s);
            let m_xyz = llr2xyz(&m);
            let b = xyz2llr(&[
                s_xyz[0] - m_xyz[0],
                s_xyz[1] - m_xyz[1], 
                s_xyz[2] - m_xyz[2]
            ]);
            
            let mut b0 = PI/2.0 + b[0];
            let b1 = PI/2.0 - b[1];
            
            if i > 0 && b0 < self.zs[k-3] {
                b0 += PI2;
            }

            // 存入贝塞尔坐标
            self.zs.extend_from_slice(&[
                b0, 
                b1,
                p_gst(t*36525.0 - self.dt, self.dt) + zd[0] * e.cos()
            ]);
        }

        // 计算辅助参数
        let p = self.zs.len() - 9;
        self.dyj = (self.zs[2] + self.zs[p+2] - self.zs[5] - self.zs[p+5])/ 
            (2.0 * CS_R_EAR);
        self.tanf1 = (CS_K0 + CS_K)/self.dyj;
        self.tanf2 = (CS_K0 - CS_K2)/self.dyj;
        self.srad = CS_K0/((self.zs[2] + self.zs[p+2])/(2.0 * CS_R_EAR));
        
        let sin_b = ((self.zs[1] + self.zs[p+1])/2.0).sin();
        self.bba = CS_BA * (1.0 + (1.0 - CS_BA2) * sin_b * sin_b/2.0);
        
        self.bhc = -(e.tan() * ((self.zs[6] + self.zs[p+6])/2.0).sin()).atan();
    }


    /// 日月坐标快速计算(贝塞尔插值法)
    /// # Arguments 
    /// * `jd` - 儒略日
    /// * `xt` - 要计算的元素类型(0日,1月,2贝塞尔坐标)
    /// # Returns
    /// 计算结果数组[3]
    pub fn chazhi(&self, jd: f64, xt: usize) -> Vec<f64> {
        let p = (xt * 3) as usize; // 计算第p个根数开始
        let m = 3; // 计算m个根数
        
        // 插值要素准备
        let n = self.zs.len() / 9; // 插值点数
        let w = self.zs.len() / n; // 每节点个数
        
        // 相对于第一点的时间距离
        let mut t = (jd - self.zjd) / self.zdt + n as f64 / 2.0 - 0.5;
        
        // 存储计算结果
        let mut z = vec![0.0; m];

        // 二点插值
        if n == 2 {
            for i in 0..m {
                z[i] = self.zs[p+i] + (self.zs[p+ i  + w ] - self.zs[p+i]) * t;
            }
            return z;
        }

        // 确定插值中心c,并限制范围
        let mut c = (t + 0.5).floor() as usize;
        if c <= 0 {
            c = 1;
        }
        if c > n - 2 {
            c = n - 2;
        }

        // 计算插值因子
        t -= c as f64;
        let p = p + c * w;

        // 三点插值计算
        for i in 0..m {
            let idx = p + i;
            z[i] = self.zs[idx] + 
                (self.zs[idx+w] - self.zs[idx-w] + 
                 (self.zs[idx+w] + self.zs[idx-w] - self.zs[idx]*2.0) * t) * t/2.0;
        }

        z
    }

    /// 获取太阳坐标
    pub fn sun(&self, jd: f64) -> Vec<f64> {
        self.chazhi(jd, 0) // 返回值可能超过360度
    }

    /// 获取月亮坐标
    pub fn moon(&self, jd: f64) -> Vec<f64> {
        self.chazhi(jd, 1)
    }

    /// 获取贝塞尔坐标
    pub fn bse(&self, jd: f64) -> Vec<f64> {
        self.chazhi(jd, 2)
    }

    /// 赤道转贝塞尔坐标
    pub fn cd2bse(&self, z: &[f64], i: &[f64]) -> [f64; 3] {
        let mut r = [z[0] - i[0], z[1], z[2]];
        r = llr_conv(&r, -i[1]);
        llr2xyz(&r)
    }

    /// 贝塞尔转赤道坐标
    pub fn bse2cd(&self, z: &[f64], i: &[f64]) -> [f64; 3] {
        let mut r = xyz2llr(z);
        r = llr_conv(&r, i[1]); 
        r[0] = rad2mrad(r[0] + i[0]);
        r
    }

    /// 贝塞尔转地标
    /// p点到原点连线与地球的交点
    /// f=1时把地球看成椭球
    pub fn bse2db(&self, z: &[f64], i: &[f64], f: bool) -> [f64; 3] {
        let mut r = xyz2llr(z);
        r = llr_conv(&r, i[1]);
        r[0] = rad2rrad(r[0] + i[0] - i[2]);
        if f {
            r[1] = (r[1].tan() / CS_BA2).atan();
        }
        r
    }

    /// 贝塞尔转地标
    /// 过p点垂直于基面的线与地球的交点
    /// p坐标为(x,y,任意z),f=1时把地球看成椭球
    pub fn bse_xy2db(&self, x: f64, y: f64, i: &[f64], f: bool) -> Vec<f64> {
        let b = if f { CS_BA } else { 1.0 };
        let p = line_ear2(x, y, 2.0, x, y, 0.0, b, 1.0, i);
        vec![p.j, p.w]
    }

    /// 计算月亮的贝塞尔坐标
    /// # Arguments
    /// * `jd` - 儒略日
    /// # Returns
    /// 月亮的贝塞尔坐标[x,y,z]
    pub fn bse_m(&self, jd: f64) -> [f64; 3] {
        let mut a = self.cd2bse(
            &self.chazhi(jd, 1), 
            &self.chazhi(jd, 2)
        );
        // 转换单位到地球半径
        a[0] /= CS_R_EAR;
        a[1] /= CS_R_EAR; 
        a[2] /= CS_R_EAR;
        a
    }

    /// 计算地球上一点的速度(贝塞尔坐标表达)
    /// # Arguments
    /// * `x` - x坐标
    /// * `y` - y坐标  
    /// * `s` - 贝赤交角
    /// * `vx` - x方向速度
    /// * `vy` - y方向速度
    /// # Returns
    /// 速度参数
    pub fn vxy(&self, x: f64, y: f64, s: f64, vx: f64, vy: f64) -> VelocityParams {
        let mut h = 1.0 - x*x - y*y;
        // 越界置0,使速度场连续,有助于迭代单向收敛
        if h < 0.0 {
            h = 0.0;
        } else {
            h = h.sqrt();
        }

        let vx_local = PI2 * (s.sin()*h - s.cos()*y);
        let vy_local = PI2 * x * s.cos();
        
        let vx_rel = vx - vx_local;
        let vy_rel = vy - vy_local;
        
        VelocityParams {
            vx: vx_local,
            vy: vy_local, 
            vx_rel,
            vy_rel,
            v: (vx_rel*vx_rel + vy_rel*vy_rel).sqrt()
        }
    }

    /// 计算半影本影半径和食分
    /// # Arguments
    /// * `mr` - 月球距离(千米) 
    /// # Returns
    /// 半影本影参数
    pub fn rsm(&self, mr: f64) -> ShadowParams {
        ShadowParams {
            r1: CS_K + self.tanf1 * mr,     // 半影半径
            r2: CS_K2 - self.tanf2 * mr,    // 本影半径
            ar2: (CS_K2 - self.tanf2 * mr).abs(),
            sf: CS_K2 / mr / CS_K0 * (self.dyj + mr)  // 食分
        }
    }

    /// 计算切入点
    /// # Arguments
    /// * `jd` - 儒略日
    /// * `dx` - x方向速度
    /// * `dy` - y方向速度
    /// * `fs` - 标志(1表示半影,0表示本影)
    /// # Returns
    /// [经度,纬度,时间]
    pub fn qrd(&self, mut jd: f64, dx: f64, dy: f64, fs: i32) -> [f64; 3] {
        let ba2 = self.bba * self.bba;
        let m = self.bse_m(jd);
        let mut x = m[0];
        let mut y = m[1];
        
        let b = self.rsm(m[2]);
        let r = if fs == 1 { b.r1 } else { 0.0 };
        
        // 计算椭圆面上的投影距离
        let d = 1.0 - (1.0/ba2 - 1.0)*y*y/(x*x + y*y)/2.0 + r;
        let t = (d*d - x*x - y*y)/(dx*x + dy*y)/2.0;
        
        // 更新位置
        x += t*dx;
        y += t*dy;
        jd += t;

        // 椭圆修正
        let c = (1.0 - ba2)*r*x*y/(d*d*d);
        x += c*y;
        y -= c*x;

        // 转换为地标系
        let mut re = self.bse2db(
            &[x/d, y/d, 0.0],
            &self.bse(jd),
            true
        );
        
        // 添加时间
        re[2] = jd;
        re
    }

    /// 计算日食基本特征
    pub fn feature(&self, jd: f64) -> EclipseFeature {
        // 使用低精度朔(误差10分钟)
        // jd = self.zjd;

        let tg = 0.04;
        
        // 计算速度等参数
        let a = self.bse_m(jd - tg);
        let b = self.bse_m(jd);
        let c = self.bse_m(jd + tg);
        
        let vx = (c[0] - a[0])/(tg * 2.0); // 影速x 
        let vy = (c[1] - a[1])/(tg * 2.0); // 影速y
        let vz = (c[2] - a[2])/(tg * 2.0);
        let ax = (c[0] + a[0] - 2.0*b[0])/(tg * tg);
        let ay = (c[1] + a[1] - 2.0*b[1])/(tg * tg);
        let v = (vx*vx + vy*vy).sqrt();
        let v2 = v*v;

        // 计算中点参数  
        let t0 = -(b[0]*vx + b[1]*vy)/v2;
        
        let mut ef = EclipseFeature {
            jd_suo: jd,                     // 朔
            jd: jd + t0,                    // 中点时间 
            d_t: self.dt,                   // deltaT
            ds: self.bhc,                   // 黄交线与赤交线夹角
            vx,                             // 影速x  
            vy,                             // 影速y
            ax,                             // x方向加速度
            ay,                             // y方向加速度 
            v,                              // 总速度
            k: vy/vx,                       // 斜率
            xc: b[0] + vx*t0,              // 中点x坐标 
            yc: b[1] + vy*t0,              // 中点y坐标
            zc: b[2] + vz*t0 - 1.37*t0*t0, // 中点z坐标
            d: ((vx*b[1] - vy*b[0])/v).abs(), // 到圆心距离
            i: self.bse(jd + t0),          // 贝塞尔坐标参数
            ..Default::default()
        };

        // 计算影轴交点
        let f = line_ear2(
            ef.xc, ef.yc, 2.0,
            ef.xc, ef.yc, 0.0, 
            CS_BA, 1.0, &ef.i
        );

        // 计算四个关键点的影子半径
        let bc = self.rsm(ef.zc);
        let mut bp = bc.clone();
        let mut b2 = bc.clone();
        let mut b3 = bc.clone();

        if f.w != 100.0 {
            bp = self.rsm(ef.zc - f.r2);
        }

        // 计算中心线参数
        if ef.d < 1.0 {
            let dt = (1.0 - ef.d * ef.d).sqrt()/v;
            let t2 = t0 - dt;
            let t3 = t0 + dt;
            b2 = self.rsm(t2*vz + b[2] - 1.37*t2*t2);
            b3 = self.rsm(t3*vz + b[2] - 1.37*t3*t3);
            
            // 中心始终
            ef.gk1 = self.qrd(t2 + jd, vx, vy, 0);
            ef.gk2 = self.qrd(t3 + jd, vx, vy, 0);
        }

        // 偏食始终参数
        let mut dt = 0.0;
        // if ef.d < 1.0 {
        //     dt = (1.0 - ef.d*ef.d).sqrt()/v;
        // }
        // let t2 = t0 - dt;
        // let t3 = t0 + dt;

        let ls = 1.0 + bc.r1;
        // dt = 0.0;
        if ef.d < ls {
            dt = (ls*ls - ef.d*ef.d).sqrt()/v;
        }
        let t4 = t0 - dt;
        let t5 = t0 + dt;
        let t6 = -b[0]/vx;

        // 计算关键时刻的位置
        ef.gk3 = self.qrd(t4 + jd, vx, vy, 1);
        ef.gk4 = self.qrd(t5 + jd, vx, vy, 1);
        ef.gk5 = self.bse_xy2db(t6*vx + b[0], t6*vy + b[1], &self.bse(t6 + jd), true);
        ef.gk5[2] = t6 + jd;

        // 根据有无中心线分类处理
        if f.w == 100.0 {
            // 无中心线情况
            let ls = self.bse2db(&[ef.xc, ef.yc, 0.0], &ef.i, false);
            ef.zx_j = ls[0];
            ef.zx_w = ls[1];
            ef.sf = (bc.r1 - (ef.d - 0.9972))/(bc.r1 - bc.r2);

            // 类型判断
            ef.lx = if ef.d > 0.9972 + bc.r1 {
                "N".to_string()  // 无食
            } else if ef.d > 0.9972 + bc.ar2 {
                "P".to_string()  // 偏食
            } else if bc.sf < 1.0 {
                "A0".to_string() // 环食
            } else {
                "T0".to_string() // 全食  
            };
        } else {
            // 有中心线情况 
            ef.zx_j = f.j;
            ef.zx_w = f.w;
            ef.sf = bp.sf;

            // 类型判断
            ef.lx = if ef.d > 0.9966 - bp.ar2 {
                if bp.sf < 1.0 {
                    "A1".to_string()
                } else { 
                    "T1".to_string()
                }
            } else if bp.sf >= 1.0 {
                if b2.sf > 1.0 && b3.sf > 1.0 {
                    "T".to_string()
                } else if b3.sf > 1.0 {
                    "H3".to_string()
                } else if b2.sf > 1.0 {
                    "H2".to_string() 
                } else {
                    "H".to_string()
                }
            } else {
                "A".to_string()
            };
        }

        // 计算日食带和时延
        if f.w != 100.0 {
            ef.dw = (2.0 * bp.r2 * CS_R_EAR).abs() / ef.sdp[1].sin();
            let ls = self.vxy(ef.xc, ef.yc, ef.i[1], ef.vx, ef.vy);
            ef.tt = 2.0 * bp.r2.abs() / ls.v;
        }

        ef
    }


    /// 添加点到界线数组中
    /// # Arguments 
    /// * `z` - 点的坐标[经度,纬度]
    /// * `p` - 存储界线的数组
    fn push_point(&self, z: &[f64], p: &mut Vec<f64>) {
        p.push(z[0]);
        p.push(z[1]);
    }

    /// 数据元素复制
    /// # Arguments
    /// * `a` - 目标数组
    /// * `n` - 目标位置(-2表示末尾,-1表示倒数第二位)
    /// * `b` - 源数组 
    /// * `m` - 源位置(-2表示末尾,-1表示倒数第二位)
    fn elm_copy(&self, a: &mut Vec<f64>, mut n: i32, b: &[f64], mut m: i32) {
        if b.is_empty() {
            return;
        }

        // 处理特殊位置参数
        if n == -2 {
            n = a.len() as i32;
        }
        if m == -2 {
            m = b.len() as i32;
        }
        if n == -1 {
            n = (a.len() - 2) as i32;
        }
        if m == -1 {
            m = (b.len() - 2) as i32;
        }

        // 复制数据
        a[n as usize] = b[m as usize];
        a[n as usize + 1] = b[m as usize + 1];
    }


    /// 计算南北界点
    /// # Arguments
    /// * `m` - 月球贝塞尔坐标[x,y,z]
    /// * `vx0` - 影足x向速度
    /// * `vy0` - 影足y向速度
    /// * `h` - 1计算北界,-1计算南界
    /// * `r` - 半影半径
    /// * `i` - 贝塞尔坐标参数
    /// # Returns
    /// [经度,纬度,x,y]
    fn nan_bei(&self, m: &[f64], vx0: f64, vy0: f64, h: f64, r: f64, i: &[f64]) -> Vec<f64> {
        let mut x = m[0] - vy0/vx0 * r * h;
        let mut y = m[1] + h * r;
        let mut js = 0;

        let mut sin_a = 0.0;
        let mut cos_a = 0.0;
        // 迭代计算交点
        for _ in 0..3 {
            let mut z = 1.0 - x*x - y*y;
            if z < 0.0 {
                if js > 0 {
                    break;
                }
                z = 0.0;
                js += 1;
            }
            z = z.sqrt();

            // 修正坐标
            x -= (x - m[0]) * z / m[2];
            y -= (y - m[1]) * z / m[2];

            // 计算速度
            let vx = vx0 - PI2 * (i[1].sin() * z - i[1].cos() * y);
            let vy = vy0 - PI2 * i[1].cos() * x;
            let v = (vx*vx + vy*vy).sqrt();

            // 计算正余弦
            sin_a = h * vy / v;
            cos_a = h * vx / v;

            // 更新坐标
            x = m[0] - r * sin_a;
            y = m[1] + r * cos_a;
        }

        // 计算交点地标
        let x2 = m[0] - CS_K * sin_a;
        let y2 = m[1] + CS_K * cos_a;
        let p = line_ear2(x2, y2, m[2], x, y, 0.0, CS_BA, 1.0, i);

        vec![p.j, p.w, x, y]
    }

    /// 界线切割计算
    /// # Arguments
    /// * `m` - 月球贝塞尔坐标
    /// * `vx0` - 影足x向速度
    /// * `vy0` - 影足y向速度  
    /// * `h` - 1计算北界,-1计算南界
    /// * `r` - 半影半径
    /// * `i` - 贝塞尔坐标参数
    /// * `a` - 切割线数据
    pub fn m_qie(&self, m: &[f64], vx0: f64, vy0: f64, h: f64, r: f64, i: &[f64], 
                 a: &mut BoundaryData) {
        // 计算界线点
        let p = self.nan_bei(m, vx0, vy0, h, r, i);

        // 初始化标记 
        if a.f2.is_none() {
            a.f2 = Some(0);
        }
        a.f = if p[1] == 100.0 { 0 } else { 1 };

        // 处理解的变化
        if a.f2.unwrap() != a.f {
            // 计算椭圆与直线交点
            let g = line_ovl(p[2], p[3], vx0, vy0, 1.0, self.bba);
            
            if g.n > 0 {
                // 确定交点和距离
                let (dj, f) = if a.f == 1 {
                    (g.r2, g.b.clone())
                } else {
                    (g.r1, g.a.clone()) 
                };

                // 计算新的贝塞尔参数
                // let mut f = f;
                // f.push(0.0);
                let i2 = vec![
                    i[0],
                    i[1], 
                    i[2] - dj/(vx0*vx0 + vy0*vy0).sqrt() * 6.28
                ];

                // 添加线头
                self.push_point(&self.bse2db(&[f[0],f[1],0.0], &i2, true), &mut a.points);
            }
        }

        // 记录本次解
        a.f2 = Some(a.f);

        // 添加界线点
        if p[1] != 100.0 {
            self.push_point(&p, &mut a.points);
        }
    }


    /// 计算日出日没食甚
    /// # Arguments
    /// * `m` - 月球贝塞尔坐标[x,y,z]
    /// * `vx0` - 影足x向速度
    /// * `vy0` - 影足y向速度
    /// * `ab` - true表示选择第一个交点,false表示选择第二个交点  
    /// * `r` - 半影半径
    /// * `i` - 贝塞尔坐标参数
    /// * `a` - 存储食甚线点的数组
    /// # Returns
    /// * true - 找到有效交点
    /// * false - 无交点或交点无效
    pub fn m_dian(&self, m: &[f64], vx0: f64, vy0: f64, ab: bool, r: f64, 
            i: &[f64], a: &mut Vec<f64>) -> bool {
        // 从月球贝塞尔坐标开始迭代
        let mut curr_a = vec![m[0], m[1]];
        let mut curr_r: f64;

        // 迭代求交点(最多2次)
        for _ in 0..2 {
        // 计算地表速度场
        let c = self.vxy(curr_a[0], curr_a[1], i[1], vx0, vy0);
        
        // 求椭圆与速度场垂线的交点
        let p = line_ovl(m[0], m[1], c.vy, -c.vx, 1.0, self.bba);
        
        // 无交点则退出
        if p.n == 0 {
            return false;
        }

        // 根据ab标志选择交点
        if ab {
            curr_a = p.a.to_vec();
            curr_r = p.r1;
        } else {
            curr_a = p.b.to_vec(); 
            curr_r = p.r2;
        }

        // 如果交点在半影范围内且有效
        if curr_r <= r {
            // 转换为地标坐标
            curr_a.push(0.0);
            let db = self.bse2db(&curr_a, i, true);
            
            // 保存到食甚线数组
            self.push_point(&db, a);
            
            return true;
        }
        }

        false
    }


    /// 计算日出日没的初亏食甚复圆线，南北界线等
    /// # Arguments
    /// * `jd` - 儒略日数
    /// # Returns
    /// 日食界线参数
    pub fn jie_x(&self, jd: f64) -> EclipseLines {
        // 计算特征参数
        let re = self.feature(jd);

        // 初始化界线数组
        let mut el = EclipseLines::default();
        
        // 计算时间参数
        let t = if 1.7*1.7 - re.d*re.d < 0.0 {
            0.0 
        } else {
            (1.7*1.7 - re.d*re.d).sqrt()/re.v + 0.01
        };
        
        let t_start = re.jd - t;
        let n = 400;
        let dt = 2.0 * t / n as f64;

        // 切入计数器
        let mut n1 = 0;
        let mut n4 = 0;

        // 日出日没食甚线预置点
        let mut ua = &mut el.q1;
        let mut ub = &mut el.q2;
        self.push_point(&[0.0, 0.0], &mut ub);
        self.push_point(&[0.0, 0.0], &mut el.q3);
        self.push_point(&[0.0, 0.0], &mut el.q4);

        // 主循环计算
        let mut t = t_start;
        for _ in 0..=n {
            // 计算速度
            let vx = re.vx + re.ax*(t - re.jd_suo);
            let vy = re.vy + re.ay*(t - re.jd_suo);
            
            // 计算月亮位置
            let m = self.bse_m(t);
            let b = self.rsm(m[2]);
            let r = b.r1;
            let i = self.bse(t);

            // 计算椭圆与圆交点
            let p = cir_ovl(1.0, self.bba, r, m[0], m[1]);
            if n1 % 2 == 1 {
                if p.n == 0 { n1 += 1; }
            } else if p.n > 0 { 
                n1 += 1;
            }

            // 处理交点
            if p.n > 0 {
                // 转换交点坐标
                let mut pa = [p.a[0], p.a[1], 0.0];
                let mut pb = [p.b[0], p.b[1], 0.0];
                pa = self.bse2db(&pa, &i, true);
                pb = self.bse2db(&pb, &i, true);

                // 保存交点
                match n1 {
                    1 => {
                        self.push_point(&pa, &mut el.p1);
                        self.push_point(&pb, &mut el.p2);
                    },
                    3 => {
                        self.push_point(&pa, &mut el.p3);
                        self.push_point(&pb, &mut el.p4);
                    },
                    _ => {}
                }
            }

            // 处理日出日没食甚线
            if !self.m_dian(&m, vx, vy, false, r, &i, ua) {
                if ua.len() > 0 {
                    ua = &mut el.q3;
                }
            }
            if !self.m_dian(&m, vx, vy, true, r, &i, ub) {
                if ub.len() > 2 {
                    ub = &mut el.q4;  
                }
            }
            if t > re.jd {
                if ua.len() == 0 { ua = &mut el.q3; }
                if ub.len() == 2 { ub = &mut el.q4; }
            }

            // 计算中心线
            let p = self.bse_xy2db(m[0], m[1], &i, true);
            if (p[1] != 100.0 && n4 == 0) || (p[1] == 100.0 && n4 == 1) {
                // 计算椭圆与直线交点
                let ls = line_ovl(m[0], m[1], vx, vy, 1.0, self.bba);
                
                let (dj, mut ls_p) = if n4 == 0 {
                    (ls.r2, ls.b.to_vec())
                } else {
                    (ls.r1, ls.a.to_vec())
                };

                ls_p.push(0.0);

                // 计算新的贝塞尔参数
                let i2 = vec![
                    i[0], 
                    i[1],
                    i[2] - dj/(vx*vx + vy*vy).sqrt() * 6.28
                ];

                self.push_point(&self.bse2db(&ls_p, &i2, true), &mut el.l0);
                n4 += 1;
            }

            if p[1] != 100.0 {
                self.push_point(&p, &mut el.l0);
            }

            // 计算南北界
            self.m_qie(&m, vx, vy,  1.0, r,           &i, &mut el.l1); // 半影北界
            self.m_qie(&m, vx, vy, -1.0, r,           &i, &mut el.l2); // 半影南界  
            self.m_qie(&m, vx, vy,  1.0, b.r2,        &i, &mut el.l3); // 本影北界
            self.m_qie(&m, vx, vy, -1.0, b.r2,        &i, &mut el.l4); // 本影南界
            self.m_qie(&m, vx, vy,  1.0, (r+b.r2)/2.0,&i, &mut el.l5); // 0.5半影北界
            self.m_qie(&m, vx, vy, -1.0, (r+b.r2)/2.0,&i, &mut el.l6); // 0.5半影南界

            t += dt;
        }

        // 连接日出日没食甚线
        self.elm_copy(&mut el.q3, 0,    &el.q1, -1); // 连接q1和q3
        self.elm_copy(&mut el.q4, 0,    &el.q2, -1); // 连接q2和q4
        self.elm_copy(&mut el.q1, -2,   &el.l1.points, 0);  // 半影北界线西端
        self.elm_copy(&mut el.q2, -2,   &el.l2.points, 0);  // 半影南界线西端
        self.elm_copy(&mut el.q3, 0,    &el.l1.points, -1); // 半影北界线东端
        self.elm_copy(&mut el.q4, 0,    &el.l2.points, -1); // 半影南界线东端
        self.elm_copy(&mut el.q2, 0,    &el.q1, 0);
        self.elm_copy(&mut el.q3, -2,   &el.q4, -1);

        el
    }


    /// 计算日食界线(力学时)
    /// # Arguments
    /// * `jd` - 儒略日数(力学时)
    /// # Returns
    /// 日食界线数据
    pub fn jie_x2(&self, jd: f64) -> EclipseLines2 {
        let mut re = EclipseLines2::default();

        // 如果时间差太大,直接返回
        if (jd - self.zjd).abs() > 0.5 {
            return re;
        }

        // 获取日月坐标
        let ss = self.sun(jd);    // 太阳赤道坐标
        let m = self.bse_m(jd);  // 月亮贝塞尔坐标
        let b = self.rsm(m[2]);  // 本半影等
        let ii = self.bse(jd);    // 贝塞尔坐标参数
        let z = m[2];            // 月球z坐标

        // 计算辅助变量
        let a0 = m[0] * m[0] + m[1] * m[1];
        let a1 = a0 - b.r2 * b.r2;
        let a2 = a0 - b.r1 * b.r1;

        // 环形计算
        let n = 200;
        for i in 0..n { //第0和第N点是同一点，可形成一个环，但不必计算，因为第0点可能在界外而无效
            // 计算角度参数
            let s = i as f64 / n as f64 * PI2;
            let (sin_s, cos_s) = s.sin_cos();

            // 计算日月位置
            let x_sun = m[0] + CS_K * cos_s; 
            let y_sun = m[1] + CS_K * sin_s;

            // 计算本影界线点
            let (x, y) = (
                m[0] + b.r2 * cos_s,
                m[1] + b.r2 * sin_s
            );
            let mut p = line_ear2(
                x_sun, y_sun, z,
                x, y, 0.0,
                CS_BA, 1.0, &ii
            );

            // 保存本影界线点
            if p.w != 100.0 {
                self.push_point(&[p.j, p.w], &mut re.p1);
            } else if (x*x + y*y).sqrt() > a1 {
                let db = self.bse2db(&[x, y, 0.0], &ii, true);
                self.push_point(&db, &mut re.p1);
            }

            // 计算半影界线点  
            let (x, y) = (
                m[0] + b.r1 * cos_s,
                m[1] + b.r1 * sin_s  
            );
            p = line_ear2(
                x_sun, y_sun, z,
                x, y, 0.0,
                CS_BA, 1.0, &ii
            );

            // 保存半影界线点
            if p.w != 100.0 {
                self.push_point(&[p.j, p.w], &mut re.p2);
            } else if (x*x + y*y).sqrt() > a2 {
                let db = self.bse2db(&[x, y, 0.0], &ii, true);
                self.push_point(&db, &mut re.p2);
            }

            // 计算晨昏圈点
            let mut p = llr_conv(&[s, 0.0, 0.0], PI2 - ss[1]);
            p[0] = rad2rrad(p[0] + ss[0] + PI2 - ii[2]);
            self.push_point(&p, &mut re.p3);
        }

        // 闭合界线
        let l1 = re.p1.len();
        let l2 = re.p2.len();
        let l3 = re.p3.len();
        if l1 >= 2 {
            re.p1.push(re.p1[0]);
            re.p1.push(re.p1[1]);
        }
        if l2 >= 2 {
            re.p2.push(re.p2[0]);
            re.p2.push(re.p2[1]);
        }
        if l3 >= 2 {
            re.p3.push(re.p3[0]);
            re.p3.push(re.p3[1]);
        }

        re
    }


    /// 计算日食界线数据表
    /// # Arguments
    /// * `jd` - 儒略日数
    /// # Returns 
    /// 界线数据数组,每个元素包含[时间,北界点,本影北界,中心线,本影南界,南界点]
    pub fn jie_x3(&self, jd: f64) -> Vec<EclipseLineData> {
        // 计算特征参数
        let re = self.feature(jd);
        
        // 设置时间参数
        let mut t = (re.jd * 1440.0).floor() / 1440.0 - 3.0/24.0;
        let n = 360;
        let dt = 1.0/1440.0;
        
        let mut lines = Vec::new();

        // 主循环计算
        for _ in 0..n {
            // 计算速度
            let vx = re.vx + re.ax * (t - re.jd_suo);
            let vy = re.vy + re.ay * (t - re.jd_suo);
            
            // 计算月亮位置
            let m = self.bse_m(t);     // 月球贝塞尔坐标
            let b = self.rsm(m[2]);    // 本半影等参数  
            let r = b.r1;              // 半影半径
            let i = self.bse(t);       // 贝塞尔坐标参数

            // 存储当前时刻的界线数据
            let mut line = EclipseLineData {
                time: t + J2000,
                points: Vec::with_capacity(5)
            };

            // 计算五条界线
            let points = [
                // 半影北界
                self.nan_bei(&m, vx, vy, 1.0, r, &i),
                // 本影北界  
                self.nan_bei(&m, vx, vy, 1.0, b.r2, &i),
                // 中心线
                self.bse_xy2db(m[0], m[1], &i, true),
                // 本影南界
                self.nan_bei(&m, vx, vy, -1.0, b.r2, &i),
                // 半影南界
                self.nan_bei(&m, vx, vy, -1.0, r, &i),
            ];

            // 检查点是否有效并保存
            let mut valid = false;
            for p in points {
                if p[1] != 100.0 {
                    line.points.push([p[0], p[1]]);
                    valid = true;
                } else {
                    line.points.push([0.0, 0.0]); // 无效点用0替代
                }
            }

            // 如果有有效点则保存该时刻数据
            if valid {
                lines.push(line);
            }

            t += dt;
        }

        lines
    }
}

/// 时刻日月坐标结构
#[derive(Default, Clone)]
struct TimeCoord {
    /// 太阳坐标
    s: Vec<f64>,
    /// 月亮坐标 
    m: Vec<f64>,
    /// 恒星时
    g: f64,
    ///视半径
    mr: f64,
    sr: f64,
    ///日面中心直角坐标(用于日食)
    x: f64,
    y: f64,
}

/// 日食批量快速计算器
pub struct RsPL {
    /// 是否采用NASA的视径比(1表示采用) 
    nasa_r: i32,
    /// 地方日食时间表
    s_t: Vec<f64>,
    a: [f64;3],  // 本影锥顶点坐标
    b: [f64;3],  // 半影锥顶点坐标
    p: TimeCoord, // t1时刻参数
    q: TimeCoord, // t2时刻参数
    v: Vec<f64>,  // 食界表
    vc: String, // 食中心类型
    vb: String, // 本影南北距离

    eclipse_data: EclipseData,
}

/// 日月坐标计算结果
#[derive(Default)]
pub struct SecXYResult {
    // 月亮参数
    pub m_cj: f64,   // 月亮视赤经
    pub m_cw: f64,   // 月亮赤纬
    pub m_r: f64,    // 地月距离
    pub m_cj2: f64,  // 修正视差后赤经 
    pub m_cw2: f64,  // 修正视差后赤纬
    pub m_r2: f64,   // 修正视差后距离
    
    // 太阳参数
    pub s_cj: f64,   // 太阳视赤经
    pub s_cw: f64,   // 太阳赤纬  
    pub s_r: f64,    // 日地距离
    pub s_cj2: f64,  // 修正视差后赤经
    pub s_cw2: f64,  // 修正视差后赤纬
    pub s_r2: f64,   // 修正视差后距离

    // 日食参数
    pub mr: f64,     // 月亮视半径
    pub sr: f64,     // 太阳视半径
    pub x: f64,      // 日面中心x坐标 
    pub y: f64,      // 日面中心y坐标
    pub t: f64,      // 计算时间
}

/// 日食数据结构
#[derive(Default)]
struct EclipseData {
    lx: String,    // 食类型
    sf: f64,       // 食分
    sf2: f64,      // 日出食分
    sf3: f64,      // 日没食分
    sflx: String,  // 食分类型
    b1: f64,       // 月日半径比
    dur: f64,      // 持续时间
    p1: f64,       // 初亏北向方位角
    v1: f64,       // 初亏顶向方位角
    p2: f64,       // 复圆北向方位角
    v2: f64,       // 复圆顶向方位角
    sun_s: f64,    // 日出时刻
    sun_j: f64,    // 日没时刻
}

/// 食界点坐标结构体
#[derive(Debug, Default, Clone)]
pub struct BoundaryPoint {
    /// 地理经度（弧度）
    pub j: f64,
    
    /// 地理纬度（弧度）
    pub w: f64,
    
    /// 食相类型（"全"或"环"）
    pub c: String,
    
    /// 标记是否有解
    pub f: i32,
    
    /// 上一次是否有解
    pub f2: i32,
}


//日食批量快速计算器
impl RsPL {
     /// 创建日食计算器实例
     pub fn new() -> Self {
        Self {
            nasa_r: 0,
            s_t: Vec::new(),
            a: [0.0,0.0,0.0],  // 本影锥顶点坐标
            b: [0.0,0.0,0.0],  // 半影锥顶点坐标
            p: TimeCoord::default(), // t1时刻参数
            q: TimeCoord::default(), // t2时刻参数
            v: Vec::new(),  // 食界表
            vc: String::new(), // 食中心类型
            vb: String::new(), // 本影南北距离
            eclipse_data: EclipseData::default(),
        }
    }

    /// 日月xy坐标计算
    /// # Arguments
    /// * `jd` - 力学时
    /// * `l` - 地理经度
    /// * `fa` - 地理纬度
    /// * `high` - 海拔(千米)
    /// * `re` - 计算结果存储
    pub fn sec_xy(&self, jd: f64, l: f64, fa: f64, high: f64, re: &mut SecXYResult) {
        // 基本参数计算
        let deltat = dt_t(jd); // TD-UT
        let zd = Nutation::nutation2(jd/36525.0);
        
        // 真恒星时(不考虑非多项式部分)
        let gst = p_gst(jd - deltat, deltat) + 
                 zd[0] * (Prece::hcjj(jd/36525.0) + zd[1]).cos();

        let rs_gs = RsGS::new();
        // =======月亮========
        let mut z = rs_gs.moon(jd);
        re.m_cj = z[0];
        re.m_cw = z[1]; 
        re.m_r = z[2]; // 月亮视赤经,月球赤纬
        
        // 得到此刻月亮时角
        let m_shi_j = rad2rrad(gst + l - z[0]); 
        
        // 修正了视差的赤道坐标
        parallax(&mut z, m_shi_j, fa, high);
        re.m_cj2 = z[0];
        re.m_cw2 = z[1];
        re.m_r2 = z[2];

        // =======太阳========
        z = rs_gs.sun(jd);
        re.s_cj = z[0];
        re.s_cw = z[1];
        re.s_r = z[2]; // 太阳视赤经,太阳赤纬

        // 得到此刻太阳时角
        let s_shi_j = rad2rrad(gst + l - z[0]);
        
        // 修正了视差的赤道坐标
        parallax(&mut z, s_shi_j, fa, high);
        re.s_cj2 = z[0];
        re.s_cw2 = z[1];
        re.s_r2 = z[2];

        // =======视半径========
        re.mr = CS_S_MOON/(re.m_r2 * RAD);
        re.sr = 959.63/(re.s_r2 * RAD) * CS_AU;
        
        if self.nasa_r == 1 {
            re.mr *= CS_S_MOON2/CS_S_MOON; // 0.99925
        }

        // =======日月赤经纬差转为日面中心直角坐标(用于日食)==============
        re.x = rad2rrad(re.m_cj2 - re.s_cj2) * 
               ((re.m_cw2 + re.s_cw2)/2.0).cos();
        re.y = re.m_cw2 - re.s_cw2;
        re.t = jd;
    }



    /// 日食的食甚计算
    /// # Arguments
    /// * `jd` - 近朔的力学时(误差几天不要紧)
    /// * `l` - 地理经度
    /// * `fa` - 地理纬度
    /// * `high` - 海拔高度(千米)
    pub fn sec_max(&mut self, mut jd: f64, l: f64, fa: f64, high: f64) {
        // 初始化日食数据
        for i in 0..5 {
            self.s_t[i] = 0.0; // 分别是:食甚,初亏,复圆,食既,生光
        }
        
        // 清空基本数据
        let mut eclipse_data = EclipseData {
            lx: String::new(), // 类型
            sf: 0.0,  // 食分
            sf2: 0.0, // 日出食分
            sf3: 0.0, // 日没食分
            sflx: " ".to_string(), // 食分类型
            b1: 1.0,  // 月日半径比(食甚时刻)
            dur: 0.0, // 持续时间
            p1: 0.0,  // 初亏方位(P北点起算)
            v1: 0.0,  // 初亏方位(V顶点起算)
            p2: 0.0,  // 复圆方位(P北点起算)
            v2: 0.0,  // 复圆方位(V顶点起算)
            sun_s: 0.0, // 日出
            sun_j: 0.0, // 日没
        };

        // 初始化插值表
        let mut rsgs = RsGS::new();
        rsgs.init(jd, 7);
        jd = rsgs.zjd; // 食甚初始估值为插值表中心时刻(粗朔)

        // 计算食甚时刻
        let mut g = SecXYResult::default();
        let mut g2 = SecXYResult::default();
        
        self.sec_xy(jd, l, fa, high, &mut g);
        jd -= g.x/0.2128; // 与食甚的误差在20分钟以内

        // 精确迭代食甚时刻
        let dt = 60.0/86400.0;
        let mut v = 0.0;
        let mut u = 0.0;
        for _ in 0..2 {
            self.sec_xy(jd, l, fa, high, &mut g);
            self.sec_xy(jd+dt, l, fa, high, &mut g2);
            // if let Err(_) = self.sec_xy(jd, l, fa, high, &mut g) {
            //     return;
            // }
            // if let Err(_) = self.sec_xy(jd+dt, l, fa, high, &mut g2) {
            //     return;
            // }

            u = (g2.y - g.y)/dt;
            v = (g2.x - g.x)/dt; 
            let dt2 = -(g.y*u + g.x*v)/(u*u + v*v);
            jd += dt2; // 极值时间
        }

        // 求直线到太阳中心的最小值
        let mut max_sf = 0.0;
        let mut max_jd = jd;
        // let r_min;
        
        // 粗算最小值
        for i in (-30..30).step_by(6) {
            let tt = jd + i as f64/86400.0;
            let _ = self.sec_xy(tt, l, fa, high, &mut g2);
            let ls = (g2.mr + g2.sr - (g2.x*g2.x + g2.y*g2.y).sqrt()) / g2.sr / 2.0;
            if ls > max_sf {
                max_sf = ls;
                max_jd = tt;
            }
        }

        // 精算最小值
        jd = max_jd;
        for i in -5..5 {
            let tt = jd + i as f64/86400.0;
            let _ = self.sec_xy(tt, l, fa, high, &mut g2);
            let ls = (g2.mr + g2.sr - (g2.x*g2.x + g2.y*g2.y).sqrt()) / g2.sr / 2.0;
            if ls > max_sf {
                max_sf = ls;
                max_jd = tt; 
            }
        }

        // 最终食甚时刻
        jd = max_jd;
        let _ = self.sec_xy(jd, l, fa, high, &mut g);
        let r_min = (g.x*g.x + g.y*g.y).sqrt();

        // 计算日出日没时刻
        let dt = dt_t(jd);
        eclipse_data.sun_s = sun_sheng_j(jd-dt+l/PI2, l, fa, -1) + dt;
        eclipse_data.sun_j = sun_sheng_j(jd-dt+l/PI2, l, fa,  1) + dt;

        let ggg = PosVel::new(g.x,g.y,g.t);
        let ggg2 = PosVel::new(g2.x,g2.y,g2.t);

        // 根据食分计算类型
        if r_min <= g.mr + g.sr {
            // 偏食计算
            self.s_t[1] = jd;
            eclipse_data.lx = String::from("偏");
            eclipse_data.sf = (g.mr + g.sr - r_min)/g.sr/2.0;
            eclipse_data.b1 = g.mr/g.sr;

            // 计算日出食分
            let _ = self.sec_xy(eclipse_data.sun_s, l, fa, high, &mut g2);
            eclipse_data.sf2 = (g2.mr + g2.sr - 
                (g2.x*g2.x + g2.y*g2.y).sqrt())/g2.sr/2.0;
            if eclipse_data.sf2 < 0.0 {
                eclipse_data.sf2 = 0.0;
            }

            // 计算日没食分
            let _ = self.sec_xy(eclipse_data.sun_j, l, fa, high, &mut g2);
            eclipse_data.sf3 = (g2.mr + g2.sr - 
                (g2.x*g2.x + g2.y*g2.y).sqrt())/g2.sr/2.0;
            if eclipse_data.sf3 < 0.0 {
                eclipse_data.sf3 = 0.0;
            }

            // 计算初亏
            self.s_t[0] = line_t(&ggg, v, u, g.mr+g.sr, false);
            for _ in 0..3 {
                let _ = self.sec_xy(self.s_t[0], l, fa, high, &mut g2);
                self.s_t[0] = line_t(&ggg2, v, u, g2.mr+g2.sr, false);
            }

            // 计算初亏方位角
            eclipse_data.p1 = rad2mrad((g2.x/g2.y).atan());
            eclipse_data.v1 = rad2mrad(eclipse_data.p1 - 
                shi_cha_j(p_gst2(self.s_t[0]), l, fa, g2.s_cj, g2.s_cw));

            // 计算复圆
            self.s_t[2] = line_t(&ggg, v, u, g.mr+g.sr, true);
            for _ in 0..3 {
                let _ = self.sec_xy(self.s_t[2], l, fa, high, &mut g2);
                self.s_t[2] = line_t(&ggg2, v, u, g2.mr+g2.sr, true);
            }

            eclipse_data.p2 = rad2mrad((g2.x/g2.y).atan());
            eclipse_data.v2 = rad2mrad(eclipse_data.p2 - 
                shi_cha_j(p_gst2(self.s_t[2]), l, fa, g2.s_cj, g2.s_cw));
        }

        // 全食计算
        if r_min <= g.mr - g.sr {
            eclipse_data.lx = String::from("全");
            self.s_t[3] = line_t(&ggg, v, u, g.mr-g.sr, false);
            let _ = self.sec_xy(self.s_t[3], l, fa, high, &mut g2);
            self.s_t[3] = line_t(&ggg2, v, u, g2.mr-g2.sr, false);

            self.s_t[4] = line_t(&ggg, v, u, g.mr-g.sr, true);
            let _ = self.sec_xy(self.s_t[4], l, fa, high, &mut g2);
            self.s_t[4] = line_t(&ggg2, v, u, g2.mr-g2.sr, true);
            eclipse_data.dur = self.s_t[4] - self.s_t[3];
        }

        // 环食计算
        if r_min <= g.sr - g.mr {
            eclipse_data.lx = String::from("环");
            self.s_t[3] = line_t(&ggg, v, u, g.sr-g.mr, false);
            let _ = self.sec_xy(self.s_t[3], l, fa, high, &mut g2);
            self.s_t[3] = line_t(&ggg2, v, u, g2.sr-g2.mr, false);

            self.s_t[4] = line_t(&ggg, v, u, g.sr-g.mr, true);
            let _ = self.sec_xy(self.s_t[4], l, fa, high, &mut g2);
            self.s_t[4] = line_t(&ggg2, v, u, g2.sr-g2.mr, true);
            eclipse_data.dur = self.s_t[4] - self.s_t[3];
        }

        // 处理日出日没情况
        if self.s_t[1] < eclipse_data.sun_s && eclipse_data.sf2 > 0.0 {
            eclipse_data.sf = eclipse_data.sf2;
            eclipse_data.sflx = "#".to_string();
        }
        if self.s_t[1] > eclipse_data.sun_j && eclipse_data.sf3 > 0.0 {
            eclipse_data.sf = eclipse_data.sf3;
            eclipse_data.sflx = "*".to_string();
        }

        // 无效时刻置0
        for i in 0..5 {
            if self.s_t[i] < eclipse_data.sun_s || 
               self.s_t[i] > eclipse_data.sun_j {
                self.s_t[i] = 0.0;
            }
        }

        // 转换时间
        eclipse_data.sun_s -= dt_t(jd);
        eclipse_data.sun_j -= dt_t(jd);

        // 保存结果
        self.eclipse_data = eclipse_data;
    }


    /// 计算本影和半影锥顶点坐标
    /// # Arguments
    /// * `jd` - 儒略日数
    fn zb0(&mut self, jd: f64) {
        // 基本参数计算
        let deltat = dt_t(jd); // TD-UT
        let e = Prece::hcjj(jd/36525.0);
        let zd = Nutation::nutation2(jd/36525.0);

        let mut rsgs = RsGS::new();
        rsgs.init(jd,7);

        // 计算t1时刻的参数
        let p = &mut self.p;
        p.g = p_gst(jd-deltat, deltat) + 
              zd[0] * (e + zd[1]).cos(); // 真恒星时(不考虑非多项式部分)
        p.s = rsgs.sun(jd);
        p.m = rsgs.moon(jd);

        // 计算t2时刻的参数
        let t2 = jd + 60.0/86400.0;
        let q = &mut self.q;
        q.g = p_gst(t2-deltat, deltat) + 
              zd[0] * (e + zd[1]).cos();
        q.s = rsgs.sun(t2);
        q.m =rsgs.moon(t2); 

        // 转为直角坐标
        let z1 = llr2xyz(&p.s);
        let z2 = llr2xyz(&p.m);

        // k为日月半径比
        let k = 959.63/CS_S_MOON * CS_AU;

        // 计算本影锥顶点坐标
        let f = [
            (z1[0] - z2[0])/(1.0 - k) + z2[0],
            (z1[1] - z2[1])/(1.0 - k) + z2[1], 
            (z1[2] - z2[2])/(1.0 - k) + z2[2]
        ];
        self.a = xyz2llr(&f);

        // 计算半影锥顶点坐标
        let f = [
            (z1[0] - z2[0])/(1.0 + k) + z2[0],
            (z1[1] - z2[1])/(1.0 + k) + z2[1],
            (z1[2] - z2[2])/(1.0 + k) + z2[2]  
        ];
        self.b = xyz2llr(&f);
    }

    /// 计算日月视差修正和视半径参数
    /// # Arguments
    /// * `p` - 时刻坐标参数
    /// * `l` - 地理经度
    /// * `fa` - 地理纬度
    fn zb_xy(&self, p: &TimeCoord, l: f64, fa: f64) -> TimeCoord {
        // 创建太阳月亮坐标数组
        let mut s = vec![p.s[0], p.s[1], p.s[2]];
        let mut m = vec![p.m[0], p.m[1], p.m[2]];

        // 视差修正
        parallax(&mut s, p.g + l - p.s[0], fa, 0.0);
        parallax(&mut m, p.g + l - p.m[0], fa, 0.0);

        // // 视半径计算
        // p.mr = CS_S_MOON/(m[2] * RAD);
        // p.sr = 959.63/(s[2] * RAD) * CS_AU;

        // // 日月赤经纬差转日面中心直角坐标
        // p.x = rad2rrad(m[0] - s[0]) * ((m[1] + s[1])/2.0).cos();
        // p.y = m[1] - s[1];

        // 创建结果
        let mut result = p.clone();
        
        // 视半径计算
        result.mr = CS_S_MOON/(m[2] * RAD);
        result.sr = 959.63/(s[2] * RAD) * CS_AU;

        // 日月赤经纬差转日面中心直角坐标
        result.x = rad2rrad(m[0] - s[0]) * ((m[1] + s[1])/2.0).cos();
        result.y = m[1] - s[1];
        
        result
    }

    /// 计算食界点坐标
    /// # Arguments  
    /// * `l` - 地理经度
    /// * `fa` - 地理纬度
    /// * `re` - 返回结果存储
    /// * `f_ab` - true选择A点,false选择B点
    /// * `f` - 取+-1
    fn p2p(&mut self, l: f64, fa: f64, re: &mut BoundaryPoint, f_ab: bool, f: f64) {
        // 计算两个时刻的日月坐标
        let p_mod = self.zb_xy(&self.p, l, fa);
        let q_mod = self.zb_xy(&self.q, l, fa); 

        // 计算速度和半径参数
        let u = q_mod.y - p_mod.y;
        let v = q_mod.x - p_mod.x;
        let a = (u*u + v*v).sqrt();
        let r = 959.63/p_mod.s[2]/RAD * CS_AU;

        // 计算交点坐标
        let w = p_mod.s[1] + f*r*v/a;
        let j = p_mod.s[0] - f*r*u/a/((w + p_mod.s[1])/2.0).cos();
        let r = p_mod.s[2];

        // 选择锥顶点并计算交点
        let a = if f_ab { &self.a } else { &self.b };
        let pp = line_ear(&[j,w,r], a, p_mod.g);
        
        // 保存结果
        re.j = pp.j;
        re.w = pp.w;
    }

    /// 计算食中心点 
    /// # Arguments
    /// * `re` - 返回结果存储
    fn pp0(&mut self, re: &mut BoundaryPoint) {
        // 计算日月连线与地球的交点
        let pp = line_ear(&self.p.m, &self.p.s, self.p.g);
        re.j = pp.j;
        re.w = pp.w;

        // 无解则返回
        if re.w == 100.0 {
            re.c = String::new();
            return;
        }

        // 确定食相类型
        re.c = String::from("全");
        let p_mod = self.zb_xy(&self.p, re.j, re.w);
        if p_mod.sr > p_mod.mr {
            re.c = String::from("环"); 
        }
    }

    /// 计算日食南北界
    /// # Arguments
    /// * `jd` - 儒略日数
    pub fn nbj(&mut self, jd: f64) {
        // 初始化插值表
        let mut rs_gs = RsGS::new();
        rs_gs.init(jd, 7);

        // 初始化食界数组和结果
        for i in 0..10 {
            self.v[i] = 100.0;
        }
        self.vc.clear();
        self.vb.clear();

        // 计算本影和半影锥顶点坐标
        self.zb0(jd);

        // 计算食中心点
        let mut g = BoundaryPoint::default();
        self.pp0(&mut g);
        self.v[0] = g.j;
        self.v[1] = g.w;
        self.vc = g.c.clone();

        // 计算本影北界(环食为南界)
        // 本影区之内,变差u,v基本不变,所以计算两次足够
        g.j = 0.0;
        g.w = 0.0;
        for _ in 0..2 {
            self.p2p(g.j, g.w, &mut g, true, 1.0);
        }
        self.v[2] = g.j;
        self.v[3] = g.w;

        // 计算本影南界(环食为北界)
        g.j = 0.0;
        g.w = 0.0;
        for _ in 0..2 {
            self.p2p(g.j, g.w, &mut g, true, -1.0);
        }
        self.v[4] = g.j;
        self.v[5] = g.w;

        // 计算半影北界
        g.j = 0.0;
        g.w = 0.0;
        for _ in 0..3 {
            self.p2p(g.j, g.w, &mut g, false, -1.0);
        }
        self.v[6] = g.j;
        self.v[7] = g.w;

        // 计算半影南界
        g.j = 0.0;
        g.w = 0.0;
        for _ in 0..3 {
            self.p2p(g.j, g.w, &mut g, false, 1.0);
        }
        self.v[8] = g.j;
        self.v[9] = g.w;

        // 粗算本影南北距离
        if self.v[3] != 100.0 && self.v[5] != 100.0 {
            let x = (self.v[2] - self.v[4]) * 
                   ((self.v[3] + self.v[5])/2.0).cos();
            let y = self.v[3] - self.v[5];
            self.vb = format!("{}千米", 
                (CS_R_EAR_A * (x*x + y*y).sqrt()).round());
        }
    }
}