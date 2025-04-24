use crate::eph::delta_t::dt_t;
use crate::eph::xl::XL;
use crate::lunar::lunar::YUE_MING;

/// 实朔实气计算器
/// 适用范围 -722年2月22日——1959年12月
/// 平气平朔计算使用古历参数进行计算
/// 定朔、定气计算使用开普勒椭圆轨道计算,同时考虑了光行差和力学时与UT1的时间差
#[derive(Debug, Default)]
pub struct SolarTerms {
    /// 朔修正表
    sb: String,
    /// 气修正表
    qb: String,
    // /// 朔直线拟合参数
    // suo_kb: [(f64, f64); 13],
    // /// 气直线拟合参数
    // qi_kb: [(f64, f64); 30],
    /// 闰月位置
    pub leap: i32,
    /// 各月名称
    pub month_names: Vec<String>,

    pub month_index: Vec<usize>,
    
    /// 中气表,其中.liqiu是节气立秋的儒略日,计算三伏时用到
    pub solar_terms: Vec<f64>,
    /// 合朔表
    pub new_moons: Vec<f64>,
    /// 各月大小
    pub month_days: Vec<i32>,
    /// 年计数
    pub year_counts: Vec<i32>,
}

impl SolarTerms {
    /// 朔直线拟合参数
    pub const SUO_KB: [(f64, f64); 13] = [
        // JD, 合朔周期(天)                       // 说明
        (1457698.231017, 29.53067166), // -721-12-17 h=0.00032 古历·春秋
        (1546082.512234, 29.53085106), // -479-12-11 h=0.00053 古历·战国
        (1640640.735300, 29.53060000), // -221-10-31 h=0.01010 古历·秦汉
        (1642472.151543, 29.53085439), // -216-11-04 h=0.00040 古历·秦汉
        (1683430.509300, 29.53086148), // -104-12-25 h=0.00313 汉书·律历志(太初历)平气平朔
        (1752148.041079, 29.53085097), //   85-02-13 h=0.00049 后汉书·律历志(四分历)
        (1807724.481520, 29.53059851), //  237-04-12 h=0.00033 晋书·律历志(景初历)
        (1883618.114100, 29.53060000), //  445-01-24 h=0.00030 宋书·律历志(何承天元嘉历)
        (1907360.704700, 29.53060000), //  510-01-26 h=0.00030 宋书·律历志(祖冲之大明历)
        (1936596.224900, 29.53060000), //  590-02-10 h=0.01010 随书·律历志(开皇历)
        (1939135.675300, 29.53060000), //  597-01-24 h=0.00890 随书·律历志(大业历)
        (1947168.000000, 29.53060000), //  619-01-21
        (2436935.000000, 29.53060000), // 边界值
    ];

    /// 气直线拟合参数
    pub const QI_KB: [(f64, f64); 36] = [
        // 节气时间点(儒略日), 合朔周期(天),        对应历法及误差说明
        (1640650.479938, 15.21842500), // -221-11-09 h=0.01709 古历·秦汉
        (1642476.703182, 15.21874996), // -216-11-09 h=0.01557 古历·秦汉
        (1683430.515601, 15.218750011), // -104-12-25 h=0.01560 汉书·律历志(太初历)平气平朔 回归年=365.25000
        (1752157.640664, 15.218749978), //   85-02-23 h=0.01559 后汉书·律历志(四分历) 回归年=365.25000
        (1807675.003759, 15.218620279), //  237-02-22 h=0.00010 晋书·律历志(景初历) 回归年=365.24689
        (1883627.765182, 15.218612292), //  445-02-03 h=0.00026 宋书·律历志(何承天元嘉历) 回归年=365.24670
        (1907369.128100, 15.218449176), //  510-02-03 h=0.00027 宋书·律历志(祖冲之大明历) 回归年=365.24278
        (1936603.140413, 15.218425000), //  590-02-17 h=0.00149 随书·律历志(开皇历) 回归年=365.24220
        (1939145.524180, 15.218466998), //  597-02-03 h=0.00121 随书·律历志(大业历) 回归年=365.24321
        (1947180.798300, 15.218524844), //  619-02-03 h=0.00052 新唐书·历志(戊寅元历)平气定朔 回归年=365.24460
        (1964362.041824, 15.218533526), //  666-02-17 h=0.00059 新唐书·历志(麟德历) 回归年=365.24480
        (1987372.340971, 15.218513908), //  729-02-16 h=0.00096 新唐书·历志(大衍历,至德历) 回归年=365.24433
        (1999653.819126, 15.218530782), //  762-10-03 h=0.00093 新唐书·历志(五纪历) 回归年=365.24474
        (2007445.469786, 15.218535181), //  784-02-01 h=0.00059 新唐书·历志(正元历,观象历) 回归年=365.24484
        (2021324.917146, 15.218526248), //  822-02-01 h=0.00022 新唐书·历志(宣明历) 回归年=365.24463
        (2047257.232342, 15.218519654), //  893-01-31 h=0.00015 新唐书·历志(崇玄历) 回归年=365.24447
        (2070282.898213, 15.218425000), //  956-02-16 h=0.00149 旧五代·历志(钦天历) 回归年=365.24220
        (2073204.872850, 15.218515221), //  964-02-16 h=0.00166 宋史·律历志(应天历) 回归年=365.24437
        (2080144.500926, 15.218530782), //  983-02-16 h=0.00093 宋史·律历志(乾元历) 回归年=365.24474
        (2086703.688963, 15.218523776), // 1001-01-31 h=0.00067 宋史·律历志(仪天历,崇天历) 回归年=365.24457
        (2110033.182763, 15.218425000), // 1064-12-15 h=0.00669 宋史·律历志(明天历) 回归年=365.24220
        (2111190.300888, 15.218425000), // 1068-02-15 h=0.00149 宋史·律历志(崇天历) 回归年=365.24220
        (2113731.271005, 15.218515671), // 1075-01-30 h=0.00038 李锐补修(奉元历) 回归年=365.24438
        (2120670.840263, 15.218425000), // 1094-01-30 h=0.00149 宋史·律历志 回归年=365.24220
        (2123973.309063, 15.218425000), // 1103-02-14 h=0.00669 李锐补修(占天历) 回归年=365.24220
        (2125068.997336, 15.218477932), // 1106-02-14 h=0.00056 宋史·律历志(纪元历) 回归年=365.24347
        (2136026.312633, 15.218472436), // 1136-02-14 h=0.00088 宋史·律历志(统元历,乾道历,淳熙历) 回归年=365.24334
        (2156099.495538, 15.218425000), // 1191-01-29 h=0.00149 宋史·律历志(会元历) 回归年=365.24220
        (2159021.324663, 15.218425000), // 1199-01-29 h=0.00149 宋史·律历志(统天历) 回归年=365.24220
        (2162308.575254, 15.218461742), // 1208-01-30 h=0.00146 宋史·律历志(开禧历) 回归年=365.24308
        (2178485.706538, 15.218425000), // 1252-05-15 h=0.04606 淳祐历 回归年=365.24220
        (2178759.662849, 15.218445786), // 1253-02-13 h=0.00231 会天历 回归年=365.24270
        (2185334.020800, 15.218425000), // 1271-02-13 h=0.00520 宋史·律历志(成天历) 回归年=365.24220
        (2187525.481425, 15.218425000), // 1277-02-12 h=0.00520 本天历 回归年=365.24220
        (2188621.191481, 15.218437494), // 1280-02-13 h=0.00015 元史·历志(郭守敬授时历) 回归年=365.24250
        (2322147.760000, 15.218425000), // 1645-09-21 边界值
    ];

    /// 低精度定朔计算
    /// 在2000年至600年, 误差在2小时以内(仍比古代历法精准很多)
    pub fn so_low(&self, w: f64) -> f64 {
        let v = 7771.37714500204;
        let mut t = (w + 1.08472) / v;

        // 修正值计算
        t -= (-0.0000331 * t * t
            + 0.10976 * (0.785 + 8328.6914 * t).cos()
            + 0.02224 * (0.187 + 7214.0629 * t).cos()
            - 0.03342 * (4.669 + 628.3076 * t).cos())
            / v
            + (32.0 * (t + 1.8) * (t + 1.8) - 20.0) / 86400.0 / 36525.0;

        t * 36525.0 + 8.0 / 24.0
    }

    /// 低精度定气计算
    /// 最大误差小于30分钟，平均5分
    pub fn qi_low(&self, w: f64) -> f64 {
        let v = 628.3319653318;
        let mut t = (w - 4.895062166) / v;

        // 第二次估算,误差2小时以内
        t -= (53.0 * t * t
            + 334116.0 * (4.67 + 628.307585 * t).cos()
            + 2061.0 * (2.678 + 628.3076 * t).cos() * t)
            / v
            / 10000000.0;

        let l = 48950621.66 + 6283319653.318 * t + 53.0 * t * t
            + 334166.0 * (4.669257 + 628.307585 * t).cos() // 地球椭圆轨道级数展开
            + 3489.0 * (4.6261 + 1256.61517 * t).cos()
            + 2060.6 * (2.67823 + 628.307585 * t).cos() * t  // 一次泊松项
            - 994.0 - 834.0 * (2.1824 - 33.75705 * t).sin(); // 光行差与章动修正

        t -= (l / 10000000.0 - w) / 628.332
            + (32.0 * (t + 1.8) * (t + 1.8) - 20.0) / 86400.0 / 36525.0;

        t * 36525.0 + 8.0 / 24.0
    }

    /// 较高精度气计算
    pub fn qi_high(&self, w: f64) -> f64 {
        // 使用天文算法计算节气时刻
        let mut t = XL::s_alon_t2(w) * 36525.0;
        t = t - dt_t(t) + 8.0 / 24.0;

        let v = ((t + 0.5) % 1.0) * 86400.0;
        if !(1200.0..=86400.0 - 1200.0).contains(&v) {
            t = XL::s_alon_t(w) * 36525.0 - dt_t(t) + 8.0 / 24.0;
        }

        t
    }

    /// 较高精度朔计算
    pub fn so_high(&self, w: f64) -> f64 {
        // 使用天文算法计算月相时刻
        let mut t = XL::ms_alon_t2(w) * 36525.0;
        t = t - dt_t(t) + 8.0 / 24.0;

        let v = ((t + 0.5) % 1.0) * 86400.0;
        if !(1800.0..=86400.0 - 1800.0).contains(&v) {
            t = XL::ms_alon_t(w) * 36525.0 - dt_t(t) + 8.0 / 24.0;
        }

        t
    }

    // /// 朔日计算
    // pub fn calc_new_moon(&self, jd: f64) -> f64 {
    //     so_low(jd);
    // }

    // /// 节气计算
    // pub fn calc_solar_term(&self, jd: f64) -> f64 {
    //     qi_low(jd);
    // }

    /// 计算节气或朔日
    /// # Arguments
    /// * `jd` - 儒略日数
    /// * `is_solar_term` - true表示计算节气,false表示计算朔日
    pub fn calc(&self, mut jd: f64, is_solar_term: bool) -> f64 {
        // 转换为完整儒略日
        jd += 2451545.0;

        // 获取拟合参数
        let (kb, pc): (&[(f64, f64)], f64) = if is_solar_term {
            (&Self::QI_KB[..], 7.0)
        } else {
            (&Self::SUO_KB[..], 14.0)
        };

        // 计算边界值
        let f1 = kb[0].0 - pc;
        let f2 = kb[kb.len() - 1].0 - pc;
        let f3 = 2436935.0;

        // 超出平气朔表范围,使用现代天文算法
        if jd < f1 || jd >= f3 {
            if is_solar_term {
                // 2451259是1999.3.21,太阳视黄经为0,春分
                let n = ((jd + pc - 2451259.0) / 365.2422 * 24.0).floor();
                let w = n * std::f64::consts::PI / 12.0;
                (self.qi_high(w) + 0.5).floor()
            } else {
                // 2451551是2000.1.7的朔日,黄经差为0
                let n = ((jd + pc - 2451551.0) / 29.5306).floor();
                let w = n * 2.0 * std::f64::consts::PI;
                (self.so_high(w) + 0.5).floor()
            }
        }
        // 平气或平朔
        else if jd >= f1 && jd < f2 {
            let mut i = 0;
            while i < kb.len() - 1 {
                if jd + pc < kb[i + 2].0 {
                    break;
                }
                i += 2;
            }
            let mut d = kb[i].0 + kb[i].1 * ((jd + pc - kb[i].0) / kb[i].1).floor();
            d = d.floor() + 0.5;

            // 修正-103年1月24日的朔日
            if d == 1683460.0 {
                d += 1.0;
            }
            d - 2451545.0
        }
        // 定气或定朔
        else if jd >= f2 && jd < f3 {
            let d = if is_solar_term {
                // 2451259是1999.3.21,太阳视黄经为0,春分
                let n = ((jd + pc - 2451259.0) / 365.2422 * 24.0).floor();
                let w = n * std::f64::consts::PI / 12.0;
                (self.qi_low(w) + 0.5).floor()
            } else {
                // 2451551是2000.1.7的朔日,黄经差为0
                let n = ((jd + pc - 2451551.0) / 29.5306).floor();
                let w = n * 2.0 * std::f64::consts::PI;
                (self.so_low(w) + 0.5).floor()
            };

            // 查找修正值
            let idx = if is_solar_term {
                ((jd - f2) / 365.2422 * 24.0) as usize
            } else {
                ((jd - f2) / 29.5306) as usize
            };

            let n = if is_solar_term {
                &self.qb[idx..=idx]
            } else {
                &self.sb[idx..=idx]
            };

            match n {
                "1" => d + 1.0,
                "2" => d - 1.0,
                _ => d,
            }
        } else {
            jd
        }
    }

    /// 农历排月序计算,定出农历年历
    /// 可计算的范围：两个冬至之间(冬至一 <= d < 冬至二)
    pub fn calc_year(&mut self, jd: f64) {
        // let mut year = LunarYear::default();

        // 该年的节气计算
        // 355是2000.12冬至,得到较靠近jd的冬至估计值
        let mut w = ((jd - 355.0 + 183.0) / 365.2422).floor() * 365.2422 + 355.0;

        // 确保冬至在jd之前
        if self.calc(w, true) > jd {
            w -= 365.2422;
        }

        // 25个节气时刻(北京时),从冬至开始到下一个冬至以后
        for i in 0..25 {
            self.solar_terms
                .push(self.calc(w + 15.2184 * i as f64, true));
        }

        // 补算前两个节气,确保一年中所有月份的"气"都被计算
        self.solar_terms.push(self.calc(w - 15.2, true)); // pe1
        self.solar_terms.push(self.calc(w - 30.4, true)); // pe2

        // 计算较靠近冬至的朔日
        let mut new_moon = self.calc(self.solar_terms[0], false);
        if new_moon > self.solar_terms[0] {
            new_moon -= 29.53;
        }

        // 计算该年所有朔日(共15个月)
        for i in 0..15 {
            self.new_moons
                .push(self.calc(new_moon + 29.5306 * i as f64, false));
        }

        // 月大小计算
        self.leap = 0;
        for i in 0..14 {
            // 月大小
            self.month_days
                .push((self.new_moons[i + 1] - self.new_moons[i]) as i32);

            // 月序初始化
            self.month_index.push(i as usize);
            self.month_names.push(i.to_string());
        }

        // -721年至-104年的后九月及月建问题处理
        let year_num = ((self.solar_terms[0] + 10.0 + 180.0) / 365.2422).floor() as i32 + 2000;
        if (-721..=-104).contains(&year_num) {
            self.fix_ancient_months(year_num);
            return;
        }

        // 无中气置闰法确定闰月
        // 第13月的月末没超过冬至,说明该年有13个月
        if self.new_moons[13] <= self.solar_terms[24] {
            // 在13个月中找第1个没有中气的月份
            let mut i = 1;
            while self.new_moons[i + 1] > self.solar_terms[2 * i] && i < 13 {
                i += 1;
            }
            self.leap = i as i32;

            // 更新月名
            for j in i..14 {
                // let month_num = self.month_names[j].parse::<i32>().unwrap();
                self.month_names[j] = (self.month_index[j] - 1 ).to_string();
            }
        }

        // for j in 0..14 { println!("ym:{}",self.month_names[j]); }

        // 月建名称转换
        self.convert_month_names(); //year_num
    }

    /// 修正古代历法月名(-721年至-104年)
    fn fix_ancient_months(&mut self, year: i32) {
        let mut ns = Vec::new();

        // 计算三年的历元
        for i in 0..3 {
            let y = year + i - 1;

            if y >= -721 {
                // 春秋历,ly为-722.12.17
                ns.push(CalendarParam {
                    jd: self.calc(
                        1457698.0 - 2451545.0
                            + ((0.342 + (y + 721) as f64 * 12.368422).floor() * 29.5306),
                        false,
                    ),
                    month_name: "十三".to_string(),
                    month_index: 2,
                });
            }

            if y >= -479 {
                // 战国历,ly为-480.12.11
                ns.push(CalendarParam {
                    jd: self.calc(
                        1546083.0 - 2451545.0
                            + ((0.500 + (y + 479) as f64 * 12.368422).floor() * 29.5306),
                        false,
                    ),
                    month_name: "十三".to_string(),
                    month_index: 2,
                });
            }

            if y >= -220 {
                // 秦汉历,ly为-221.10.31
                ns.push(CalendarParam {
                    jd: self.calc(
                        1640641.0 - 2451545.0
                            + ((0.866 + (y + 220) as f64 * 12.369000).floor() * 29.5306),
                        false,
                    ),
                    month_name: "后九".to_string(),
                    month_index: 11,
                });
            }
        }

        // 更新月名
        for i in 0..14 {
            // 找到对应的历法参数
            let mut idx = 2;
            while self.new_moons[i] < ns[idx].jd {
                idx -= 1;
            }

            // 计算月积数
            let acc = ((self.new_moons[i] - ns[idx].jd + 15.0) / 29.5306).floor();

            if acc < 12.0 {
                // 正常月
                let month_num = ((acc as i32 + ns[idx].month_index) % 12) as usize;
                self.month_names[i] = YUE_MING[month_num].to_string();
            } else {
                // 闰月
                self.month_names[i] = ns[idx].month_name.clone();
            }
        }
    }

    /// 月建名称转换
    fn convert_month_names(&mut self) {
        
        // 转换各月名称
        for i in 0..14 {
            // Dm为初一的儒略日,v2为月建序号
            let dm = self.new_moons[i] + 2451545.0;
            let v2 = self.month_names[i].parse::<i32>().unwrap();

            // 根据月建序号获取默认月名称
            let mut mc = YUE_MING[v2 as usize % 12].to_string();

            // 根据不同历法时期调整月名称
            if (1724360.0..=1729794.0).contains(&dm) {
                // 8.01.15至23.12.02 建子为十二,其它顺推
                mc = YUE_MING[((v2 + 1) % 12) as usize].to_string();
            } else if (1807724.0..=1808699.0).contains(&dm) {
                // 237.04.12至239.12.13 建子为十二,其它顺推
                mc = YUE_MING[((v2 + 1) % 12) as usize].to_string();
            } else if (1999349.0..=1999467.0).contains(&dm) {
                // 761.12.02至762.03.30 建子为正月,其它顺推
                mc = YUE_MING[((v2 + 2) % 12) as usize].to_string();
            } else if (1973067.0..=1977052.0).contains(&dm) {
                // 689.12.18至700.11.15 建子为正月,建寅为一月,其它不变
                if v2 % 12 == 0 {
                    mc = "正".to_string();
                }
                if v2 == 2 {
                    mc = "一".to_string();
                }
            }

            // 调整特殊月名称
            // 239.12.13及23.12.02均为十二月,避免两个连续十二月
            // if (dm - 1729794.0).abs() < 0.1 || (dm - 1808699.0).abs() < 0.1 {
            //     mc = "拾贰".to_string();
            // }

            self.month_names[i] = mc;
        }
    }

    // /// 朔日修正
    // pub fn correct_new_moon(&self, jd: f64) -> f64 {
    //     let pc = 14.0;
    //     let (f1, f2, f3) = (self.suo_kb[0].0 - pc,
    //                       self.suo_kb[self.suo_kb.len()-1].0 - pc,
    //                       2436935.0);

    //     // 超出平朔表范围,使用现代方法计算
    //     if jd < f1 || jd >= f3 {
    //         let w = ((jd + 8.0) / 29.5306).floor() * 2.0 * PI;
    //         return self.so_accurate(w);
    //     }

    //     // 平朔计算
    //     if jd >= f1 && jd < f2 {
    //         let mut i = 0;
    //         while i < self.suo_kb.len()-1 {
    //             if jd + pc < self.suo_kb[i+2].0 { break; }
    //             i += 2;
    //         }
    //         let mut d = self.suo_kb[i].0 + self.suo_kb[i].1
    //                  * ((jd + pc - self.suo_kb[i].0) / self.suo_kb[i].1).floor();
    //         d = d.floor() + 0.5;

    //         // 修正-103年1月24日的朔日
    //         if d == 1683460.0 { d += 1.0; }

    //         return d - 2451545.0;
    //     }

    //     // 定朔计算
    //     if jd >= f2 && jd < f3 {
    //         let d = self.calc_new_moon(jd);
    //         let n = &self.sb[(((jd-f2)/29.5306) as usize)..=((jd-f2)/29.5306) as usize];

    //         match n {
    //             "1" => d + 1.0,
    //             "2" => d - 1.0,
    //             _ => d
    //         }
    //     } else {
    //         jd
    //     }
    // }

    // /// 节气修正
    // pub fn correct_solar_term(&self, jd: f64) -> f64 {
    //     let pc = 7.0;
    //     let (f1, f2, f3) = (self.qi_kb[0].0 - pc,
    //                       self.qi_kb[self.qi_kb.len()-1].0 - pc,
    //                       2436935.0);

    //     // 超出平气表范围,使用现代方法计算
    //     if jd < f1 || jd >= f3 {
    //         let w = ((jd + 293.0)/365.2422*24.0).floor() * PI / 12.0;
    //         let d = self.qi_accurate(w);

    //         if (d - jd) > 5.0 {
    //             return self.qi_accurate(w - PI/12.0);
    //         }
    //         if (d - jd) < -5.0 {
    //             return self.qi_accurate(w + PI/12.0);
    //         }
    //         return d;
    //     }

    //     // 平气计算
    //     if jd >= f1 && jd < f2 {
    //         let mut i = 0;
    //         while i < self.qi_kb.len()-1 {
    //             if jd + pc < self.qi_kb[i+2].0 { break; }
    //             i += 2;
    //         }
    //         let d = self.qi_kb[i].0 + self.qi_kb[i].1
    //                * ((jd + pc - self.qi_kb[i].0) / self.qi_kb[i].1).floor();
    //         return d.floor() + 0.5 - 2451545.0;
    //     }

    //     // 定气计算
    //     if jd >= f2 && jd < f3 {
    //         let d = self.calc_solar_term(jd);
    //         let n = &self.qb[(((jd-f2)/365.2422*24.0) as usize)..=((jd-f2)/365.2422*24.0) as usize];

    //         match n {
    //             "1" => d + 1.0,
    //             "2" => d - 1.0,
    //             _ => d
    //         }
    //     } else {
    //         jd
    //     }
    // }

    // /// 根据朔直线拟合参数计算合朔时刻
    // pub fn calc_new_moon_accurate(&self, w: f64) -> f64 {
    //     // 寻找适用的拟合参数
    //     let mut i = 0;
    //     while i < Self::SUO_KB.len()-1 {
    //         if w < Self::SUO_KB[i+1].0 {
    //             break;
    //         }
    //         i += 1;
    //     }

    //     // 计算合朔时刻
    //     let (jd0, p) = Self::SUO_KB[i];
    //     let n = ((w - jd0) / p).floor();

    //     jd0 + n * p
    // }

    // /// 根据气直线拟合参数计算节气时刻
    // pub fn calc_solar_term_accurate(&self, w: f64) -> f64 {
    //     // 寻找适用的拟合参数
    //     let mut i = 0;
    //     while i < Self::QI_KB.len()-1 {
    //         if w < Self::QI_KB[i+1].0 {
    //             break;
    //         }
    //         i += 1;
    //     }

    //     // 计算节气时刻
    //     let (jd0, p) = Self::QI_KB[i];
    //     let n = ((w - jd0) / p).floor();

    //     jd0 + n * p
    // }

    /// 定朔修正表解压数据
    const SUO_S: &'static str = "EqoFscDcrFpmEsF2DfFideFelFpFfFfFiaipqti1ksttikptikqckstekqttgkqttgkqteksttikptikq2fjstgjqtt\
         jkqttgkqtekstfkptikq2tijstgjiFkirFsAeACoFsiDaDiADc1AFbBfgdfikijFifegF1FhaikgFag1E2btaieeibgg\
         iffdeigFfqDfaiBkF1kEaikhkigeidhhdiegcFfakF1ggkidbiaedksaFffckekidhhdhdikcikiakicjF1deedFhFcc\
         gicdekgiFbiaikcfi1kbFibefgEgFdcFkFeFkdcfkF1kfkcickEiFkDacFiEfbiaejcFfffkhkdgkaiei1ehigikhdFi\
         kfckF1dhhdikcfgjikhfjicjicgiehdikGcikggcifgiejF1jkieFhegikggcikFegiegkfjebhigikggcikdgkaFkij\
         cfkcikfkcifikiggkaeeigefkcdfcfkhkdgkegieidhijcFfakhfgeidieidiegikhfkfckfcjbdehdikggikgkfkicji\
         cjF1dbidikFiggcifgiejkiegkigcdiegfggcikdbgfgefjF1kfegikggcikdgFkeeijcfkcikfkekcikdgkabhkFika\
         ffcfkhkdgkegbiaekfkiakicjhfgqdq2fkiakgkfkhfkfcjiekgFebicggbedF1jikejbbbiakgbgkacgiejkijjgigf\
         iakggfggcibFifjefjF1kfekdgjcibFeFkijcfkfhkfkeaieigekgbhkfikidfcjeaibgekgdkiffiffkiakF1jhbakg\
         dki1dj1ikfkicjicjieeFkgdkicggkighdF1jfgkgfgbdkicggfggkidFkiekgijkeigfiskiggfaidheigF1jekijci\
         kickiggkidhhdbgcfkFikikhkigeidieFikggikhkffaffijhidhhakgdkhkijF1kiakF1kfheakgdkifiggkigicjiej\
         kieedikgdfcggkigieeiejfgkgkigbgikicggkiaideeijkefjeijikhkiggkiaidheigcikaikffikijgkiahi1hhdik\
         gjfifaakekighie1hiaikggikhkffakicjhiahaikggikhkijF1kfejfeFhidikggiffiggkigicjiekgieeigikggiff\
         iggkidheigkgfjkeigiegikifiggkidhedeijcfkFikikhkiggkidhh1ehigcikaffkhkiggkidhh1hhigikekfiFkFik\
         cidhh1hitcikggikhkfkicjicghiediaikggikhkijbjfejfeFhaikggifikiggkigiejkikgkgieeigikggiffiggkig\
         ieeigekijcijikggifikiggkideedeijkefkfckikhkiggkidhh1ehijcikaffkhkiggkidhh1hhigikhkikFikfckctd\
         hh1hiaikgjikhfjicjicgiehdikcikggifikigiejfejkieFhegikggifikiggfghigkfjeijkhigikggifikiggkigie\
         eijcijcikfksikifikiggkidehdeijcfdckikhkiggkhghh1ehijikifffffkhsFngErD1pAfBoDd1BlEtFqA2AqoEpDq\
         ElAEsEeB2BmADlDkqBtC1FnEpDqnEmFsFsAFnllBbFmDsDiCtDmAB2BmtCgpEplCpAEiBiEoFqFtEqsDcCnFtADnFlE\
         gdkEgmEtEsCtDmADqFtAFrAtEcCqAE1BoFqC1F1DrFtBmFtAC2ACnFaoCgADcADcCcFfoFtDlAFgmFqBq2bpEoAEmkq\
         nEeCtAE1bAEqgDfFfCrgEcBrACfAAABqAAB1AAClEnFeCtCgAADqDoBmtAAACbFiAAADsEtBqAB2FsDqpFqEmFsCeDt\
         FlCeDtoEpClEqAAFrAFoCgFmFsFqEnAEcCqFeCtFtEnAEeFtAAEkFnErAABbFkADnAAeCtFeAfBoAEpFtAABtFqAApD\
         cCGJ";

    /// 定气修正表解压数据
    const QI_S: &'static str = "FrcFs22AFsckF2tsDtFqEtF1posFdFgiFseFtmelpsEfhkF2anmelpFlF1ikrotcnEqEq2FfqmcDsrFor22FgFrcgDsc\
         Fs22FgEeFtE2sfFs22sCoEsaF2tsD1FpeE2eFsssEciFsFnmelpFcFhkF2tcnEqEpFgkrotcnEqrEtFermcDsrE222Fg\
         Bmcmr22DaEfnaF222sD1FpeForeF2tssEfiFpEoeFssD1iFstEqFppDgFstcnEqEpFg11FscnEqrAoAF2ClAEsDmDtCt\
         BaDlAFbAEpAAAAAD2FgBiBqoBbnBaBoAAAAAAAEgDqAdBqAFrBaBoACdAAf1AACgAAAeBbCamDgEifAE2AABa1C1BgF\
         diAAACoCeE1ADiEifDaAEqAAFe1AcFbcAAAAAF1iFaAAACpACmFmAAAAAAAACrDaAAADG0";

    /// 解压缩定气朔修正表数据
    fn decompress(&self, s: &str) -> String {
        let mut result = String::new();
        let chars = s.chars();

        for c in chars {
            match c {
                'J' => result.push_str("00"),
                'I' => result.push_str("000"),
                'H' => result.push_str("0000"),
                'G' => result.push_str("00000"),
                't' => result.push_str("02"),
                's' => result.push_str("002"),
                'r' => result.push_str("0002"),
                'q' => result.push_str("00002"),
                'p' => result.push_str("000002"),
                'o' => result.push_str("0000002"),
                'n' => result.push_str("00000002"),
                'm' => result.push_str("000000002"),
                'l' => result.push_str("0000000002"),
                'k' => result.push_str("01"),
                'j' => result.push_str("0101"),
                'i' => result.push_str("001"),
                'h' => result.push_str("001001"),
                'g' => result.push_str("0001"),
                'f' => result.push_str("00001"),
                'e' => result.push_str("000001"),
                'd' => result.push_str("0000001"),
                'c' => result.push_str("00000001"),
                'b' => result.push_str("000000001"),
                'a' => result.push_str("0000000001"),
                'A' => result.push_str("000000000000000000000000000000"),
                'B' => result.push_str("00000000000000000000"),
                'C' => result.push_str("0000000000000000"),
                'D' => result.push_str("000000000000"),
                'E' => result.push_str("00000000"),
                'F' => result.push_str("0000"),
                _ => result.push(c),
            }
        }

        result
    }

    /// 初始化朔气修正表
    pub fn new() -> Self {
        //                   /// 闰月位置
        // pub leap: i32,
        // /// 各月名称
        // pub month_names: Vec<String>,
        // /// 中气表,其中.liqiu是节气立秋的儒略日,计算三伏时用到
        // pub solar_terms: Vec<f64>,
        // /// 合朔表
        // pub new_moons: Vec<f64>,
        // /// 各月大小
        // pub month_days: Vec<i32>,
        // /// 年计数
        // pub year_counts: Vec<i32>,
        let mut st = SolarTerms {
            sb: String::new(),
            qb: String::new(),
            // suo_kb: [(0.0, 0.0); 13],
            // qi_kb: [(0.0, 0.0); 30],
            leap: 0,
            month_names: Vec::new(),
            month_index: Vec::new(),
            solar_terms: Vec::new(),
            new_moons: Vec::new(),
            month_days: Vec::new(),
            year_counts: Vec::new(),
        };

        // 解压朔日修正表
        st.sb = st.decompress(Self::SUO_S);
        // 解压节气修正表
        st.qb = st.decompress(Self::QI_S);

        st
    }
}

// 定义一个结构体来存储历法参数
struct CalendarParam {
    jd: f64,            // 儒略日
    month_name: String, // 月名
    month_index: i32,   // 月索引
}
