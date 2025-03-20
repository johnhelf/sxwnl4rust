// use chrono::{Datelike};
use crate::lunar::lunar::LunarDay;
use lazy_static::lazy_static;

/// 节日信息结构
#[derive(Debug, Default, Clone)]
pub struct FestivalInfo {
    /// 重要节日(A级)
    pub major: Vec<String>,
    /// 一般节日(B级)
    pub minor: Vec<String>,
    /// 其他纪念日(C级)
    pub other: Vec<String>,
    /// 是否放假
    pub is_holiday: bool,
}

impl FestivalInfo {
    /// 添加节日
    fn add_festival(&mut self, festival_type: &str, name: &str) {
        match festival_type {
            "#" => {
                self.major.push(name.to_string());
                self.is_holiday = true;
            }
            "I" => self.minor.push(name.to_string()),
            "." => self.other.push(name.to_string()),
            _ => {}
        }
    }
}

lazy_static! {
    /// 周节日数据
    static ref WEEK_FESTIVALS: &'static str = concat!(
        "0150I世界麻风日|", // 一月最后一个星期日
        "0520.国际母亲节|",
        "0530I全国助残日|",
        "0630.父亲节|",
        "0730.被奴役国家周|",
        "0932I国际和平日|",
        "0940.国际聋人节 世界儿童日|",
        "0950I世界海事日|",
        "1011.国际住房日|",
        "1013I国际减轻自然灾害日|",
        "1144I感恩节"
    );

    /// 公历节日数据(按月份存储)
    static ref SOLAR_FESTIVALS: [&'static str; 12] = [
        // 1月
        "01#元旦|",
        // 2月
        "02I世界湿地日,10.国际气象节,14I情人节|",
        // 3月
        concat!(
            "01.国际海豹日,03.全国爱耳日,05.1963-9999学雷锋纪念日,08I妇女节,12I植树节,",
            "12.1925-9999孙中山逝世纪念日,14.国际警察日,15I1983-9999消费者权益日,",
            "17.中国国医节,17.国际航海日,21.世界森林日,21.消除种族歧视国际日,",
            "21.世界儿歌日,22I世界水日,23I世界气象日,24.1982-9999世界防治结核病日,",
            "25.全国中小学生安全教育日,30.巴勒斯坦国土日|"
        ),
        // 4月
        concat!(
            "01I1564-9999愚人节,01.全国爱国卫生运动月,01.税收宣传月,07I世界卫生日,",
            "22I世界地球日,23.世界图书和版权日,24.亚非新闻工作者日|"
        ),
        // 5月
        concat!(
            "01#1889-9999劳动节,04I青年节,05.碘缺乏病防治日,08.世界红十字日,",
            "12I国际护士节,15I国际家庭日,17.国际电信日,18.国际博物馆日,",
            "20.全国学生营养日,23.国际牛奶日,31I世界无烟日|"
        ),
        // 6月
        concat!(
            "01I1925-9999国际儿童节,05.世界环境保护日,06.全国爱眼日,",
            "17.防治荒漠化和干旱日,23.国际奥林匹克日,25.全国土地日,26I国际禁毒日|"
        ),
        // 7月
        concat!(
            "01I1997-9999香港回归纪念日,01I1921-9999中共诞辰,01.世界建筑日,",
            "02.国际体育记者日,07I1937-9999抗日战争纪念日,11I世界人口日,30.非洲妇女日|"
        ),
        // 8月
        "01I1927-9999建军节,08.中国男子节|",
        // 9月
        concat!(
            "03I1945-9999抗日战争胜利纪念,08.国际扫盲日,08.国际新闻工作者日,",
            "09.毛泽东逝世纪念,10I中国教师节,14.世界清洁地球日,16.国际臭氧层保护日,",
            "18I九一八事变纪念日,20.国际爱牙日,27.世界旅游日,28I孔子诞辰|"
        ),
        // 10月
        concat!(
            "01#1949-9999国庆节,01.国际老人节,02#1949-9999国庆节假日,",
            "03#1949-9999国庆节假日,04.世界动物日,06.老人节,08.全国高血压日,",
            "08.世界视觉日,09.世界邮政日,09.万国邮联日,10I辛亥革命纪念日,",
            "10.世界精神卫生日,13.世界保健日,14.世界标准日,15.国际盲人节,",
            "16.世界粮食日,17.世界消除贫困日,22.世界传统医药日,24.联合国日,31.世界勤俭日|"
        ),
        // 11月
        concat!(
            "07.十月社会主义革命纪念日,08.中国记者日,09.全国消防安全宣传教育日,",
            "10.世界青年节,11.国际科学与和平周,12.孙中山诞辰纪念日,14.世界糖尿病日,",
            "17.国际大学生节,17.世界学生节,20.彝族年,21.世界问候日,21.世界电视日,",
            "22.彝族年,29.国际声援巴勒斯坦人民国际日|"
        ),
        // 12月
        concat!(
            "01I1988-9999世界艾滋病日,03.世界残疾人日,05.国际经济和社会发展志愿人员日,",
            "08.国际儿童电视日,09.世界足球日,10.世界人权日,12I西安事变纪念日,",
            "13I南京大屠杀纪念日,20.澳门回归纪念,21.国际篮球日,24I平安夜,25I圣诞节,",
            "26.毛泽东诞辰纪念|"
        )
    ];
}

impl FestivalInfo {
    /// 计算节日(包含公历、农历、周历节日)
    pub fn calc_festivals(&mut self, lunar_day: &LunarDay) -> FestivalInfo {
        // let mut festivals = FestivalInfo::default();

        // 1. 处理公历节日
        self.add_solar_festivals(lunar_day);

        // 2. 处理农历节日
        self.add_lunar_festivals(lunar_day);

        // 3. 处理周历节日
        self.add_week_festivals(lunar_day);

        self.clone()
    }

    /// 添加公历节日
    fn add_solar_festivals(&mut self, date: &LunarDay) {
        // 获取月日字符串,如"0102"表示1月2日
        let md = format!("{:02}{:02}", date.day.m, date.day.d);

        // 获取当月节日表
        let month_festivals = SOLAR_FESTIVALS[(date.day.m as usize) - 1]
            .split(',')
            .collect::<Vec<_>>();

        // 遍历当月节日
        for fest in month_festivals {
            // 跳过空条目
            if fest.is_empty() {
                continue;
            }

            // 获取日期部分
            let day = &fest[0..2];
            if day != &md[2..] {
                continue;
            }

            // 获取节日类型和名称
            let fest_type = &fest[2..3];
            let name = if fest[3..].contains('-') {
                // 处理有年份范围的节日
                let year_range: Vec<_> = fest[3..].split('-').collect();
                let start_year = year_range[0].parse::<i32>().unwrap_or(0);
                let end_year = year_range[1]
                    .split(' ')
                    .next()
                    .unwrap()
                    .parse::<i32>()
                    .unwrap_or(9999);

                if date.day.y >= start_year && date.day.y <= end_year {
                    fest.split(' ').last().unwrap()
                } else {
                    continue;
                }
            } else {
                // 1850年以后的节日
                if date.day.y >= 1850 {
                    &fest[3..]
                } else {
                    continue;
                }
            };

            self.add_festival(fest_type, name);
        }

        // 周末放假
        if date.day.week == 0 || date.day.week == 6 {
            self.is_holiday = true;
        }
    }

    /// 添加周历节日
    fn add_week_festivals(&mut self, lunar: &LunarDay) {
        // 计算周序号
        let mut w = lunar.day.weeki;
        if lunar.day.week >= lunar.day.week0 {
            w += 1;
        }

        // w2用于处理月末周
        let mut w2 = w;
        if lunar.day.weeki == lunar.day.weekn - 1 {
            w2 = 5;
        }

        // 生成查找键
        // 格式:"0523"表示5月第2周周3
        let m0 = format!("{:02}", lunar.day.m); // 补零的月份
        let w_key = format!("{}{}{}", m0, w, lunar.day.week);
        let w2_key = format!("{}{}{}", m0, w2, lunar.day.week);

        // 遍历周历节日表
        for fest in WEEK_FESTIVALS.split('|') {
            // 获取节日键值(前4位)
            let fest_key = &fest[0..4];

            // 不匹配则跳过
            if fest_key != w_key && fest_key != w2_key {
                continue;
            }

            // 获取节日类型和名称
            let fest_type = &fest[4..5];
            let name = &fest[5..];

            self.add_festival(fest_type, name);
        }
    }

    /// 添加农历节日
    fn add_lunar_festivals(&mut self, lunar_day: &LunarDay) {
        // 农历日期字符串,如"正月初一"
        let lunar_date = format!(
            "{}月{}",
            lunar_day.lunar_month_name, lunar_day.lunar_day_name
        );

        // 非闰月的农历节日
        if !lunar_day.is_leap_month {
            match lunar_date.as_str() {
                "正月初一" => {
                    self.add_festival("#", "春节");
                }
                "正月初二" => {
                    self.add_festival(".", "大年初二");
                }
                "五月初五" => {
                    self.add_festival("#", "端午节");
                }
                "八月十五" => {
                    self.add_festival("#", "中秋节");
                }
                "正月十五" => {
                    self.add_festival(".", "元宵节");
                    self.add_festival(".", "上元节");
                    self.add_festival(".", "壮族歌墟节");
                    self.add_festival(".", "苗族踩山节");
                    self.add_festival(".", "达斡尔族卡钦");
                }
                "正月十六" => {
                    self.add_festival(".", "侗族芦笙节(至正月二十)");
                    // festivals.other.push("侗族芦笙节(至正月二十)".to_string());
                }
                "正月廿五" => {
                    self.add_festival(".", "填仓节");
                }
                "正月廿九" => {
                    self.add_festival(".", "送穷日");
                }
                "二月初一" => {
                    self.add_festival(".", "瑶族忌鸟节");
                }
                "二月初二" => {
                    self.add_festival(".", "春龙节(龙抬头)");
                    self.add_festival(".", "畲族会亲节");
                }
                "二月初八" => {
                    self.add_festival(".", "傈傈族刀杆节");
                }
                "三月初三" => {
                    self.add_festival("I", "北帝诞");
                    self.add_festival(".", "苗族黎族歌墟节");
                }
                "三月十五" => {
                    self.add_festival(".", "白族三月街(至三月二十)");
                }
                "三月廿三" => {
                    self.add_festival("I", "天后诞 妈祖诞");
                }
                "四月初八" => {
                    self.add_festival("I", "牛王诞");
                }
                "四月十八" => {
                    self.add_festival(".", "锡伯族西迁节");
                }
                "五月十三" => {
                    self.add_festival("I", "关帝诞");
                    self.add_festival(".", "阿昌族泼水节");
                }
                "五月廿二" => {
                    self.add_festival(".", "鄂温克族米阔鲁节");
                }
                "五月廿九" => {
                    self.add_festival(".", "瑶族达努节");
                }
                "六月初六" => {
                    self.add_festival("I", "姑姑节 天贶节");
                    self.add_festival(".", "壮族祭田节");
                    self.add_festival(".", "瑶族尝新节");
                }
                "六月廿四" => {
                    self.add_festival(".", "火把节、星回节(彝、白、佤、阿昌、纳西、基诺族)");
                }
                "七月初七" => {
                    self.add_festival("I", "七夕(中国情人节,乞巧节,女儿节)");
                }
                "七月十三" => {
                    self.add_festival(".", "侗族吃新节");
                }
                "七月十五" => {
                    self.add_festival("I", "中元节 鬼节");
                }
                "九月初九" => {
                    self.add_festival("I", "重阳节");
                }
                "十月初一" => {
                    self.add_festival("I", "祭祖节(十月朝)");
                }
                "十月十五" => {
                    self.add_festival("I", "下元节");
                }
                "十月十六" => {
                    self.add_festival(".", "瑶族盘王节");
                }
                "十二初八" => {
                    self.add_festival("I", "腊八节");
                }
                _ => {}
            }
        }

        // 判断除夕
        if lunar_day.next_month_name == "正" {
            if (lunar_day.lunar_day_name == "三十" && lunar_day.lunar_month_days == 30)
                || (lunar_day.lunar_day_name == "廿九" && lunar_day.lunar_month_days == 29)
            {
                self.add_festival("#", "除夕");
            }
            if lunar_day.lunar_day_name == "廿三" {
                self.add_festival("I", "小年");
            }
        }

        // 添加节气节日
        if !lunar_day.jie_qi.is_empty() {
            if lunar_day.jie_qi == "清明" {
                self.add_festival("#", "清明");
            } else {
                self.add_festival("I", lunar_day.jie_qi.as_str());
            }
        }

        // 添加数九等农历杂节
        if lunar_day.cur_dz >= 0 && lunar_day.cur_dz < 81 {
            let num = ["零", "一", "二", "三", "四", "五", "六", "七", "八", "九"];
            let n = lunar_day.cur_dz / 9 + 1;

            if lunar_day.cur_dz % 9 == 0 {
                self.add_festival("I", format!("『{}九』", num[n as usize]).as_str());
            } else {
                self.add_festival(
                    ".",
                    format!("{}九第{}天", num[n as usize], lunar_day.cur_dz % 9 + 1).as_str(),
                );
            }
        }

        // 计算三伏
        // let day_gz = self.calc_day_gan_zhi(solar_date.timestamp() as f64);
        let mut chars = lunar_day.gan_zhi_day.chars();
        let gz1 = chars.next().unwrap().to_string(); // 天干
        let gz2 = chars.next().unwrap().to_string(); // 地支

        if *gz1 == *"庚" {
            if lunar_day.cur_xz >= 20 && lunar_day.cur_xz < 30 {
                self.add_festival("I", "初伏");
            } else if lunar_day.cur_xz >= 30 && lunar_day.cur_xz < 40 {
                self.add_festival("I", "中伏");
            } else if lunar_day.cur_lq >= 0 && lunar_day.cur_lq < 10 {
                self.add_festival("I", "末伏");
            }
        }

        // 入梅、出梅
        if *gz1 == *"丙" && lunar_day.cur_mz >= 0 && lunar_day.cur_mz < 10 {
            self.add_festival("I", "入梅");
        }
        if *gz2 == *"未" && lunar_day.cur_xs >= 0 && lunar_day.cur_xs < 12 {
            self.add_festival("I", "出梅");
        }
    }
}
