use crate::eph::eph_base::{PI, RAD, rad2mrad};
// use crate::eph::prece::Prece;
use lazy_static::lazy_static;

//=================================章动计算=========================================
pub struct Nutation {}

lazy_static! {
    /// 章动计算表
    static ref NU_TAB: [(i32, i32, i32, i32, i32, i32, i32, i32, i32, i32, i32); 77] = [
        (0, 0, 0, 0, 1, -172064161, -174666, 33386, 92052331, 9086, 15377),
        (0, 0, 2,-2, 2, -13170906, -1675, -13696, 5730336, -3015, -4587),
        (0, 0, 2, 0, 2, -2276413, -234, 2796, 978459, -485, 1374),
        (0, 0, 0, 0, 2, 2074554, 207, -698, -897492, 470, -291),
        (0, 1, 0, 0, 0, 1475877, -3633, 11817, 73871, -184, -1924),
        (0, 1, 2,-2, 2, -516821, 1226, -524, 224386, -677, -174),
        (1, 0, 0, 0, 0, 711159, 73, -872, -6750, 0, 358),
        (0, 0, 2, 0, 1, -387298, -367, 380, 200728, 18, 318),
        (1, 0, 2, 0, 2, -301461, -36, 816, 129025, -63, 367),
        (0,-1, 2,-2, 2, 215829, -494, 111, -95929, 299, 132),

        (0, 0, 2,-2, 1, 128227, 137, 181, -68982, -9, 39),
        (-1, 0, 2, 0, 2, 123457, 11, 19, -53311, 32, -4),
        (-1, 0, 0, 2, 0, 156994, 10, -168, -1235, 0, 82),
        (1, 0, 0, 0, 1, 63110, 63, 27, -33228, 0, -9),
        (-1, 0, 0, 0, 1, -57976, -63, -189, 31429, 0, -75),
        (-1, 0, 2, 2, 2, -59641, -11, 149, 25543, -11, 66),
        (1, 0, 2, 0, 1, -51613, -42, 129, 26366, 0, 78),
        (-2, 0, 2, 0, 1, 45893, 50, 31, -24236, -10, 20),
        (0, 0, 0, 2, 0, 63384, 11, -150, -1220, 0, 29),
        (0, 0, 2, 2, 2, -38571, -1, 158, 16452, -11, 68),

        (0,-2, 2,-2, 2, 32481, 0, 0, -13870, 0, 0),
        (-2, 0, 0, 2, 0, -47722, 0, -18, 477, 0, -25),
        (2, 0, 2, 0, 2, -31046, -1, 131, 13238, -11, 59),
        (1, 0, 2,-2, 2, 28593, 0, -1, -12338, 10, -3),
        (-1, 0, 2, 0, 1, 20441, 21, 10, -10758, 0, -3),
        (2, 0, 0, 0, 0, 29243, 0, -74, -609, 0, 13),
        (0, 0, 2, 0, 0, 25887, 0, -66, -550, 0, 11),
        (0, 1, 0, 0, 1, -14053, -25, 79, 8551, -2, -45),
        (-1, 0, 0, 2, 1, 15164, 10, 11, -8001, 0, -1),
        (0, 2, 2,-2, 2, -15794, 72, -16, 6850, -42, -5),

        (0, 0,-2, 2, 0, 21783, 0, 13, -167, 0, 13),
        (1, 0, 0,-2, 1, -12873, -10, -37, 6953, 0, -14),
        (0,-1, 0, 0, 1, -12654, 11, 63, 6415, 0, 26),
        (-1, 0, 2, 2, 1, -10204, 0, 25, 5222, 0, 15),
        (0, 2, 0, 0, 0, 16707, -85, -10, 168, -1, 10),
        (1, 0, 2, 2, 2, -7691, 0,  44, 3268, 0, 19),
        (-2, 0, 2, 0, 0, -11024, 0, -14, 104, 0, 2),
        (0, 1, 2, 0, 2,  7566, -21, -11, -3250, 0, -5),
        (0, 0, 2, 2, 1, -6637, -11, 25, 3353, 0, 14),
        (0,-1, 2, 0, 2, -7141, 21, 8, 3070, 0, 4),

        (0, 0, 0, 2, 1,      -6302,     -11,      2,     3272,     0,     4),
        (1, 0, 2,-2, 1,       5800,      10,      2,    -3045,     0,    -1),
        (2, 0, 2,-2, 2,       6443,       0,     -7,    -2768,     0,    -4),
        (-2, 0, 0, 2, 1,      -5774,     -11,    -15,     3041,     0,    -5),
        (2, 0, 2, 0, 1,      -5350,       0,     21,     2695,     0,    12),
        (0,-1, 2,-2, 1,      -4752,     -11,     -3,     2719,     0,    -3),
        (0, 0, 0,-2, 1,      -4940,     -11,    -21,     2720,     0,    -9),
        (-1,-1, 0, 2, 0,       7350,       0,     -8,      -51,     0,     4),
        (2, 0, 0,-2, 1,       4065,       0,      6,    -2206,     0,     1),
        (1, 0, 0, 2, 0,       6579,       0,    -24,     -199,     0,     2),

        (0, 1, 2,-2, 1,       3579,       0,      5,    -1900,     0,     1),
        (1,-1, 0, 0, 0,       4725,       0,     -6,      -41,     0,     3),
        (-2, 0, 2, 0, 2,      -3075,       0,     -2,     1313,     0,    -1),
        (3, 0, 2, 0, 2,      -2904,       0,     15,     1233,     0,     7),
        (0,-1, 0, 2, 0,       4348,       0,    -10,      -81,     0,     2),
        (1,-1, 2, 0, 2,      -2878,       0,      8,     1232,     0,     4),
        (0, 0, 0, 1, 0,      -4230,       0,      5,      -20,     0,    -2),
        (-1,-1, 2, 2, 2,      -2819,       0,      7,     1207,     0,     3),
        (-1, 0, 2, 0, 0,      -4056,       0,      5,       40,     0,    -2),
        (0,-1, 2, 2, 2,      -2647,       0,     11,     1129,     0,     5),

        (-2, 0, 0, 0, 1,      -2294,       0,    -10,     1266,     0,    -4),
        (1, 1, 2, 0, 2,       2481,       0,     -7,    -1062,     0,    -3),
        (2, 0, 0, 0, 1,       2179,       0,     -2,    -1129,     0,    -2),
        (-1, 1, 0, 1, 0,       3276,       0,      1,       -9,     0,     0),
        (1, 1, 0, 0, 0,      -3389,       0,      5,       35,     0,    -2),
        (1, 0, 2, 0, 0,       3339,       0,    -13,     -107,     0,     1),
        (-1, 0, 2,-2, 1,      -1987,       0,     -6,     1073,     0,    -2),
        (1, 0, 0, 0, 2,      -1981,       0,      0,      854,     0,     0),
        (-1, 0, 0, 1, 0,       4026,       0,   -353,     -553,     0,  -139),
        (0, 0, 2, 1, 2,       1660,       0,     -5,     -710,     0,    -2),

        (-1, 0, 2, 4, 2,      -1521,       0,      9,      647,     0,     4),
        (-1, 1, 0, 1, 1,       1314,       0,      0,     -700,     0,     0),
        (0,-2, 2,-2, 1,      -1283,       0,      0,      672,     0,     0),
        (1, 0, 2, 2, 1,      -1331,       0,      8,      663,     0,     4),
        (-2, 0, 2, 2, 2,       1383,       0,     -2,     -594,     0,    -2),
        (-1, 0, 0, 0, 2,       1405,       0,      4,     -610,     0,     2),
        (1, 1, 2,-2, 2,       1290,       0,      0,     -556,     0,     0)

    ];


}

/// 章动计算表B
const NUT_B: [(f64, f64, f64, f64, f64); 10] = [
    (2.1824, -33.75705, 36e-6, -1720.0, 920.0),
    (3.5069, 1256.66393, 11e-6, -132.0, 57.0),
    (1.3375, 16799.4182, -51e-6, -23.0, 10.0),
    (4.3649, -67.5141, 72e-6, 21.0, -9.0),
    (0.04, -628.302, 0.0, -14.0, 0.0),
    (2.36, 8328.691, 0.0, 7.0, 0.0),
    (3.46, 1884.966, 0.0, -5.0, 2.0),
    (5.44, 16833.175, 0.0, -4.0, 2.0),
    (3.69, 25128.110, 0.0, -3.0, 0.0),
    (3.55, 628.362, 0.0, 2.0, 0.0),
];

impl Nutation {
    /// 计算地球章动效应
    ///
    /// # Arguments
    ///
    /// * `t` - J2000.0起算的儒略世纪数
    /// * `zq` - 只计算周期大于zq(天)的项,如果为0则计算所有项
    ///
    /// # Returns
    ///
    /// 返回包含两个元素的数组 [黄经章动, 交角章动]
    pub fn nutation(t: f64, zq: i32) -> [f64; 2] {
        //t的二、三、四次方
        let t2 = t * t;
        let t3 = t2 * t;
        let t4 = t3 * t;

        // 五个基本参数
        let l =
            485868.249036 + 1717915923.2178 * t + 31.8792 * t2 + 0.051635 * t3 - 0.00024470 * t4;
        let l1 = 1287104.79305 + 129596581.0481 * t - 0.5532 * t2 + 0.000136 * t3 - 0.00001149 * t4;
        let f =
            335779.526232 + 1739527262.8478 * t - 12.7512 * t2 - 0.001037 * t3 + 0.00000417 * t4;
        let d = 1072260.70369 + 1602961601.2090 * t - 6.3706 * t2 + 0.006593 * t3 - 0.00003169 * t4;
        let om = 450160.398036 - 6962890.5431 * t + 7.4722 * t2 + 0.007702 * t3 - 0.00005939 * t4;

        let mut dl = 0.0;
        let mut de = 0.0;

        // 遍历章动表计算周期项
        for i in 0..77 {
            let nu_tab_cell = NU_TAB[i];
            let c = (nu_tab_cell.0 as f64 * l
                + nu_tab_cell.1 as f64 * l1
                + nu_tab_cell.2 as f64 * f
                + nu_tab_cell.3 as f64 * d
                + nu_tab_cell.4 as f64 * om)
                / RAD;

            // 只计算周期大于zq天的项
            if zq > 0 {
                let q = 36526.0 * 2.0 * PI * RAD
                    / (1717915923.2178 * nu_tab_cell.0 as f64
                        + 129596581.0481 * nu_tab_cell.1 as f64
                        + 1739527262.8478 * nu_tab_cell.2 as f64
                        + 1602961601.2090 * nu_tab_cell.3 as f64
                        + 6962890.5431 * nu_tab_cell.4 as f64);
                if q < zq as f64 {
                    continue;
                }
            }

            dl += (nu_tab_cell.5 as f64 + nu_tab_cell.6 as f64 * t) * c.sin()
                + nu_tab_cell.7 as f64 * c.cos();
            de += (nu_tab_cell.8 as f64 + nu_tab_cell.9 as f64 * t) * c.cos()
                + nu_tab_cell.10 as f64 * c.sin();
        }

        // 返回IAU2000B章动值
        [dl / (10000000.0 * RAD), de / (10000000.0 * RAD)]
    }

    /// 计算赤经章动及赤纬章动
    ///
    /// # Arguments
    ///
    /// * `z` - 天体坐标数组 '\['赤经,赤纬,距离'\]'
    /// * `e` - 黄赤交角
    /// * `dl` - 黄经章动
    /// * `de` - 交角章动
    ///
    /// # Returns
    ///
    /// 返回修正后的坐标数组 '\['赤经,赤纬,距离'\]'
    pub fn cd_nutation(z: &[f64], e: f64, dl: f64, de: f64) -> [f64; 3] {
        // 创建新数组存储结果
        let mut r = [z[0], z[1], z[2]];

        // 计算赤经章动
        r[0] += (e.cos() + e.sin() * z[0].sin() * z[1].tan()) * dl - z[0].cos() * z[1].tan() * de;

        // 计算赤纬章动
        r[1] += e.sin() * z[0].cos() * dl + z[0].sin() * de;

        // 限制赤经在0-2π之间
        r[0] = rad2mrad(r[0]);

        r
    }

    /// 中精度章动计算
    ///
    /// # Arguments
    ///
    /// * `t` - 世纪数
    ///
    /// # Returns
    ///
    /// 返回'\['黄经章动,交角章动'\]'数组
    pub fn nutation2(t: f64) -> [f64; 2] {
        let t2 = t * t;
        let mut dl = 0.0;
        let mut de = 0.0;

        // 遍历章动表B进行计算
        for i in 0..NUT_B.len() {
            let nut_b = NUT_B[i];
            let c = nut_b.0 + nut_b.1 * t + nut_b.2 * t2;

            // 第一项特殊处理
            let a = if i == 0 { -1.742 * t } else { 0.0 };

            dl += (nut_b.3 + a) * c.sin();
            de += nut_b.4 * c.cos();
        }

        // 返回黄经章动和交角章动
        [dl / (100.0 * RAD), de / (100.0 * RAD)]
    }

    /// 只计算黄经章动
    ///
    /// # Arguments
    ///
    /// * `t` - 世纪数
    ///
    /// # Returns
    ///
    /// 返回黄经章动值(弧度)
    pub fn nutation_lon2(t: f64) -> f64 {
        let t2 = t * t;
        let mut dl = 0.0;

        // 遍历章动表B只计算黄经章动
        for i in 0..NUT_B.len() {
            // 第一项特殊处理
            let a = if i == 0 { -1.742 * t } else { 0.0 };

            let nut_b = NUT_B[i];

            // 计算正弦项的参数
            let c = nut_b.0 + nut_b.1 * t + nut_b.2 * t2;

            // 累加章动项
            dl += (nut_b.3 + a) * c.sin();
        }

        // 返回换算成弧度的黄经章动值
        dl / (100.0 * RAD)
    }
}
