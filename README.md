#  **EnergySyst1** 
<br>

## 一、 许可说明

项目脚本基于：
1. Python 环境 @ PSFL 许可 <br>
1. Pyomo 平台 @ BSD 许可<br>
2. Pypower 库 @ BSD 许可 <br>
3. Ipopt 库 @ EPL 许可<br>
<br>
<br>


## 二、功能介绍

通过解析 [Pypower](https://github.com/rwl/PYPOWER/tree/master)中集成的 IEEE cases，得到系统参数，基于 [Pyomo](https://www.pyomo.org/documentation) 框架实现 IEEE 测试系统的经济最优（网损最优）. 电力系统潮流计算方式为[AC-OPF](https://www.youtube.com/watch?v=5MwNL2SuEaI&t=1238s&ab_channel=GurobiOptimization)，最优化问题求解选用内点法求解器 [IPOPT](https://pypi.org/project/ipopt/).<br>
<br>
<br>

## 三、使用说明

项目含多个测试算例，脚本可直接运行.
1. IEEE 9 节点系统<br>
2. IEEE 118 节点系统<br>
3. IEEE 300 节点系统<br>
<br>
<br>

## 四、作者 & 维护
  
Author: HiFery    :bowtie:  <br>
Email：howardyzy@outlook.com <br>
Last update：2025-03-17
<br>
<br>

欢迎交流、讨论。【邮箱】 常年【在线】 :punch: <br> 
欢迎任何形式的。【转载】 :heartpulse: 【复用】 :clap:  <br>
 **但需『提前告知』** 


