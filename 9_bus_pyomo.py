# -*- coding: utf-8 -*-
"""
This project is to test Pyomo-based optimaztion platform using IEEE 9-bus system
Designed by Zhiyuan @ 10th, April, 2023

updated by Zhiyuan @ 9th, Sept., 2024


Using BMF formulation

Important notes:

    # Gij = model.Y[i, j].real
    # Bij = model.Y[i, j].imag

    Branch_[i,j] = Bij*(delta[i]-delta[j]) = model.Y[i, j].imag*(delta[i]-delta[j])

    model.pprint()
    
 添加支路功率限制

"""

import multiprocessing

from pyomo.environ import *
from pyomo.opt import SolverFactory
import math
import cmath
import numpy as np
from pypower.api import case9, case118, ppoption, runpf, loadcase
from numpy import r_, c_, ix_, zeros, pi, ones, exp, argmax
from pypower.makeYbus import makeYbus
from pypower.ext2int import ext2int
from pypower.idx_brch import PF, PT, QF, QT
from pypower.makeBdc import makeBdc
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
path = 'C:\\Users\\Zhiyuan\\Desktop\\Zhiyuan_Testing\\IPOPT\\Ipopt-3.10.1-win64-intel11.1\\Ipopt-3.10.1-win64-intel11.1\\bin\\ipopt.exe'

from pyomo.util.infeasible import log_infeasible_constraints
import logging


# get the structured data from pypower
ppc = loadcase(case9())
ppc = ext2int(ppc)
baseMVA, bus, gen, branch = ppc["baseMVA"], ppc["bus"], ppc["gen"], ppc["branch"]
if ppc["branch"].shape[1] < QT:
    ppc["branch"] = c_[ppc["branch"], zeros((ppc["branch"].shape[0], QT - ppc["branch"].shape[1] + 1))]

Ybus, Yf, Yt = makeYbus(baseMVA, bus, branch)
d = Ybus.todok()
d = dict(d.items()) # very useful command by transefering sparse matrix to dictionary

# get the Zbus and theta matrix (seems not necessary)

# Y_bus = Ybus.todense()
# Zbus = np.reciprocal(Y_bus)
# theta =np.angle(Zbus)
# theta[np.isnan(theta)] = 0

# theta_dic = {}
# for i in range(len(theta)):
#     length, width = theta[i].shape
#     for j in range(width):
#         key = (i, j)
#         theta_value = theta[i,j]
#         theta_dic[key] = theta_value


# get the from_to list in the branch
row_branch, column = branch.shape
line_list = []
for i in range(row_branch):
    from_bus = branch[i,0].astype(int)
    to_bus = branch[i,1].astype(int)
    temp_tuple = (from_bus, to_bus)
    line_list.append(temp_tuple)

# generate bus-gen_matrix for parameters initializing
row_bus, column_bus = bus.shape
row1, column1 = gen.shape

pg_max_list = np.zeros(row_bus)
pg_min_list = np.zeros(row_bus)
qg_max_list = np.zeros(row_bus)
qg_min_list = np.zeros(row_bus)
pg_list = np.zeros(row_bus)
qg_list = np.zeros(row_bus)

for i in range(row1):
    gen_order = gen[i,0]
    pg_max_list[gen_order] = gen[i, 8]/baseMVA
    pg_min_list[gen_order] = gen[i, 9]/baseMVA
    qg_max_list[gen_order] = gen[i, 3]/baseMVA
    qg_min_list[gen_order] = gen[i, 4]/baseMVA
    pg_list[gen_order] = gen[i, 1]/baseMVA
    qg_list[gen_order] = gen[i, 2]/baseMVA

type_list = bus[:,1]

# get the bus-gencost matrix
gencost_mtrx = ppc['gencost']
# gencost_mtrx = np.array(gencost_mtrx)
row_cost, column_cost = gencost_mtrx.shape
gen_order_list = gen[:,0]
bus_gen_mtrx = np.zeros((row_bus,column_cost))

for i in range(len(gen_order_list)):
    bus_gen_mtrx[gen_order_list[i], :] = gencost_mtrx[i, :]

## the gencost normalization process
norm_gencost = 0 # 用来定义是否需要 normalize gen cost

if norm_gencost == 1:
    # column_4 是二次项； column_5 是一次项， column_6 是常数项
    l2_norm_4 = np.linalg.norm(bus_gen_mtrx[:,4])
    l2_norm_5 = np.linalg.norm(bus_gen_mtrx[:,5])
    l2_norm_6 = np.linalg.norm(bus_gen_mtrx[:,6])
    bus_gen_mtrx[:, 4] /= l2_norm_4
    bus_gen_mtrx[:, 5] /= l2_norm_5
    bus_gen_mtrx[:, 6] /= l2_norm_6


# transefer the matrix to dic
cost_dic = {}
for i in range(len(bus_gen_mtrx)):
    for j in range(len(bus_gen_mtrx[i])):
        key = (i, j)
        cost_value = bus_gen_mtrx[i][j]
        cost_dic[key] = cost_value


# Get the carbon Profile
csv_file_path = '9_bus_carbon_intensity.csv'
data_carbon = np.loadtxt(csv_file_path, delimiter=',')

test_carbon_ = data_carbon[1]

carbon_cost = [0]*row_bus

for i in range(len(gen_order_list)):
    carbon_cost[gen_order_list[i]] = test_carbon_[i]

## the carbon_intensity normalization process
norm_carbon = 0 # 用来定义是否需要 normalize gen cost

if norm_carbon == 1:
    l2_norm_carbon = np.linalg.norm(carbon_cost)
    carbon_cost /= l2_norm_carbon


# get the branch list
branch_list = ppc['branch']
MVA_list = branch_list[:,5]



# Define the model
model = ConcreteModel()

# Define the sets
model.buses = Set(initialize = bus[:,0].astype(int))
model.lines = Set(initialize = line_list)
model.cost_dims = Set(initialize = range(0, column_cost))

model.ratio_type = Set(initialize = range(0,4))     # explanation， 0: wind;      1: pv;       2: ES;       4: other

# Define the parameters (change the bus-> dict)
model.PD = Param(model.buses, initialize = bus[:, 2]/baseMVA, default = 0.0) # default meaning?
model.QD = Param(model.buses, initialize = bus[:, 3]/baseMVA, default = 0.0)
model.PG_MAX = Param(model.buses, initialize = pg_max_list, default = 0.0)
model.PG_MIN = Param(model.buses, initialize = pg_min_list, default = 0.0)
model.QG_MAX = Param(model.buses, initialize = qg_max_list, default = 0.0)
model.QG_MIN = Param(model.buses, initialize = qg_min_list, default = 0.0)
model.V_MAX = Param(model.buses, initialize = bus[:, 11], default = 0.0)
model.V_MIN = Param(model.buses, initialize = bus[:, 12], default = 0.0)
model.Y = Param(model.buses*model.buses, initialize = d, default = 0.0)
# model.theta = Param(model.buses*model.buses, initialize = theta_dic, default=0.0)
model.C = Param(model.buses*model.cost_dims, initialize = cost_dic, default = 0.0)
model.Carbon = Param(model.buses, initialize = carbon_cost, default = 0.0)


# print(model)

# Define the variables
model.PG = Var(model.buses, domain = Reals)
model.QG = Var(model.buses, domain = Reals)
model.V = Var(model.buses, initialize=1.0, bounds=(0.95, 1.05))
model.delta = Var(model.buses, initialize=0, bounds=(-math.pi/8, math.pi/8))
model.ratio = Var(model.ratio_type, initialize = 0.02, bounds = (0.01, 0.99))

for i in range(row_bus):
    if type_list[i] == 3 or type_list[i] == 2:
        model.PG[i].lb = pg_min_list[i]
        model.PG[i].ub = pg_max_list[i]
        model.QG[i].lb = qg_min_list[i]
        model.QG[i].ub = qg_max_list[i]
        if type_list[i] == 3:
            model.delta[i].fix(0)
            model.V[i].fix(1.0)

    else:
        model.PG[i].fix(0)
        model.QG[i].fix(0)


# set stale and fixed
# model.V[:].stale = True
# model.delta[:].stale = True
#
# # # Define the Y variable (admittance)
# # model.Y = Var(model.LINES, within=Reals)

#
# # Define the objective function
# model.obj = Objective(expr=sum(model.PG[i]**3+model.PG[i]**2*model.C[i,4]+ model.PG[i]*model.C[i,5]+model.C[i,6]+model.C[i,1] for i in model.buses), sense=minimize)



# Define the power balance constraints
model.power_balance_constraints = ConstraintList()
for i in model.buses:
    # Gij = model.Y[i, j].real
    # Bij = model.Y[i, j].imag
    real_power_inj = sum(model.V[i] * model.V[j] * (model.Y[i, j].real * cos(model.delta[i] - model.delta[j]) + (model.Y[i, j].imag) * sin(model.delta[i] - model.delta[j])) for j in model.buses)
    reactive_power_inj = sum(model.V[i] * model.V[j] * ((model.Y[i, j].real) * sin(model.delta[i] - model.delta[j]) - (model.Y[i, j].imag) * cos(model.delta[i] - model.delta[j])) for j in model.buses)


    if type_list[i] == 3 or type_list[i] == 2:
        model.power_balance_constraints.add(expr = model.PG[i] - model.PD[i] == real_power_inj)
        model.power_balance_constraints.add(expr = model.QG[i] - model.QD[i] ==  reactive_power_inj)

    if type_list[i] == 1:
        model.power_balance_constraints.add(expr= - model.PD[i] == real_power_inj)
        model.power_balance_constraints.add(expr= - model.QD[i] == reactive_power_inj)



# power supply-demand
model.power_supply_demand_balance = ConstraintList()

model.power_supply_demand_balance.add(expr = sum(model.PG[i] for i in model.buses) - sum(model.PD[i] for i in model.buses) <= 0.05) #     存在网损 0.05 的网损
# model.power_supply_demand_balance.add(expr = sum(model.PG[i] for i in model.buses) - sum(model.PD[i] for i in model.buses) >= 0)
# model.power_supply_demand_balance.add(expr = sum(model.QG[i] for i in model.buses) - sum(model.QD[i] for i in model.buses) == 0)

model.generator_limits = ConstraintList()
for i in model.buses:
    if type_list[i] == 3 or type_list[i] == 2:

        model.generator_limits.add(expr = model.PG_MIN[i] - model.PG[i] <= 0)
        model.generator_limits.add(expr = model.PG[i] - model.PG_MAX[i] <= 0)
        model.generator_limits.add(expr = model.QG_MIN[i] - model.QG[i] <= 0)
        model.generator_limits.add(expr = model.QG[i] - model.QG_MAX[i] <= 0)

# power ratio constraint
model.ratio_limits = Constraint(expr = sum(model.ratio[i] for i in model.ratio_type) - 1 == 0)

# branch limits
model.branch_limits = ConstraintList()
for i in range(len(MVA_list)):
    line_limit = (MVA_list[i]/baseMVA)**2
    f_bus = branch_list[i][0]
    t_bus = branch_list[i][1]
    Z_f_t = -1 / model.Y[f_bus, t_bus]

    R_f_t = Z_f_t.real
    X_f_t = Z_f_t.imag

    U_f_real = model.V[f_bus] * cos(model.delta[f_bus])
    U_f_imag = model.V[f_bus] * sin(model.delta[f_bus])
    U_t_real = model.V[t_bus] * cos(model.delta[t_bus])
    U_t_imag = model.V[t_bus] * sin(model.delta[t_bus])

    dif_U_real = U_f_real - U_t_real
    dif_U_imag = U_f_imag - U_t_imag

    I_f_t_real = (dif_U_real * R_f_t + X_f_t * dif_U_imag)/(R_f_t**2 + X_f_t**2)
    I_f_t_imag = (dif_U_imag * R_f_t - X_f_t * dif_U_real)/(R_f_t**2 + X_f_t**2)

    ### Sf_t = U_f*I_f_t.conjugate() 由此，复功率的实部虚部如下：
    P_f_t = U_f_real * I_f_t_real + U_f_imag * I_f_t_imag
    Q_f_t = U_f_imag * I_f_t_real - U_f_real * I_f_t_imag

    S_f_t_2 = P_f_t**2+Q_f_t**2
    model.branch_limits.add(expr = S_f_t_2 - line_limit <= 0 )
    model.branch_limits.add(expr = S_f_t_2 + line_limit >= 0 )


# # power balance ratio on the bus-2
# model.ratio_balance =


# 网损最优
# model.obj = Objective(expr = (sum(model.PG[i] for i in model.buses) - sum(model.PD[i] for i in model.buses))**2, sense=minimize)

# 潮流经济最优
# model.obj = Objective(expr = sum((model.PG[i]*baseMVA)**2*model.C[i,4]+ (model.PG[i]*baseMVA)*model.C[i,5]+model.C[i,6]+model.C[i,1] for i in model.buses), sense=minimize)

# 碳排最优 + 经济最优
model.obj_cost = Objective(expr =   sum(  (model.PG[i]*baseMVA)**2*model.C[i,4]+ (model.PG[i]*baseMVA)*model.C[i,5]+ model.C[i,6]for i in model.buses), sense = minimize)
model.obj_carbon = Objective(expr = sum( model.PG[i]*baseMVA*model.Carbon[i]  for i in model.buses), sense = minimize)

model.obj_carbon.deactivate()
solver = SolverFactory('ipopt', executable=path)
results = solver.solve(model, tee = True)

model.obj_carbon.activate()
model.obj_cost.deactivate()
solver = SolverFactory('ipopt', executable=path)
results = solver.solve(model, tee = True)



### model print after solving
print(model.obj_cost())
print(model.obj_carbon())
# PG = [model.PG[i].value  for i in range(row_bus)]
# QG = [model.QG[i].value  for i in range(row_bus)]
# V =  [model.V[i].value  for i in range(row_bus)]
# delta =  [model.delta[i].value  for i in range(row_bus)]
# PD = sum(model.PD[i] for i in model.buses)
#
# print('PG', PG)
# print('QG', QG)
# print('Voltage', V)
# print('Delta', delta)
# print('loss', sum(PG) - PD )


# check the power flow of each branch
# quick test for line 3-4 and 4-3
# test
# branch_list = ppc['branch']
# MVA_list = branch_list[:,5]
#
# for i in range(len(MVA_list)):
#     f_bus = branch_list[i][0]
#     t_bus = branch_list[i][1]
#     P_f_t = model.V[f_bus] * model.V[t_bus] * (model.Y[f_bus, t_bus].real * cos(model.delta[f_bus] - model.delta[t_bus]) + model.Y[f_bus, t_bus].imag * sin(model.delta[f_bus] - model.delta[t_bus]))
#     Q_f_t = model.V[f_bus] * model.V[t_bus] * (model.Y[f_bus, t_bus].real * sin(model.delta[f_bus] - model.delta[t_bus]) - model.Y[f_bus, t_bus].imag * cos(model.delta[f_bus] - model.delta[t_bus]))
#     S_f_t = sqrt(P_f_t**2+Q_f_t**2)
#     print(P_f_t(),Q_f_t(), S_f_t())
# #
# for i in range(len(MVA_list)):
#     f_bus = branch_list[i][0]
#     t_bus = branch_list[i][1]
#     Z_f_t = -1 / model.Y[f_bus, t_bus]
#     U_f = model.V[f_bus]() * cmath.cos(model.delta[f_bus]()) + model.V[f_bus]() * cmath.sin(model.delta[f_bus]()) * 1j
#     U_t = model.V[t_bus]() * cmath.cos(model.delta[t_bus]()) + model.V[t_bus]() * cmath.sin(model.delta[t_bus]()) * 1j
#     I_f_t = (U_f - U_t) / Z_f_t
#     S_f_t = U_f * I_f_t.conjugate()
#     print(value(I_f_t), value(S_f_t))
#
# for i in range(len(MVA_list)):
#     f_bus = branch_list[i][0]
#     t_bus = branch_list[i][1]
#     Z_f_t = -1 / model.Y[f_bus, t_bus]
#
#     R_f_t = Z_f_t.real
#     X_f_t = Z_f_t.imag
#
#     U_f_real = model.V[f_bus] * cos(model.delta[f_bus])
#     U_f_imag = model.V[f_bus] * sin(model.delta[f_bus])
#     U_t_real = model.V[t_bus] * cos(model.delta[t_bus])
#     U_t_imag = model.V[t_bus] * sin(model.delta[t_bus])
#
#     dif_U_real = U_f_real - U_t_real
#     dif_U_imag = U_f_imag - U_t_imag
#
#     I_f_t_real = (dif_U_real * R_f_t + X_f_t * dif_U_imag)/(R_f_t**2 + X_f_t**2)
#     I_f_t_imag = (dif_U_imag * R_f_t - X_f_t * dif_U_real)/(R_f_t**2 + X_f_t**2)
#
#     P_f_t = U_f_real * I_f_t_real + U_f_imag * I_f_t_imag
#     Q_f_t = U_f_imag * I_f_t_real - U_f_real * I_f_t_imag
#
#     S_f_t = sqrt(P_f_t**2+Q_f_t**2)
#
#
#
#     print(value(P_f_t),value(Q_f_t))





