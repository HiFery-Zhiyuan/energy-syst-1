# -*- coding: utf-8 -*-
"""
This is the main script to generate the cabon and OPF calculations

updated by Zhiyuan @ 29th, March, 2024


Using BMF formulation

Important notes:

    # Gij = model.Y[i, j].real
    # Bij = model.Y[i, j].imag

    Branch_[i,j] = Bij*(delta[i]-delta[j]) = model.Y[i, j].imag*(delta[i]-delta[j])

    model.pprint()


THREE multiple-objective solutions are incorporated:

1) BASIC compromise programming
2) epsilon-constraint programming
3) Weighted minmax programming


"""

import multiprocessing

from pyomo.environ import *
from pyomo.opt import SolverFactory
import math
import cmath
import numpy as np
from pypower.api import case9,case57, case118, ppoption, runpf, loadcase
from numpy import r_, c_, ix_, zeros, pi, ones, exp, argmax
from pypower.makeYbus import makeYbus
from pypower.ext2int import ext2int
from pypower.idx_brch import PF, PT, QF, QT
from pypower.makeBdc import makeBdc
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from min_max_9_bus import min_max_9bus

if __name__ == "__main__":

    path = 'C:\\Users\\Zhiyuan\\Desktop\\Zhiyuan_Testing\\IPOPT\\Ipopt-3.10.1-win64-intel11.1\\Ipopt-3.10.1-win64-intel11.1\\bin\\ipopt.exe'
    # get the structured data from pypower
    ppc = loadcase(case9())
    ppc = ext2int(ppc)
    baseMVA, bus, gen, branch = ppc["baseMVA"], ppc["bus"], ppc["gen"], ppc["branch"]
    if ppc["branch"].shape[1] < QT:
        ppc["branch"] = c_[ppc["branch"], zeros((ppc["branch"].shape[0], QT - ppc["branch"].shape[1] + 1))]

    Ybus, Yf, Yt = makeYbus(baseMVA, bus, branch)
    d = Ybus.todok()
    d = dict(d.items())  # very useful command by transefering sparse matrix to dictionary

    # get the from_to list in the branch
    row_branch, column = branch.shape
    line_list = []
    for i in range(row_branch):
        from_bus = branch[i, 0].astype(int)
        to_bus = branch[i, 1].astype(int)
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
        gen_order = gen[i, 0]
        pg_max_list[gen_order] = gen[i, 8] / baseMVA
        pg_min_list[gen_order] = gen[i, 9] / baseMVA
        qg_max_list[gen_order] = gen[i, 3] / baseMVA
        qg_min_list[gen_order] = gen[i, 4] / baseMVA
        pg_list[gen_order] = gen[i, 1] / baseMVA
        qg_list[gen_order] = gen[i, 2] / baseMVA

    type_list = bus[:, 1]

    # get the bus-gencost matrix
    gencost_mtrx = ppc['gencost']
    # gencost_mtrx = np.array(gencost_mtrx)
    row_cost, column_cost = gencost_mtrx.shape
    gen_order_list = gen[:, 0]
    bus_gen_mtrx = np.zeros((row_bus, column_cost))

    for i in range(len(gen_order_list)):
        bus_gen_mtrx[gen_order_list[i], :] = gencost_mtrx[i, :]

    ## the gencost normalization process
    norm_gencost = 0  # 用来定义是否需要 normalize gen cost

    if norm_gencost == 1:
        # column_4 是二次项； column_5 是一次项， column_6 是常数项
        l2_norm_4 = np.linalg.norm(bus_gen_mtrx[:, 4])
        l2_norm_5 = np.linalg.norm(bus_gen_mtrx[:, 5])
        l2_norm_6 = np.linalg.norm(bus_gen_mtrx[:, 6])
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

    carbon_cost = [0] * row_bus

    for i in range(len(gen_order_list)):
        carbon_cost[gen_order_list[i]] = test_carbon_[i]

    ## the carbon_intensity normalization process
    norm_carbon = 0  # 用来定义是否需要 normalize gen cost

    if norm_carbon == 1:
        l2_norm_carbon = np.linalg.norm(carbon_cost)
        carbon_cost /= l2_norm_carbon

    # get the branch list
    branch_list = ppc['branch']
    MVA_list = branch_list[:, 5]

    # Define the model
    main_model = ConcreteModel()

    # Define the sets
    main_model.buses = Set(initialize=bus[:, 0].astype(int))
    main_model.lines = Set(initialize=line_list)
    main_model.cost_dims = Set(initialize=range(0, column_cost))

    # Define the parameters (change the bus-> dict)
    main_model.PD = Param(main_model.buses, initialize=bus[:, 2] / baseMVA, default=0.0)  # default meaning?
    main_model.QD = Param(main_model.buses, initialize=bus[:, 3] / baseMVA, default=0.0)
    main_model.PG_MAX = Param(main_model.buses, initialize=pg_max_list, default=0.0)
    main_model.PG_MIN = Param(main_model.buses, initialize=pg_min_list, default=0.0)
    main_model.QG_MAX = Param(main_model.buses, initialize=qg_max_list, default=0.0)
    main_model.QG_MIN = Param(main_model.buses, initialize=qg_min_list, default=0.0)
    main_model.V_MAX = Param(main_model.buses, initialize=bus[:, 11], default=0.0)
    main_model.V_MIN = Param(main_model.buses, initialize=bus[:, 12], default=0.0)
    main_model.Y = Param(main_model.buses * main_model.buses, initialize=d, default=0.0)
    # model.theta = Param(model.buses*model.buses, initialize = theta_dic, default=0.0)
    main_model.C = Param(main_model.buses * main_model.cost_dims, initialize=cost_dic, default=0.0)
    main_model.Carbon = Param(main_model.buses, initialize=carbon_cost, default=0.0)

    # print(model)

    # Define the variables
    main_model.PG = Var(main_model.buses, domain=Reals)
    main_model.QG = Var(main_model.buses, domain=Reals)
    main_model.V = Var(main_model.buses, initialize=1.0, bounds=(0.95, 1.05))
    main_model.delta = Var(main_model.buses, initialize=0, bounds=(-math.pi / 8, math.pi / 8))

    for i in range(row_bus):
        if type_list[i] == 3 or type_list[i] == 2:
            main_model.PG[i].lb = pg_min_list[i]
            main_model.PG[i].ub = pg_max_list[i]
            main_model.QG[i].lb = qg_min_list[i]
            main_model.QG[i].ub = qg_max_list[i]
            if type_list[i] == 3:
                main_model.delta[i].fix(0)
                main_model.V[i].fix(1.0)

        else:
            main_model.PG[i].fix(0)
            main_model.QG[i].fix(0)

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
    main_model.power_balance_constraints = ConstraintList()
    for i in main_model.buses:
        # Gij = model.Y[i, j].real
        # Bij = model.Y[i, j].imag
        real_power_inj = sum(main_model.V[i] * main_model.V[j] * (
                main_model.Y[i, j].real * cos(main_model.delta[i] - main_model.delta[j]) + (main_model.Y[i, j].imag) * sin(
            main_model.delta[i] - main_model.delta[j])) for j in main_model.buses)
        reactive_power_inj = sum(main_model.V[i] * main_model.V[j] * (
                (main_model.Y[i, j].real) * sin(main_model.delta[i] - main_model.delta[j]) - (main_model.Y[i, j].imag) * cos(
            main_model.delta[i] - main_model.delta[j])) for j in main_model.buses)

        if type_list[i] == 3 or type_list[i] == 2:
            main_model.power_balance_constraints.add(expr=main_model.PG[i] - main_model.PD[i] == real_power_inj)
            main_model.power_balance_constraints.add(expr=main_model.QG[i] - main_model.QD[i] == reactive_power_inj)

        if type_list[i] == 1:
            main_model.power_balance_constraints.add(expr=- main_model.PD[i] == real_power_inj)
            main_model.power_balance_constraints.add(expr=- main_model.QD[i] == reactive_power_inj)

    # power supply-demand
    main_model.power_supply_demand_balance = ConstraintList()

    main_model.power_supply_demand_balance.add(
        expr=sum(main_model.PG[i] for i in main_model.buses) - sum(main_model.PD[i] for i in main_model.buses) <= 0.03)  # 存在网损 0.05 的网损
    # model.power_supply_demand_balance.add(expr = sum(model.PG[i] for i in model.buses) - sum(model.PD[i] for i in model.buses) >= 0)
    # model.power_supply_demand_balance.add(expr = sum(model.QG[i] for i in model.buses) - sum(model.QD[i] for i in model.buses) == 0)

    main_model.generator_limits = ConstraintList()
    for i in main_model.buses:
        if type_list[i] == 3 or type_list[i] == 2:
            main_model.generator_limits.add(expr=main_model.PG_MIN[i] - main_model.PG[i] <= 0)
            main_model.generator_limits.add(expr=main_model.PG[i] - main_model.PG_MAX[i] <= 0)
            main_model.generator_limits.add(expr=main_model.QG_MIN[i] - main_model.QG[i] <= 0)
            main_model.generator_limits.add(expr=main_model.QG[i] - main_model.QG_MAX[i] <= 0)

   # get the min_max value for diverse object

    min_carbon, min_cost = min_max_9bus(path)
    # branch limits
    main_model.branch_limits = ConstraintList()
    for i in range(len(MVA_list)):
        line_limit = (MVA_list[i] / baseMVA) ** 2
        f_bus = branch_list[i][0]
        t_bus = branch_list[i][1]
        Z_f_t = -1 / main_model.Y[f_bus, t_bus]

        R_f_t = Z_f_t.real
        X_f_t = Z_f_t.imag

        U_f_real = main_model.V[f_bus] * cos(main_model.delta[f_bus])
        U_f_imag = main_model.V[f_bus] * sin(main_model.delta[f_bus])
        U_t_real = main_model.V[t_bus] * cos(main_model.delta[t_bus])
        U_t_imag = main_model.V[t_bus] * sin(main_model.delta[t_bus])

        dif_U_real = U_f_real - U_t_real
        dif_U_imag = U_f_imag - U_t_imag

        I_f_t_real = (dif_U_real * R_f_t + X_f_t * dif_U_imag) / (R_f_t ** 2 + X_f_t ** 2)
        I_f_t_imag = (dif_U_imag * R_f_t - X_f_t * dif_U_real) / (R_f_t ** 2 + X_f_t ** 2)

        ### Sf_t = U_f*I_f_t.conjugate() 由此，复功率的实部虚部如下：
        P_f_t = U_f_real * I_f_t_real + U_f_imag * I_f_t_imag
        Q_f_t = U_f_imag * I_f_t_real - U_f_real * I_f_t_imag

        S_f_t_2 = P_f_t ** 2 + Q_f_t ** 2
        main_model.branch_limits.add(expr=S_f_t_2 - line_limit <= 0)
        main_model.branch_limits.add(expr=S_f_t_2 + line_limit >= 0)

    # 碳排最优 + 经济最优 多目标建模和求解

    # （1）Basic Compromise programming (BCP) with diverse settings of weights w1*|f1 - f1_min| + w2*|f2-f2_min|
    BCP = 0
    if BCP == 1:
        main_model.mutli_obj = Objective(expr=
                                     0.5*(sum((main_model.PG[i]*baseMVA) ** 2 * main_model.C[i, 4] + (main_model.PG[i]*baseMVA) * main_model.C[i, 5] + main_model.C[i, 6] for i in main_model.buses) / min_cost)  +
                                     0.5*(sum(main_model.PG[i]*baseMVA*main_model.Carbon[i] for i in main_model.buses) / min_carbon),
                               sense=minimize)

    # （2）Epsilon-Constraint (EC)
    EC = 1
    if EC == 1:
        # add Epsilon constraint
        epsilon = 1.02
        main_model.epsilon_constraint = ConstraintList()
        main_model.epsilon_constraint.add(
            expr = sum(main_model.PG[i]*baseMVA*main_model.Carbon[i] for i in main_model.buses) <= epsilon * min_carbon
        )

        # add objective function
        main_model.obj_cost = Objective(expr=
                                        sum((main_model.PG[i]*baseMVA) ** 2 * main_model.C[i, 4] + (main_model.PG[i]*baseMVA) * main_model.C[i, 5] + main_model.C[i, 6] for i in main_model.buses),
                                        sense=minimize
        )

    # (3) weighted min max (MM) formulation
    MM = 0
    if MM == 1:
        # add Variables
        main_model.number_obj = Set(initialize=range(1,3))
        main_model.y = Var(main_model.number_obj, domain = NonNegativeReals)
        main_model.z = Var( domain = NonNegativeReals)

        # add constraint of two objectives
        main_model.objective_constraint_01 = Constraint(expr=
                                                        main_model.y[1] == sum((main_model.PG[i]*baseMVA) ** 2 * main_model.C[i, 4] + (main_model.PG[i]*baseMVA) * main_model.C[i, 5] + main_model.C[i, 6] for i in main_model.buses) / min_cost
        )
        main_model.objective_constraint_02 = Constraint(expr=
                                                        main_model.y[2] == sum(main_model.PG[i]*baseMVA*main_model.Carbon[i] for i in main_model.buses) / min_carbon
        )

        main_model.objective_constraint_03 = Constraint(expr = 0.5*main_model.y[1] <= main_model.z)
        main_model.objective_constraint_04 = Constraint(expr = 0.5*main_model.y[2] <= main_model.z)

        # add objective
        main_model.mutli_obj = Objective(expr=main_model.z, sense=minimize)



    solver = SolverFactory('ipopt', executable=path)
    results = solver.solve(main_model, tee=True)

    # verify purpose
    carbon_ = sum(main_model.PG[i]*baseMVA*main_model.Carbon[i] for i in main_model.buses)
    cost_ = sum((main_model.PG[i]*baseMVA) ** 2 * main_model.C[i, 4] + (main_model.PG[i]*baseMVA) * main_model.C[i, 5] + main_model.C[i, 6] for i in main_model.buses)


    print('The ORIGINAL carbon optimal result: ', min_carbon)
    print('The ORIGINAL cost optimal result: ', min_cost)
    print('***************** AFTER multi-objective formulation*****************')
    print('carbon-related related optimal result: ', carbon_())
    print('cost-related optimal result: ', cost_())

    PG = [main_model.PG[i].value  for i in range(row_bus)]
    QG = [main_model.QG[i].value  for i in range(row_bus)]
    V =  [main_model.V[i].value  for i in range(row_bus)]
    delta =  [main_model.delta[i].value  for i in range(row_bus)]
    PD = sum(main_model.PD[i] for i in main_model.buses)

    print('PG', PG)
    print('QG', QG)
    print('Voltage', V)
    print('Delta', delta)
    print('loss', (sum(PG) - PD)*baseMVA )








