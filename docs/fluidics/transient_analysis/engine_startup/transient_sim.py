import numpy as np
from turborocket.transient.start_up import (
    GasGenerator,
    MainEngine,
    Pump,
    Turbine,
    LiquidValve,
    Cavity,
)
from turborocket.tools.loaders import TransientLoader
from turborocket.fluids.fluids import IncompressibleFluid
import pathlib
import pandas as pd
import matplotlib.pyplot as plt

######################################### LOADING YAML FILES ##########################################
target_file = "config.yaml"

local_dir = pathlib.Path(__file__).parent.resolve().joinpath(target_file)

loader = TransientLoader(config_file=local_dir)


######################################### Gas Generator Sizing ########################################

gg_inj_ox = IncompressibleFluid(
    rho=loader._propellants["Oxidiser"]._density,
    P=loader._engines["GAS_GENERATOR"]._P_inj_ox,
)
gg_inj_fu = IncompressibleFluid(
    rho=loader._propellants["Fuel"]._density,
    P=loader._engines["GAS_GENERATOR"]._P_inj_fu,
)

GG = GasGenerator(
    Ox=loader._propellants["Oxidiser"]._name,
    Fu=loader._propellants["Fuel"]._name,
    Pcc=loader._engines["GAS_GENERATOR"]._P_cc,
    MR=loader._engines["GAS_GENERATOR"]._MR,
)

GG.injector_cond(
    ox_in=gg_inj_ox,
    fu_in=gg_inj_fu,
    cdo=loader._engines["GAS_GENERATOR"]._Cd_o,
    cdf=loader._engines["GAS_GENERATOR"]._Cd_f,
)
_ = GG.size_system(
    m_dot=loader._engines["GAS_GENERATOR"]._m_dot,
    eta_c=loader._engines["GAS_GENERATOR"]._eta_c,
)
_ = GG.get_geometry()

GG.set_l_star(L_star=loader._engines["GAS_GENERATOR"]._L_star)
GG.set_pcc_transient(P_cc_transient=loader._engines["GAS_GENERATOR"]._P_init)

############################################ Engine Sizing ############################################

me_inj_ox = IncompressibleFluid(
    rho=loader._propellants["Oxidiser"]._density,
    P=loader._engines["MAIN_ENGINE"]._P_inj_ox,
)
me_inj_fu = IncompressibleFluid(
    rho=loader._propellants["Fuel"]._density,
    P=loader._engines["MAIN_ENGINE"]._P_inj_fu,
)

ME = MainEngine(
    Ox=loader._propellants["Oxidiser"]._name,
    Fu=loader._propellants["Fuel"]._name,
    Pcc=loader._engines["MAIN_ENGINE"]._P_cc,
    MR=loader._engines["MAIN_ENGINE"]._MR,
)

ME.injector_cond(
    ox_in=me_inj_ox,
    fu_in=me_inj_fu,
    cdo=loader._engines["MAIN_ENGINE"]._Cd_o,
    cdf=loader._engines["MAIN_ENGINE"]._Cd_f,
)
_ = ME.size_system(
    m_dot=loader._engines["MAIN_ENGINE"]._m_dot,
    eta_c=loader._engines["MAIN_ENGINE"]._eta_c,
)
_ = ME.get_geometry()

ME.set_l_star(L_star=loader._engines["MAIN_ENGINE"]._L_star)
ME.set_pcc_transient(P_cc_transient=loader._engines["MAIN_ENGINE"]._P_init)


######################################### Valve Instantiation #########################################

GG_MOV = LiquidValve(
    cda=loader._valves["GG_MOV"]._cda,
    tau=loader._valves["GG_MOV"]._tau,
    s_pos_init=loader._valves["GG_MOV"]._s_pos_init,
)

GG_MFV = LiquidValve(
    cda=loader._valves["GG_MFV"]._cda,
    tau=loader._valves["GG_MFV"]._tau,
    s_pos_init=loader._valves["GG_MFV"]._s_pos_init,
)

ME_MOV = LiquidValve(
    cda=loader._valves["ME_MOV"]._cda,
    tau=loader._valves["ME_MOV"]._tau,
    s_pos_init=loader._valves["ME_MOV"]._s_pos_init,
)

ME_MFV = LiquidValve(
    cda=loader._valves["ME_MFV"]._cda,
    tau=loader._valves["ME_MFV"]._tau,
    s_pos_init=loader._valves["ME_MFV"]._s_pos_init,
)

######################################### Cavity Instantiation #########################################

ox_gg_injector_fluid = IncompressibleFluid(
    rho=loader._propellants["Oxidiser"]._density,
    P=loader._engines["GAS_GENERATOR"]._P_inj_ox_init,
    B=loader._propellants["Oxidiser"]._B,
)
fu_gg_injector_fluid = IncompressibleFluid(
    rho=loader._propellants["Fuel"]._density,
    P=loader._engines["GAS_GENERATOR"]._P_inj_fu_init,
    B=loader._propellants["Fuel"]._B,
)

inj_ox_gg = Cavity(
    fluid=ox_gg_injector_fluid, V=loader._engines["GAS_GENERATOR"]._V_ox_inj
)
inj_fu_gg = Cavity(
    fluid=fu_gg_injector_fluid, V=loader._engines["GAS_GENERATOR"]._V_fu_inj
)

# Main Engine
ox_me_injector_fluid = IncompressibleFluid(
    rho=loader._propellants["Oxidiser"]._density,
    P=loader._engines["MAIN_ENGINE"]._P_inj_ox_init,
    B=loader._propellants["Oxidiser"]._B,
)
fu_me_injector_fluid = IncompressibleFluid(
    rho=loader._propellants["Fuel"]._density,
    P=loader._engines["MAIN_ENGINE"]._P_inj_fu_init,
    B=loader._propellants["Fuel"]._B,
)

inj_ox_me = Cavity(
    fluid=ox_me_injector_fluid, V=loader._engines["MAIN_ENGINE"]._V_ox_inj
)

inj_fu_me = Cavity(
    fluid=fu_me_injector_fluid, V=loader._engines["MAIN_ENGINE"]._V_fu_inj
)


#################################### Pump and Turbine Instantiaton #####################################

pump = Pump(
    D=loader._turbo_pumps["MERMAID_MAN"]["Pumps"]["FUEL_PUMP"]._D_nom,
    Q_nom=loader._turbo_pumps["MERMAID_MAN"]["Pumps"]["FUEL_PUMP"]._Q_nom,
    eta_nom=loader._turbo_pumps["MERMAID_MAN"]["Pumps"]["FUEL_PUMP"]._eta_nom,
    N_nom=loader._turbo_pumps["MERMAID_MAN"]["N_nom"],
)

turbine = Turbine(
    delta_b=loader._turbo_pumps["MERMAID_MAN"]["Turbines"]["FUEL_TURBINE"]._delta_b,
    a_rat=loader._turbo_pumps["MERMAID_MAN"]["Turbines"]["FUEL_TURBINE"]._a_rat,
    D_m=loader._turbo_pumps["MERMAID_MAN"]["Turbines"]["FUEL_TURBINE"]._D_m,
    eta=loader._turbo_pumps["MERMAID_MAN"]["Turbines"]["FUEL_TURBINE"]._eta,
)


############################################ Simulation Setup ##########################################

# Start and stop time, and relaxation parameter
t_init = loader._t_start
t_stop = loader._t_stop

alpha = loader._alpa
dp = loader._max_dp
dt_fix = loader._dt_fix
dt_init = loader._dt_init

# Initialise our first time step
dt = dt_init

# Defining the pump inlet conditions and oxidiser delivery conditions
fu_deliver = IncompressibleFluid(
    rho=loader._propellants["Fuel"]._density, P=loader._propellants["Fuel"]._P_deliver
)

ox_deliver = IncompressibleFluid(
    rho=loader._propellants["Oxidiser"]._density,
    P=loader._propellants["Oxidiser"]._P_deliver,
)

# We initialise our system parameters
N = loader._turbo_pumps["MERMAID_MAN"]["N_init"]
I = loader._turbo_pumps["MERMAID_MAN"]["I"]

m_dot_gg_inj_fu = 0  # TODO: Fix this magic Number
m_dot_gg_inj_ox = 0  # TODO: Fix this magic Number

m_dot_me_inj_fu = 0  # TODO: Fix this magic Number
m_dot_me_inj_ox = 0  # TODO: Fix this magic Number

m_dot_gg_mov = 0  # TODO: Fix this magic Number
m_dot_gg_mfv = 0  # TODO: Fix this magic Number

m_dot_me_mov = 0  # TODO: Fix this magic Number
m_dot_me_mfv = 0  # TODO: Fix this magic Number

tau_turbine = 0  # TODO: Fix this magic Number
tau_pump = 0  # TODO: Fix this magic Number

# Initialising our arrays for data storage
t_array = [t_init]

# ---------- PUMP ----------
N_pump_array = [N]
p_pump_array = [pump.get_exit_condition(inlet=fu_deliver, N=N, m_dot=0).get_pressure()]
eta_pump_array = [pump.get_eta(Q=0, N=N)]
tau_pump_array = [tau_pump]

# ---------- Turbine --------
tau_turbine_array = [tau_turbine]

# ---------- Gas Generator --------
T_gg_comb_array = [295]
Pcc_gg_array = [loader._engines["GAS_GENERATOR"]._P_init]
MR_gg_array = [0]

dp_dt_gg_array = [0]

m_dot_gg_inj_fu_array = [m_dot_gg_inj_fu]
m_dot_gg_inj_ox_array = [m_dot_gg_inj_ox]

p_gg_inj_ox_array = [inj_ox_gg.get_fluid().get_pressure()]
p_gg_inj_fu_array = [inj_fu_gg.get_fluid().get_pressure()]

# ---------- Main Engine -----------
T_me_comb_array = [295]
Pcc_me_array = [loader._engines["MAIN_ENGINE"]._P_init]
MR_me_array = [0]

dp_dt_me_array = [0]

m_dot_me_inj_fu_array = [m_dot_me_inj_fu]
m_dot_me_inj_ox_array = [m_dot_me_inj_ox]

p_me_inj_ox_array = [inj_ox_me.get_fluid().get_pressure()]
p_me_inj_fu_array = [inj_fu_me.get_fluid().get_pressure()]

# ---------- Valves --------
m_dot_gg_mov_array = [m_dot_gg_mov]
m_dot_gg_mfv_array = [m_dot_gg_mfv]

m_dot_me_mov_array = [m_dot_me_mov]
m_dot_me_mfv_array = [m_dot_me_mfv]

# ----------- Ambeint Conditions --------------
P_exit = 1e5  # TODO: Fix this Magic Number


# ----------- Sequences -----------
seq_gg_mov = loader._sequences["GG_MOV"]
seq_gg_mfv = loader._sequences["GG_MFV"]

seq_me_mov = loader._sequences["ME_MOV"]
seq_me_mfv = loader._sequences["ME_MFV"]


############################################ System Integration ##########################################

while t_array[-1] < t_stop:

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Valve Position Updating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    GG_MOV.actuate(position=seq_gg_mov(t_array[-1]))
    GG_MFV.actuate(position=seq_gg_mfv(t_array[-1]))

    ME_MOV.actuate(position=seq_me_mov(t_array[-1]))
    ME_MFV.actuate(position=seq_me_mfv(t_array[-1]))

    GG_MOV.update_pos(dt=dt)
    GG_MFV.update_pos(dt=dt)

    ME_MOV.update_pos(dt=dt)
    ME_MFV.update_pos(dt=dt)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Updating Pump Conditions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    pump_exit = pump.get_exit_condition(
        inlet=fu_deliver, N=N, m_dot=(m_dot_gg_mfv + m_dot_me_mfv)
    )
    pump_t = pump.get_torque(inlet=fu_deliver, N=N, m_dot=(m_dot_gg_mfv + m_dot_me_mfv))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Evaluating for Valve M_dot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    m_dot_gg_mov = (
        alpha * GG_MOV.get_mdot(upstr=ox_deliver, downstr=inj_ox_gg.get_fluid(), dt=dt)
    ) + (1 - alpha) * m_dot_gg_mov_array[-1]

    m_dot_gg_mfv = (
        alpha * GG_MFV.get_mdot(upstr=pump_exit, downstr=inj_fu_gg.get_fluid(), dt=dt)
        + (1 - alpha) * m_dot_gg_mfv_array[-1]
    )

    m_dot_me_mov = (
        alpha * ME_MOV.get_mdot(upstr=ox_deliver, downstr=inj_ox_me.get_fluid(), dt=dt)
    ) + (1 - alpha) * m_dot_me_mov_array[-1]

    m_dot_me_mfv = (
        alpha * ME_MFV.get_mdot(upstr=pump_exit, downstr=inj_fu_me.get_fluid(), dt=dt)
        + (1 - alpha) * m_dot_me_mfv_array[-1]
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Updating Engines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    # Gas Generator
    gg_dic = GG.transient_time_step(
        ox_in=inj_ox_gg.get_fluid(),
        fu_in=inj_fu_gg.get_fluid(),
        eta_c=loader._engines["GAS_GENERATOR"]._eta_c,
    )

    # Main Engine
    me_dic = ME.transient_time_step(
        ox_in=inj_ox_me.get_fluid(),
        fu_in=inj_fu_me.get_fluid(),
        eta_c=loader._engines["MAIN_ENGINE"]._eta_c,
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Updating Injector Cavitities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    m_ox_gg_cav = (m_dot_gg_mov - gg_dic["m_dot_o"]) * dt
    m_fu_gg_cav = (m_dot_gg_mfv - gg_dic["m_dot_f"]) * dt

    m_ox_me_cav = (m_dot_me_mov - me_dic["m_dot_o"]) * dt
    m_fu_me_cav = (m_dot_me_mfv - me_dic["m_dot_f"]) * dt

    inj_ox_gg.update_pressure(m_dot=m_ox_gg_cav)
    inj_fu_gg.update_pressure(m_dot=m_fu_gg_cav)

    inj_ox_me.update_pressure(m_dot=m_ox_me_cav)
    inj_fu_me.update_pressure(m_dot=m_fu_me_cav)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Updating Turbine Performance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    R_gas = gg_dic["Cp"] * (gg_dic["gamma"] - 1) / gg_dic["gamma"]

    turbine_t = gg_dic["m_dot_t"] * turbine.get_torque(
        T_o=gg_dic["T_o"],
        P_o=gg_dic["P_cc"],
        gamma=gg_dic["gamma"],
        R=R_gas,
        P_exit=P_exit,
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Logging Data from the Curent Timestep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    # ---------- PUMP ----------
    N_pump_array.append(N)
    p_pump_array.append(pump_exit.get_pressure())
    eta_pump_array.append(
        pump.get_eta(Q=(m_dot_gg_mfv + m_dot_me_mfv) / fu_deliver.get_density(), N=N)
    )
    tau_pump_array.append(pump_t)

    # ---------- Turbine --------
    tau_turbine_array.append(turbine_t)

    # ---------- Gas Generator --------
    T_gg_comb_array.append(gg_dic["T_o"])
    Pcc_gg_array.append(gg_dic["P_cc"])
    MR_gg_array.append(gg_dic["MR"])

    dp_dt_gg_array.append(gg_dic["dp_dt"])

    m_dot_gg_inj_fu_array.append(gg_dic["m_dot_f"])
    m_dot_gg_inj_ox_array.append(gg_dic["m_dot_o"])

    # ---------- Main Engine --------
    T_me_comb_array.append(me_dic["T_o"])
    Pcc_me_array.append(me_dic["P_cc"])
    MR_me_array.append(me_dic["MR"])

    dp_dt_me_array.append(me_dic["dp_dt"])

    m_dot_me_inj_fu_array.append(me_dic["m_dot_f"])
    m_dot_me_inj_ox_array.append(me_dic["m_dot_o"])

    # ---------- Valves --------
    m_dot_gg_mov_array.append(m_dot_gg_mov)
    m_dot_gg_mfv_array.append(m_dot_gg_mfv)

    m_dot_me_mov_array.append(m_dot_me_mov)
    m_dot_me_mfv_array.append(m_dot_me_mfv)

    # ---------- Injector Heads -----------
    p_gg_inj_ox_array.append(inj_ox_gg.get_fluid().get_pressure())
    p_gg_inj_fu_array.append(inj_fu_gg.get_fluid().get_pressure())

    p_me_inj_ox_array.append(inj_ox_me.get_fluid().get_pressure())
    p_me_inj_fu_array.append(inj_fu_me.get_fluid().get_pressure())

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Calculating the rate of Change of Pressures for the Cavity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    dp_dt_ox_gg = (p_gg_inj_ox_array[-2] - p_gg_inj_ox_array[-1]) / dt
    dp_dt_fu_gg = (p_gg_inj_fu_array[-2] - p_gg_inj_fu_array[-1]) / dt

    dp_dt_ox_me = (p_me_inj_ox_array[-2] - p_me_inj_ox_array[-1]) / dt
    dp_dt_fu_me = (p_me_inj_fu_array[-2] - p_me_inj_fu_array[-1]) / dt

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Solving for our Next Time Step ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    dp_dt_time = max(
        abs(dp_dt_ox_gg),
        abs(dp_dt_fu_gg),
        abs(dp_dt_ox_me),
        abs(dp_dt_fu_me),
        abs(me_dic["dp_dt"]),
        abs(gg_dic["dp_dt"]),
    )

    if dp_dt_time == 0:
        dt = dt_fix
    else:
        dt = dp / dp_dt_time

    if abs(dt) > dt_fix:
        dt = dt_fix

    t_array.append(t_array[-1] + dt)

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Updating our chamber pressures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    GG.set_pcc_transient(P_cc_transient=(Pcc_gg_array[-1] + dp_dt_gg_array[-1] * dt))
    ME.set_pcc_transient(P_cc_transient=(Pcc_me_array[-1] + dp_dt_me_array[-1] * dt))

    print(f"Current Timestep = {t_array[-1]}")

    # We then perform our acceleration based on the delivered torque, and suplied torque
    dn_dt = (turbine_t - pump_t) / I

    N += dn_dt * dt


############################################ DATA STORAGE ############################################

# Creating a pandas array for all the key pieces of information

data = {
    "t": t_array,
    # Pump Arrays
    "N_pump": N_pump_array,
    "p_pump": p_pump_array,
    "eta_pump": eta_pump_array,
    "tau_pump": tau_pump_array,
    # Turbine Arrays
    "tau_turbine": tau_turbine_array,
    # Gas Generator Arrays
    "T_gg_comb": T_gg_comb_array,
    "Pcc_gg": Pcc_gg_array,
    "MR_gg": MR_gg_array,
    "dp_dt_gg": dp_dt_gg_array,
    "m_dot_gg_inj_fu": m_dot_gg_inj_fu_array,
    "m_dot_gg_inj_ox": m_dot_gg_inj_ox_array,
    "p_gg_inj_ox": p_gg_inj_ox_array,
    "p_gg_inj_fu": p_gg_inj_fu_array,
    # Main Engine
    "T_me_comb": T_me_comb_array,
    "Pcc_me": Pcc_me_array,
    "MR_me": MR_me_array,
    "dp_dt_me": dp_dt_me_array,
    "m_dot_me_inj_fu": m_dot_me_inj_fu_array,
    "m_dot_me_inj_ox": m_dot_me_inj_ox_array,
    "p_me_inj_ox": p_me_inj_ox_array,
    "p_me_inj_fu": p_me_inj_fu_array,
    # Valves
    "m_dot_gg_mov": m_dot_gg_mov_array,
    "m_dot_gg_mfv": m_dot_gg_mfv_array,
    "m_dot_me_mov": m_dot_me_mov_array,
    "m_dot_me_mfv": m_dot_me_mfv_array,
}

df = pd.DataFrame(data=data)

file_dir = pathlib.Path(__file__).parent.resolve().joinpath("simulation.csv")
df.to_csv(file_dir)

print(df["p_pump"])

############################################ Plotting ############################################

# ||||||||||||||||||||||||||||||||| FIG 1 |||||||||||||||||||||||||||||||||||

# First Plot will mainly focus on valve conditions and gas generator conditions
fig, ax = plt.subplots(4, 2, sharex=True)

fig.set_size_inches(18.5, 10.5)

fig.suptitle("TPA + Valve + Injector")

# Pump and Turbine Torques
ax[0][0].plot(df["t"], df["tau_turbine"], label="Turbine Torque")
ax[0][0].plot(df["t"], df["tau_pump"], label="Pump Torque")
ax[0][0].set_ylabel("Torque (Nm)")
ax[0][0].legend()
ax[0][0].set_title("Turbopump Properties")

# Pump Shaft Speed
ax[1][0].plot(df["t"], df["N_pump"] * 60 / (2 * np.pi))
ax[1][0].set_ylabel("Shaft Speed (rpm)")

# Pump Increase in Head Pressure
ax[2][0].plot(
    df["t"],
    (df["p_pump"] - loader._propellants["Fuel"]._P_deliver) / 1e5,
)
ax[2][0].set_ylabel("Pump $\Delta P$ (Bar)")

# Pump Efficiency
ax[3][0].plot(df["t"], df["eta_pump"] * 100, label="Pump Efficiency")
ax[3][0].set_ylabel("Pump $\eta$ (%)")
ax[3][0].set_xlabel("Time (s)")

# GG Valve Mass Flows
ax[0][1].plot(df["t"], df["m_dot_gg_mov"] * 1e3, label="GG_MOV")
ax[0][1].plot(df["t"], df["m_dot_gg_mfv"] * 1e3, label="GG_MFV")

ax[0][1].legend()
ax[0][1].set_ylabel("Mass Flow (g/s)")
ax[0][1].set_title("Valve and Injector")

# ME Valve Mass Flows
ax[1][1].plot(df["t"], df["m_dot_me_mov"] * 1e3, label="ME_MOV")
ax[1][1].plot(df["t"], df["m_dot_me_mfv"] * 1e3, label="ME_MFV")

ax[1][1].legend()
ax[1][1].set_ylabel("Mass Flow (g/s)")

# GG Cavity Pressures
ax[2][1].plot(df["t"], df["p_gg_inj_ox"] / 1e5, label="GG Oxidiser Dome")
ax[2][1].plot(df["t"], df["p_gg_inj_fu"] / 1e5, label="GG Fuel Dome")
ax[2][1].legend()
ax[2][1].set_ylabel("Injector Pressure (Bar)")

# ME Cavity Pressures
ax[3][1].plot(df["t"], df["p_me_inj_ox"] / 1e5, label="ME Oxidiser Dome")
ax[3][1].plot(df["t"], df["p_me_inj_fu"] / 1e5, label="ME Fuel Dome")
ax[3][1].legend()
ax[3][1].set_ylabel("Injector Pressure (Bar)")


# ||||||||||||||||||||||||||||||||| FIG 2 |||||||||||||||||||||||||||||||||||

# First Plot will mainly focus on valve conditions and gas generator conditions
fig2, ax2 = plt.subplots(4, 2, sharex=True)

fig2.set_size_inches(18.5, 10.5)

fig2.suptitle("Main Engine and Gas Generator TCA")

# Gas Generator Combustion Temperature
ax2[0][0].plot(df["t"], df["T_gg_comb"])
ax2[0][0].set_ylabel("Combustion Temperature (K)")
ax2[0][0].set_title("Gas Generator Properties")

# Gas Generator Mixture Ratio
ax2[1][0].plot(df["t"], df["MR_gg"])
ax2[1][0].set_ylabel("MR (n.d)")

# Gas Generator Chamber Pressure (Bar)
ax2[2][0].plot(
    df["t"],
    df["Pcc_gg"] / 1e5,
)
ax2[2][0].set_ylabel("Chamber Pressure (Bar)")

# Engine Mass Flow Rate
ax2[3][0].plot(
    df["t"],
    (df["m_dot_gg_inj_fu"] + df["m_dot_gg_inj_ox"]) * 1e3,
)
ax2[3][0].set_ylabel("Mass Flow Rate (g/s)")
ax2[3][0].set_xlabel("Time (s)")


# Gas Generator Combustion Temperature
ax2[0][1].plot(df["t"], df["T_me_comb"], color="r")
ax2[0][1].set_ylabel("Combustion Temperature (K)")
ax2[0][1].set_title("Main Engine Properties")

# Gas Generator Mixture Ratio
ax2[1][1].plot(df["t"], df["MR_me"], color="r")
ax2[1][1].set_ylabel("MR (n.d)")

# Gas Generator Chamber Pressure (Bar)
ax2[2][1].plot(df["t"], df["Pcc_me"] / 1e5, color="r")
ax2[2][1].set_ylabel("Chamber Pressure (Bar)")

# Engine Mass Flow Rate
ax2[3][1].plot(
    df["t"], (df["m_dot_me_inj_fu"] + df["m_dot_me_inj_ox"]) * 1e3, color="r"
)
ax2[3][1].set_ylabel("Mass Flow Rate (g/s)")
ax2[3][1].set_xlabel("Time (s)")


plt.show()
