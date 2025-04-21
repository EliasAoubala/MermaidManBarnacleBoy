from turborocket.meanline.meanline_relations import TurbineStageDesign
from turborocket.fluids.ideal_gas import IdealFluid

import json

fluid = IdealFluid(p=37.9e5, t=889, R=258, gamma=1.108, cp=2640)

stage = TurbineStageDesign(gas=fluid, m_dot=3.24, omega=20000, alpha=20)


stage.set_operating_point(u_cis=0.21, Rt=15, b=15e-3, t=5e-3, delta_r=0.5e-3, N=13)

result = stage.solve_performance()

print("Pressures: ")
print(json.dumps(result["pressure"], indent=1))
print("\n")

print("Temperatures: ")
print(json.dumps(result["temperature"], indent=1))
print("\n")

print("Velocity: ")
print(json.dumps(result["velocity"], indent=1))
print("\n")

print("Performance: ")
print(json.dumps(result["performance"], indent=1))
print("\n")

print("Geometry: ")
print(json.dumps(result["geometry"], indent=1))
print("\n")

print("Mach: ")
print(json.dumps(result["mach"], indent=1))
print("\n")

print("Angles: ")
print(json.dumps(result["angles"], indent=1))
print("\n")

# Validated up to equation 46
#


# Testing the critical velocity theory, as presented in the paper

c = 1000  # m/s

t_o = 800  # K

R = 200

gamma = 1.1

a_star = (((2 * gamma) / (gamma + 1)) * R * t_o) ** (1 / 2)

m_star = c / a_star

t_test = t_o * (1 - ((gamma - 1) / (gamma + 1)) * m_star**2)

cp = R * (gamma) / (gamma - 1)

t_act = t_o - (c**2) / (2 * cp)

print(f"Their method: {t_test}")
print(f"Our method: {t_act}")
