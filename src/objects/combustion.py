"""
This file contains all the combustion objects for modelling the off-nominal performance

"""

from rocketcea.cea_obj import CEA_Obj
import CoolProp.CoolProp as CP
from scipy.optimize import minimize


class GasGenerator:

    def __init__(
        self,
        m_dot: float,
        o_f: float,
        P_c: float,
        Oxid: str = "Nitrous",
        Fuel: str = "IPA",
    ):

        self._m_dot = m_dot
        self._o_f = o_f
        self._P_c = P_c
        self._oxid = Oxid
        self._Fuel = Fuel

        self._cea = CEA_Obj(
            oxName=Oxid,
            fuelName=Fuel,
            pressure_units="Bar",
            temperature_units="K",
            density_units="kg/m^3",
            specific_heat_units="J/kg-K",
        )

        return

    def evaluate_geometry(self, P_f: float, P_o: float, c_star_eta: float) -> None:
        """
        Function that solves for the injector geometry of the gas generator
        This equatino assumes the SPI model

        Args:
            P_f (float): Nominal Fuel Supply Pressure
            P_o (float): Nominal Oxidiser Supply Pressure
            c_star_eta (float): C* efficiency
        """

        # Fuel Orifice

        self._rho_f = 786

        dp_f = (P_f - self._P_c) * 1e5

        m_dot_f = self._m_dot / (1 + self._o_f)

        self._CdAf = m_dot_f / (2 * self._rho_f * dp_f) ** (1 / 2)

        # Oxidiser Orifice

        self._rho_o = CP.PropsSI("D", "P", P_f * 1e5, "Q", 0, self._oxid)

        dp_o = (P_o - self._P_c) * 1e5

        m_dot_o = self._m_dot - m_dot_f

        self._CdAo = m_dot_o / (2 * self._rho_o * dp_o) ** (1 / 2)

        # Combustion Chamber

        c_star_t = self._cea.get_Cstar(Pc=self._P_c, MR=self._o_f)

        c_star_a = c_star_t * c_star_eta

        self._CdAc = (c_star_a * self._m_dot) / (self._P_c * 1e5)

        return

    def injector_flow(self, P_d: float, P_c: float, prop: str) -> float:
        """This equation solves for the mass flow rate through a given injector

        Args:
            P_d (float): Injector Delivery Pressure
            P_c (float): Combustion Chamber Pressure

        Returns:
            float: Mass Flow rate through injector
        """

        if prop == "Oxidiser":
            cd_a = self._CdAo
            rho_a = CP.PropsSI("D", "P", P_d * 1e5, "Q", 0, self._oxid)

        elif prop == "Fuel":
            cd_a = self._CdAf
            rho_a = self._rho_f

        m_dot = cd_a * (2 * rho_a * (P_d - P_c) * 1e5) ** (1 / 2)

        return m_dot

    def error_func(
        self: float,
        P_c_g: float,
        P_f: float,
        P_o: float,
        c_star_eta: float,
    ) -> float:

        # We evaluate for our injector flow rates

        m_dot_f = self.injector_flow(P_d=P_f, P_c=P_c_g, prop="Fuel")

        m_dot_o = self.injector_flow(P_d=P_o, P_c=P_c_g, prop="Oxidiser")

        m_dot_inj = m_dot_o + m_dot_f

        # We now evaluate for the combustion chamber cnoditions

        o_f = m_dot_o / m_dot_f

        C_s_t = self.get_Cstar(Pc=P_c_g, MR=o_f)

        C_s_a = C_s_t * c_star_eta

        m_dot_c = P_c_g * 1e5 * self._CdAc / C_s_a

        error = abs(m_dot_c - m_dot_inj)

        return error

    def performance_evaluation(
        self,
        P_f: float,
        P_o: float,
        c_star_eta: float,
    ) -> list[float]:
        """
        This equation solves for the updated performance conditions of the gas generator given a change of inlet pressures.
        This process

        This is inherintly an iterative process to figure out what the new chamber pressure and flow rates of the gas generator will be, hence we use an auxilary function which we can then use the scipy library to figure out what the requried mass flow rates are.

        Args:
            P_f (float): Fuel Delivery Pressure
            P_o (float): Oxidiser Delivery Pressure
            c_star (float): C* efficiency

        Returns:
            list[float]:
        """

        res = minimize(
            fun=self.error_func,
            x0=self._P_c,
            tol=1e-5,
            args=[P_f, P_o, c_star_eta],
            bounds=(0.1, 20),
        )

        P_c = res.x[0]

        # We can now derive performance accordingly

        m_dot_f = self.injector_flow(P_d=P_f, P_c=P_c, prop="Fuel")

        m_dot_o = self.injector_flow(P_d=P_o, P_c=P_c, prop="Oxidiser")

        m_dot_inj = m_dot_o + m_dot_f

        # We now evaluate for the combustion chamber cnoditions

        o_f = m_dot_o / m_dot_f

        C_s_t = self.get_Cstar(Pc=P_c, MR=o_f)

        C_s_a = C_s_t * c_star_eta

        m_dot_c = P_c * 1e5 * self._CdAc / C_s_a

        print(f"Iteration Error {(m_dot_c-m_dot_inj)*1e2/m_dot_inj} %")

        # What key information do we still need

        Cp = self._cea.get_Chamber_Cp(Pc=P_c, MR=o_f, frozen=0)

        t = self._cea.get_Temperatures(Pc=P_c, MR=o_f)

        return [P_c, t, h, Cp]
