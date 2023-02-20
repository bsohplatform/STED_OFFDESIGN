from dataclasses import dataclass, field
from CoolProp.CoolProp import PropsSI

@dataclass
class Fluid_flow:
    Y: str
    m: float=0.0
    T: float=0.0
    Ts: float=0.0
    p: float=0.0
    h: float=0.0
    hl: float=0.0
    hg: float=0.0
    s: float=0.0
    d: float=0.0
    dl: float=0.0
    dg: float=0.0
    c: float=0.0
    l: float=0.0
    v: float=0.0
    pr: float=0.0
    Tcr: float=0.0
    Pcr: float=0.0
    
    
@dataclass
class Cycle_Inputs:
    cycle_layout: str = 'bas'
    cycle_DSH: float = 0.0    
    cycle_tol: float = 1.0e-3

@dataclass
class Comp_Inputs:
    comp_V_dis: float = 0.0 # Displacement volume
    comp_C_gap: float = 0.0 # Clearance factor (clearance/displacement)
    comp_n_poly: float = 0.0 # Polytropic number
    comp_eff_mech: float = 1.0 # mechanical_efficiency
    comp_frequency: float = 60.0 # Compressor frequency [Hz]

@dataclass
class PHE_Inputs:
    type = 'phe'
    phx_N_element: int = 30
    phx_N_plate: int = 0 # Number of Plates
    phx_phi: float = 0.0 # Ratio of developed length to projected length
    phx_thk_plate: float = 0.0 # Single plate thickness
    phx_thk_tot: float = 0.0 # Total thickness of PHE
    phx_L_vert: float = 0.0 # Vertical length of PHE (Center to center of inlet and outlet ports)
    phx_L_hor: float = 0.0 # Horizontal length of PHE (Center to center of inlet and outlet ports)
    phx_L_width: float = 0.0 # total horizontal length of PHE
    phx_D_p: float = 0.0 # Port Diameter
    phx_beta: float = 0.0 # chevron angle

@dataclass
class Outputs:
    comp_W: float = 0.0
    evap_Q: float = 0.0
    cond_Q: float = 0.0
    
    comp_eff_isen: float = 0.0


class Aux_fn:
    @staticmethod
    def PropCal(Fluid_flow, outpar: str, par1: str, par2: str):
        prop = {'T':Fluid_flow.T, 'P':Fluid_flow.p, 'H':Fluid_flow.h, 'S': Fluid_flow.s, 
                         'D':Fluid_flow.d, 'C': Fluid_flow.c, 'L': Fluid_flow.l}
                
        outval = PropsSI(outpar,par1,prop[par1],par2,prop[par2],Fluid_flow.Y)
        
        return outval
    

if __name__ == '__main__':
    cond = Fluid_flow(Y='Water', T=300, p = 101300)
    
    cond.d = Aux_fn.PropCal(cond, 'D', 'T', 'P')
    print(cond.d)