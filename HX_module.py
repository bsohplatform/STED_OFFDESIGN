import numpy as np
import HP_dataclass
from CoolProp.CoolProp import PropsSI

class PHE_module:
    def __init__(self):
        a = 1
        
    def PHE_cond(self, primary_in, primary_out, secondary_in, secondary_out, PHE_Inputs, noHX):
        x_primary = np.zeros(shape=(PHE_Inputs.phx_N_element+1))
        d_primary = np.zeros(shape=(PHE_Inputs.phx_N_element+1))
        T_primary = np.zeros(shape=(PHE_Inputs.phx_N_element+1))
        htc_primary = np.zeros(shape=(PHE_Inputs.phx_N_element+1))
        Re_primary = np.zeros(shape=(PHE_Inputs.phx_N_element+1))
        h_primary = np.zeros(shape=(PHE_Inputs.phx_N_element+1))
        T_secondary = np.zeros(shape=(PHE_Inputs.phx_N_element+1))
        htc_secondary = np.zeros(shape=(PHE_Inputs.phx_N_element+1))
        Re_secondary = np.zeros(shape=(PHE_Inputs.phx_N_element+1))
        a_cond_PHE = 1 
        if noHX == 0:
            T_sec_out_lb = secondary_out.T
            T_sec_out_ub = secondary_out.p
            Cp_sec = secondary_out.c
        elif noHX == 1:
            T_sec_out_lb = secondary_in.T
            T_sec_out_ub = PropsSI('T','P',secondary_out.p,'Q',0.0,secondary_out.Y)
            Cp_sec = secondary_in.c
        else:
            T_sec_out_lb = secondary_out.T
            T_sec_out_ub = secondary_out.p
            mdot_sec_lb = 0;
            mdot_sec_ub = 2*(primary_in.h - primary_in.hl)*primary_in.m/(secondary_out.h - secondary_in.h)
            Cp_sec = 0.5*(secondary_in.c + secondary_out.c)
            
        G_cond_ref = (primary_in.m/PHE_Inputs.)
        
        while a_cond_PHE:
            T_sec_out = 0.5*(T_sec_out_lb + T_sec_out_ub)
            if noHX == 2:
                mdot_sec = 0.5*(mdot_sec_lb + mdot_sec_ub)
            
            T_primary[0] = primary_in.T
            h_primary[0] = primary_in.h
            T_secondary[0] = T_sec_out
            for jj in range(PHE_Inputs.phx_N_element):
                Re_secondary[jj] = 
                
            
        
    def PHE_evap(self, PHE_Inputs, noHX):
        if noHX == 0:
            
        elif noHX == 1:
        else: