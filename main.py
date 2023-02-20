from HP_dataclass import *
import numpy as np
from CoolProp.CoolProp import PropsSI
import COMP_module

class VCHP_off:
    def __init__(self):
        a = 1
        
    def Input_Processing(self, InCond, OutCond, InEvap, OutEvap):
        if InCond.p <= 0.0:
            InCond.p = 101300.0
            
        if OutCond.p <= 0.0:
            OutCond.p = 101300.0
        
        if InCond.T <= 0.0:
            noCond = 0
        
        if OutCond.T <= 0.0:
            noCond = 1
        
        if InCond.m <= 0.0 and OutCond.m <= 0:
            noCond = 2
        else:
            if InCond.m == 0:
                InCond.m = OutCond.m
            else:
                OutCond.m = InCond.m
        
        if InEvap.p <= 0.0:
            InEvap.p = 101300.0
            
        if OutEvap.p <= 0.0:
            OutEvap.p = 101300.0
        
        if InEvap.T <= 0.0:
            noEvap = 0
        
        if OutEvap.T <= 0.0:
            noEvap = 1
        
        if InEvap.m <= 0.0 and OutEvap.m <= 0:
            noEvap = 2
        else:
            if InEvap.m == 0:
                InEvap.m = OutEvap.m
            else:
                OutEvap.m = InEvap.m
        
        if noCond == 0:
            OutCond.c = Aux_fn.PropCal(OutCond, 'C', 'T', 'P')
            OutCond.v = Aux_fn.PropCal(OutCond, 'V', 'T', 'P')
            OutCond.pr = Aux_fn.PropCal(OutCond, 'Prandtl', 'T', 'P')
            OutCond.l = Aux_fn.PropCal(OutCond, 'L', 'T', 'P')
            OutCond.d = Aux_fn.PropCal(OutCond, 'D', 'T', 'P')
        elif noCond == 1:
            InCond.c = Aux_fn.PropCal(InCond, 'C', 'T', 'P')
            InCond.v = Aux_fn.PropCal(InCond, 'V', 'T', 'P')
            InCond.pr = Aux_fn.PropCal(InCond, 'Prandtl', 'T', 'P')
            InCond.l = Aux_fn.PropCal(InCond, 'L', 'T', 'P')
            InCond.d = Aux_fn.PropCal(InCond, 'D', 'T', 'P')            
        else:
            InCond.c = Aux_fn.PropCal(InCond, 'C', 'T', 'P')
            InCond.v = Aux_fn.PropCal(InCond, 'V', 'T', 'P')
            InCond.pr = Aux_fn.PropCal(InCond, 'Prandtl', 'T', 'P')
            InCond.l = Aux_fn.PropCal(InCond, 'L', 'T', 'P')
            InCond.d = Aux_fn.PropCal(InCond, 'D', 'T', 'P')
            OutCond.c = Aux_fn.PropCal(OutCond, 'C', 'T', 'P')
            OutCond.v = Aux_fn.PropCal(OutCond, 'V', 'T', 'P')
            OutCond.pr = Aux_fn.PropCal(OutCond, 'Prandtl', 'T', 'P')
            OutCond.l = Aux_fn.PropCal(OutCond, 'L', 'T', 'P')
            OutCond.d = Aux_fn.PropCal(OutCond, 'D', 'T', 'P')
        
        if noEvap == 0:
            OutEvap.c = Aux_fn.PropCal(OutEvap, 'C', 'T', 'P')
            OutEvap.v = Aux_fn.PropCal(OutEvap, 'V', 'T', 'P')
            OutEvap.pr = Aux_fn.PropCal(OutEvap, 'Prandtl', 'T', 'P')
            OutEvap.l = Aux_fn.PropCal(OutEvap, 'L', 'T', 'P')
            OutEvap.d = Aux_fn.PropCal(OutEvap, 'D', 'T', 'P')
        elif noEvap == 1:
            InEvap.c = Aux_fn.PropCal(InEvap, 'C', 'T', 'P')
            InEvap.v = Aux_fn.PropCal(InEvap, 'V', 'T', 'P')
            InEvap.pr = Aux_fn.PropCal(InEvap, 'Prandtl', 'T', 'P')
            InEvap.l = Aux_fn.PropCal(InEvap, 'L', 'T', 'P')
            InEvap.d = Aux_fn.PropCal(InEvap, 'D', 'T', 'P')
        else:
            InEvap.c = Aux_fn.PropCal(InEvap, 'C', 'T', 'P')
            InEvap.v = Aux_fn.PropCal(InEvap, 'V', 'T', 'P')
            InEvap.pr = Aux_fn.PropCal(InEvap, 'Prandtl', 'T', 'P')
            InEvap.l = Aux_fn.PropCal(InEvap, 'L', 'T', 'P')
            InEvap.d = Aux_fn.PropCal(InEvap, 'D', 'T', 'P')
            OutEvap.c = Aux_fn.PropCal(OutEvap, 'C', 'T', 'P')
            OutEvap.v = Aux_fn.PropCal(OutEvap, 'V', 'T', 'P')
            OutEvap.pr = Aux_fn.PropCal(OutEvap, 'Prandtl', 'T', 'P')
            OutEvap.l = Aux_fn.PropCal(OutEvap, 'L', 'T', 'P')
            OutEvap.d = Aux_fn.PropCal(OutEvap, 'D', 'T', 'P')
            
        return (InCond, OutCond, InEvap, OutEvap, noCond, noEvap)
    
    def HX_Geometry(self, HX_input):
        if HX_input == 'phe':
            A_cross = 
            
        
        
    
    def OffDesign_Solver(self, InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, Cycle_Inputs, Comp_Inputs, Cond_Inputs, Evap_Inputs, Outputs, noCond, noEvap):
        cond_p_ub = PropsSI('PCRIT','',0,'',0,InCond_REF.Y)
        if noCond == 0:
            cond_p_lb = PropsSI('P','T',OutCond.T, 'Q', 1.0, OutCond_REF.Y)
        else:
            cond_p_lb = PropsSI('P','T',InCond.T, 'Q', 1.0, InCond_REF.Y)
        
        cond_a = 1
        while cond_a:
            InCond_REF.p = 0.5*(cond_p_ub + cond_p_lb)
            InCond_REF.Ts = PropsSI('T','P',InCond_REF.p,'Q',1.0, InCond_REF.Y )
            InCond_REF.hl = PropsSI('H','P',InCond_REF.p,'Q',0.0, InCond_REF.Y )
            InCond_REF.hg = PropsSI('H','P',InCond_REF.p,'Q',1.0, InCond_REF.Y )
            InCond_REF.dl = PropsSI('D','P',InCond_REF.p,'Q',0.0, InCond_REF.Y )
            InCond_REF.dg = PropsSI('D','P',InCond_REF.p,'Q',1.0, InCond_REF.Y )
            
            if noEvap == 0:
                evap_p_ub = PropsSI('P','T',OutEvap.T, 'Q', 1.0, OutEvap_REF.Y)
            else:
                evap_p_ub = PropsSI('P','T',InEvap.T, 'Q', 1.0, InEvap_REF.Y)
            evap_p_lb = 101300.0
            
            
            evap_a = 1
            while evap_a:
                OutEvap_REF.p = 0.5*(evap_p_ub+evap_p_lb)
                OutEvap_REF.Ts = PropsSI('T','P',OutEvap_REF.p,'Q',1.0, OutEvap_REF.Y )
                OutEvap_REF.hl = PropsSI('H','P',OutEvap_REF.p,'Q',0.0, OutEvap_REF.Y )
                OutEvap_REF.hg = PropsSI('H','P',OutEvap_REF.p,'Q',1.0, OutEvap_REF.Y )
                OutEvap_REF.dl = PropsSI('D','P',OutEvap_REF.p,'Q',0.0, OutEvap_REF.Y )
                OutEvap_REF.dg = PropsSI('D','P',OutEvap_REF.p,'Q',1.0, OutEvap_REF.Y )
                
                comp = COMP_module()
                (OutEvap_REF, InCond_REF, Outputs.comp_W, Outputs.comp_eff_isen, Cycle_Inputs.DSH, evap_a) = comp.Off(OutEvap_REF, InCond_REF, Comp_Inputs, DSH = Cycle_Inputs.DSH)
                
if __name__ == '__main__':
    inputs = {
        'hot_fluid': 'water',
        'hot_T_in': 305.15,
        'hot_T_out': 0.0,
        'hot_p': 101300,
        'hot_m': 1.0,
        'cold_fluid': 'water',
        'cold_T_in': 285.15,
        'cold_T_out': 0.0,
        'cold_p':1013000,
        'cold_m':1.0,
        'layout':'bas',
        'DSH':5.0,
        }