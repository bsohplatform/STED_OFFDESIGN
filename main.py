from HP_dataclass import *
from COMP_module import *
from HX_module import *
from CoolProp.CoolProp import PropsSI

class VCHP_off:
    def __init__(self, evap_BC, cond_BC):
        self.evap_BC = evap_BC
        self.cond_BC = cond_BC
        
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
                
                comp = COMP_module()
                (OutEvap_REF, InCond_REF, Outputs.comp_W, Outputs.comp_eff_isen, Outputs.DSH, evap_a) = comp.Off(OutEvap_REF, InCond_REF, Comp_Inputs, DSH = Cycle_Inputs.DSH)
                OutCond_REF.m = InCond_REF.m
                InEvap_REF.m = OutEvap_REF.m
                
                cond = HX_module(hx_type='phx', cor = False, Inputs=Cond_Inputs)
                (InCond_REF, OutCond_REF, InCond, OutCond, cond_Q, cond_rho)=cond.PHX('cond',InCond_REF, OutCond_REF, InCond, OutCond, noCond)
                Outputs.DSC = OutCond_REF.Ts - OutCond_REF.T
                
                evap = HX_module(hx_type='phx', cor = False, Inputs=Evap_Inputs)
                (InEvap_REF, OutEvap_REF, InEvap, OutEvap, evap_Q, evap_rho)=evap.PHX('evap',InEvap_REF, OutEvap_REF, InEvap, OutEvap, noEvap)
                
                if InEvap_REF.h == 0:
                    evap_p_lb = OutEvap_REF.p
                    err_evap_p = 1
                else:
                    if self.evap_BC == 'q':
                        err_evap_p = (Cycle_Inputs.evap_Q - evap_Q)/Cycle_Inputs.evap_Q_cal
                    else:    
                        err_evap_p = (InEvap_REF.h - OutCond_REF.h)/OutCond_REF.h
                
                    if err_evap_p < 0:
                        evap_p_lb = OutEvap_REF.p
                    else:
                        evap_p_ub = OutEvap_REF.p
                
                if abs(err_evap_p) < 1.0e-3:
                    evap_a = 0
                elif evap_p_ub - evap_p_lb < 0.1:
                    evap_a = 0
                    
            if self.cond_BC == 'dsc':
                err_cond_p = (Outputs.DSC - Cycle_Inputs.DSC)/Cycle_Inputs.DSC 
            elif self.cond_BC == 'q':
                err_cond_p = (cond_Q - Cycle_Inputs.cond_Q)/Cycle_Inputs.cond_Q
            else:
                M_ref_cal = (cond.V*cond_rho+evap.V*evap_rho)
                err_cond_p = (M_ref_cal - Cycle_Inputs.M_ref)/Cycle_Inputs.M_ref
            
            if err_cond_p < 0:
                cond_p_lb = InCond_REF.p
            else:
                cond_p_ub = InCond_REF.p    
            
            if abs(err_cond_p) < 1.0e-3:
                cond_a = 0
            elif cond_p_ub - cond_p_lb < 0.1:
                cond_a = 0
        
        Outputs.M_ref = cond.V*cond_rho+evap.V*evap_rho
        Outputs.mean_rho = Outputs.M_ref/(cond.V + evap.V)
        Outputs.evap_Q = evap_Q
        Outputs.cond_Q = cond_Q
        
        return (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, Outputs)
        
if __name__ == '__main__':
    cycle_inputs = Cycle_Inputs()
    comp_inputs = Comp_Inputs()
    cond_inputs = PHX_Inputs()
    evap_inputs = PHX_Inputs()
    InEvap = Fluid_flow(Y='water',m=0.191, p = 101300.0)
    OutEvap = Fluid_flow(Y='water',m=0.191, T=280.15, p = 101300.0)
    InCond = Fluid_flow(Y='water',m=0.228,T=305.15, p = 101300.0)
    OutCond = Fluid_flow(Y='water',m=0.228, p = 101300.0)
    InEvap_REF = Fluid_flow(Y='R410A')
    OutEvap_REF = Fluid_flow(Y='R410A')
    InCond_REF = Fluid_flow(Y='R410A')
    OutCond_REF = Fluid_flow(Y='R410A')
    outputs = Outputs()
    
    cycle_inputs.layout = 'bas'
    cycle_inputs.DSH = 3.0
    cycle_inputs.DSC = 0.01
    
    comp_inputs.n_poly = 2.5
    comp_inputs.V_dis = 12.0e-6
    comp_inputs.frequency = 60
    comp_inputs.C_gap = 0.05
    
    cond_inputs.UA = 1600.0
    cond_inputs.dp = 0.01
    cond_inputs.mdot_nominal = 0.228
    
    evap_inputs.UA = 1350.0
    evap_inputs.dp = 0.01
    evap_inputs.mdot_nominal = 0.191
    
    
    bas_off = VCHP_off('h','dsc')
    (InCond, OutCond, InEvap, OutEvap, noCond, noEvap) = bas_off.Input_Processing(InCond, OutCond, InEvap, OutEvap)
    (InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, outputs) = bas_off.OffDesign_Solver(InCond, OutCond, InEvap, OutEvap, InCond_REF, OutCond_REF, InEvap_REF, OutEvap_REF, cycle_inputs, comp_inputs, cond_inputs, evap_inputs, outputs, noCond, noEvap)