from CoolProp.CoolProp import PropsSI
from HP_dataclass import *

class COMP_module:
    def __init__(self, primary_in, primary_out, Inputs):
        self.primary_in = primary_in
        self.primary_out = primary_out
        self.Inputs = Inputs
        
    def Off(self, mode='poly', DSH=5.0):
        if mode == 'poly':
            self.primary_in.T =  DSH + self.primary_in.Ts
            self.primary_in.d = Aux_fn.PropCal(self.primary_in, 'D', 'T', 'P')
            self.primary_in.h = Aux_fn.PropCal(self.primary_in, 'H', 'T', 'P')
            self.primary_in.s = Aux_fn.PropCal(self.primary_in, 'S', 'T', 'P')
            
            n_comp = self.Inputs.comp_n_poly
            V_comp = self.Inputs.comp_V_dis
            f_comp = self.Inputs.comp_frequency
            C_comp = self.Inputs.comp_C_gap
            
            eff_vol = 1 + C_comp - C_comp*pow(self.primary_out.p/self.primary_in.p,1/n_comp)
            self.primary_in.m = eff_vol*self.primary_in.d*V_comp*f_comp
            self.primary_out.m = self.primary_in.m
            
            w_comp = (n_comp/(n_comp-1))*(self.primary_in.p/self.primary_in.d)*(pow(self.primary_out.p/self.primary_in.p,(n_comp-1)/n_comp) - 1); 
            comp_W = self.primary_in.m*w_comp/self.Inputs.comp_eff_mech
            
            self.primary_out.h = self.primary_in.h + w_comp
            self.primary_out.T = Aux_fn.PropCal(self.primary_out, 'T', 'H', 'P')
            h_comp_out_ideal = PropsSI('H','P',self.primary_out.p,'S',self.primary_in.s,self.primary_out.Y)
            comp_eff_isen = (h_comp_out_ideal - self.primary_in.h)/(self.primary_out.h - self.primary_in.h)
            
            try:
                self.primary_out.Ts = PropsSI('T','P',self.primary_out.T, 'Q', 1.0, self.primary_out.Y)
            except:
                self.primary_out.Ts = PropsSI('TCRIT','',0,'',0,self.primary_out.Y)
            if abs(self.primary_out.T - self.primary_out.Ts) < 0.1:
                DSH = DSH + 1.0
                a = 0
            else:
                a = 1
            
            return (comp_W, comp_eff_isen, DSH, a)
            

if __name__ == '__main__':
    comp_in = Fluid_flow(Y='R410A', T=290.0, p = 1.0e6)
    comp_in.Ts = PropsSI('T','P',1.0e6,'Q',1.0,comp_in.Y)
    comp_out = Fluid_flow(Y='R410A',p = 3.0e6)
    inputs = Inputs()
    inputs.comp_n_poly = 1.5
    inputs.comp_V_dis = 21.0e-6
    inputs.comp_frequency = 60
    inputs.comp_C_gap = 0.05
    inputs.DSH = 5.0
    comp = COMP_module(comp_in, comp_out, inputs)
    (comp.comp_W, comp.comp_eff_isen, DSH, cond_a) = comp.Off(DSH = 5.0)
    
    print(comp.comp_W)