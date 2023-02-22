import numpy as np
from CoolProp.CoolProp import PropsSI

class HX_module:
    def __init__(self, hx_type, Inputs):
        if hx_type == 'phx':
            self.phx_inputs = Inputs
            b = self.phx_inputs.thk_tot/self.phx_inputs.N_plate - self.phx_inputs.thk_plate
            self.A_cx = self.phx_inputs.N_plate*b*self.phx_inputs.L_width
            self.Dh = 2*b
            self.R_plate = self.phx_inputs.thk_plate/15.0
            self.dA_flow = self.phx_inputs.L_width*self.phx_inputs.L_vert*(2*self.phx_inputs.N_plate - 1)/(self.phx_inputs.N_element - 1)
            self.V = b*self.phx_inputs.L_width*self.phx_inputs.L_vert
            
    def PHX_cond(self, primary_in, primary_out, secondary_in, secondary_out, noHX):
        
        A_cx = self.A_cx
        Dh = self.Dh
        
        G_primary = np.zeros(shape=(self.phx_inputs.N_element+1))
        x_primary = np.zeros(shape=(self.phx_inputs.N_element+1))
        T_primary = np.zeros(shape=(self.phx_inputs.N_element+1))
        h_primary = np.zeros(shape=(self.phx_inputs.N_element+1))
        htc_primary = np.zeros(shape=(self.phx_inputs.N_element+1))
        Re_primary = np.zeros(shape=(self.phx_inputs.N_element+1))
        pr_primary = np.zeros(shape=(self.phx_inputs.N_element+1))
        f_primary = np.zeros(shape=(self.phx_inputs.N_element+1))
        p_primary = np.zeros(shape=(self.phx_inputs.N_element+1))
        d_primary = np.zeros(shape=(self.phx_inputs.N_element+1))
        
        T_secondary = np.zeros(shape=(self.phx_inputs.N_element+1))    
        
        U_tot = np.zeros(shape=(self.phx_inputs.N_element+1))
        eps = np.zeros(shape=(self.phx_inputs.N_element+1))
        Q_trans = np.zeros(shape=(self.phx_inputs.N_element+1))
        
        
        a_cond_PHX = 1 
        if noHX == 0:
            T_sec_out_lb = secondary_out.T
            T_sec_out_ub = secondary_out.T
            c_sec = secondary_out.c
            v_sec = secondary_out.v
            l_sec = secondary_out.l
            pr_sec = secondary_out.pr
            try:
                mdot_sec = secondary_out.m
            except:
                mdot_sec = secondary_in.m
        elif noHX == 1:
            T_sec_out_lb = secondary_in.T
            T_sec_out_ub = PropsSI('T','P',secondary_out.p,'Q',0.0,secondary_out.Y)
            c_sec = secondary_in.c
            v_sec = secondary_in.v
            l_sec = secondary_in.l
            pr_sec = secondary_in.pr
            try:
                mdot_sec = secondary_in.m
            except:
                mdot_sec = secondary_out.m
        else:
            T_sec_out_lb = secondary_out.T
            T_sec_out_ub = secondary_out.T
            mdot_sec_lb = 0;
            mdot_sec_ub = 2*(primary_in.h - primary_in.hl)*primary_in.m/(secondary_out.h - secondary_in.h)
            c_sec = 0.5*(secondary_in.c + secondary_out.c)
            v_sec = 0.5*(secondary_in.v + secondary_out.v)
            l_sec = 0.5*(secondary_in.l + secondary_out.l)
            pr_sec = c_sec*v_sec/l_sec
            
        T_primary[0] = primary_in.T
        h_primary[0] = primary_in.h
        p_primary[0] = primary_in.p
        hl_primary = primary_in.hl
        hg_primary = primary_in.hg
        
        x_primary[0] = (h_primary[0]-hl_primary)/(hg_primary - hl_primary)
        
        G_primary_ref = primary_in.m/A_cx
        if x_primary[0] > 0 and x_primary[0] < 1:
            dl_primary = PropsSI('D','P',p_primary[0],'Q',0.0, primary_in.Y)
            dg_primary = PropsSI('D','P',p_primary[0],'Q',1.0, primary_in.Y)
            d_primary[0] = 1/(x_primary[0]/dg_primary+(1-x_primary[0])/dl_primary)
            vl_primary = PropsSI('V','P',p_primary[0],'Q',0.0, primary_in.Y)
            vg_primary = PropsSI('V','P',p_primary[0],'Q',1.0, primary_in.Y)
            cl_primary = PropsSI('C','P',p_primary[0],'Q',0.0, primary_in.Y)
            ll_primary = PropsSI('L','P',p_primary[0],'Q',0.0, primary_in.Y)
    
            Xtt = pow((1-x_primary[0])/x_primary[0],0.9)*pow(dg_primary/dl_primary, 0.5)*pow(vl_primary/vg_primary, 0.1)
            
            G_primary[0] = G_primary_ref*((1-Xtt)+Xtt*pow(dl_primary/dg_primary,0.5))
            Re_primary[0] = G_primary[0]*Dh/vl_primary
            pr_primary[0] = cl_primary*vl_primary/ll_primary
            htc_primary[0] = 4.118*pow(Re_primary[0],0.4)*pow(pr_primary[0],1/3)
            htc_primary[0] = htc_primary[0]*ll_primary/Dh
            d_primary[0] = 1/(x_primary[0]/dg_primary+(1-x_primary[0])/dl_primary)
        else:
            d_primary[0] = PropsSI('D','T',T_primary[0],'P',p_primary[0],primary_in.Y)
            v_primary = PropsSI('V','T',T_primary[0],'P',p_primary[0],primary_in.Y)
            c_primary = PropsSI('C','T',T_primary[0],'P',p_primary[0],primary_in.Y)
            l_primary = PropsSI('L','T',T_primary[0],'P',p_primary[0],primary_in.Y)
            
            G_primary[0] = G_primary_ref
            Re_primary[0] = G_primary_ref*Dh/v_primary
            pr_primary[0] = c_primary*v_primary/l_primary
            htc_primary[0] = 0.724*pow(self.phx_inputs.beta/30,0.646)*pow(Re_primary[0], 0.583)*pow(pr_primary[0], 1/3)
            htc_primary[0] = htc_primary[0]*l_primary/Dh
        
        while a_cond_PHX:
            T_sec_out = 0.5*(T_sec_out_lb + T_sec_out_ub)
            if noHX == 2:
                mdot_sec = 0.5*(mdot_sec_lb + mdot_sec_ub)
            
            T_secondary[0] = T_sec_out
            G_secondary = mdot_sec/A_cx
            Re_secondary = G_secondary*Dh/v_sec
            htc_secondary = 0.724*pow(self.phx_inputs.beta/30,0.646)*pow(Re_secondary, 0.583)*pow(pr_sec, 1/3)
            htc_secondary = htc_secondary*l_sec/Dh
            
            U_tot[0] = 1/(1/htc_primary[0]+self.R_plate+1/htc_secondary)
            if x_primary[0] > 0 and x_primary[0] < 1:
                C_min = c_sec*mdot_sec
                NTU = U_tot[0]*self.dA_flow/C_min
                eps[0] = 1 - np.exp(-NTU)
            else:
                C_min = min(c_primary*primary_in.m, c_sec*mdot_sec)
                C_max = max(c_primary*primary_in.m, c_sec*mdot_sec)
                C_r = C_min/C_max
                NTU = U_tot[0]*self.dA_flow/C_min
                eps[0] = (1-np.exp(-(1-C_r)*NTU))/(1-C_r*np.exp(-(1-C_r)*NTU))
            
            Q_trans[0] = eps[0]*C_min*(T_primary[0] - T_secondary[0])
            
            if x_primary[0] > 0 and x_primary[0] < 1:
                Bo = Q_trans[0]/self.dA_flow/(G_primary[0]*(hg_primary - hl_primary))
                f_primary[0] = 21500*pow(Re_primary[0], -1.14)*pow(Bo, -0.085)
            else:
                f_primary[0] = 32/Re_primary[0]
                
            dp_primary = 2*f_primary[0]*G_primary[0]**2*self.phx_inputs.L_vert/self.phx_inputs.N_element/(d_primary[0]*Dh)
            
            for jj in range(self.phx_inputs.N_element):
                h_primary[jj+1] = h_primary[jj] - Q_trans[jj]/primary_in.m
                p_primary[jj+1] = p_primary[jj] - dp_primary
                T_primary[jj+1] = PropsSI('T','H',h_primary[jj+1],'P',p_primary[jj+1],primary_in.Y)
                hl_primary = PropsSI('H','P',p_primary[jj+1],'Q',0.0, primary_in.Y)
                hg_primary =PropsSI('H','P',p_primary[jj+1],'Q',1.0, primary_in.Y)
                x_primary[jj+1] = (h_primary[jj+1]-hl_primary)/(hg_primary - hl_primary)
                
                if x_primary[jj+1] > 0 and x_primary[jj+1] < 1:
                    dl_primary = PropsSI('D','P',p_primary[jj+1],'Q',0.0, primary_in.Y)
                    dg_primary = PropsSI('D','P',p_primary[jj+1],'Q',1.0, primary_in.Y)
                    d_primary[jj+1] = 1/(x_primary[jj+1]/dg_primary+(1-x_primary[jj+1])/dl_primary)
                    vl_primary = PropsSI('V','P',p_primary[jj+1],'Q',0.0, primary_in.Y)
                    vg_primary = PropsSI('V','P',p_primary[jj+1],'Q',1.0, primary_in.Y)
                    cl_primary = PropsSI('C','P',p_primary[jj+1],'Q',0.0, primary_in.Y)
                    ll_primary = PropsSI('L','P',p_primary[jj+1],'Q',0.0, primary_in.Y)
                    
                    Xtt = pow((1-x_primary[jj+1])/x_primary[jj+1],0.9)*pow(dg_primary/dl_primary, 0.5)*pow(vl_primary/vg_primary, 0.1)
                    
                    G_primary[jj+1] = G_primary_ref*((1-Xtt)+Xtt*pow(dl_primary/dg_primary,0.5))
                    Re_primary[jj+1] = G_primary[jj+1]*Dh/vl_primary
                    pr_primary[jj+1] = cl_primary*vl_primary/ll_primary
                    htc_primary[jj+1] = 4.118*pow(Re_primary[jj+1],0.4)*pow(pr_primary[jj+1],1/3)
                    htc_primary[jj+1] = htc_primary[jj+1]*ll_primary/Dh
                else:
                    d_primary[jj+1] = PropsSI('D','T',T_primary[jj+1],'P',p_primary[jj+1],primary_in.Y)
                    v_primary = PropsSI('V','T',T_primary[jj+1],'P',p_primary[jj+1],primary_in.Y)
                    c_primary = PropsSI('C','T',T_primary[jj+1],'P',p_primary[jj+1],primary_in.Y)
                    l_primary = PropsSI('L','T',T_primary[jj+1],'P',p_primary[jj+1],primary_in.Y)
                    
                    G_primary[jj+1] = G_primary_ref
                    Re_primary[jj+1] = G_primary_ref*Dh/v_primary
                    pr_primary[jj+1] = c_primary*v_primary/l_primary
                    htc_primary[jj+1] = 0.724*pow(self.phx_inputs.beta/30,0.646)*pow(Re_primary[jj+1], 0.583)*pow(pr_primary[jj+1], 1/3)
                    htc_primary[jj+1] = htc_primary[jj+1]*l_primary/Dh
                    
                T_secondary[jj+1] = T_secondary[jj] - Q_trans[jj]/(c_sec*mdot_sec)
                G_secondary = mdot_sec/A_cx
                Re_secondary = G_secondary*Dh/v_sec
                htc_secondary = 0.724*pow(self.phx_inputs.beta/30,0.646)*pow(Re_secondary, 0.583)*pow(pr_sec, 1/3)
                htc_secondary = htc_secondary*l_sec/Dh
                
                U_tot[jj+1] = 1/(1/htc_primary[jj+1]+self.R_plate+1/htc_secondary)
                
                if x_primary[jj+1] > 0 and x_primary[jj+1] < 1:
                    C_min = c_sec*mdot_sec
                    NTU = U_tot[jj+1]*self.dA_flow/C_min
                    eps[jj+1] = 1 - np.exp(-NTU)
                else:
                    C_min = min(c_primary*primary_in.m, c_sec*mdot_sec)
                    C_max = max(c_primary*primary_in.m, c_sec*mdot_sec)
                    C_r = C_min/C_max
                    NTU = U_tot[jj+1]*self.dA_flow/C_min
                    eps[jj+1] = (1-np.exp(-(1-C_r)*NTU))/(1-C_r*np.exp(-(1-C_r)*NTU))
                
                Q_trans[jj+1] = eps[jj+1]*C_min*(T_primary[jj+1] - T_secondary[jj+1])
                if x_primary[jj+1] > 0 and x_primary[jj+1] < 1:
                    Bo = Q_trans[jj+1]/self.dA_flow/(G_primary[jj+1]*(hg_primary - hl_primary))
                    f_primary[jj+1] = 21500*pow(Re_primary[jj+1], -1.14)*pow(Bo, -0.085)
                else:
                    f_primary[jj+1] = 32/Re_primary[jj+1]
                    
                    
                dp_primary = 2*f_primary[jj+1]*G_primary[jj+1]**2*self.phx_inputs.L_vert/(self.phx_inputs.N_element-1)/(d_primary[jj+1]*Dh)
                
            T_secondary_in_cal = T_secondary[-1]
            
            if noHX == 0:
                a_cond_PHX = 0
            elif noHX == 1:
                err_T_sec = (secondary_in.T - T_secondary_in_cal)/secondary_in.T
                if err_T_sec > 0:
                    T_sec_out_lb = T_sec_out
                else:
                    T_sec_out_ub = T_sec_out
                
                if abs(err_T_sec) < 1.0e-4:
                    a_cond_PHX = 0
                elif T_sec_out_ub - T_sec_out_lb < 0.001:
                    a_cond_PHX = 0
            else:
                err_T_sec = (secondary_in.T - T_secondary_in_cal)/secondary_in.T
                if err_T_sec > 0:
                    T_sec_out_lb = mdot_sec
                else:
                    T_sec_out_ub = mdot_sec
                
                if abs(err_T_sec) < 1.0e-4:
                    a_cond_PHX = 0
                elif T_sec_out_ub - T_sec_out_lb < 0.001:
                    a_cond_PHX = 0
        
        mean_d = d_primary.mean()
        Q = Q_trans.sum()
        
        primary_out.T=T_primary[-1]
        primary_out.p=p_primary[-1]
        primary_out.Ts=PropsSI('T','P',primary_out.p,'Q',0.0,primary_out.Y)
        primary_out.h=h_primary[-1]
        primary_out.hl=hl_primary
        primary_out.hg=hg_primary
        
        DSC_cal = primary_out.Ts - primary_out.T
        secondary_in.T = T_secondary_in_cal
        secondary_in.m = mdot_sec
        secondary_out.T = T_sec_out
        secondary_out.m = mdot_sec
        
        return(primary_out, secondary_in, secondary_out, Q, mean_d, DSC_cal)
        
if __name__ == '__main__':
    from HP_dataclass import*
    inputs = PHX_Inputs()
    inputs.thk_tot = 0.411
    inputs.thk_plate = 0.6e-3
    inputs.beta = 60
    inputs.L_width = 0.6
    inputs.L_vert = 0.8
    inputs.N_plate = 187
    
    noCond = 1
    InCond = Fluid_flow(Y='Water', m=10.0, T=300, p = 101300.0)
    OutCond = Fluid_flow(Y='Water', m=10.0, p = 101300.0)
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
    
    
    
    InCond_REF = Fluid_flow(Y='R410A', m=1.65, T=333.15, p = 3.0e6)
    OutCond_REF = Fluid_flow(Y='R410A')
    InCond_REF.h = Aux_fn.PropCal(InCond_REF, 'H','T','P')
    InCond_REF.hl = PropsSI('H','P',InCond_REF.p,'Q',0.0, InCond_REF.Y )
    InCond_REF.hg = PropsSI('H','P',InCond_REF.p,'Q',1.0, InCond_REF.Y )
    cond = HX_module(hx_type='phx',Inputs=inputs)
    (OutCond_REF, InCond, OutCond, Q, mean_d, DSC_cal)=cond.PHX_cond(InCond_REF, OutCond_REF, InCond, OutCond, noCond)