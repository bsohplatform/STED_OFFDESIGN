from scipy.special import ellipe
from math import pi
import numpy as np
from CoolProp.CoolProp import PropsSI

class HX_module:
    def __init__(self, hx_type, cor, Inputs):
        self.cor = cor
        if hx_type == 'phx':
            self.phx_inputs = Inputs
            if cor == True:
                self.crg_depth = self.phx_inputs.thk_tot/self.phx_inputs.N_plate - self.phx_inputs.thk_plate if self.phx_inputs.crg_depth == 0 else self.phx_inputs.crg_depth
                self.crg_pitch = 4*self.crg_depth if self.phx_inputs.crg_pitch == 0 else self.phx_inputs.crg_pitch
                self.enlargement = 2/pi*pow((self.crg_depth*pi/self.crg_pitch)**2+1,0.5)*ellipe((self.crg_depth*pi/self.crg_pitch)**2/((self.crg_depth*pi/self.crg_pitch)**2+1))
                self.A_cx = self.phx_inputs.N_plate*self.crg_depth*self.phx_inputs.L_width
                self.Dh = 2*self.crg_depth
                self.R_plate = self.phx_inputs.thk_plate/15.0
                self.A_flow = self.phx_inputs.L_width*self.phx_inputs.L_vert*(self.enlargement)*(2*self.phx_inputs.N_plate - 1) if self.phx_inputs.A_flow == 0 else self.phx_inputs.A_flow
                self.V = self.crg_depth*self.phx_inputs.L_width*self.phx_inputs.L_vert
                self.beta = 60.0 if self.phx_inputs.beta == 0 else self.phx_inputs.beta
            else:
                self.V = 1.0
    def PHX(self, purpose, primary_in, primary_out, secondary_in, secondary_out, noHX):
        cor = self.cor
        n = 0
        
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
        
        UA_tot = np.zeros(shape=(self.phx_inputs.N_element+1))
        eps = np.zeros(shape=(self.phx_inputs.N_element+1))
        Q_trans = np.zeros(shape=(self.phx_inputs.N_element+1))
        
        
        a_PHX = 1 
        if noHX == 0:
            T_sec_lb = secondary_out.T
            T_sec_ub = secondary_out.T if purpose == 'cond' else PropsSI('T','P',secondary_in.p,'Q',0.0,secondary_in.Y)
            c_sec = secondary_out.c
            v_sec = secondary_out.v
            l_sec = secondary_out.l
            pr_sec = secondary_out.pr
            try:
                mdot_sec = secondary_out.m
            except:
                mdot_sec = secondary_in.m
            if cor == False:
                UA = self.phx_inputs.UA*pow(mdot_sec/self.phx_inputs.mdot_nominal,0.8)
            
        elif noHX == 1:
            T_sec_lb = secondary_in.T
            T_sec_ub = PropsSI('T','P',secondary_out.p,'Q',0.0,secondary_out.Y) if purpose == 'cond' else secondary_in.T
            
            c_sec = secondary_in.c
            v_sec = secondary_in.v
            l_sec = secondary_in.l
            pr_sec = secondary_in.pr
            try:
                mdot_sec = secondary_in.m
            except:
                mdot_sec = secondary_out.m
            if cor == False:
                UA = self.phx_inputs.UA*pow(mdot_sec/self.phx_inputs.mdot_nominal,0.8)    
                
        else:
            c_sec = 0.5*(secondary_in.c + secondary_out.c)
            v_sec = 0.5*(secondary_in.v + secondary_out.v)
            l_sec = 0.5*(secondary_in.l + secondary_out.l)
            pr_sec = c_sec*v_sec/l_sec
            
            T_sec_lb = secondary_out.T if purpose == 'cond' else secondary_in.T
            T_sec_ub = secondary_out.T if purpose == 'cond' else secondary_in.T
            mdot_sec_lb = 0;
            mdot_sec_ub = 3*(primary_in.hg - primary_in.hl)*primary_in.m/c_sec/(secondary_out.T - secondary_in.T) if purpose == 'cond' else 3*(primary_out.hg - primary_out.hl)*primary_out.m/c_sec/(secondary_in.T - secondary_out.T)
            
        
        T_primary[0] = primary_in.T if purpose == 'cond' else primary_out.T
        h_primary[0] = primary_in.h if purpose == 'cond' else primary_out.h
        p_primary[0] = primary_in.p if purpose == 'cond' else primary_out.p
        hl_primary = primary_in.hl if purpose == 'cond' else primary_out.hl
        hg_primary = primary_in.hg if purpose == 'cond' else primary_out.hg
        
        x_primary[0] = (h_primary[0]-hl_primary)/(hg_primary - hl_primary)
        
        if cor == True:
            G_primary_ref = primary_in.m/self.A_cx
            if x_primary[0] > 0 and x_primary[0] < 1:
                dl_primary = PropsSI('D','P',p_primary[0],'Q',0.0, primary_in.Y)
                dg_primary = PropsSI('D','P',p_primary[0],'Q',1.0, primary_in.Y)
                d_primary[0] = 1/(x_primary[0]/dg_primary+(1-x_primary[0])/dl_primary)
                vl_primary = PropsSI('V','P',p_primary[0],'Q',0.0, primary_in.Y)
                vg_primary = PropsSI('V','P',p_primary[0],'Q',1.0, primary_in.Y)
                cl_primary = PropsSI('C','P',p_primary[0],'Q',0.0, primary_in.Y)
                ll_primary = PropsSI('L','P',p_primary[0],'Q',0.0, primary_in.Y)
        
                Xtt = pow((1-max(min(x_primary[0],1.0),0))/max(min(x_primary[0],1.0),0),0.9)*pow(dg_primary/dl_primary, 0.5)*pow(vl_primary/vg_primary, 0.1)
                
                G_primary[0] = G_primary_ref*((1-Xtt)+Xtt*pow(dl_primary/dg_primary,0.5))
                Re_primary[0] = G_primary[0]*self.Dh/vl_primary
                pr_primary[0] = cl_primary*vl_primary/ll_primary
            else:
                d_primary[0] = PropsSI('D','T',T_primary[0],'P',p_primary[0],primary_in.Y)
                v_primary = PropsSI('V','T',T_primary[0],'P',p_primary[0],primary_in.Y)
                c_primary = PropsSI('C','T',T_primary[0],'P',p_primary[0],primary_in.Y)
                l_primary = PropsSI('L','T',T_primary[0],'P',p_primary[0],primary_in.Y)
                
                G_primary[0] = G_primary_ref
                Re_primary[0] = G_primary_ref*self.Dh/v_primary
                pr_primary[0] = c_primary*v_primary/l_primary
        
        while a_PHX:
            n = n+1
            T_sec = 0.5*(T_sec_lb + T_sec_ub)
            if noHX == 2:
                mdot_sec = 0.5*(mdot_sec_lb + mdot_sec_ub)
                if cor == False:
                    UA = self.phx_inputs.UA*pow(mdot_sec/self.phx_inputs.mdot_nominal,0.8)
            
            T_secondary[0] = T_sec
            if cor == True:
                if x_primary[0] > 0 and x_primary[0] < 1:
                    if purpose == 'cond':
                        htc_primary[0] = 4.118*pow(Re_primary[0],0.4)*pow(pr_primary[0],1/3)
                        htc_primary[0] = htc_primary[0]*ll_primary/self.Dh
                    else:
                        i_primary = PropsSI('I','H',h_primary[0],'P',p_primary[0],primary_in.Y)
                        We_primary = pow(G_primary[0],2)*self.Dh/d_primary[0]/i_primary
                        Rel_primary = G_primary[0]*self.Dh/vl_primary
                        Bd_primary = 9.806*(dl_primary-dg_primary)*self.Dh**2/i_primary
                        htc_primary[0] = 1.48e3*pow(We_primary,-3.22e-2)*pow(dl_primary/dg_primary,-3.38e-1)*pow(Rel_primary, 4.51e-1)*pow(Bd_primary,-4.69e-1)
                else:
                    htc_primary[0] = 0.724*pow(self.beta/30,0.646)*pow(Re_primary[0], 0.583)*pow(pr_primary[0], 1/3)
                    htc_primary[0] = htc_primary[0]*l_primary/self.Dh1
                    
                G_secondary = mdot_sec/self.A_cx
                Re_secondary = G_secondary*self.Dh/v_sec
                htc_secondary = 0.724*pow(self.beta/30,0.646)*pow(Re_secondary, 0.583)*pow(pr_sec, 1/3)
                htc_secondary = htc_secondary*l_sec/self.Dh
                
                UA_tot[0] = 1/(1/htc_primary[0]+self.R_plate+1/htc_secondary)*self.A_flow/(self.phx_inputs.N_element-1)
            else:
                if x_primary[0] < 0 or x_primary[0] > 1:
                    c_primary = PropsSI('C','T',T_primary[0],'P',p_primary[0],primary_in.Y)
                    
                UA_tot[0] = UA/(self.phx_inputs.N_element-1)
                
            if x_primary[0] > 0 and x_primary[0] < 1:
                C_min = c_sec*mdot_sec
                NTU = UA_tot[0]/C_min
                eps[0] = 1 - np.exp(-NTU)
            else:
                C_min = min(c_primary*primary_in.m, c_sec*mdot_sec)
                C_max = max(c_primary*primary_in.m, c_sec*mdot_sec)
                C_r = C_min/C_max
                NTU = UA_tot[0]/C_min
                eps[0] = (1-np.exp(-(1-C_r)*NTU))/(1-C_r*np.exp(-(1-C_r)*NTU))
            
            Q_trans[0] = eps[0]*C_min*((T_primary[0] - T_secondary[0])) if purpose=='cond' else eps[0]*C_min*((T_secondary[0] - T_primary[0]))
            
            if cor == True:
                if x_primary[0] > 0 and x_primary[0] < 1:
                    if purpose == 'cond':    
                        Bo = Q_trans[0]/self.dA_flow/(G_primary[0]*(hg_primary - hl_primary))
                        f_primary[0] = 21500*pow(Re_primary[0], -1.14)*pow(Bo, -0.085)
                        dp_primary = 2*f_primary[0]*G_primary[0]**2*self.phx_inputs.L_vert/(self.phx_inputs.N_element-1)/(d_primary[0]*self.Dh)
                        dp_primary = 138*G_primary[0]/2/d_primary[0]
                    else:
                        dp_primary = -138*G_primary[0]/2/d_primary[0]
                    dp_primary = -138*G_primary[0]/2/d_primary[0]
                else:
                    f_primary[0] = 32/Re_primary[0]
                    dp_primary = 2*f_primary[0]*G_primary[0]**2*self.phx_inputs.L_vert/(self.phx_inputs.N_element-1)/(d_primary[0]*self.Dh)
                    dp_primary = dp_primary if purpose == 'cond' else -dp_primary
            else:
                dp_primary = self.phx_inputs.dp/(self.phx_inputs.N_element-1)
                dp_primary = dp_primary*primary_in.p if purpose == 'cond' else -dp_primary*primary_out.p
                
            for jj in range(self.phx_inputs.N_element):
                h_primary[jj+1] = h_primary[jj] - Q_trans[jj]/primary_in.m
                T_secondary[jj+1] = T_secondary[jj] - Q_trans[jj]/(c_sec*mdot_sec)
                p_primary[jj+1] = p_primary[jj] - dp_primary
                try:
                    hl_primary = PropsSI('H','P',p_primary[jj+1],'Q',0.0, primary_in.Y)
                    hg_primary = PropsSI('H','P',p_primary[jj+1],'Q',1.0, primary_in.Y)
                except:
                    hl_primary = PropsSI('H','P',p_primary[jj+1]*0.9,'Q',0.0, primary_in.Y)
                    hg_primary = PropsSI('H','P',p_primary[jj+1]*0.9,'Q',1.0, primary_in.Y)
                
                x_primary[jj+1] = (h_primary[jj+1]-hl_primary)/(hg_primary - hl_primary)
                if x_primary[jj+1] < 0 and purpose == 'evap':
                    break
                
                try:
                    T_primary[jj+1] = PropsSI('T','H',h_primary[jj+1],'P',p_primary[jj+1],primary_in.Y)
                except:
                    break
                if cor == True:
                    if x_primary[jj+1] > 0 and x_primary[jj+1] < 1:
                        try:    
                            dl_primary = PropsSI('D','P',p_primary[jj+1],'Q',0.0, primary_in.Y)
                            dg_primary = PropsSI('D','P',p_primary[jj+1],'Q',1.0, primary_in.Y)
                            vl_primary = PropsSI('V','P',p_primary[jj+1],'Q',0.0, primary_in.Y)
                            vg_primary = PropsSI('V','P',p_primary[jj+1],'Q',1.0, primary_in.Y)
                            cl_primary = PropsSI('C','P',p_primary[jj+1],'Q',0.0, primary_in.Y)
                            ll_primary = PropsSI('L','P',p_primary[jj+1],'Q',0.0, primary_in.Y)
                        except:
                            dl_primary = PropsSI('D','P',p_primary[jj+1]*0.9,'Q',0.0, primary_in.Y)
                            dg_primary = PropsSI('D','P',p_primary[jj+1]*0.9,'Q',1.0, primary_in.Y)
                            vl_primary = PropsSI('V','P',p_primary[jj+1]*0.9,'Q',0.0, primary_in.Y)
                            vg_primary = PropsSI('V','P',p_primary[jj+1]*0.9,'Q',1.0, primary_in.Y)
                            cl_primary = PropsSI('C','P',p_primary[jj+1]*0.9,'Q',0.0, primary_in.Y)
                            ll_primary = PropsSI('L','P',p_primary[jj+1]*0.9,'Q',0.0, primary_in.Y)
                        d_primary[jj+1] = 1/(x_primary[jj+1]/dg_primary+(1-x_primary[jj+1])/dl_primary)
                        Xtt = pow((1-max(min(x_primary[jj+1],1.0),0))/max(min(x_primary[jj+1],1.0),0),0.9)*pow(dg_primary/dl_primary, 0.5)*pow(vl_primary/vg_primary, 0.1)
                        
                        G_primary[jj+1] = G_primary_ref*((1-Xtt)+Xtt*pow(dl_primary/dg_primary,0.5))
                        Re_primary[jj+1] = G_primary[jj+1]*self.Dh/vl_primary
                        pr_primary[jj+1] = cl_primary*vl_primary/ll_primary
                        
                        if purpose == 'cond':
                            htc_primary[jj+1] = 4.118*pow(Re_primary[jj+1],0.4)*pow(pr_primary[jj+1],1/3)
                            htc_primary[jj+1] = htc_primary[jj+1]*ll_primary/self.Dh
                        else:
                            i_primary = PropsSI('I','H',h_primary[jj+1],'P',p_primary[jj+1],primary_in.Y)
                            We_primary = pow(G_primary[jj+1],2)*self.Dh/d_primary[jj+1]/i_primary
                            Rel_primary = G_primary[jj+1]*self.Dh/vl_primary
                            Bd_primary = 9.806*(dl_primary-dg_primary)*self.Dh**2/i_primary
                            htc_primary[jj+1] = 1.48e3*pow(We_primary,-3.22e-2)*pow(dl_primary/dg_primary,-3.38e-1)*pow(Rel_primary, 4.51e-1)*pow(Bd_primary,-4.69e-1)
                    else:
                        d_primary[jj+1] = PropsSI('D','T',T_primary[jj+1],'P',p_primary[jj+1],primary_in.Y)
                        v_primary = PropsSI('V','T',T_primary[jj+1],'P',p_primary[jj+1],primary_in.Y)
                        c_primary = PropsSI('C','T',T_primary[jj+1],'P',p_primary[jj+1],primary_in.Y)
                        l_primary = PropsSI('L','T',T_primary[jj+1],'P',p_primary[jj+1],primary_in.Y)
                        
                        G_primary[jj+1] = G_primary_ref
                        Re_primary[jj+1] = G_primary_ref*self.Dh/v_primary
                        pr_primary[jj+1] = c_primary*v_primary/l_primary
                        htc_primary[jj+1] = 0.724*pow(self.beta/30,0.646)*pow(Re_primary[jj+1], 0.583)*pow(pr_primary[jj+1], 1/3)
                        htc_primary[jj+1] = htc_primary[jj+1]*l_primary/self.Dh
                        
                    G_secondary = mdot_sec/self.A_cx
                    Re_secondary = G_secondary*self.Dh/v_sec
                    htc_secondary = 0.724*pow(self.beta/30,0.646)*pow(Re_secondary, 0.583)*pow(pr_sec, 1/3)
                    htc_secondary = htc_secondary*l_sec/self.Dh
                    
                    UA_tot[jj+1] = 1/(1/htc_primary[jj+1]+self.R_plate+1/htc_secondary)*self.A_flow/(self.phx_inputs.N_element-1)
                else:
                    if x_primary[jj+1] < 0 or x_primary[jj+1] > 1:
                        c_primary = PropsSI('C','T',T_primary[jj+1],'P',p_primary[jj+1],primary_in.Y)
                    UA_tot[jj+1] = UA/(self.phx_inputs.N_element-1)
                
                if x_primary[jj+1] > 0 and x_primary[jj+1] < 1:
                    C_min = c_sec*mdot_sec
                    NTU = UA_tot[jj+1]/C_min
                    eps[jj+1] = 1 - np.exp(-NTU)
                else:
                    C_min = min(c_primary*primary_in.m, c_sec*mdot_sec)
                    C_max = max(c_primary*primary_in.m, c_sec*mdot_sec)
                    C_r = C_min/C_max
                    NTU = UA_tot[jj+1]/C_min
                    eps[jj+1] = (1-np.exp(-(1-C_r)*NTU))/(1-C_r*np.exp(-(1-C_r)*NTU))
                
                Q_trans[jj+1] = eps[jj+1]*C_min*(T_primary[jj+1] - T_secondary[jj+1]) if purpose == 'cond' else eps[jj+1]*C_min*(T_secondary[jj+1] - T_primary[jj+1])
                if cor == True:
                    if x_primary[jj+1] > 0 and x_primary[jj+1] < 1:
                        if purpose == 'cond':
                            '''Bo = abs(Q_trans[jj+1])/self.A_flow/(self.phx_inputs.N_element-1)/(G_primary[jj+1]*(hg_primary - hl_primary))
                            f_primary[jj+1] = 21500*pow(Re_primary[jj+1], -1.14)*pow(Bo, -0.085)
                            dp_primary = 2*f_primary[jj+1]*G_primary[jj+1]**2*self.phx_inputs.L_vert/(self.phx_inputs.N_element-1)/(d_primary[jj+1]*self.Dh)'''
                            dp_primary =  138*G_primary[jj+1]/2/d_primary[jj+1]
                        else:
                            dp_primary = -138*G_primary[jj+1]/2/d_primary[jj+1]
                    else:
                        f_primary[jj+1] = 32/Re_primary[jj+1]
                        dp_primary = 2*f_primary[jj+1]*G_primary[jj+1]**2*self.phx_inputs.L_vert/(self.phx_inputs.N_element-1)/(d_primary[jj+1]*self.Dh)
                        dp_primary = dp_primary if purpose == 'cond' else -dp_primary
                else:
                    dp_primary = self.phx_inputs.dp/(self.phx_inputs.N_element-1)
                    dp_primary = dp_primary*primary_in.p if purpose == 'cond' else -dp_primary*primary_out.p
                    
            T_secondary_cal = T_secondary[-1]
            
            if noHX == 0:
                if purpose == 'cond':
                    err_T_sec = 0.0
                else:
                    if T_secondary_cal == 0.0:
                        if n < 20:
                            T_sec_ub = T_sec
                            err_T_sec = 1
                        else:
                            err_T_sec = 0
                    else:    
                        err_T_sec = (secondary_out.T - T_secondary_cal)/secondary_out.T
                        if err_T_sec > 0:
                            T_sec_lb = T_sec
                        else:
                            T_sec_ub = T_sec
            elif noHX == 1:
                if purpose == 'cond':
                    err_T_sec = (secondary_in.T - T_secondary_cal)/secondary_in.T
                    if err_T_sec > 0:
                        T_sec_lb = T_sec
                    else:
                        T_sec_ub = T_sec
                else:
                    err_T_sec = 0.0
            else:
                if purpose == 'cond':
                    err_T_sec = (secondary_in.T - T_secondary_cal)/secondary_in.T
                    if err_T_sec > 0:
                        mdot_sec_lb = mdot_sec
                    else:
                        mdot_sec_ub = mdot_sec
                else:
                    if T_secondary_cal == 0.0:
                        if n < 20:
                            mdot_sec_ub = mdot_sec
                            err_T_sec = 1
                        else:
                            err_T_sec = 0
                    else:
                        err_T_sec = (secondary_out.T - T_secondary_cal)/secondary_out.T        
                        if err_T_sec > 0:
                            mdot_sec_lb = mdot_sec
                        else:
                            mdot_sec_ub = mdot_sec
            
            if abs(err_T_sec) < 1.0e-4:
                a_PHX = 0
            else:
                if noHX == 2:
                    if mdot_sec_ub - mdot_sec_lb < 0.001:
                        a_PHX = 0
                else:
                    if T_sec_ub - T_sec_lb < 0.1:
                        a_PHX = 0
            
        
        mean_d = d_primary.mean()
        
        if purpose == 'cond':
            primary_out.T=T_primary[-1]
            primary_out.p=p_primary[-1]
            try:
                primary_out.Ts=PropsSI('T','P',primary_out.p,'Q',0.5,primary_out.Y)
            except:
                primary_out.Ts = primary_out.T
            primary_out.h=h_primary[-1]
            primary_out.hl=hl_primary
            primary_out.hg=hg_primary
        else:
            primary_in.T=T_primary[-1]
            primary_in.p=p_primary[-1]
            try:
                primary_in.Ts=PropsSI('T','P',primary_in.p,'Q',0.5,primary_in.Y)
            except:
                primary_in.Ts=primary_in.T
            primary_in.h=h_primary[-1]
            primary_in.hl=hl_primary
            primary_in.hg=hg_primary
            
        if noHX == 0:
            secondary_in.T = T_secondary[-1] if purpose == 'cond' else T_secondary[0]
        elif noHX == 1:
            secondary_out.T = T_secondary[0] if purpose == 'cond' else T_secondary[-1]
        else:
            secondary_in.m = mdot_sec
            secondary_out.m = mdot_sec
        
        Q = (h_primary[0] - h_primary[-1])*primary_in.m
        return(primary_in, primary_out, secondary_in, secondary_out, Q, mean_d)
    
if __name__ == '__main__':
    from HP_dataclass import*
    inputs = PHX_Inputs()
    inputs.thk_tot = 0.411
    inputs.thk_plate = 0.6e-3
    inputs.beta = 60
    inputs.L_width = 0.6
    inputs.L_vert = 0.8
    inputs.N_plate = 187
    inputs.cor_pitch = 0.002
    inputs.UA = 50000
    inputs.dp = 0.005
    inputs.mdot_nominal = 9.0
    
    '''
    noCond = 1
    InCond = Fluid_flow(Y='Water', m=10.0, T=300, p = 101300.0)
    OutCond = Fluid_flow(Y='Water', m=10.0, T=0, p = 101300.0)
    if noCond == 0:
        OutCond.c = Aux_fn.PropCal(OutCond, 'C', 'T', 'P') 
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
    
    InCond_REF = Fluid_flow(Y='R410A', m=0.8, T=333.15, p = 3.0e6)
    OutCond_REF = Fluid_flow(Y='R410A')
    InCond_REF.h = Aux_fn.PropCal(InCond_REF, 'H','T','P')
    InCond_REF.hl = PropsSI('H','P',InCond_REF.p,'Q',0.0, InCond_REF.Y )
    InCond_REF.hg = PropsSI('H','P',InCond_REF.p,'Q',1.0, InCond_REF.Y )
    cond = HX_module(hx_type='phx', cor=False, Inputs=inputs)
    (InCond_REF, OutCond_REF, InCond, OutCond, Q, mean_d)=cond.PHX('cond',InCond_REF, OutCond_REF, InCond, OutCond, noCond)
    '''
    
    noEvap = 2
    InEvap = Fluid_flow(Y='Water', m=0.0, T=280.15, p = 101300.0)
    OutEvap = Fluid_flow(Y='Water', m=0.0, T=278.15, p = 101300.0)
    if noEvap == 0:
        OutEvap.c = Aux_fn.PropCal(OutEvap, 'C', 'T', 'P') 
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
    
    InEvap_REF = Fluid_flow(Y='R410A')
    OutEvap_REF = Fluid_flow(Y='R410A', m=0.8, T=275.15, p = 0.7e6)
    OutEvap_REF.h = Aux_fn.PropCal(OutEvap_REF, 'H','T','P')
    OutEvap_REF.hl = PropsSI('H','P',OutEvap_REF.p,'Q',0.0, OutEvap_REF.Y )
    OutEvap_REF.hg = PropsSI('H','P',OutEvap_REF.p,'Q',1.0, OutEvap_REF.Y )
    Evap = HX_module(hx_type='phx',  cor=False, Inputs=inputs)
    (OutEvap_REF, OutEvap_REF, InEvap, OutEvap, Q, mean_d)=Evap.PHX('evap', OutEvap_REF, OutEvap_REF, InEvap, OutEvap, noEvap)
    