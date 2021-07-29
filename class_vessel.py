import numpy as np

class vessels:
    
    gpd = 11000/1440
    ndpo = 125
    area = 400
    tcff1 = [2400,2700]
    dpfactor = [0.01,1.4]
    
    # ions data......................................
    ionsrej = [0.25,0.25,0.68]
    ionsp = [0.02,0.02,0.02]
    iontcf1 = [1200,1500,3000]
    iontcf2 = [475,475,3000]
    ionfcf = [0.006,0.006,0.044]
    dpfactor = [0.01,1.4]
    aw = [23.0,35.5,61.0]
    
    @staticmethod
    def tds(fc):
        tds = 0
        for i in range(3):
            tds+=fc[i]
        return tds
    
    @classmethod
    def osp(cls,fc,T):
        C = 0
        m = np.zeros(3)
        for i in range(3):
            m[i] = fc[i]/(cls.aw[i])
            C+= m[i]
        cw = 1.000000339
        pi = 1.1266*0.001*(T+273.15)*C/cw
        return [pi,m]
    
    @staticmethod
    def tcff(T):
        if T>=25:
            tcff1 = 2400
        else:
            tcff1 = 2700
        return tcff1
    
    def balance(diff,bal_cat,bal_ani):
        tds_balance = np.ones(3)
        if diff >0:
            tds_balance[0] = tds_balance[0]*(1-(bal_cat/2))
            tds_balance[1] = tds_balance[1]*(1+(bal_ani/2))
            tds_balance[2] = tds_balance[2]*(1+(bal_ani/2))
        if diff <0:
            tds_balance[0] = tds_balance[0]*(1+(bal_cat/2))
            tds_balance[1] = tds_balance[1]*(1-(bal_ani/2)) 
            tds_balance[2] = tds_balance[2]*(1-(bal_ani/2))
        return tds_balance
    
    
    @classmethod
    def elements(cls,fc,T,r,pfa,pbp,fp,ff,n_memb):
        tcffa = np.exp(cls.tcff(T)*((1/(T+273.15))-(1/298))) 
        permflux = pfa*1440/(cls.area*n_memb)
        
       # ftds = cls.tds(fc)

        rc = np.zeros(3)
        #fc = np.zeros(3)
        pc = np.zeros(3)
        rc = np.zeros(3)
        lma = np.zeros(3)
        spions = np.zeros(3)
        tcfions = np.zeros(3)
        tcfflux = np.zeros(3)
        ptds = 0
        rtds = 0
        ospr =  cls.osp(fc,T)[0]
    
        dp = cls.dpfactor[0]*1.5*(ff**cls.dpfactor[1])
    
        for i in range(3):
            spions[i] = 1-(cls.ionsrej[i])
            
            if T<=25:
                tcfions[i] = np.exp(cls.iontcf1[i]*((1/(T+273))-(1/298)))
            else:
                tcfions[i] = np.exp(cls.iontcf2[i]*((1/(T+273))-(1/298)))
                
            tcfflux[i] = (np.exp(-cls.ionfcf[i]*permflux))/(np.exp(-cls.ionfcf[i]*17.64706))

            #print(fc)
#t-loop.........................................................................
        t = 0
        nn = 0
        while t <=2:
        
            ndp = fp - pbp - ospr
        
            if ndp<0:
                ndp = 10
    
            pf = ((cls.gpd*ndp/(tcffa*cls.ndpo)))
            rf = ff - pf
            #print(t,ospr,ndp,pf,ptds,pc)
        
            while nn<=5:
                cation_total = pc[0]*50/23
                #cation_total = pc[0]*2.18
                anion_total = pc[1]*50/35.45 + (pc[2]*50/61)
                #anion_total = pc[1]*1.41
                diff = cation_total-anion_total
                #print(diff)
                if cation_total == 0 or anion_total== 0:
                    bal_c = 0
                    bal_a = 0
                else:
                    bal_c = diff/cation_total
                    bal_a = diff/anion_total
                tds_balance_factor = cls.balance(diff,bal_c,bal_a)
                #print(tds_balance_factor)
                #print(nn,rc,pc)
            
                for i in range(3):
                    if rc[i] == 0 or rc[i] == fc[i] or rc[i]*fc[i]<0:
                        lma[i] = 0.1
                    else:
                        lma[i] = (rc[i]-fc[i])/np.log(rc[i]/fc[i])
                    pc[i] = (lma[i]*spions[i]/tcfions[i])*(tcfflux[i])*tds_balance_factor[i]
                    rc[i] = ((fc[i]*ff)-(pc[i]*pf))/rf
                    #ptds=pc[2]
                    #rtds=rc[2]
                
                nn = nn+1
                ptds = pc[0] + pc[1] + pc[2]
                rtds = rc[0] + rc[1] + rc[2]
    
            fosp = cls.osp(fc,T)
            posp = cls.osp(pc,T)
            rosp = cls.osp(rc,T)
            ospr = ((fosp[0] + rosp[0])/2) - posp[0]
            ndp = fp - pbp - ospr
            t = t+1
            rp = fp - dp
        
        return[pf,ptds,dp,rp,rtds,rf]

    
    @classmethod
    def vessel(cls,fc,pfa,T,r,pbp,n_memb,n_vess,fp):
        
        pc = np.zeros(n_memb)
        pfx = np.zeros(n_memb)
        pftotal = 0
        ff = (pfa*100/r)/n_vess

        for en in range(n_memb):
            #print(fp)
            output = cls.elements(fc,T,r,pfa,pbp,fp,ff,n_memb)
            pc[en] = output[1] 
            pfx[en] = output[0]
            fp = output[3]
            fc = output[4]
            ff = output[5]
            pftotal += output[0]
            #print(output,pftotal)
        
        return [pftotal,pc,pfx]
    
    @classmethod
    def Isolver(cls,fc,pfa,T,r,pbp,n_memb,n_vess):
        it = 50000
        n = 0
        I = 0.5
        e = 0
        X = np.zeros(it)
        Y = np.zeros(it)
        if pfa <=0.5:
            X[0] = 500
        else:
            X[0] = 50
        
        #print(n,X[0],Y[0])
        for i in range(it):
            e = e + (pfa-Y[i-1])
            X[i] = X[0] if i == 0 else X[0] +  I*e 
            Y[i] = cls.vessel(fc,pfa,T,r,pbp,n_memb,n_vess,X[i])[0]
            n = n+1
            print(n,X[i],Y[i])
            if n == it:
                output = "No solution"
                #print(n,X[i],Y[i])
            if abs(pfa-Y[i])<= 0.0001*pfa:
                output = X[i]
                break
            if n>10000:
                if abs(Y[i-1]-Y[i])<= 0.00000001:
                    output = X[i]
                    break
        return output
    
    
obj = vessels()
fc = [590.23,909.67,0.5]
pfa = 10
T = 30
r = 15
pbp = 0
n_memb = 1
n_vess = 1
fpa = obj.Isolver(fc, pfa, T, r, pbp, n_memb, n_vess)
print(obj.vessel(fc,pfa,T,r,pbp,n_memb,n_vess,fpa))
    
