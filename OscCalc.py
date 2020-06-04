# This is a slightly sleazy steal of Josh Boehm's C++ OscCalc class
# Although all bugs are almost certainly mine
# Ryan Nichol 3rd June 2020
# The original C++ header is below.
##*********************************************************************
## Welcome the the OscCalc class - hopefully your one stop shop for 
##   computing three flavor neutrino oscillations.
##
## This code works under the assumption that you want a quick solution
##  and so it stores and holds values for the sins and cosines, and 
##  other temporary values.  If you want a high precision answer you should find 
##  Mark Messier's code which was living in AtNu last I looked.
##  This code combines expansions in alpha and sinth13 to produce a pretty robust 
##  and rapid calculation.  For an evaluation of the difference between these 
##  formula and the exact solutions please see Appendix A in my thesis.  You will 
##  also find a derivation for these formula there - Josh Boehm 04-27-2009
## **************************************************************************

import math

class OscCalc:
    #Numbers pulled from 2006 PDG pg 97
    Z_A = 0.5 #average Z/A
    N_A = 6.0221415e23 #avogadro's number
    INV_CM_TO_EV = 1.97326968e-5 #convert 1/cm to eV
    INV_CM_TO_GEV = 1.97326968e-14 #convert 1/cm to GeV
    INV_KM_TO_EV = 1.97326968e-10 #convert 1/km to eV
    GF = 1.166371e-5 #Fermi constant (GeV^{-2})
    SQRT2 = math.sqrt(2) #Sqrt 2 ... why??

    def __init__(self,sinSqTheta12=0.307, sinSqTheta13=0.0218, sinSqTheta23=0.536,deltamSq21=7.53e-5,deltamSq32=2.444e-3,dcp=1.37*math.pi,density=2.7,L=811,isAntiNu=1):
        self.sinSqTheta12 = sinSqTheta12
        self.sinSqTheta13 = sinSqTheta13
        self.sinSqTheta23 = sinSqTheta23 
        self.deltamSq21 = deltamSq21 #Units??
        self.deltamSq32 = deltamSq32 #Units??
        self.deltamsq13 = deltamSq21+deltamSq32 #So this isn't necessarily correctly named?
        self.dcp = dcp
        self.density = density #Units??
        self.L = L  #Units??
        self.isAntiNu = isAntiNu #-1 for antineutrino, +1 for neutrino
        ne=self.Z_A*self.N_A*self.density
        self.elecDensity=ne*self.INV_CM_TO_EV*self.INV_CM_TO_GEV*self.INV_CM_TO_GEV
        self.fV=self.SQRT2*self.GF*self.elecDensity
        self.calcSinCos()
    
    def calcSinCos(self):
        self.sin12=math.sqrt(self.sinSqTheta12)
        self.cosSqTheta12=1-self.sinSqTheta12
        self.cos12=math.sqrt(self.cosSqTheta12)
        self.sin212=2*self.sin12*self.cos12
        self.cos212=self.cosSqTheta12-self.sinSqTheta12
        self.sinSq2Theta12=self.sin212*self.sin212
        self.sin23=math.sqrt(self.sinSqTheta23)
        self.cosSqTheta23=1-self.sinSqTheta23
        self.cos23=math.sqrt(self.cosSqTheta23)
        self.sin223=2*self.sin23*self.cos23
        self.cos223=self.cosSqTheta23-self.sinSqTheta23
        self.sinSq2Theta23=self.sin223*self.sin223
        self.sin13=math.sqrt(self.sinSqTheta13)
        self.cosSqTheta13=1-self.sinSqTheta13
        self.cos13=math.sqrt(self.cosSqTheta13)
        self.sin213=2*self.sin13*self.cos13
        self.cos213=self.cosSqTheta13-self.sinSqTheta13
        self.sinSq2Theta13=self.sin213*self.sin213

        self.sin_dcp = math.sin(self.dcp)
        self.cos_dcp = math.cos(self.dcp)
    

    def MuToElec(self,E):
        sinsq_2th12 = self.sinSq2Theta12
        sinsq_2th13 = self.sinSq2Theta13                                                                   
        cos_th23 = self.cos23
        cos_th12 = self.cos12
        sin_th13 = self.sin13
        cos_th13 = self.cos13  
                                                                              
        sin_2th23 = self.sin223
        sin_2th12 = self.sin212
        cos_2th13 = self.cos213
        cos_2th12 = self.cos212

        sinsq_th23 = self.sinSqTheta23
        sinsq_th12 = self.sinSqTheta12
        cos_dcp = self.cos_dcp
        

 
        #Building the more complicated terms
        Delta = self.deltamsq13*self.L/(4*E*1e9*self.INV_KM_TO_EV)
        A = 2*self.fV*E*1e9/(self.deltamsq13)
        alpha = self.deltamSq21/self.deltamsq13

        # A and d_cp both change sign for antineutrinos
        plusminus = int(self.isAntiNu)
        A *= plusminus
        d_cp = self.dcp*plusminus
        sin_dcp = self.sin_dcp* plusminus

        #Now calculate the resonance terms for alpha expansion (C13) and s13 expansion (C12)
        C13 = math.sqrt(sinsq_2th13+(A-cos_2th13)*(A-cos_2th13))
        C12 = 1  #really C12 -> infinity when alpha = 0 but not an option really
        if(math.fabs(alpha) > 1e-10):  #want to be careful here
            temp = cos_2th12 - A/alpha
            C12 = math.sqrt(sinsq_2th12+(temp*temp))


        #More complicated sin and cosine terms
        cosC13Delta = math.cos(C13*Delta)
        sinC13Delta = math.sin(C13*Delta)
  
        sin1pADelta = math.sin((A+1)*Delta)
        cos1pADelta = math.cos((A+1)*Delta)

        sinADelta = math.sin(A*Delta)
        sinAm1Delta = math.sin((A-1)*Delta)
        cosdpDelta = math.cos(d_cp+Delta)
        sinApam2Delta = math.sin((A+alpha-2)*Delta)
        cosApam2Delta = math.cos((A+alpha-2)*Delta)

        cosaC12Delta = 0
        sinaC12Delta = 0 
  
        if(math.fabs(alpha) > 1e-10):
            cosaC12Delta = math.cos(alpha*C12*Delta)
            sinaC12Delta = math.sin(alpha*C12*Delta)

        #First we calculate the terms for the alpha expansion (good to all orders in th13)
        # this is the equivalent of Eq 47 & 48 corrected for Mu to E instead of E to Mu

        # Leading order term 
        p1 = sinsq_th23*sinsq_2th13*sinC13Delta*sinC13Delta/(C13*C13)

        # terms that appear at order alpha
        #first work out the vacuum case since we get 0/0 otherwise.......
        p2Inner = Delta*cosC13Delta

        if(math.fabs(A) > 1e-9):
            p2Inner = Delta*cosC13Delta*(1-A*cos_2th13)/C13 -A*sinC13Delta*(cos_2th13-A)/(C13*C13)
        p2 = -2*sinsq_th12*sinsq_th23*sinsq_2th13*sinC13Delta/(C13*C13)*p2Inner*alpha


        #again working out vacuum first.....
        p3Inner = Delta* cos_th13* cos_th13*(-2*sin_dcp*sinC13Delta*sinC13Delta+2*cos_dcp*sinC13Delta*cosC13Delta)

        if(math.fabs(A) > 1e-9):
            p3Inner = (sinC13Delta/(A*C13*C13))*(- sin_dcp*(cosC13Delta - cos1pADelta)*C13+ cos_dcp*(C13*sin1pADelta - (1-A*cos_2th13)*sinC13Delta))
            p3 = sin_2th12*sin_2th23*sin_th13*p3Inner*alpha

        #  p1 + p2 + p3 is the complete contribution for this expansion
        # Now for the expansion in orders of math.sin(th13) (good to all order alpha) 
        #  this is the equivalent of Eq 65 and 66 

        # leading order term
        pa1 = 0.0
        pa2 = 0.0

        # no problems here when A -> 0
        if(math.fabs(alpha) > 1e-10):
            # leading order term
            pa1 = cos_th23*cos_th23*sinsq_2th12*sinaC12Delta*sinaC12Delta/(C12*C12)

            # and now to calculate the first order in s13 term
            t1 = (cos_2th12 - A/alpha)/C12 - alpha*A*C12*sinsq_2th12/(2*(1-alpha)*C12*C12)
            t2 = -cos_dcp*(sinApam2Delta-sinaC12Delta*t1)
            t3 = -(cosaC12Delta-cosApam2Delta)*sin_dcp
 
            denom = (1-A-alpha+A*alpha*cos_th12*cos_th12)*C12
            t4 = sin_2th12*sin_2th23*(1-alpha)*sinaC12Delta/denom

            pa2 = t4*(t3+t2)*sin_th13
  
        #pa1+pa2 is the complete contribution from this expansion

        # In order to combine the information correctly we need to add the two
        #  expansions and subtract off the terms that are in both (alpha^1, s13^1) 
        #  these may be taken from the expansion to second order in both parameters
        #  Equation 31 

        t1 = Delta*sinC13Delta*cosdpDelta
        if(math.fabs(A) > 1e-9): 
            t1 = sinADelta*cosdpDelta*sinAm1Delta/(A*(A-1))

        repeated = 2*alpha*sin_2th12*sin_2th23*sin_th13*t1

        #  Calculate the total probability
        totalP = p1+p2+p3 + (pa1+pa2) - repeated
        return totalP


    def MuToTau(self,E):
        sinsq_2th12 = self.sinSq2Theta12
        sinsq_2th13 = self.sin213*self.sin213   
        sinsq_2th23 = self.sin223*self.sin223                                                                        
        cos_th23 = self.cos23
        cos_th12 = self.cos12
        sin_th12 = self.sin12
        sin_th13 = self.sin13
        cos_th13 = self.cos13  
                                                                              
        sin_2th23 = self.sin223
        sin_2th12 = self.sin212
        cos_2th13 = self.cos213
        cos_2th23 = self.cos223
        cos_2th12 = self.cos212

        sinsq_th23 = self.sin23*self.sin23
        sinsq_th12 = self.sinSqTheta12
        #Building the more complicated terms                                                                              
        Delta = self.deltamsq13*self.L/(4*E*1e9*self.INV_KM_TO_EV)
        A = 2*self.fV*E*1e9/(self.deltamsq13)
        alpha = self.deltamSq21/self.deltamsq13
  
        # A and d_cp both change sign for antineutrinos
        plusminus = self.isAntiNu
        A *= plusminus
        d_cp = self.dcp*plusminus
        sin_dcp = self.sin_dcp* plusminus 
        cos_dcp = self.cos_dcp

        #Now calculate the resonance terms for alpha expansion (C13) and s13 expansion (C12)
        C13 = math.sqrt(sinsq_2th13+(A-cos_2th13)*(A-cos_2th13))                                                                                                                   
        C12 = 1  #really C12 -> infinity when alpha = 0 but not an option really
        if(math.fabs(alpha) > 1e-10):  #want to be careful here
            temp = cos_2th12 - A/alpha
            C12 = math.sqrt(sinsq_2th12+(temp*temp))
  
                                                                                                                       
        #More complicated sin and cosine terms
        cosC13Delta = math.cos(C13*Delta)
        sinC13Delta = math.sin(C13*Delta)

        sin1pADelta = math.sin((A+1)*Delta)
        cos1pADelta = math.cos((A+1)*Delta)

        sinADelta = math.sin((A)*Delta)

        sin1pAmCDelta = math.sin(0.5*(A+1-C13)*Delta)
        sin1pApCDelta = math.sin(0.5*(A+1+C13)*Delta)

        cosaC12Delta = 0
        sinaC12Delta = 0
                                                                                                                       
        if(math.fabs(alpha) > 1e-10):
            cosaC12Delta = math.cos(alpha*C12*Delta)
            sinaC12Delta = math.sin(alpha*C12*Delta)

        sinApam2Delta = math.sin((A+alpha-2)*Delta)
        cosApam2Delta = math.cos((A+alpha-2)*Delta)

        sinAm1Delta = math.sin((A-1)*Delta)
        cosAm1Delta = math.cos((A-1)*Delta)

        sinDelta = math.sin(Delta)
        sin2Delta = math.sin(2*Delta)

        cosaC12pApam2Delta = 0                                                                                                                       
        if(math.fabs(alpha) > 1e-10):
            cosaC12pApam2Delta = math.cos((alpha*C12+A+alpha-2)*Delta)

        #First we calculate the terms for the alpha expansion (good to all orders in th13)
        # this is the equivalent of Eq 49 & 50 corrected for Mu to E instead of E to Mu

        # Leading order term
        pmt_0 = 0.5*sinsq_2th23
        pmt_0 *= (1 - (cos_2th13-A)/C13)*sin1pAmCDelta*sin1pAmCDelta +  (1 + (cos_2th13-A)/C13)*sin1pApCDelta*sin1pApCDelta - 0.5*sinsq_2th13*sinC13Delta*sinC13Delta/(C13*C13)

        # terms that appear at order alpha
        t0 = (cos_th12*cos_th12-sin_th12*sin_th12*sin_th13*sin_th13*(1+2*sin_th13*sin_th13*A+A*A)/(C13*C13))*cosC13Delta*sin1pADelta*2
        t1 = 2*(cos_th12*cos_th12*cos_th13*cos_th13-cos_th12*cos_th12*sin_th13*sin_th13+sin_th12*sin_th12*sin_th13*sin_th13+(sin_th12*sin_th12*sin_th13*sin_th13-cos_th12*cos_th12)*A)
        t1 *= sinC13Delta*cos1pADelta/C13

        t2 =  sin_th12*sin_th12*sinsq_2th13*sinC13Delta/(C13*C13*C13)
        t2 *= A/Delta*sin1pADelta+A/Delta*(cos_2th13-A)/C13*sinC13Delta- (1-A*cos_2th13)*cosC13Delta

        pmt_1 = -0.5*sinsq_2th23*Delta*(t0+t1+t2)   

        t0 = t1 = t2 = t3 = 0.0

        t0 = cosC13Delta-cos1pADelta
        t1 = 2*cos_th13*cos_th13*math.sin(d_cp)*sinC13Delta/C13*t0
        t2 = -cos_2th23*cos_dcp*(1+A)*t0*t0

        t3  = cos_2th23*cos_dcp*(sin1pADelta+(cos_2th13-A)/C13*sinC13Delta)
        t3 *= (1+2*sin_th13*sin_th13*A + A*A)*sinC13Delta/C13 - (1+A)*sin1pADelta


        if(math.fabs(A) > 1e-9): 
            pmt_1 = pmt_1 + (t1+t2+t3)*sin_th13*sin_2th12*sin_2th23/(2*A*cos_th13*cos_th13)
        else:
            pmt_1 = pmt_1 + sin_th13*sin_2th12*sin_2th23*cos_th13*cos_th13*Delta*(2*sin_dcp*sinC13Delta*sinC13Delta+cos_dcp*cos_2th23*2*sinC13Delta*cosC13Delta)

        pmt_1 *= alpha

        #  pmt_0 + pmt_1 is the complete contribution for this expansion
                                                                                                                       
        # Now for the expansion in orders of math.sin(th13) (good to all order alpha)
        #  this is the equivalent of Eq 67 and 68
                                                                                                                       
        # leading order term
        pmt_a0 =  0.5*sinsq_2th23

        pmt_a0 *= 1 - 0.5*sinsq_2th12*sinaC12Delta*sinaC12Delta/(C12*C12)- cosaC12pApam2Delta- (1 - (cos_2th12 - A/alpha)/C12)*sinaC12Delta*sinApam2Delta    
        denom = (1-A-alpha+A*alpha*cos_th12*cos_th12)*C12

        t0 = (cosaC12Delta-cosApam2Delta)*(cosaC12Delta-cosApam2Delta)          
        t1 = (cos_2th12 - A/alpha)/C12*sinaC12Delta+sinApam2Delta          
        t2 = ((cos_2th12 - A/alpha)/C12+2*(1-alpha)/(alpha*A*C12))*sinaC12Delta+ sinApam2Delta
        t3 = (alpha*A*C12)/2.0*cos_2th23*cos_dcp*(t0 + t1*t2)
        t3 += sin_dcp*(1-alpha)*(cosaC12Delta-cosApam2Delta)*sinaC12Delta
        pmt_a1 = sin_th13*sin_2th12*sin_2th23/denom*t3

        # pmt_a1+pmt_a2 is the complete contribution from this expansion
                                                                                                                       
        # In order to combine the information correctly we need to add the two
        #  expansions and subtract off the terms that are in both (alpha^1, s13^1)
        #  and lower order terms
        #  these may be taken from the expansion to second order in both parameters
        #  Equation 34


        # Now for the term of order alpha * s13 or lower order!
        t0 = t1 = t2 = t3 = 0.0

        t1 = +sin_dcp*sinDelta*sinADelta*sinAm1Delta/(A*(A-1))
        t2 = -1/(A-1)*cos_dcp*sinDelta*(A*sinDelta-sinADelta*cosAm1Delta/A)*cos_2th23
        t0 =  2*alpha*sin_2th12*sin_2th23*sin_th13*(t1+t2)

        t1 = sinsq_2th23*sinDelta*sinDelta - alpha*sinsq_2th23*cos_th12*cos_th12*Delta*sin2Delta

        repeated = t0+t1
        #  Calculate the total probability
        totalP = pmt_0 + pmt_1 + pmt_a0 + pmt_a1 - repeated
        return totalP


    def MuToMu(self,E):
        sinsq_2th12 = self.sinSq2Theta12
        sinsq_2th13 = self.sin213*self.sin213
        sinsq_2th23 = self.sin223*self.sin223
        #  std::cout << E << "\t" << L  << "\t" << self.deltamSq32 << "\t" << self.deltamSq21 << "\t" << sinsq_2th23 << "\t" << sinsq_2th13 << "\n" 
                                          
        cos_th23 = self.cos23                                                                            
        cos_th13 = self.cos13
        cos_th12 = self.cos12
        sin_th12 = self.sin12
        sin_th13 = self.sin13

        #  std::cout << "th23: " <<  cos_th23 << "\t" << sinsq_2th23 << "\n"

  
        sin_2th23 = self.sin223
        #  sin_2th13 = self.sin213
        sin_2th12 = self.sin212
                                                                                
        cos_2th23 = self.cos223
        cos_2th13 = self.cos213
        cos_2th12 = self.cos212


        sinsq_th23 = self.sin23*self.sin23
        sinsq_th12 = self.sinSqTheta12

        d_cp = self.dcp
        cos_dcp = self.cos_dcp
        sin_dcp = self.sin_dcp

        #Building the more complicated terms                                                                              
        Delta = self.deltamsq13*self.L/(4*E*1e9*self.INV_KM_TO_EV)
        A = 2*self.fV*E*1e9/(self.deltamsq13)
        alpha = self.deltamSq21/self.deltamsq13
  
        # A and d_cp both change sign for antineutrinos  
        plusminus = self.isAntiNu
        A *= plusminus
        d_cp *= plusminus
        sin_dcp *= plusminus
  
        #Now calculate the resonance terms for alpha expansion (C13) and s13 expansion (C12)  
        C13 = math.sqrt(sinsq_2th13+(A-cos_2th13)*(A-cos_2th13))
        C12 = 1  #really C12 -> infinity when alpha = 0 but not an option really
        if(math.fabs(alpha) > 1e-10):  #want to be careful here
            temp = cos_2th12 - A/alpha
            C12 = math.sqrt(sinsq_2th12+(temp*temp))
  
                                                                                                                       
        #More complicated sin and cosine terms
        cosC13Delta = math.cos(C13*Delta)
        sinC13Delta = math.sin(C13*Delta)

  
        sin1pADelta = math.sin((A+1)*Delta)
        cos1pADelta = math.cos((A+1)*Delta)

        sinADelta = math.sin((A)*Delta)

        sin1pAmCDelta = math.sin(0.5*(A+1-C13)*Delta)
        sin1pApCDelta = math.sin(0.5*(A+1+C13)*Delta)

  
        cosaC12Delta = 0
        sinaC12Delta = 0
                                                                                                                       
        if(math.fabs(alpha) > 1e-10):
            cosaC12Delta = math.cos(alpha*C12*Delta)
            sinaC12Delta = math.sin(alpha*C12*Delta)
        # otherwise not relevant

        sinApam2Delta = math.sin((A+alpha-2)*Delta)
        cosApam2Delta = math.cos((A+alpha-2)*Delta)

        sinAm1Delta = math.sin((A-1)*Delta)
        cosAm1Delta = math.cos((A-1)*Delta)

        sinDelta = math.sin(Delta)
        sin2Delta = math.sin(2*Delta)
        cosdpDelta = math.cos(d_cp+Delta)

        cosaC12pApam2Delta = 0                                                                                                                       
        if(math.fabs(alpha) > 1e-10):
            cosaC12pApam2Delta = math.cos((alpha*C12+A+alpha-2)*Delta)
  

        #This bit is the mu-to-tau part
        #First we calculate the terms for the alpha expansion (good to all orders in th13)
        # this is the equivalent of Eq 49 & 50 corrected for Mu to E instead of E to Mu  
        pMuToTau=0
  
        # Leading order term
        pmt_0 = 0.5*sinsq_2th23
        pmt_0 *= (1 - (cos_2th13-A)/C13)*sin1pAmCDelta*sin1pAmCDelta +  (1 + (cos_2th13-A)/C13)*sin1pApCDelta*sin1pApCDelta- 0.5*sinsq_2th13*sinC13Delta*sinC13Delta/(C13*C13)

        #    if(E>1.6 && E<1.606)
    
        # terms that appear at order alpha
        t0 = (cos_th12*cos_th12-sin_th12*sin_th12*sin_th13*sin_th13*(1+2*sin_th13*sin_th13*A+A*A)/(C13*C13))*cosC13Delta*sin1pADelta*2
        t1 = 2*(cos_th12*cos_th12*cos_th13*cos_th13-cos_th12*cos_th12*sin_th13*sin_th13+sin_th12*sin_th12*sin_th13*sin_th13+(sin_th12*sin_th12*sin_th13*sin_th13-cos_th12*cos_th12)*A)
        t1 *= sinC13Delta*cos1pADelta/C13

        t2 =  sin_th12*sin_th12*sinsq_2th13*sinC13Delta/(C13*C13*C13)
        t2 *= A/Delta*sin1pADelta+A/Delta*(cos_2th13-A)/C13*sinC13Delta- (1-A*cos_2th13)*cosC13Delta

        pmt_1 = -0.5*sinsq_2th23*Delta*(t0+t1+t2)   

        t0 = t1 = t2 = t3 = 0.0

        t0 = cosC13Delta-cos1pADelta
        t1 = 2*cos_th13*cos_th13*sin_dcp*sinC13Delta/C13*t0
        t2 = -cos_2th23*cos_dcp*(1+A)*t0*t0

        t3  = cos_2th23*cos_dcp*(sin1pADelta+(cos_2th13-A)/C13*sinC13Delta)
        t3 *= (1+2*sin_th13*sin_th13*A + A*A)*sinC13Delta/C13 - (1+A)*sin1pADelta


        if(math.fabs(A) > 1e-9): 
            pmt_1 = pmt_1 + (t1+t2+t3)*sin_th13*sin_2th12*sin_2th23/(2*A*cos_th13*cos_th13)
        else:
            pmt_1 = pmt_1 + sin_th13*sin_2th12*sin_2th23*cos_th13*cos_th13*Delta*(2*sin_dcp*sinC13Delta*sinC13Delta+cos_dcp*cos_2th23*2*sinC13Delta*cosC13Delta)

        pmt_1 *= alpha

        #  pmt_0 + pmt_1 is the complete contribution for this expansion
                                                                                                                       
        # Now for the expansion in orders of math.sin(th13) (good to all order alpha)
        #  this is the equivalent of Eq 67 and 68
                                                                                                                       
        # leading order term
        pmt_a0 =  0.5*sinsq_2th23

        pmt_a0 *= 1 - 0.5*sinsq_2th12*sinaC12Delta*sinaC12Delta/(C12*C12)- cosaC12pApam2Delta- (1 - (cos_2th12 - A/alpha)/C12)*sinaC12Delta*sinApam2Delta
            
        denom = (1-A-alpha+A*alpha*cos_th12*cos_th12)*C12

        t0 = (cosaC12Delta-cosApam2Delta)*(cosaC12Delta-cosApam2Delta)         
        t1 = (cos_2th12 - A/alpha)/C12*sinaC12Delta+sinApam2Delta   
        t2 = ((cos_2th12 - A/alpha)/C12+2*(1-alpha)/(alpha*A*C12))*sinaC12Delta+ sinApam2Delta
        t3 = (alpha*A*C12)/2.0*cos_2th23*cos_dcp*(t0 + t1*t2)
        t3 += sin_dcp*(1-alpha)*(cosaC12Delta-cosApam2Delta)*sinaC12Delta

        pmt_a1 = sin_th13*sin_2th12*sin_2th23/denom*t3

        # pmt_a1+pmt_a2 is the complete contribution from this expansion
                                                                                                                       
        # In order to combine the information correctly we need to add the two
        #  expansions and subtract off the terms that are in both (alpha^1, s13^1)
        #  and lower order terms
        #  these may be taken from the expansion to second order in both parameters
        #  Equation 34


        # Now for the term of order alpha * s13 or lower order!
        t0 = t1 = t2 = t3 = 0.0

        t1 = +sin_dcp*sinDelta*sinADelta*sinAm1Delta/(A*(A-1))
        t2 = -1/(A-1)*cos_dcp*sinDelta*(A*sinDelta-sinADelta*cosAm1Delta/A)*cos_2th23
        t0 =  2*alpha*sin_2th12*sin_2th23*sin_th13*(t1+t2)

        t1 = sinsq_2th23*sinDelta*sinDelta - alpha*sinsq_2th23*cos_th12*cos_th12*Delta*sin2Delta

        repeated = t0+t1

        #  Calculate the total probability
        pMuToTau = pmt_0 + pmt_1 + pmt_a0 + pmt_a1 - repeated
   


        pMuToElec=0
        #Now for the MuToElec part
  
        # this is the equivalent of Eq 47 & 48 corrected for Mu to E instead of E to Mu

        # Leading order term 
        p1 = sinsq_th23*sinsq_2th13*sinC13Delta*sinC13Delta/(C13*C13)

        # terms that appear at order alpha
        #first work out the vacuum case since we get 0/0 otherwise.......
        p2Inner = Delta*cosC13Delta

        if(math.fabs(A) > 1e-9):
            p2Inner = Delta*cosC13Delta*(1-A*cos_2th13)/C13 -A*sinC13Delta*(cos_2th13-A)/(C13*C13)

        p2 = -2*sinsq_th12*sinsq_th23*sinsq_2th13*sinC13Delta/(C13*C13)*p2Inner*alpha


        #again working out vacuum first.....
        p3Inner = Delta* cos_th13* cos_th13*(-2*sin_dcp*sinC13Delta*sinC13Delta+2*cos_dcp*sinC13Delta*cosC13Delta)

        if(math.fabs(A) > 1e-9):
            p3Inner = (sinC13Delta/(A*C13*C13))*(- sin_dcp*(cosC13Delta - cos1pADelta)*C13+ cos_dcp*(C13*sin1pADelta - (1-A*cos_2th13)*sinC13Delta))

        p3 = sin_2th12*sin_2th23*sin_th13*p3Inner*alpha

        #  p1 + p2 + p3 is the complete contribution for this expansion
  
        # Now for the expansion in orders of math.sin(th13) (good to all order alpha) 
        #  this is the equivalent of Eq 65 and 66

        # leading order term
        pa1 = 0.0 
        pa2 = 0.0

        # no problems here when A -> 0
        if(math.fabs(alpha) > 1e-10):
            # leading order term
            pa1 = cos_th23*cos_th23*sinsq_2th12*sinaC12Delta*sinaC12Delta/(C12*C12)

            # and now to calculate the first order in s13 term
            t1 = (cos_2th12 - A/alpha)/C12- alpha*A*C12*sinsq_2th12/(2*(1-alpha)*C12*C12)
            t2 = -cos_dcp*(sinApam2Delta-sinaC12Delta*t1)
            t3 = -(cosaC12Delta-cosApam2Delta)*sin_dcp
            denom = (1-A-alpha+A*alpha*cos_th12*cos_th12)*C12
            t4 = sin_2th12*sin_2th23*(1-alpha)*sinaC12Delta/denom
            pa2 = t4*(t3+t2)*sin_th13
    
        #pa1+pa2 is the complete contribution from this expansion

        # In order to combine the information correctly we need to add the two
        #  expansions and subtract off the terms that are in both (alpha^1, s13^1) 
        #  these may be taken from the expansion to second order in both parameters
        #  Equation 31 

        t1 = Delta*sinC13Delta*cosdpDelta
        if(math.fabs(A) > 1e-9):
            t1 = sinADelta*cosdpDelta*sinAm1Delta/(A*(A-1))

        repeated = 2*alpha*sin_2th12*sin_2th23*sin_th13*t1

        #  Calculate the total probability
        pMuToElec = p1+p2+p3 + (pa1+pa2) - repeated

        p1 = 1. - pMuToTau - pMuToElec
        if(p1 < 1e-6): 
            #  cout<<"P(mu->mu) less than zero Damnation! "<<x<<" "<<p1<<endl 
            p1 = 0
        return p1


    def ElecToTau(self,E):
        #  EtoTau is the same as E->Mu wich sinsq_th23 <-> cossq_th23, math.sin(2th23) <->-sin(2th23)
        origCos = self.cos23
        origSin = self.sin23
        orig2Sin = self.sin223

        self.cos23 = origSin
        self.sin23 = origCos
        self.sin223 = -orig2Sin

        prob = ElecToMu(E)

        #restore the world
        self.cos23 = origCos
        self.sin23 = origSin
        self.sin223 = orig2Sin
        return prob


    def ElecToMu(self, E):
        # Flip delta to reverse direction
        oldSinDelta = self.sin_dcp
        oldDelta =  self.dcp

        self.dcp = -oldDelta
        self.sin_dcp = -oldSinDelta
                                                                                                                       
        prob = self.MuToElec(E)
                                                                                                                       
        #restore the world
        self.dcp = oldDelta
        self.sin_dcp = oldSinDelta
        return prob   


    def ElecToElec(self,E):
        sinsq_2th12 = self.sinSq2Theta12
        sinsq_2th13 = self.sin213*self.sin213
                                                                                                                       
        sin_th12 = self.sin12
 
        #  cos_2th23 = self.cos223
        cos_2th13 = self.cos213
        cos_2th12 = self.cos212
                                                                                                                       
        d_cp = self.dcp
        sin_dcp = self.sin_dcp
                                                                                                                       
        #Building the more complicated terms
        Delta = self.deltamsq13*self.L/(4*E*1e9*self.INV_KM_TO_EV)
        A = 2*self.fV*E*1e9/(self.deltamsq13)
        alpha = self.deltamSq21/self.deltamsq13
                                                                                                                       
        # A and d_cp both change sign for antineutrinos
        plusminus = self.isAntiNu
        A *= plusminus
        d_cp *= plusminus
        sin_dcp *= plusminus
                                                                                                                       
        #Now calculate the resonance terms for alpha expansion (C13) and s13 expansion (C12)
        C13 = math.sqrt(sinsq_2th13+(A-cos_2th13)*(A-cos_2th13))
                                                                                                                       
        C12 = 1  #really C12 -> infinity when alpha = 0 but not an option really
        if(math.fabs(alpha) > 1e-10):  #want to be careful here
            temp = cos_2th12 - A/alpha
            C12 = math.sqrt(sinsq_2th12+(temp*temp))
                                                                                                                       
        #More complicated sin and cosine terms
        cosC13Delta=math.cos(C13*Delta)
        sinC13Delta=math.sin(C13*Delta)


        cosaC12Delta = 0
        sinaC12Delta = 0
                                                                                                                       
        if(math.fabs(alpha) > 1e-10):
            cosaC12Delta = math.cos(alpha*C12*Delta)
            sinaC12Delta = math.sin(alpha*C12*Delta)
    
                                                                                                                       
        #First we calculate the terms for the alpha expansion (good to all orders in th13)
        # this is the equivalent of Eq 45 & 46 corrected for Mu to E instead of E to Mu
                                                                                                                       
        # Leading order term
        p1 = 1 - sinsq_2th13*sinC13Delta*sinC13Delta/(C13*C13)
                                                                                                                       
        # terms that appear at order alpha
        p2Inner = Delta*cosC13Delta*(1-A*cos_2th13)/C13 -A*sinC13Delta*(cos_2th13-A)/(C13*C13)
                                                                                                                       
        p2 = +2*sin_th12*sin_th12*sinsq_2th13*sinC13Delta/(C13*C13)*p2Inner*alpha
        #  p1 + p2 is the complete contribution for this expansion
                                                                                                                       
        # Now for the expansion in orders of math.sin(th13) (good to all order alpha)
        #  this is the equivalent of Eq 63 and 64
                                                                                                                       
        # leading order term
        pa1 = 1.0
        pa2 = 0.0
                                                                                                                       
        if(math.fabs(alpha) > 1e-10):
            # leading order term
            pa1 = 1.0 - sinsq_2th12*sinaC12Delta*sinaC12Delta/(C12*C12)
  
        #pa1 is the complete contribution from this expansion, there is no order s13^1 term
                                                                                                                       
        # In order to combine the information correctly we need to add the two
        #  expansions and subtract off the terms that are in both (alpha^1, s13^1)
        #  these may be taken from the expansion to second order in both parameters
        #  Equation 30
                                                                                                                       
        repeated = 1
                                                                                                                       
        #  Calculate the total probability
        totalP = p1+p2 + (pa1+pa2) - repeated
                                                                                                                       
        return totalP


    def TauToTau(self,E):
        #  TautoTau is the same as Mu->Mu wich sinsq_th23 <-> cossq_th23, sin(2th23) <->-sin(2th23)
        origCos = self.cos23
        origSin = self.sin23
        orig2Sin = self.sin223
        origCosSq=self.cosSqTheta23
        origSinSq=self.sinSqTheta23
        

        self.cos23 = origSin
        self.sin23 = origCos
        self.cosSqTheta23=origSinSq
        self.sinSqTheta23= origCosSq
        self.sin223 = -orig2Sin

        prob = self.MuToMu(E)
  
        #restore the world
        self.cos23 = origCos
        self.sin23 = origSin
        self.sin223 = orig2Sin
        self.cosSqTheta23=origCosSq
        self.sinSqTheta23= origSinSq
  
        return prob

  
    def TauToMu(self,E):
        # Flip delta to reverse direction
        oldSinDelta = self.sin_dcp
        oldDelta =  self.dcp
  
        self.dcp = -oldDelta
        self.sin_dcp = -oldSinDelta
                              
        prob = MuToTau(E)
                                                                                      
        #restore the world
        self.dcp = oldDelta
        self.sin_dcp = oldSinDelta
  
        return prob

  
    def TauToElec(self, E):
        # Flip delta to reverse direction
        oldSinDelta = self.sin_dcp
        oldDelta =  self.dcp
   
        self.dcp = -oldDelta
        self.sin_dcp = -oldSinDelta 
                               
        prob = ElecToTau(E)
                                 
        #restore the world 
        self.dcp = oldDelta
        self.sin_dcp = oldSinDelta 
        return prob 

if __name__ == '__main__':
    calcy = OscCalc()
    print(calcy.MuToMu(2))    
