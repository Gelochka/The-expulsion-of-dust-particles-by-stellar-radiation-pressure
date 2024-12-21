#!/usr/bin/python
import numpy as np 
import miepython

def reading_file(file_dust):
    lambda_m = []   
    m = []       
    with open(file_dust, "r") as file1:
        for line in file1:
            parts = line.split()
            if len(parts) < 3:
                 continue
            
            try:
                lambda_val = float(parts[0]) 
                n_val = float(parts[1])
                k_val = float(parts[2])
                 
                if file_dust == "AC_RM.ri":
                    lambda_val *= 10000  # Convert to micrometers
                
                # Specific case for SilicateDraine03.RI (adjust n)
                if file_dust == "SilicateDraine03.RI":
                    n_val += 1  # Increase n by 1
                
                # Append complex refractive index to 'm' list
                m.append(complex(n_val, -k_val))
                
                # Append wavelength to 'lambda_m' list
                lambda_m.append(lambda_val)

            except ValueError:
                # If there is an issue with conversion, log an error and continue
                #print(f"Could not parse line: {line.strip()}")
                continue  # Skip to the next line
    
    # Return the wavelength list and refractive index list
    return lambda_m, m  
 
def Qpr(m,a,lambda1):
    lambda1=np.array(lambda1)
    qpr=[]
    x=2*np.pi *a/lambda1
    qext, qsca, qback, g = miepython.mie(m,x)
    qpr_value  = qext - g*qsca 
    qpr.append(qpr_value)
    return(qpr)
     
def Qbar(qpr,T,lambda1):
    h=6.26;c=2.96; k=1.38; Tk=T/1000
    qpbar=[]
    qpr = np.array(qpr)
    lambda1 = np.array(lambda1)
 #   print(len(qpr))
 #   print(len(lambda1))
    for tt in Tk:
        integral = np.array(qpr / ((np.exp(h*c/(lambda1 *k*tt)) - 1) *  lambda1**5))
        integral = np.trapezoid(integral, lambda1)
        qprbar_val = 15* h**4 *c**4 /(np.pi**4 * k**4 *tt**4)  * integral
        qpbar.append(qprbar_val)
    return(qpbar)

 
def writea_Qbar(substance,m,lambda1):
    T = np.arange(1000,50010,10);
    for a in [0.01, 0.05, 0.1,0.5,0.8,1,5,10]:
        q_pr = Qpr(m, a, lambda1) 
        qpbar = Qbar(q_pr,  T, lambda1)
        file2 = open (str(substance)+ str(a)+".txt", "w" )
        for t, q in zip(T, qpbar):
            file2.write(f"{t} {q}\n")
    file2.close();
 
lambda_ice,m_ice =  reading_file("ICE_WAR.RI")
lambda_ac,m_ac =  reading_file("AC_RM.ri") 
print(len(lambda_ac))
print(len(m_ac))
lambda_si,m_si =  reading_file("SilicateDraine03.RI") 
#print(len(lambda_si))
writea_Qbar("ICE",m_ice,lambda_ice)
writea_Qbar("Amor",m_ac,lambda_ac)
writea_Qbar("Silicat",m_si,lambda_si)


def betta_for_star(star_type,substanse,m,lambda1):
    pho_dictionary= {'ice': 0.9, 'amor': 2,'si': 3}
    pho=pho_dictionary[substanse]
    print(pho)
    h=6.26;c=2.96; k=1.38; 
    if (star_type=="mainsequence"):
        class_spec = ["O5" ,"B0","B5","A0","A5","F0","G0","G5","K0","K5","M0","M5","M8"]
        T = np.array([40000 , 28000, 15500, 9900, 8500, 7400, 6030, 5220, 4900, 4130, 3480, 2800, 2400])
        M = np.array([1.6  , 1.25, 0.81, 0.51, 0.32, 0.23, 0.004, -0.03, -0.11,-0.16, -0.33, -0.67, -1.0])
        R = np.array([5.8 , 4.3, 2.9, 1.9, 1.3, 0.8, 0.1, -0.1, -0.4,-0.8, -1.2, -2.1, -3.1])
    if (star_type=="giant"): 
        class_spec = ["G0","G5", "K0", "K5", "M0"]
        T = np.array([5600, 5000, 4500, 3800, 3200])
        M = np.array([0.4,0.5, 0.6, 0.7, 0.8])
        R = np.array([0.8, 1, 1.2, 1.4,1.6]) #последнее из головы, его нет в книжке
    if (star_type=="supergiant"): 
        class_spec = ["B0","A0", "F0", "G0", "G5", "K0", "K5"]
        T = np.array([30000, 12000,7000, 5700, 4850, 4100, 3500]) 
        M = np.array([1.7, 1.2, 1.1, 1.0, 1.1, 1.1,1.2])
        R = np.array([1.3, 1.6, 1.8, 2.0, 2.1, 2.3, 2.6])
    R = np.exp(R*np.log(10))
    M = np.exp(M*np.log(10))
    a = np.arange(0.01,10,0.1)
    betta =[]
    
    with open(str(star_type) + "_" +str(substanse) +".txt", "w") as file_result:
        for (spec, t, mM, rR) in zip(class_spec, T, M, R):
            for a1 in a:
                tt = t/1000
                qpr = Qpr(m, a1, lambda1)
                qpr = np.array(qpr)
                lambda1 = np.array(lambda1)
                integral = np.array(qpr / ((np.exp(h*c/(lambda1 *k*tt)) - 1) *  lambda1**5))
                integral = np.trapezoid(integral, lambda1)
                qprbar_val = 15* h**4 *c**4 /(np.pi**4 * k**4 *tt**4)  * integral
               
                betta_val = 51.9 * 1e-19 * rR**2 * t**4 / mM * qprbar_val/(a1*pho)
                betta_val = np.array(betta_val)
                file_result.write(f"{a1} {betta_val[0]} {spec}\n")
                betta.append(betta_val)
    return(betta)



betta1 = betta_for_star("supergiant", 'ice',m_ice, lambda_ice)
betta1 = betta_for_star("supergiant", 'amor',m_ac, lambda_ac)
betta1 = betta_for_star("supergiant", 'si',m_si, lambda_si)


betta1 = betta_for_star("mainsequence", 'ice',m_ice, lambda_ice)
betta1 = betta_for_star("mainsequence", 'amor',m_ac, lambda_ac)
betta1 = betta_for_star("mainsequence", 'si',m_si, lambda_si)

betta1 = betta_for_star("giant", 'ice',m_ice, lambda_ice)
betta1 = betta_for_star("giant", 'amor',m_ac, lambda_ac)
betta1 = betta_for_star("giant", 'si',m_si, lambda_si)





  
