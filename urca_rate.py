import numpy as np

#Calculate net urca rate over chempo_delta.
#with electron contribution only!!!

G_Fermi = 1.16637e-11  # MeV^-2
cabbibo_angle = 13.04
G = G_Fermi * np.cos(np.radians(cabbibo_angle))
fN_pi = 1
cA = 1.26
mpi0 = 134.977  # MeV
def mU_spec_net_rate(mp,mn,kfp,kfn,kfe,kfmu,T):
    
    theta_n=np.where(kfn >= kfp + kfe,1,1 - 3/8 * (kfp + kfe - kfn)**2 / (kfp * kfe))
    
    epsilon_Fn = np.sqrt(kfn**2 + mn**2)
    epsilon_Fp = np.sqrt(kfp**2 + mp**2)
    
    rate_n = (
        1 / (5760 * np.pi**9) * G**2 * cA**2 * fN_pi**4 * epsilon_Fn**3 * 
        epsilon_Fp / mpi0**4 * kfn**4 * kfp / (kfn**2 + mpi0**2)**2 * 
        theta_n * 1835 * np.pi**6 * T**6
    )
    
    condition0=kfn > 3 * kfp + kfe
    theta1=(3 * kfp + kfe - kfn)**2 / (kfn * kfe)
    condition1=np.logical_and(3 * kfp + kfe > kfn,kfn > 3 * kfp - kfe)
    theta2=4 * (3 * kfp - kfn) / kfn
    condition2=np.logical_and(3 * kfp - kfe > kfn,kfn > kfp + kfe)
    theta3=2 + 3 * (2 * kfp - kfn) / kfe - 3 * (kfp - kfe)**2 / (kfn * kfe)
    condition3=kfn < kfp + kfe
    
    theta_p=np.where(condition0,0,np.where(condition1,theta1,np.where(condition2,theta2,np.where(condition3,theta3,1))))
    
    rate_p = (
        1 / (40320 * np.pi**9) * G**2 * cA**2 * fN_pi**4 *
        epsilon_Fp**3 * epsilon_Fn / mpi0**4 *
        (kfn - kfp)**4 * kfn / ((kfn - kfp)**2 + mpi0**2)**2 *
        theta_p  * 1835 * np.pi**6 * T**6
    )
    
    theta_dUrca=np.where(kfn<kfp + kfe,1,0)
    rate_dUrca = 1 / (240 * np.pi**5) * G**2 * (1+3*cA**2) * epsilon_Fn*epsilon_Fp*kfe*theta_dUrca*17*np.pi**4 * T**4

    return rate_n+rate_p+rate_dUrca