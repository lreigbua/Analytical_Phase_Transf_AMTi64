# Microstructure evolution model that calculates phase fractions for a given thermal history:
import numpy as np
import matplotlib.pyplot as plt
import math as m
from pathlib import Path
import os
import sys
import json

def Run(input_csv, data_folder, create_output_json = True, create_plots = True, create_output_CSVs = True):

    user_path = Path.cwd()

    data = np.genfromtxt(input_csv, delimiter=',')

    module_path = Path(__file__).resolve().parent

    os.chdir(module_path)



    # read temperature history:
    # data = np.genfromtxt('../input/amp_temp_amprin_and_HT850.csv', delimiter=',')
    # data = np.genfromtxt('../input/Abaqus_one_elem_tests/amp_temp_amprint.txt', delimiter=',')
    # data = np.genfromtxt('../input/Abaqus_one_elem_tests/one_elem.csv', delimiter=',')
    # data = np.genfromtxt('../input/cooling.csv', delimiter=',')
    # data = np.genfromtxt('../input/XRD_test.csv', delimiter=',')
    # data = np.genfromtxt('../input/half_hatch_rehatched/Temperature_Element_at_heigt_1.98.csv', delimiter=',')
    # data = np.genfromtxt('../input/output_6_scans_half_power/Temperature_Element_at_heigt_0.06.csv', delimiter=',')
    # data = np.genfromtxt('../input/240_6_rescans/Temperature_Element_at_heigt_8.28.csv', delimiter=',')
    # data = np.genfromtxt('../input/amp_temp_amprin_and_HT900.csv', delimiter=',')

    #Final Inputs
    # data = np.genfromtxt('../Final_Inputs/C1-B1/Temperature_Element_at_heigt_8.28.csv', delimiter=',')
    # data = np.genfromtxt('../Final_Inputs/C2-B2/Temperature_Element_at_heigt_8.28.csv', delimiter=',')
    # data = np.genfromtxt('../Final_Inputs/C3-B2/Temperature_Element_at_heigt_8.28.csv', delimiter=',')
    # data = np.genfromtxt('../Final_Inputs/C4-B3/Temperature_Element_at_heigt_8.28.csv', delimiter=',')
    # data = np.genfromtxt('../Final_Inputs/C5-B3/Temperature_Element_at_heigt_1.98.csv', delimiter=',')
    # data = np.genfromtxt('../Final_Inputs/Cyl_t8s/Temperature_Element_at_heigt_53.04.csv', delimiter=',')
    # data = np.genfromtxt('../Final_Inputs/C2-HT800/HT1.csv', delimiter=';')
    # data = np.genfromtxt('../Final_Inputs/C5-HT850/HT2.csv', delimiter=';')




    t_input = data[:,0] # time
    T_input = data[:,1] # temperature

    t = np.array(t_input[0])

    # if a given time step is larger than min_t_step, we interpolate between the two points with min_t_step length time steps,
    # otherwise we interpolate with a time step that is 1/3 of the given time step:
    min_t_step = 10
    for i in range(1, len(t_input)):
        if t_input[i] - t_input[i-1] > min_t_step:
            t = np.append(t, np.arange(t_input[i-1] + min_t_step, t_input[i], min_t_step))
        else:
            min_t_step_c = (t_input[i]-t_input[i-1])/50 + 0.0001
            t = np.append(t, np.arange(t_input[i-1] + min_t_step_c, t_input[i], min_t_step_c))

        t = np.append(t, t_input[i])

    T = np.interp(t, t_input, T_input)

        
    # print('Time:', t)
    # print('Temperature:', T)
    # sys.exit()

    #Discretize the given temperature history with smaller time steps using interpolation:
    # time_step = 0.001
    # t = np.arange(t[0], t[-1], time_step)
    # T = np.interp(t, data[:,0], data[:,1])

    # initialize arrays:
    beta_f = np.zeros(len(t))
    alpha_f = np.zeros(len(t))
    mart_f = np.zeros(len(t))
    RLS = np.zeros(len(t))
    t_lath = np.zeros(len(t))

    #Initial values:
    # beta_f[0] = 0.01
    # alpha_f[0] = 0.0
    # mart_f[0] = 0.99
    # RLS[0] = 1
    # t_lath[0] = 1

    beta_f[0] = 0.0
    alpha_f[0] = 0.0
    mart_f[0] = 0.0
    RLS[0] = -1
    t_lath[0] = 0


    #for martensitic transf
    mart_trans_started_flag = False
    f_p_T0 = 0
    f_c_T0 = 1


    #lath thickness parameters:
    k = 0.45
    R=294
    # R = 350

    Q = 97000
    k0 = 19
    c = 3
    t_lath_0 = 0.25
    Rg = 8.314


    ### Melting temperature:
    T_melt = 1610

    #### Beta to Alpha transformation data:
    # JMAK parameters:
    data = np.genfromtxt('./JMAK_params_BtoA.csv', delimiter=',', skip_header=1)
    k_BtoA = data[:,0]
    n_BtoA = data[:,1]
    temp_jmak_BtoA = data[:,2]

    # Alpha equilibrium phase fraction
    data = np.genfromtxt('./Equilibrium_alpha_fraction.csv', delimiter=',', skip_header=1)
    T_alpha_eq = data[:,1]
    alpha_eq = data[:,0]

    # KM parameters:
    b_km = 0.015
    T_mart = 851

    # Mart to alpha + beta transformation data:
    T_mart_dis = 400
    data = np.genfromtxt('./JMAK_params_mart_to_AandB.csv', delimiter=',', skip_header=1)
    k_mart_to_AandB = data[:,0]
    n_mart_to_AandB = data[:,1]
    temp_jmak_mart_to_AandB = data[:,2]

    data = np.genfromtxt('./Equilibrium_mart_fraction.csv', delimiter=',', skip_header=1)
    T_mart_eq = data[:,1]
    mart_eq = data[:,0]

    for i in range(1, len(t)): #iterate over all time steps

        dt = t[i]-t[i-1] # time increment
        dT = T[i]-T[i-1] # temperature increment

        if dt == 0:
                alpha_f[i] = alpha_f[i-1]
                beta_f[i] = beta_f[i-1]
                mart_f[i] = mart_f[i-1]
                t_lath[i] = t_lath[i-1]
                RLS[i] = RLS[i-1]
                continue


        # print('Time:', t[i])
        # print('Temperature:', T[i])


    #First asses RLS state:  
        # if it is powder:
        if RLS[i-1] == -1:
            if T[i] <= T_melt: #if next temperature is below melting temperature:
                RLS[i] = -1 # keep powder state
                continue # still powder, nothing happens, go to next iteration
            else: # if temperature is above melting temperature:
                RLS[i] = 0 # change to liquid state
                continue # no transformation happens, go to next iteration
        
        # if it is liquid:
        if RLS[i-1] == 0:
            if T[i] >= T_melt: #if above T_melt:
                RLS[i] = 0 # keep liquid state
                continue # still liquid, nothing happens, go to next iteration
            else: # if below T_melt:
                RLS[i] = 1 # change to solid state
                beta_f[i] = 1.0 # set beta phase fraction to 1
                alpha_f[i] = 0.0 # set alpha phase fraction to 0
                mart_f[i] = 0.0 # set martensite phase fraction to 0
                t_lath[i] = t_lath_0 # set lath thickness to t_lath_0
                continue
        
        if RLS[i-1] == 1: # if it is solid:
            if T[i] >= T_melt:
                RLS[i] = 0 # change to liquid state
                beta_f[i] = 0.0 # set phase fractions to zero
                alpha_f[i] = 0.0
                mart_f[i] = 0.0
                t_lath[i] = 0.0
                continue # no transformation happens, go to next iteration

            else:
                RLS[i] = 1 # keep solid state

                # if it is solid, we calculate phase transformations depending on heating:

                #data for some transformations:
                f_mart_eq = np.interp(T[i], T_mart_eq, mart_eq)
                # f_mart_eq = ( 1 - m.exp( -b_km*(T_mart - T[i]) ) ) * (beta_f[i-1]+mart_f[i-1])
                beta_eq_i = 1 - np.interp(T[i], T_alpha_eq, alpha_eq)
                alpha_eq_i = np.interp(T[i], T_alpha_eq, alpha_eq)


                #in case there is no phase transformation, we keep the phase fractions the same:
                alpha_f[i] = alpha_f[i-1]
                beta_f[i] = beta_f[i-1]
                mart_f[i] = mart_f[i-1]
                t_lath[i] = t_lath[i-1]

                

                if dT/dt > -20.0 and dT/dt < 0 and beta_f[i-1] > (beta_eq_i*(beta_f[i-1]+alpha_f[i-1])) + 0.00001: # JMAK model beta to alpha

                    # if T[i] > 480 and T[i] < 500:
                    #     print('Time:', t[i-1])
                    #     print('Temperature:', T[i-1])
                    #     print('beta_f[i-1]:', alpha_f[i-1])                
                    # if T[i] < 300:
                    #     continue

                    # f_tot = beta_f[i-1] + alpha_f[i-1]

                    alpha_eq_i = np.interp(T[i], T_alpha_eq, alpha_eq)

                    # linearly interpolate JMAK parameters k_i and n_i according to current temperature:
                    k_i = np.interp(T[i], temp_jmak_BtoA, k_BtoA)
                    n_i = np.interp(T[i], temp_jmak_BtoA, n_BtoA)

                    #Calculate Tau:
                    tau = ( -1/k_i*m.log(1 - (alpha_f[i-1]/alpha_eq_i)/(beta_f[i-1]+alpha_f[i-1]) ))**(1/n_i)
                    # tau = ( -1/k * m.log( (beta_f[i-1] - beta_eq_i * f_tot) / (f_tot * (1 - beta_eq_i)) ) ) ** 1/n_i

                    # print('tau:',tau)

                    # caculate change in alpha phase fraction:
                    alpha_f[i] = (1 - m.exp(-k_i * (tau+dt)**n_i)) * alpha_eq_i * (beta_f[i-1]+alpha_f[i-1])
                    # beta_f[i] = f_tot*(1 - (1 - m.exp(-k_i * (tau+dt)**n_i)) * (1 - beta_eq_i))

                    # calculate change in beta phase fraction:
                    
                    # alpha_f[i] = f_tot - beta_f[i]
                    beta_f[i] = beta_f[i-1] - (alpha_f[i] - alpha_f[i-1])
                    mart_f[i] = mart_f[i-1]

                ## Beta to Martensite Mode 1
                # elif dT/dt <= -410.0 and beta_f[i-1] > beta_eq_i and T[i] < T_mart: # KM model beta to martensite fast cooling:
                                
                #     mart_f[i] = ( 1 - m.exp( -b_km*(T_mart - T[i]) ) ) * (beta_f[i-1]+mart_f[i-1])
                #     beta_f[i] = beta_f[i-1] - (mart_f[i] - mart_f[i-1])
                #     alpha_f[i] = alpha_f[i-1]

                # elif dT/dt > -410.0 and dT/dt <= -20.0 and beta_f[i-1] > beta_eq_i and T[i] < T_mart: # KM model beta to martensite slow cooling:

                #     mart_f[i] = ( 1 - m.exp( -b_km*(T_mart - T[i]) ) ) * (beta_f[i-1]+mart_f[i-1]-beta_eq_i)
                #     beta_f[i] = beta_f[i-1] - (mart_f[i] - mart_f[i-1])
                #     alpha_f[i] = alpha_f[i-1]


                ## Beta to Martensite Zhang:
                if dT/dt <=-20 and T[i] < T_mart :

                    if mart_trans_started_flag == False: #if start of cooling cycle
                        mart_trans_started_flag = True

                        f_p_T0 = beta_f[i-1]
                        f_c_T0 = mart_f[i-1]

                        if f_p_T0 < 0.1:
                            f_beta_r = f_p_T0
                        else:
                            f_beta_r = 0.1 * ( 1 - f_p_T0)
                    
                    # f_beta_r = 0
                    if f_p_T0 > f_beta_r:

                        mart_f[i] = ( mart_f[i-1] - b_km * (T[i]-T[i-1]) * ( f_p_T0 - f_beta_r + f_c_T0) )/ ( 1 - b_km * (T[i]-T[i-1]) )
                        beta_f[i] = beta_f[i-1] - (mart_f[i] - mart_f[i-1])
                        alpha_f[i] = alpha_f[i-1]

                else:
                        mart_trans_started_flag = False
                        f_p_T0 = 0
                        f_c_T0 = 1

                if (dT/dt >= 0):
                    
                    if mart_f[i-1] > f_mart_eq and T[i] > T_mart_dis: # JMAK model of mart to alpha + beta transformation:
                                            
                        # linearly interpolate JMAK parameters k_i and n_i according to current temperature:
                        k_mart_i = np.interp(T[i], temp_jmak_mart_to_AandB, k_mart_to_AandB)
                        n_mart_i = np.interp(T[i], temp_jmak_mart_to_AandB, n_mart_to_AandB)

                        #Calculate Tau:
                        tau_mart_i = ( -1/k_mart_i * m.log( (mart_f[i-1] - f_mart_eq) / (1-f_mart_eq) )) ** (1/n_mart_i)
                        # tau_mart_i = (-1/k_mart_i * m.log(1 - (mart_f[i-1] - f_mart_eq)/(1- alpha_f[i-1] - f_mart_eq)) )**(1/n_mart_i)
                
                        mart_f_inter = 1 - (1 - m.exp(-k_mart_i*(tau_mart_i+dt)**n_mart_i)) * (1 - f_mart_eq)
                        # mart_f_inter = f_mart_eq - ( m.exp( -k_mart_i * ( tau_mart_i + dt )**n_mart_i) ) * (mart_f[i-1] + beta_f[i-1] - f_mart_eq)

                        alpha_f_inter = alpha_f[i-1] + (mart_f[i-1] - mart_f_inter) * alpha_eq_i
                        beta_f_inter = beta_f[i-1] + (mart_f[i-1] - mart_f_inter) * beta_eq_i
                    
                    else:
                        mart_f_inter = mart_f[i-1]
                        alpha_f_inter = alpha_f[i-1]
                        beta_f_inter = beta_f[i-1]


                    if beta_f_inter < (beta_eq_i*(beta_f_inter+alpha_f_inter)) + 0.00001 and (alpha_f_inter + mart_f_inter) >= 0.001: #Parabolic growth of alpha to beta transformation
                    # if beta_f[i-1] < beta_eq_i and (alpha_f_inter + mart_f_inter) >= 0.001: #Parabolic growth of alpha to beta transformation
                        

                        f_diss = 2.2 * (10**-31) * ((T[i] + 273) **9.89)
                        t_star = (beta_f_inter/ (beta_eq_i*f_diss))**2
                        t_crit = f_diss ** -2

                        if 0 < dt + t_star and dt + t_star < t_crit:
                            alpha_final = 1 - beta_eq_i * f_diss * m.sqrt(dt + t_star)
                        elif dt + t_star >= t_crit:
                            alpha_final = 1 - beta_eq_i

                        alpha_inc = alpha_final - (alpha_f_inter + mart_f_inter)
                        beta_inc = -alpha_inc

                        alpha_f[i] = alpha_inc * alpha_f_inter / (alpha_f_inter + mart_f_inter) + alpha_f_inter
                        mart_f[i] = alpha_inc * mart_f_inter / (alpha_f_inter + mart_f_inter) + mart_f_inter
                        beta_f[i] = beta_f_inter + beta_inc

                    else:
                        alpha_f[i] = alpha_f_inter
                        beta_f[i] = beta_f_inter
                        mart_f[i] = mart_f_inter

                #Calculate lath thickness:
                # if alpha_f[i] + mart_f[i] > 0.001 and T[i] < 1000:
                if alpha_f[i] > 0.001 and T[i] < 1000:
                # if alpha_f[i] > 0.001:
                    # alpha_tot_i = alpha_f[i] + mart_f[i]
                    # alpha_tot_i_m1 = alpha_f[i-1] + mart_f[i-1]

                    alpha_tot_i = alpha_f[i]
                    alpha_tot_i_m1 = alpha_f[i-1]

                    teq = k * m.exp(-R/(T[i] + 273))

                    # t_lath[i] = 1/alpha_tot_i * ( t_lath[i-1] * alpha_tot_i_m1 + teq * (alpha_tot_i - alpha_tot_i_m1))
                    t_lath_int = (t_lath[i-1] - teq) * alpha_tot_i_m1 / alpha_tot_i + teq
                    # t_lath[i] = t_lath[i-1] * alpha_tot_i_m1/alpha_tot_i + teq * (alpha_tot_i - alpha_tot_i_m1)/alpha_tot_i

                    if T[i] < 1000:
                        t_lath[i] = ( t_lath_int**c + k0*m.exp(-Q/(Rg*(T[i]+273))) * dt)**(1/c)
                        # t_lath[i] = (t_lath_int**c)**(1/c)
                    else:
                        t_lath[i] = t_lath_int
                
                else :
                    t_lath[i] = t_lath_0

    # End of time-series iteration

    output_dict = {}
    output_dict["final_alpha_f"] = alpha_f[-1]
    output_dict["final_beta_f"] = beta_f[-1]
    output_dict["final_mart_f"] = mart_f[-1]
    output_dict["t_lath"] = t_lath[-1]



    os.chdir(user_path)


    # Plots:

    output_dir = Path(data_folder)

    #if output_dir does not exist create it:
    if not output_dir.exists():
        output_dir.mkdir()

    if create_output_json == True:

        output_dict_path = output_dir / "output.json"
        with open(output_dict_path, 'w') as json_file:
            json.dump(output_dict, json_file)

    if create_plots:



        fontsize = 20

        # Print temperature:
        plt.figure()
        plt.plot(t, T)
        plt.xlabel('Time (s)')
        plt.ylabel('Temperature (C)')
        plt.title('Temperature vs Time')
        plt.savefig('temp_vs_time.png')


        #plot RSL vs temperature:
        plt.figure()
        plt.plot(t, RLS)
        plt.xlabel('Time (s)')
        plt.ylabel('RLS')
        plt.title('RLS vs Time')
        plt.savefig('RLS_vs_time.png')

        #plot graph with two lines, with one of them having their values represented on the right axis and the other on the left axis:
        plt.figure()
        plt.rcParams['font.size'] = fontsize    

        fig, ax1 = plt.subplots()
        ax1.plot(t, T, 'b-')
        ax1.set_xlabel('time (s)')
        ax1.set_ylabel('Temperature (°C)', color='b')
        ax1.set_ylim(0, 1500)
        # ax1.set_xlim(0, 300)
        ax1.tick_params('y', colors='b')
        ax1.legend(["T"], loc='upper left')

        ax2 = ax1.twinx()
        ax2.plot(t, alpha_f, 'r--')
        ax2.plot(t, beta_f, 'k-.')
        ax2.plot(t, mart_f, 'g:')
        ax2.set_ylabel('Phase Fraction')
        ax2.set_ylim(0, 1)
        ax2.tick_params('y')


        ax2.legend(['$f_{\\alpha}$', '$f_{\\beta}$', "$f_{\\alpha '} $"], loc='upper right')

        fig.tight_layout()
        plt.savefig(output_dir / 'Phase_Fractions.png')


        #plot lath thickness
        plt.figure()
        #increase size of fonts:
        plt.rcParams['font.size'] = fontsize

        fig, ax1 = plt.subplots()
        ax1.plot(t, T, 'b')
        # ax1.set_xlim(0, 100)
        ax1.set_ylim(0, 1800)
        ax1.set_xlabel('time (s)')
        ax1.set_ylabel('Temperature (°C)', color='b')
        ax1.tick_params('y', colors='b')

        ax2 = ax1.twinx()
        ax2.plot(t, t_lath, 'r-')
        # ax2.plot(t, alpha_f + mart_f, 'r--')
        ax2.set_ylabel('Lath thickness ($\\mu m$)', color='r')
        ax2.tick_params('y', colors='r')
        ax2.set_ylim(0, 1.5)

        fig.tight_layout()
        plt.savefig(output_dir / 'Lath_thickness.png')

        


    
    if create_output_CSVs:
        np.savetxt(output_dir / "Alpha_Fraction.csv", np.column_stack((t, alpha_f)), delimiter=',', header='Time (s), Alpha Fraction')
        np.savetxt(output_dir / "Beta_Fraction.csv", np.column_stack((t, beta_f)), delimiter=',', header='Time (s), Beta Fraction')
        np.savetxt(output_dir / "Martensite_Fraction.csv", np.column_stack((t, beta_f)), delimiter=',', header='Time (s), Martensite Fraction')
        np.savetxt(output_dir / "Lath_thickness.csv", np.column_stack((t, t_lath)), delimiter=',', header='Time (s), Lath thickness (um)')

           

    #Compare with XRD data:
    
    # real_alpha = np.genfromtxt('./Predicted_total_alpha_XRD.csv', delimiter=',')
    # plt.figure()
    # plt.rcParams['font.size'] = fontsize
    # fig, ax1 = plt.subplots()
    # ax1.plot(t, T, 'b--')
    # ax1.set_xlabel('time (s)')
    # ax1.set_ylabel('Temperature (°C)', color='b')
    # ax1.tick_params('y', colors='b')

    # ax2 = ax1.twinx()
    # ax2.plot(t, alpha_f + mart_f, 'r-')
    # ax2.plot(real_alpha[:,0], real_alpha[:,1], 'r--')
    # ax2.set_ylabel('Phase Fraction', color='r')
    # ax2.tick_params('y', colors='r')
    # ax2.set_ylim(0, 1)
    # ax1.legend(['T'], loc='upper left')
    # ax2.legend(["$(\\alpha + \\alpha')$ predicted", "$(\\alpha + \\alpha')$ measured XRD"], loc='lower right')

    # fig.tight_layout()
    # plt.savefig('only_alpha.png')


    # # plot this function: alpha_eq_i += P[n]*(T[i]/1000)**(8-n)
    # plt.figure()
    # T_plot = np.linspace(600, 1100, 100)
    # alpha_eq = np.zeros(len(T_plot))

    # for i in range(0, len(T_plot)):
    #     for n in range(0, 9):
    #         alpha_eq[i] += P[n]*(T_plot[i]/1000)**(8-n)

    # plt.plot(T_plot, alpha_eq)
    # #save plot to file:
    # plt.savefig('alpha_eq.png')





    # Plot this function: teq = k * m.exp(-R/(T[i] + 273))

    # plt.figure()
    # T_plot = np.linspace(25, 1100, 50)
    # teq = np.zeros(len(T_plot))

    # for i in range(0, len(T_plot)):
    #     teq[i] = k * m.exp(-R/(T_plot[i] + 273))

    # plt.plot(T_plot, teq)

    # #save plot to file:
    # plt.savefig('teq.png')


    # save t, alpha_f, beta_f, mart_f, RLS, t_lath to separate file with time:
    # np.savetxt('alpha.csv', np.column_stack((t, T, alpha_f)), delimiter=',')


    return output_dict

    # # Print final phase fractions:
    # print('Final phase fractions:')
    # print('Alpha:', alpha_f[-1])
    # print('Beta:', beta_f[-1])
    # print('Martensite:', mart_f[-1])
    # print('Lath thickness (um):', t_lath[-1])





