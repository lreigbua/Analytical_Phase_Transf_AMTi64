from pathlib import Path
import os
import numpy as np
import math
import json



def Run(input_csv, data_folder, output_json = ""):

    print(f'Running Structure Model for {input_csv}...')

    # Save paths
    temperature_history_file = Path(input_csv)
    module_path = Path(__file__).parent
    output_json_path = Path(output_json)
    data_folder_path = Path(data_folder)
    user_folder = Path.cwd()

    data = np.genfromtxt(temperature_history_file, delimiter=',')

    # Change to the data folder
    if not os.path.exists(data_folder_path):
        os.mkdir(data_folder_path)
    os.chdir(data_folder_path)

    # Remove duplicate time points
    indices_to_remove = []  
    for i in range(1, len(data)):
        if abs(data[i,0] - data[i-1,0]) < 1e-6:
            indices_to_remove.append(i)

    # Remove specified rows
    data = np.delete(data, indices_to_remove, axis=0)

    time = data[:,0] - data[0,0] # start time from 0
    temperature = data[:,1]


    # Write time marks to inp file
    time_mark_text=f"*TIME POINTS, NAME=LASERON"
    for i in range(len(time)):
        time_mark_text += f'\n{time[i]},'

    output_file = Path('time_marks.inp')
    output_file.write_text(time_mark_text)

    # Write temperature history to amplitude inp file
    amplitude_text = "*AMPLITUDE, NAME=ASSIGN_LASER, TIME = TOTAL TIME"
    for i in range(len(time)):
        amplitude_text += f'\n{time[i]}, {temperature[i]},'

    output_file = Path('amp_temp_amprint.inp')
    output_file.write_text(amplitude_text)

    step_time = time[-1]

    
    #Modify and save main INP file
    text = Path(module_path / 'INP_default.inp').read_text()

    text = text.replace('step_time', str(step_time))
    text = text.replace('increment', '5.0')

    path = module_path / 'ABQ_phase_trans_types.inp'
    text = text.replace('ABQ_phase_trans_types.inp', str(path.resolve()))

    path = module_path / 'mat_phase_trans_ti64.inp'
    text = text.replace('mat_phase_trans_ti64.inp', str(path.resolve()))

    Path('structure_model.inp').write_text(text)

    #Run Abaqus
    os.system(f'module load abaqus/2020 && abaqus job=structure_model ask_delete=OFF interactive')

    #Read history output
    path = module_path / 'read_history_outputs.py'
    os.system(f'module load abaqus/2020 && abaqus cae noGUI={str(path.resolve())}')



    

    alpha = np.genfromtxt('alpha.csv', delimiter=',')
    beta = np.genfromtxt('Beta.csv', delimiter=',')
    temp = np.genfromtxt('TEMP.csv', delimiter=',')
    mart = np.genfromtxt('mart.csv', delimiter=',')

    # Calculate lath thickness
    time = temp[:,0]
    temp = temp[:,1]
    alpha = alpha[:,1]
    
    k = 0.5
    R=294

    Q = 97000;
    k0 = 20;
    t0 = 0.2;
    Rg = 8.314;

    t_lath = np.zeros(len(alpha))

    t_lath[0] = t0
    for i in range(1, len(alpha)):

        if alpha[i] > 0.001 and temp[i] < 1610:
            teq = k * math.exp(-R/temp[i]);
            Tavg = (temp[i] + temp[i-1]) / 2;
            t_lath[i] = ( (t_lath[i-1] - teq) * alpha[i-1] / alpha[i] + teq )
            # t_lath[i] = t_lath[i-1]

            # if Tavg > 500 and Tavg < 1000:
            #     t_lath[i]+= k0*math.exp(-Q/(Rg*Tavg))*(time[i]-time[i-1])
            # # t_lath[i] = ( (t_lath[i-1] - teq) * alpha[i-1] / alpha[i] + teq )
            # # if Tavg > 400 and time :
            # #     t_lath[i] = k0*math.exp(-Q/(Rg*Tavg))*(time[i]-time[i-1]) + t_lath[i-1]
        
        else:
            t_lath[i] = t0

    
    # add time and save to csv:
    t_lath = np.vstack((time, t_lath)).T
    np.savetxt("t_lath.csv", t_lath, delimiter=",")

    
    os.chdir(user_folder)

    if output_json != "":
        # Save final microstructure data as dict:
        data_dict = {}
        data_dict["alpha"] = alpha[-1]
        data_dict["beta"] = beta[-1,1]
        data_dict["mart"] = mart[-1,1]
        data_dict["t_lath"] = t_lath[-1,1]

        with open(output_json, 'w') as fp:
            json.dump(data_dict, fp)

        

    print(f'Finished Structure Model!')