from pathlib import Path
import numpy as np
import shutil
import os 
import sys

# read arg

temperature_file = sys.argv[1]

module_path = Path('../src/Abaqus_PBF_Structure_model')


temperature_history_file = '../input' / Path(temperature_file)

data = np.genfromtxt(temperature_history_file, delimiter=',')

# Remove duplicate time points
indices_to_remove = []  
for i in range(1, len(data)):
    if abs(data[i,0] - data[i-1,0]) < 1e-6:
        indices_to_remove.append(i)

# Remove specified rows
data = np.delete(data, indices_to_remove, axis=0)

time = data[:,0] - data[0,0] # start time from 0
temperature = data[:,1]


# Write time marks to file
time_mark_text=f"*TIME POINTS, NAME=LASERON"
for i in range(len(time)):
    time_mark_text += f'\n{time[i]},'

output_file = Path('time_marks.inp')
output_file.write_text(time_mark_text)

# Write temperature history to file
amplitude_text = "*AMPLITUDE, NAME=ASSIGN_LASER, TIME = TOTAL TIME"
for i in range(len(time)):
    amplitude_text += f'\n{time[i]}, {temperature[i]},'

output_file = Path('amp_temp_amprint.inp')
output_file.write_text(amplitude_text)

step_time = time[-1]



text = Path(module_path / 'INP_default.inp').read_text()

text = text.replace('step_time', str(step_time))
# text = text.replace('increment', str(time.max()+0.001))
text = text.replace('increment', '5.0')

Path('structure_model.inp').write_text(text)