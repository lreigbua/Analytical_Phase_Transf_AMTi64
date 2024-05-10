import sys
sys.path.insert(0, '../../../')
import Abaqus_PBF_Structure_model

# read input file form sys args:
input_file = sys.argv[1]


Abaqus_PBF_Structure_model.Run(input_csv = input_file,
                                data_folder = "./data",
                                 output_json = f"./plots_output.json")