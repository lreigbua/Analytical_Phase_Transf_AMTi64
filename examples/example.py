import sys
sys.path.insert(0, '../')
import Analytical_Phase_Transf_AMTi64


Analytical_Phase_Transf_AMTi64.Run(input_csv = "./Inputs/amp_temp_amprint.txt",
                                data_folder = "./output",
                                 output_json = f"output.json")