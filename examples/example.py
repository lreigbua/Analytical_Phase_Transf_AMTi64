import sys
sys.path.insert(0, '../')
import Analytical_Phase_Transf_AMTi64


final_microstructure_dict = Analytical_Phase_Transf_AMTi64.Run(input_csv = "./Inputs/amp_temp_amprint.txt",
                                 data_folder = "./output",
                                 Output_Json = "output.json",
                                 create_output_CSVs = True,
                                 create_plots = True
                                 )

print(final_microstructure_dict)