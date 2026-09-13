# Analytical Microstructure Evolution Model of Ti64 during Additive Manufacturing

## Introduction

This code reads a temperature history csv file and calculates the microstructure evolution in terms of phase fractions and lath thicknesses for each time increment.
![image](https://github.com/user-attachments/assets/4ef3d9f3-72fa-4753-8d80-ff54927f5101)


### Flow-Chart of Microstructure Evolution Model calculations in each time increment:

![FLow-chart_MEM](https://github.com/lreigbua/Structure_Model_Abaqus/assets/93150422/3ef10158-58e6-496d-acfa-3bead2cd4fb7)

## How to run example

```bash
# Download repository and build
git clone https://github.com/lreigbua/Analytical_Phase_Transf_AMTi64.git
cd Analytical_Phase_Transf_AMTi64
pip install -e .

#Go to example folder
cd examples

#Run python script
python example.py
```

## Citation & Documentation

If you use this code in your research, please cite the foundational doctoral thesis. The full text also serves as the primary documentation for this repository, containing comprehensive details on the computational model's theoretical background, architecture, and validation.

**DOI:** [10.48730/pkpm-t966](https://doi.org/10.48730/pkpm-t966)

```bibtex
@phdthesis{reigbuades2026,
  author       = {Reig Buades, Luis Miguel},
  title        = {An integrated process-structure-property-performance modelling framework for additive layer manufacturing of Ti-6Al-4V},
  school       = {University of Strathclyde},
  year         = {2026},
  doi          = {10.48730/pkpm-t966},
  url          = {[https://doi.org/10.48730/pkpm-t966](https://doi.org/10.48730/pkpm-t966)}
}
```
