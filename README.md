# READ ME

TPMS GENERATOR & ANALYSIS TOOL is a script that was developed with the objective of improving the efficiency of developing optimised TPMS scaffolds for orthopaedic implants. This tool requires and is best operated with nTop Automate version 5.16.2 and Ansys Mechanical 2024 R2 via PyMAPDL version 0.68.2. 

This particular version of PyMAPDL was used due to the compatibility with the digital licensing system used. Instability have been observed in more recent versions of this library and should be used with caution.

This tool allows design iterations to be completed seamlessly and most importantly, quickly with accuracy. Features include:
- list processing of stepped geometry
- secant method root finding of optimised geometry
- modular scripts to allow manual execution if necessary

Run times should vary between 15 minutes up to 2 hours.

A list of libraries used in this script have been provided in *requirements.txt*.

## Documentation
This tool works in the CLI with Python 3.10.2 or newer. Make sure to fork or download ALL the files and save them in an accessible location.

### Setup
Prerequisites:
1. nTop Automate version 5.16.2
2. Ansys Mechanical 2024 R2
3. PyMAPDL version 0.68.2

Recommendations: 
1. Update the execution path (*exec_file*) parameter in *digital_lab.py* to the Ansys Mechanical version on your device
2. Ensure parameters in *main.py* are adjusted to the desired values
3. Update the *tpms_no_porosity_dynamic.inp* file and ensure that files are written to the write directory. Ensure that the 'data' folder directory is specified correctly so that files are written correctly

### Guidance
Once the script has been downloaded, open the directory where you have saved the files in CLI.

For modularity, an argument parser has been used, requiring you to clarify which operation to execute. This requires you to input the *{mode_name}*, which is either *secant*, *converge*, or *custom_converge*. You can also add an optional parameter to change the *{tpms_type}* if necessary. TPMS types are specified by numbers followed by their name: *0-Gyroid, 1-Schwarz, 2-Diamond, 3-Lidinoid, 4-SplitP, 5-Neovius*.

Input line is as follows:

```
python main.py {mode_name} -l {tpms_type}
```

If necessary, direct execution of other modules can be done. Modification of the *main* section should provide access to the specific functions within each module.

## Disclosure
This project was supported by the University of Sheffield's Undergraduate Research Scheme (SURE) 2024 for its initial development between July 2024 and August 2024. Final development was conducted from September 2024 to May 2025 as part of a dissertation by the author. The author shall not be liable for any misuse or accidents from any further development of this script.

## Contact
For further enquiries, please do not hesitate to contact:
Zhe-Wei Kho - [zhewei.kho@gmail.com](mailto:zhewei.kho@gmail.com)
Vee San Cheong - [v.cheong@sheffield.ac.uk](mailto:v.cheong@sheffield.ac.uk)

Last Updated: 12th June 2025
