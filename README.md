# ERT_permafrost_pingos
A python code for multidimensional processing and visualization of ERT data including options for model sensitivity plots and the creation of standardized zones of interest of mound landforms such as pingos in permafrost environments.

## Installation
Set up pyGIMLi according to the pyGIMLi documentation (https://www.pygimli.org/installation.html, access: 30 Okt 2025). All Python scripts were tested and run in Jupter Lab.

## Running pyGIMLi_proc_2D
- Change `infolder` to the input and workspace directory that includes the dat-files to be processed. All dat-files will be targeted, therefore, make sure that the folder only contains dat-files that have a proper BERT / pyGIMLi structure and are already filtered based on resistivity or stack deviation. The script does not include data preprocessing and filtering.
- Based on the infile names, results will be stored in _results folders, error statistics including chi2 and NRMSE (%) in _errors.txt files, and information on the inversion process in _info.log files. An examplary visualization plot is given as png. However, it is highly recommended to run `pyGIMLi_visu_2D` for more advanced visualization including sensitivity.

## Running pyGIMLi_visu_2D

## Running pyGIMLi_proc_3D

## Citation
The codes partly built on examples from the pyGIMLi documentation (https://www.pygimli.org, access: 30 Oct 2025), T. Herring, 2024 (https://github.com/teddiherring/CPERS, access: 30 Oct 2025) and AI chatbots.
Therefore, please, additionally cite following resources:
- Rücker, C., Günther, T., Wagner, F. M., 2017. pyGIMLi: An open-source library for modelling and inversion in geophysics. Computers & Geosciences 109, pp. 106-123.
- Herring, T., Lewkowicz, A. G., Chiasson, A., Wang, Y., Way, R. G., Young, J. M., Froese, D., Smith, S. L., Andersen, B., Bellehumeur-Génier, O., Bevington, A. R., Bonnaventure, P. P., Duguay, M. A., Etzelmüller, B., Gooseff, M. N., Godsey, S. E., Miceli, C. M., 2024. The Canadian Permafrost Electrical Resistivity Survey (CPERS) database: 15 years of permafrost resistivity data. Arctic Science 10(4), pp. 850-856.
