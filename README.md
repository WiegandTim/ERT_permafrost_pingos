# ERT_permafrost_pingos
Python codes for multidimensional processing and visualization of ERT data including options for model sensitivity plots and the creation of standardized zones of interest below mound landforms such as pingos in permafrost environments.

## Installation
Set up pyGIMLi according to the pyGIMLi documentation (https://www.pygimli.org/installation.html, access: 30 Oct 2025). All Python scripts were tested and run in Jupter Lab.

## Running pyGIMLi_proc_2D
- Change `infolder` to the input and workspace directory that includes the dat-files to be processed. All dat-files will be targeted, therefore, make sure that the folder only contains dat-files that have a proper BERT/pyGIMLi structure and are already filtered based on resistivity or stack deviation. The script does not include data preprocessing and filtering.
- Based on the infile names, results will be stored in _results folders, error statistics including chi2 and NRMSE (%) in _errors.txt files, and information on the inversion process in _info.log files. An examplary visualization plot is given as png. However, it is highly recommended to run `pyGIMLi_visu_2D` for more advanced visualization including sensitivity.

## Running pyGIMLi_visu_2D
Specify the parameters below. You can do so in the list `inparams` to enable multiprocessing. For some parameters, different values can be provided as the idea of the script is to plot different versions of a profile together like, e.g., WSl, DipDip, combined WSl-DipDip.

### Parameters
- `infolder`: Folder with input files and workspace; necessary parameter
- `inlist`: List of input files; necessary parameter
- `ymin`: Minimum y-axis value of output plot; let more space when plotting additional scatter plot and violin plots; necessary parameter
- `ymax`: Maximum y-axis value of output plot; necessary parameter
- `y_steps`: Distance between y axis ticks/labels in m; default: 10
- `fig_width`: Width of figure in cm; there is no full control on `fig_width` and `fig_height` due to equal aspect ratio of x and y axis; default: 16
- `fig_height`: Height of figure in cm; there is no full control on `fig_width` and `fig_height` due to equal aspect ratio of x and y axis; default: 16
- `apply_scatter`: Whether to apply a scatter plot of resistivity and sensitivity; default: True
- `apply_ZoI`: Whether to apply a Zone of Interest (ZOI) below the mound feature(s); default: True
- `apply_violin`: Whether to apply violin plots of resistivity and sensitivity of ZOI and TES (full trapezoidal/triangular ERT section that is covered by data points); default: True
- `pingo_xmin`: Distance along the profile where the pingo/permafrost mound starts ("foot" of the pingo); relevant for ZOI creation; give different start values when more than one pingo/mound in the transect should be considered; default: []
- `pingo_xmax`: Distance along the profile where the pingo/permafrost mound ends ("foot" of the pingo); relevant for ZOI creation; give different end values when more than one pingo/mound in the transect should be considered; default: []
- `spacing`: Electrode spacing in m; relevant for ZOI creation, as ZOI is located half the spacing below the surface of the pingo/mound; default: 2
- `zoi_scaling`: Scaling factor for the lower, mirrored boundary of ZOI; default: 1 (no scaling)
- `alayer_factor`: Additional factor for moving the whole ZOI; e.g., 5 moves the ZOI 5 * spacing downwards (e.g., to exclude the active layer); default: 0
- `clabel_levels`: Factor to individually specify sensitivity contour label levels; default: None
- `roll_along`: Whether to exclude a triangle in between of two ERT measurements that form a roll-along measurement from TES; if applying, specify x and y coordinates below; default: False
- `empty_layout`: Whether to generate the plot without sensitivity contour lines, electrode positions and only showing TES; default: False
- `ncol`: Number of columns in layout; default: 1
- `cmin`: Minimum resistivity value for colourbar (logarithmic); default: 100
- `cmax`: Maximum resistivity value for colourbar (logarithmic); default: 100000
- `levels`: Number of logarithmic value steps on colourbar; default: 4 (here with `cmin` and `cmax`: 100, 1000, 10000, 100000)
- `cross_ele`: Possibility to indicate a crossing ERT profile in the plot; default: None
- `res_contours`: Possibility to draw specific resistivity contour lines (e.g., assumed frozen-unfrozen boundary); default: None
- `mirror_profile`: Whether to mirror the profile horizontally; default: False
- `outside_DOI_alpha`: Transparency of white bleaching of the area around TES; default: 0.5
- `roll_along_triangle_x1`: X-coordinates for the lower left corner of the roll-along triangle *); default: []
- `roll_along_triangle_y1`: Y-coordinates for the lower left corner of the roll-along triangle *); default: []
- `roll_along_triangle_x2`: X-coordinates for the upper central corner of the roll-along triangle *); default: []
- `roll_along_triangle_y2`: Y-coordinates for the upper central corner of the roll-along triangle *); default: []
- `roll_along_triangle_x3`: X-coordinates for the lower right corner of the roll-along triangle *); default: []
- `roll_along_triangle_y3`: Y-coordinates for the lower right corner of the roll-along triangle *); default: []

*) Provide different values when plotting different versions of a profile together (e.g., WSl, DipDip, combined WSl-DipDip)

## Running pyGIMLi_proc_3D

## Citation
The codes partly built on examples from the pyGIMLi documentation (https://www.pygimli.org, access: 30 Oct 2025), T. Herring, 2024 (https://github.com/teddiherring/CPERS, access: 30 Oct 2025) and AI chatbots.
Therefore, please, additionally cite following resources:
- Rücker, C., Günther, T., Wagner, F. M., 2017. pyGIMLi: An open-source library for modelling and inversion in geophysics. Computers & Geosciences 109, pp. 106-123.
- Herring, T., Lewkowicz, A. G., Chiasson, A., Wang, Y., Way, R. G., Young, J. M., Froese, D., Smith, S. L., Andersen, B., Bellehumeur-Génier, O., Bevington, A. R., Bonnaventure, P. P., Duguay, M. A., Etzelmüller, B., Gooseff, M. N., Godsey, S. E., Miceli, C. M., 2024. The Canadian Permafrost Electrical Resistivity Survey (CPERS) database: 15 years of permafrost resistivity data. Arctic Science 10(4), pp. 850-856.

Sensitivity (Coverage) acc. to:
- Herring, T., Lewkowicz, A. G., 2022. A systematic evaluation of electrical resistivity tomography for permafrost interface detection using forward modeling. Permafrost and Periglacial Processes 33(2), pp. 134-146.
