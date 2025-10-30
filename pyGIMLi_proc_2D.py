# ---------------
# pyGIMLi_proc_2D
# ---------------

# Processing (Inversion) of 2D BERT/pyGIMLi ERT dat-files with pyGIMLi
# Important remark: This script is partly based on code examples from the pyGIMLi documentation (https://www.pygimli.org, access: 30 Oct 2025), T. Herring, 2024 (https://github.com/teddiherring/CPERS, access: 30 Oct 2025) and AI chatbots.

infolder = '/Users/YourName/Projects/data/HuskyLakes' # Change folder with BERT/pyGIMli dat-files here

import os
os.chdir(infolder)

def pyGIMLi_proc_2D(infile):
    import sys
    import pygimli as pg
    from pygimli.physics import ert
    import numpy as np  
    import os

    infolder = os.path.dirname(os.path.abspath(infile))

    if infile.endswith('.dat'):

        # Redirect stdout to a file
        class LoggerWriter:
            def __init__(self, filename):
                self.file = open(filename, 'w')
            def write(self, message):
                self.file.write(message)
            def flush(self):
                self.file.flush()

        sys.stdout = LoggerWriter(f'{infile[:-4]}_info.log')
        sys.stderr = LoggerWriter(f'{infile[:-4]}_info.log')

        # Read in data
        data = pg.load(os.path.abspath(infile))

        # Prepare topo -> along ground distance to horizontal distance # alternative: pygimli.meshtools.interpolateAlongCurve(), methods: linear, spline, harmonic
        topo = np.array(data.sensors())
        delta_x = np.diff(topo[:, 0], prepend = topo[0, 0])
        delta_y = np.diff(topo[:, 2], prepend = topo[0, 2])
        delta_x_hori = np.sqrt(delta_x**2 - delta_y**2)
        x_hori = np.cumsum(delta_x_hori)
        topo[:, 0] = x_hori
        data.setSensorPositions(topo)
        spacing = round(delta_x_hori[1] - delta_x_hori[0])

        # Prepare data
        data['k'] = ert.createGeometricFactors(data, numerical = True)
        mgr = ert.ERTManager(sr = False) # sr = False is different from Herring, 2024
        k0 = ert.createGeometricFactors(data)
#        ert.showData(data, vals=k0/data['k'], label='Topography effect') #Showing data not possible in multiprocessing
        mgr.checkData(data)
        print(data)
        data['err'] = ert.estimateError(data, relativeError = 0.03, absoluteUError = 5e-5)
#        ert.show(data, data['err'] * 100)

        # Inversion
        mgr.inv.inv.setOptimizeLambda(True) # From Herring, 2024, however, not working!
        mgr.inv.inv.setBlockyModel(True) # From Herring, 2024, however, probably not working!
        invdata = mgr.invert(data, verbose = True, paraDX = 0.1 * spacing, paraMaxCellSize = spacing) #paraDepth = 60 ## paraDZ = 0.2 with little influence
#        mgr.showResult()

        # Error statistics
        lam = mgr.inv.inv.getLambda()
        chi2 = mgr.inv.inv.getChi2()
        rmse = np.sqrt(np.mean((data['rhoa'] - mgr.inv.response) ** 2))
        nrmse = np.sqrt(np.mean(((data['rhoa'] - mgr.inv.response) / data['rhoa']) ** 2)) * 100
        mae = np.mean(((data['rhoa'] - mgr.inv.response) / data['rhoa'])) * 100
        
        err_file = f"{infolder}/{infile[:-4]}_errors.txt"
        with open(err_file, 'w') as f:
            f.write(f"lam: {round(lam, 2)}\n")
            f.write(f"chi2: {round(chi2, 2)}\n")
            f.write(f"RMSE: {round(rmse, 2)}\n")
            f.write(f"NRMSE: {round(nrmse, 2)}\n")
            f.write(f"MAE: {round(mae, 2)}\n")
            f.write(f"mean_rho: {round(np.mean(data['rhoa']), 2)}\n")
            f.write(f"sd_rho: {round(np.std(data['rhoa']), 2)}\n")
            f.write(f"mean_calc: {round(np.mean(mgr.inv.response), 2)}\n")
            f.write(f"sd_calc: {round(np.std(mgr.inv.response), 2)}\n")
            f.write(f"spacing: {spacing}\n")

        import matplotlib.pyplot as plt

        # Output plot
        def output_plot(invdata, infolder, infile, cmin, cmax):
            fig, ax = plt.subplots(1, 1, figsize = [12, 10])
            ax, cBar = mgr.showResult( # shows last inverted result
               invdata,
                ax = ax,
                cMap = 'turbo_r',
                coverage = mgr.coverage(),
                cMin = cmin,
                cMax = cmax,
                logScale = True,
                orientation = 'vertical',
                label = 'Resistivity ($\Omega$m)'
            )
            def number_form(x):
                if x >= 1000:
                    return f"{int(x / 1000)}k"
                else:
                    return str(int(x))
            cmin_exp = number_form(cmin)
            cmax_exp = number_form(cmax)
            fig.savefig(infolder + '/' + infile[:-4] + '_' + cmin_exp + '_' + cmax_exp + '.png', format = 'png', bbox_inches = 'tight', dpi = 500)
#            print(f'Processing of {infile} with cmin = {cmin} | cmax = {cmax} finished')

        cmin=np.round(np.sort(mgr.paraModel())[int(len(mgr.paraModel())/100)])
        cmax=np.round(np.sort(mgr.paraModel())[-int(len(mgr.paraModel())/100)]) ## Example: array has 2418 values (len) --> 24.18 (/100) --> -24 (-int) --> 24.-last value of the sorted array = 142222.2258 --> 142222 (round)
        output_plot(invdata, infolder, infile, cmin, cmax)

        # Output files
        outfolder = f"{infolder}/{infile[:-4]}_results"
        mgr.saveResult(folder = outfolder)
        for files in os.listdir(f"{outfolder}/ERTManager"): # move files from ERTManager to main folder as renaming of folder ERTManager due to multiprocessing not possible
            orig = os.path.join(f"{outfolder}/ERTManager", files)
            dest = os.path.join(outfolder, files)
            os.rename(orig, dest)
        os.rmdir(f"{outfolder}/ERTManager")

# Multiprocessing of input dat-files
import multiprocess # Important: when working in Jupyter, multiprocess instead of multiprocessing necessary
from multiprocess import Pool
max_proc = os.cpu_count() - 1
with Pool(max_proc) as pool:
    pool.map(pyGIMLi_proc_2D, os.listdir(infolder))