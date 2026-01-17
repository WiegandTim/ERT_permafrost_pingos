# ---------------
# pyGIMLi_proc_3D
# ---------------

# Processing (Inversion) of 3D BERT/pyGIMLi ERT dat-files with pyGIMLi
# Important remark: This script is partly based on code examples from the pyGIMLi documentation (https://www.pygimli.org, access: 30 Oct 2025) and AI chatbots.

insite = 'HuskyLakes'

# Example: HuskyLakes
if insite == 'HuskyLakes':
    infolder = '/Users/YourName/Projects/data/HuskyLakes' # Change folder with BERT/pyGIMli dat-files here
    infiles = ['WSl5_3D.dat', 'Dip5_3D.dat', 'WSl5_Dip5_3D.dat']
    inparaDepth = 60
    inparaMaxCellFactor = 15 # factor will be multiplied with spacing for paraMaxCellSize; iteratively choosen
    inparaDX = 0

def pyGIMLI_3D_proc(
    infolder,
    infile,
    inparaMaxCellFactor,
    inparaDepth,
    inparaDX
):
    import sys
    import numpy as np
    import time
    import pygimli as pg
    import pygimli.meshtools as mt
    from pygimli.physics import ert

    start = time.time()

    # Read in data
    data = pg.load(f'{infolder}/{infile}')
    
    # Prepare topo
    topo = np.array(data.sensors())
    spacing = round(((topo[1, 0] - topo[0, 0])**2 + (topo[1, 1] - topo[0, 1])**2)**0.5)

    # Redirect stdout to a file
    class LoggerWriter:
        def __init__(self, filename):
            self.file = open(filename, 'w')
        def write(self, message):
            self.file.write(message)
        def flush(self):
            self.file.flush()

    sys.stdout = LoggerWriter(f'{infolder}/{infile[:-4]}_cell{spacing * inparaMaxCellFactor}_info.log')
    sys.stderr = LoggerWriter(f'{infolder}/{infile[:-4]}_cell{spacing * inparaMaxCellFactor}_info.log')

    # Prepare data
    data["k"] = ert.geometricFactors(data, dim = 3)
    print(data)
    mgr = ert.ERTManager(data)
    mgr.checkData(data)
    data.remove(data['rhoa'] < 0)
    print(data)
    plc = mt.createParaMeshPLC3D(data, paraDepth = inparaDepth, paraMaxCellSize = spacing * inparaMaxCellFactor, surfaceMeshQuality = 34, paraDX = inparaDX)
    mesh = mt.createMesh(plc, quality = 1.3)
    print(mesh)
    data["err"] = ert.estimateError(data, relativeError = 0.02)

    # Invert data
    mgr.inv.inv.setOptimizeLambda(True) # not working!
    mgr.inv.inv.setBlockyModel(True) # not working!
    mgr.invert(mesh = mesh, verbose = True)

    # Error statistics
    lam = mgr.inv.inv.getLambda()
    chi2 = mgr.inv.inv.getChi2()
    rmse = np.sqrt(np.mean((data['rhoa'] - mgr.inv.response) ** 2))
    nrmse = np.sqrt(np.mean(((data['rhoa'] - mgr.inv.response) / data['rhoa']) ** 2)) * 100
    mae = np.mean(((data['rhoa'] - mgr.inv.response) / data['rhoa'])) * 100

    err_file = f"{infolder}/{infile[:-4]}_cell{spacing * inparaMaxCellFactor}_errors.txt"
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

    # Output files
    outfolder = f"{infolder}/{infile[:-4]}_cell{spacing * inparaMaxCellFactor}_results"
    mgr.saveResult(folder = outfolder)

    # Processing time
    end = time.time()
    print(f'Processing time: {round(end - start)} sec')

# Processing of input dat-files
for i in infiles:
    pyGIMLI_3D_proc(infolder, i, inparaMaxCellFactor, inparaDepth, inparaDX)