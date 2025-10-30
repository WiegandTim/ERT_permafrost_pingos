# ---------------
# pyGIMLi_visu_2D
# ---------------

# Visualization of 2D pyGIMLi ERT profiles including options for model sensitivity plots and standardized zones of interest below mound landforms such as pingos.
# Important remark: This script is partly based on code examples from the pyGIMLi documentation (https://www.pygimli.org, https://github.com/gimli-org/pyGIMLi, access: 30 Oct 2025), T. Herring, 2024 (https://github.com/teddiherring/CPERS, access: 30 Oct 2025) and AI chatbots.

insite = 'HuskyLakes'

# Example: Husky Lakes
if insite == 'HuskyLakes':
    infolder = '/Users/YourName/Projects/data/HuskyLakes' # Change infolder here
    inprofiles = ['Q2', 'L13']
    ymin_vals = [-95, -40]
    ymax_vals = [35, 17.5]
    y_steps_vals = [20, 10] 
    fig_width = 16
    fig_height = 16
    apply_scatter = True
    apply_ZoI = True
    apply_violin = True
    spacing = 5
    zoi_scaling = 1
    alayer_factor = 0
    clabel_levels = None
    roll_along = False
    empty_layout = False
    ncol = 1
    cmin = 100
    cmax = 100000
    levels = 4
    cross_ele = None
    res_contours = None
    mirror_profile = False
    outside_DOI_alpha = 0.5

    inparams = []
    for i, profile in enumerate(inprofiles):
        inlist = [f'WSl5_{profile}', f'Dip5_{profile}', f'WSl5_Dip5_{profile}']

        if profile == 'Q2':
            pingo_xmin = [70]
            pingo_xmax = [270]
        elif profile == 'L13':
            pingo_xmin = [0]
            pingo_xmax = [70]
        else:
            pingo_xmin = []
            pingo_xmax = []

        ymin = ymin_vals[i]
        ymax = ymax_vals[i]
        y_steps = y_steps_vals[i]

        x = (infolder, inlist, ymin, ymax, y_steps, fig_width, fig_height, apply_scatter, apply_ZoI, apply_violin, pingo_xmin, pingo_xmax, spacing, zoi_scaling, alayer_factor, clabel_levels, roll_along, empty_layout, ncol, cmin, cmax, levels, cross_ele, res_contours)
        inparams.append(x) # inparams for multiprocess

def pyGIMLi_visu_2D(
    infolder,
    inlist,
    ymin,
    ymax,
    y_steps = 10,
    fig_width = 16,
    fig_height = 16,
    apply_scatter = True,
    apply_ZoI = True,
    apply_violin = True,
    pingo_xmin = [],
    pingo_xmax = [],
    spacing = 2,
    zoi_scaling = 1,
    alayer_factor = 0,
    clabel_levels = None,
    roll_along = False,
    empty_layout = False,
    ncol = 1,
    cmin = 100,
    cmax = 100000,
    levels = 4,
    cross_ele = None,
    res_contours = None,
    mirror_profile = False,
    outside_DOI_alpha = 0.5,
    roll_along_triangle_x1 = [],
    roll_along_triangle_y1 = [],
    roll_along_triangle_x2 = [],
    roll_along_triangle_y2 = [],
    roll_along_triangle_x3 = [],
    roll_along_triangle_y3 = []
):
    import pygimli as pg
    from pygimli.physics import ert
    import math
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.interpolate import griddata
    from scipy.spatial import ConvexHull
    from shapely.geometry import Polygon, box, Point
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.ticker import LogLocator

    cm = 1/2.54
    nrow = math.ceil(len(inlist) / ncol)

    cbar_params = {
        'cMin': cmin,
        'cMax': cmax,
        'logScale': True,
        'cMap': 'Spectral',
        'nLevs': levels,
        'label': 'Resistivity ($\Omega$m)',
        'orientation': 'horizontal'
    }

    for _ in range(2): # running the script two times is necessary to remove white space between subplots -> needs smarter solution!

        xmin_all = []
        xmax_all = []
        for i in range(len(inlist)):
            data_folder = infolder + '/' + inlist[i] + '_results'
            vtk = pg.load(data_folder + '/resistivity.vtk')
            positions = vtk.positions()
            pos_x = positions[:, 0]
            xmin = np.min(pos_x)
            xmax = np.max(pos_x)
            xmin_all.append(xmin)
            xmax_all.append(xmax)
        xmin = min(xmin_all)
        xmax = max(xmax_all)

        x_label = (xmin + ((xmax - xmin) / 2)) / (xmax - xmin)

        scat_DOI_list = []
        zoi_res_list = []
        zoi_cov_list = []
        doi_res_list = []
        doi_cov_list = []
        scat_xmin_all = []
        scat_xmax_all = []
        scat_ymin_all = []
        scat_ymax_all = []

        subplot_height = ((ymax - ymin) / (xmax - xmin)) * fig_height * cm
        total_height = nrow * subplot_height

        fig, ax = plt.subplots(nrow, ncol, figsize = [fig_width * cm, total_height])
        ax = np.atleast_1d(ax) # for nrow = 1, ncol = 1
        ax = ax.flatten()

        mpl.rcParams['axes.linewidth'] = 0.6

        for i in range(len(inlist)):
            data_folder = infolder + '/' + inlist[i] + '_results'
            vtk = pg.load(data_folder + '/resistivity.vtk')
            res = pg.load(data_folder + '/resistivity.vector')
            cover = vtk.data('Coverage')
            dat = pg.load(data_folder + '/data.dat')
            topo = np.array(dat.sensors())
            cell_centers = np.array(vtk.cellCenters()) # use dir(vtk) to show all attributes of vtk, e.g. cellCenters
            cell_nodes = np.array([node.pos() for node in vtk.nodes()])
            
            # Prepare coverage acc. to Herring/Lewkowicz 2022
            cover = np.array(cover)
            cover = 10 ** cover # resolve log10
            cover = cover / max(cover) # normalization acc. to maximum value

            cover_centers = []
            for j in range(vtk.cellCount()):
                cover_centers.append(cell_centers[j])
            cover_centers = np.unique(np.array(cover_centers), axis = 0)
            cover_centers = cover_centers[:, :2]
            x = cell_centers[:, 0]
            y = cell_centers[:, 1]
            z = cover

            interp_x, interp_y = np.meshgrid(
                np.linspace(x.min(), x.max(), 100),
                np.linspace(y.min(), y.max(), 100)
            )
            interp_z = griddata((x, y), z, (interp_x, interp_y), method = 'linear')
            interp_res = griddata((x, y), res, (interp_x, interp_y), method = 'linear')

            def round_up(x):
                a = 2.5
                b = a * (x // a + 1)
                if b - x < 1:
                    return b + a + 1
                return b + 1

            positions = vtk.positions()
            pos_y = positions[:, 1]

            points = np.vstack([
                [xmin, topo[np.argmin(topo[:, 0]), 1]],  # First point
                topo[:, :2],                            # Topo points
                [xmax, topo[np.argmax(topo[:, 0]), 1]], # Last point on topo
                [xmax, ymax],                           # Upper right corner
                [xmin, ymax]                            # Upper left corner
            ])
            polygon = MplPolygon(points, closed = True, color = 'white')
            
            # Calculate apparent res points and penetration depth
            a_xy = topo[dat['a']][:, :2]
            b_xy = topo[dat['b']][:, :2]
            m_xy = topo[dat['m']][:, :2]
            n_xy = topo[dat['n']][:, :2]

            # A (0), B (1), M (2), N (3), A_x (4), A_y (5), B_x (6), B_y (7), C_x (8), C_y(9), D_x (10), D_y (11)
            quadru_xy = np.hstack([dat['a'][:, None], dat['b'][:, None], dat['m'][:, None], dat['n'][:, None], a_xy, b_xy, m_xy, n_xy])

            x_coords = []
            y_coords = []

            for j in range(0, len(quadru_xy)):
            
                # Check max dist between electrodes of quadrupole (depending on array)
                a = np.array([quadru_xy[j, 0], quadru_xy[j, 1], quadru_xy[j, 2], quadru_xy[j, 3]])
                L1 = np.argmin(a)
                L2 = np.argmax(a)

                # Depth of Investigation indices after Barker 1989
                if L1 == 0 and L2 == 1: # Wenner
                    f = 0.17
                elif L1 == 0 and L2 == 3: # Dipole
                    f = 0.25
            
                x1 = quadru_xy[j, 3 + ((L1 + 1) * 2) - 1]
                y1 = quadru_xy[j, 3 + ((L1 + 1) * 2)]
                x2 = quadru_xy[j, 3 + ((L2 + 1) * 2) - 1]
                y2 = quadru_xy[j, 3 + ((L2 + 1) * 2)]

                # Calculate distance between A and B (Pythagoras)
                a = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

                # Calculate depth d
                d = f * a

                # Midpoint of A and B
                x_mid = (x1 + x2) / 2
                y_mid = (y1 + y2) / 2

                # Calculate vector AB and its perpendicular unit vector
                v_x = x2 - x1
                v_y = y2 - y1
                norm_v = np.sqrt(v_x**2 + v_y**2)
                
                # Perpendicular unit vector (-v_y/norm_v, v_x/norm_v)
                u_x = v_y / norm_v
                u_y = -v_x / norm_v

                # Scale perpendicular unit vector by depth d
                perp_x = u_x * d
                perp_y = u_y * d

                # Data point coordinates (perpendicular to surface)
                x_data_point = x_mid + perp_x
                y_data_point = y_mid + perp_y

                # Append results to lists
                x_coords.append(x_data_point)
                y_coords.append(y_data_point)

            chull_xy = np.column_stack((x_coords, y_coords))
            np.savetxt(infolder + '/' + inlist[i] + '_datapoints.txt', chull_xy, delimiter = ' ', fmt = '%.6f')

            def round_down(x):
                a = 2.5
                b = a * (x // a)
                if x - b < 1:
                    return b - a - 1
                return b - 1

            chull = ConvexHull(chull_xy)
            chull = chull_xy[chull.vertices]
            chull = np.array(chull)
            
            col_0 = chull[:, 0]
            min_val = np.min(col_0)
            max_val = np.max(col_0)
            min_idx = np.argmin(col_0)
            max_idx = np.argmax(col_0)
            
            if min_idx < max_idx:
                core = chull[min_idx:max_idx+1]
            else:
                core = np.vstack((chull[min_idx:], chull[:max_idx+1]))

            # DOI area
            elements_before_min = np.array([
                [xmin - 1, ymin - 1],
                [xmin - 1, ymax + 1],
                [min_val, ymax + 1]
            ])
            elements_after_max = np.array([
                [max_val, ymax + 1],
                [xmax + 1, ymax + 1],
                [xmax + 1, ymin - 1]
            ])
            DOI_area = np.vstack((
                elements_before_min,
                core,
                elements_after_max
            ))
            DOI_area = np.vstack((DOI_area, DOI_area[0:1]))  # Add copy of first point at end

            # Data area
            elements_before_min = np.array([
                [min_val, ymax + 1]
            ]) 
            elements_after_max = np.array([
                [max_val, ymax + 1],
            ])           
            data_area = np.vstack((
                elements_before_min,
                core,
                elements_after_max
            ))
            data_area = np.vstack((data_area, data_area[0:1]))  # Add copy of first point at end

            # Prepare cover and res for scatterplot
            scat = np.column_stack((cell_centers, cover, res))

            # Limit scat points to DOI_area
            DOI_poly = Polygon(DOI_area)
            data_poly = Polygon(data_area)
            
            if roll_along == True:
                points_to_find = [(roll_along_triangle_x1[i], roll_along_triangle_y1[i]), (roll_along_triangle_x2[i], roll_along_triangle_y2[i]), (roll_along_triangle_x3[i], roll_along_triangle_y3[i])]
                chull_xy_rounded = np.round(chull_xy, decimals=6)
                indices = []
                for point in points_to_find:
                    matches = np.where((chull_xy_rounded[:, 0] == point[0]) & (chull_xy_rounded[:, 1] == point[1]))
                    indices.append(matches[0][0])
                roll_along_triangle = chull_xy[indices]
                roll_along_poly = Polygon(roll_along_triangle)
                
                from shapely.ops import unary_union
                DOI_poly = unary_union([roll_along_poly, DOI_poly])
                data_poly = data_poly.difference(roll_along_poly)
                
                coords = np.array(DOI_poly.exterior.coords)
                DOI_area = np.column_stack((coords[:, 0], coords[:, 1]))
            
            scat_DOI = scat[np.array([Point(x, y).within(data_poly) for x, y in scat[:, :2]])]
            scat_DOI_list.append(scat_DOI)
            
            scat_xmin = min(scat_DOI[:, 3])
            scat_xmax = max(scat_DOI[:, 3])
            scat_ymin = min(scat_DOI[:, 4])
            scat_ymax = max(scat_DOI[:, 4])
            
            scat_xmin_all.append(scat_xmin)
            scat_xmax_all.append(scat_xmax)
            scat_ymin_all.append(scat_ymin)
            scat_ymax_all.append(scat_ymax)

            plot = pg.show( # ax[i], plot = pg.show( --> only pg.show() necessary, probably due to ax.flatten()
                vtk,
                res,
                ax = ax[i],
                colorBar = False,
                **cbar_params
            )

            ax[i].set_aspect('equal', adjustable = 'box')
            if not empty_layout == True:
                ax[i].scatter(topo[:, 0], topo[:, 1], color = 'black', label = 'Electrode Positions', s = 2.5, zorder = 3)
                ax[i].scatter(x_coords, y_coords, facecolors = 'black', marker = '.', linewidth = 0, s = 1, zorder = 2)
            if cross_ele is not None:
                ax[i].scatter(topo[np.array(cross_ele) - 1, 0], topo[np.array(cross_ele) - 1, 1], color = 'black', s = 2.5, marker = 'x', zorder = 3)
            polygon.set_zorder(2)
            ax[i].add_patch(polygon)
            if empty_layout == True:
                ax[i].fill(DOI_area[:, 0], DOI_area[:, 1], color = 'white', edgecolor = 'none', zorder = 2)
            else:
                ax[i].fill(DOI_area[:, 0], DOI_area[:, 1], color = 'white', edgecolor = 'none', alpha = outside_DOI_alpha, zorder = 2)
                cover_cont = ax[i].contour(interp_x, interp_y, interp_z, levels = [0.001, 0.003, 0.005, 0.01], colors = 'black', linewidths = 0.3, linestyles = 'solid', zorder = 1)
                if mirror_profile == False:
                    ax[i].clabel(cover_cont, inline = True, fontsize = 5, inline_spacing = 5, zorder = 1, levels = clabel_levels)
            ticks = np.arange(0, xmax, 20)
            ax[i].set_xticks(ticks)
            ax[i].tick_params(axis = 'both', which = 'major', labelsize = 8, length = 2, width = 0.7)
            ax[i].set_xlim(xmin, xmax)
            ax[i].set_ylim(ymin, ymax)
            ax[i].set_xlabel('m', fontsize = 8, labelpad = 2)
            ax[i].set_ylabel('m a.s.l.', fontsize = 8, labelpad = 2)
            ax[i].tick_params(axis = 'y', labelrotation = 90)
            for label in ax[i].get_yticklabels():
                label.set_verticalalignment('center')

            lower_tick = np.ceil(ymin / y_steps) * y_steps
            upper_tick = np.floor(ymax / y_steps) * y_steps
            major_ticks = np.arange(lower_tick, upper_tick, y_steps)
            major_ticks = np.append(major_ticks, upper_tick) # upper tick is not considered by np.arange
            upper_tick_diff = ymax - upper_tick
            if upper_tick_diff < (0.25 * y_steps):
                major_ticks = major_ticks[:-1]  # remove uppermost tick if not enough space
            lower_tick_diff = lower_tick - ymin
            if lower_tick_diff < (0.25 * y_steps):
                major_ticks = major_ticks[1:]
            major_labels = [str(int(tick)) for tick in major_ticks]
            ax[i].set_yticks(major_ticks)
            ax[i].set_yticklabels(major_labels)

            # Resistivity contour
            if res_contours is not None:
                res_cont = ax[i].contour(interp_x, interp_y, interp_res, levels = res_contours, zorder = 1, colors = 'darkred', linewidths=1)
                if mirror_profile == False:
                    ax[i].clabel(res_cont, inline = True, fontsize = 5, inline_spacing = 5, zorder = 1)

            if math.ceil((i + 1) / ncol) < nrow: # remove except of last plot
                ax[i].set_xticks([])
                ax[i].set_xticklabels([])
                ax[i].set_xlabel('')
            
            if ncol > 1 and i % 2 != 0: # remove at right plot
                ax[i].set_yticks([])
                ax[i].set_yticklabels([])
                ax[i].set_ylabel('')
            
            with open(f"{infolder}/{inlist[i]}_errors.txt", 'r') as errfile: # 'r' -> read only
                rmse_text = errfile.readlines()[3].strip()
                rmse = float(rmse_text.split(': ')[1])

            ax[i].text(0.025, 0.9, inlist[i], fontsize = 8, transform = ax[i].transAxes)
            if not empty_layout == True:
                ax[i].text(0.025, 0.825, f'RMSE: {rmse:.1f}%', fontsize = 8, transform = ax[i].transAxes)
                letter = chr(ord('a') + i)
                if len(inlist) > 1:
                    ax_label = f"{letter})"
                    ax[i].text(0.975, 0.9, ax_label, fontweight = 'bold', fontsize = 9, ha = 'right', transform = ax[i].transAxes)

            # Add zone of interest
            if apply_ZoI == True:
                from scipy.interpolate import interp1d

                def mirror_point(point, line_point1, line_point2):
                    P = np.array(point)
                    A = np.array(line_point1)
                    B = np.array(line_point2)

                    AB = B - A # direction vector

                    t = np.dot(P - A, AB) / np.dot(AB, AB)
                    proj = A + t * AB

                    P_mirrored = proj + zoi_scaling * (proj - P)

                    return P_mirrored

                def find_intersections(f1, f2, x_range, num_points = 1000):
                    x = np.linspace(x_range[0], x_range[1], num_points)
                    y1 = f1(x)
                    y2 = f2(x)
                    diff = y1 - y2
                    sign_changes = np.where(np.diff(np.sign(diff)))[0]
                    intersections = []
                    for i in sign_changes:
                        x_intersect = x[i] + (x[i+1] - x[i]) * diff[i] / (diff[i] - diff[i+1])
                        intersections.append(x_intersect)
                    return np.array(intersections)

                mask = np.zeros(len(cell_centers), dtype=bool)
                for j in range(len(pingo_xmin)):
                    pingo_a = pingo_xmin[j]
                    pingo_b = pingo_xmax[j]

                    surface = interp1d(tuple(topo[:, 0]), tuple(topo[:, 1]), kind = 'cubic', fill_value = 'extrapolate')
                    surface_min1 = interp1d(tuple(topo[:, 0]), tuple(topo[:, 1] - (spacing / 2)), kind = 'cubic', fill_value = 'extrapolate')

                    line_point1 = [pingo_a, surface(pingo_a)]
                    line_point2 = [pingo_b, surface(pingo_b)]

                    points_mirror = []
                    for k in topo:
                        point = [k[0], k[1] - (spacing / 2)] # points to be mirrored are half the spacing below surface
                        mirrored_point = mirror_point(point, line_point1, line_point2)
                        points_mirror.append(mirrored_point)
                    points_mirror = np.vstack(points_mirror)
                    points_mirror[:, 1] = points_mirror[:, 1] - alayer_factor
                    
                    surface_mirror = interp1d(points_mirror[:, 0], points_mirror[:, 1], kind = 'cubic', fill_value = 'extrapolate')
                    
                    # move upper surface after the ZOI has created by alayer_factor
                    surface_min1 = interp1d(tuple(topo[:, 0]), tuple(topo[:, 1] - (spacing / 2) - alayer_factor), kind = 'cubic', fill_value = 'extrapolate')

                    x = cell_centers[:, 0]
                    y = cell_centers[:, 1]

                    mask_cells = (x >= pingo_a) & (x <= pingo_b) & (y <= surface_min1(x)) & (y >= surface_mirror(x))
                    mask = mask | mask_cells

                    # Prepare line
                    x_range = (pingo_a, pingo_b)
                    x_intersect = find_intersections(surface_min1, surface_mirror, x_range)

                    x_between = np.linspace(x_intersect[0], x_intersect[1], 100)
                    ax[i].plot(x_between, surface_min1(x_between), color='black', linewidth = 0.7)
                    ax[i].plot(x_between, surface_mirror(x_between), color='black', linewidth = 0.7)
                
                points_select_x = x[mask]
                points_select_y = y[mask]
                
                zoi_res = res[mask]
                np.savetxt(infolder + '/' + inlist[i] + '_zoi_res.txt', zoi_res)
                zoi_cov = cover[mask]
                
                doi_cov = scat_DOI[:, 3]
                doi_res = scat_DOI[:, 4]
                np.savetxt(infolder + '/' + inlist[i] + '_doi_res.txt', doi_res)
                
                zoi_res_list.append(zoi_res)
                zoi_cov_list.append(zoi_cov)
                doi_res_list.append(doi_res)
                doi_cov_list.append(doi_cov)

            # mirror profile
            if mirror_profile == True:
                ax[i].invert_xaxis()
                if res_contours is not None:
                    ax[i].clabel(res_cont, inline = True, fontsize = 5, inline_spacing = 5, zorder = 1)
                if empty_layout == False:
                    ax[i].clabel(cover_cont, inline = True, fontsize = 5, inline_spacing = 5, zorder = 1, levels = clabel_levels)

        plt.subplots_adjust(wspace = 0.025, hspace = 0)

    scatplot_xmin = min(scat_xmin_all) / 2
    scatplot_xmax = 2
    scatplot_ymin = 10
    scatplot_ymax = (math.ceil(max(scat_ymax_all) / 100000) * 100000) * 2

    log_ymin = np.log10(scatplot_ymin)
    log_ymax = np.log10(scatplot_ymax)
    log_y_ticks = np.arange(np.ceil(log_ymin), np.floor(log_ymax) + 1)
    y_ticks = 10 ** log_y_ticks
    y_labels = []
    for j, tick in enumerate(y_ticks):
        if j % 2 == 1:  # Every second tick
            y_labels.append(f'$10^{int(np.log10(tick))}$')
        else:
            y_labels.append('')

    log_xmin = np.log10(scatplot_xmin)
    log_xmax = np.log10(scatplot_xmax)
    log_x_ticks = np.arange(np.ceil(log_xmin), np.floor(log_xmax) + 1)
    x_ticks = 10 ** log_x_ticks
    x_labels = []
    for j, tick in enumerate(x_ticks):
        if j % 2 == 1:
            x_labels.append(f'$10^{{{int(np.log10(tick))}}}$')
        else:
            x_labels.append('')

    if apply_scatter == True:
        for i in range(len(inlist)):

            # Polyline
            x = scat_DOI_list[i][:, 3]
            y = scat_DOI_list[i][:, 4]
            log_x = np.log10(x)
            log_y = np.log10(y)

            z = np.polyfit(log_x, log_y, deg = 1)   # slope and intercept
            poly = np.poly1d(z)

            log_y_trend = poly(log_x)
            y_trend = 10**log_y_trend

            # R2
            ss_tot = np.sum((log_y - np.mean(log_y))**2)
            ss_res = np.sum((log_y - log_y_trend)**2)
            r_squared = 1 - (ss_res / ss_tot)
            
            x_y_ratio = (xmax - xmin) / (ymax - ymin)  
            ax_inset = ax[i].inset_axes([0, 0, 0.45 / x_y_ratio, 0.45])
            ax_inset.set_xlim(scatplot_xmin, scatplot_xmax)
            ax_inset.set_ylim(scatplot_ymin, scatplot_ymax * 1.3) # * 1.3 to make upper end corresponding to the violin plot
            ax_inset.set_xscale('log')
            ax_inset.set_yscale('log')
            ax_inset.yaxis.tick_right()
            ax_inset.yaxis.set_label_position("right")
            ax_inset.scatter(x, y, color = 'dimgray', marker='.', linewidth = 0, s = 3, alpha = 0.4)
            ax_inset.plot(x, y_trend, color = 'black', linewidth = 0.5, linestyle = 'solid')
            ax_inset.text(0.05, 0.05, f'R²={r_squared:.2f}', transform = ax_inset.transAxes, fontsize = 6)
            ax_inset.xaxis.set_label_position('top')
            ax_inset.xaxis.tick_top()
            ax_inset.set_xlabel('Sensitivity', fontsize = 6, labelpad = 2, ha = 'center', va = 'bottom')
            ax_inset.set_ylabel('Res. [$\Omega$m]', fontsize = 6, labelpad = 2, ha = 'center', va = 'top')
            ax_inset.set_xticks(x_ticks)
            ax_inset.set_xticklabels(x_labels)
            ax_inset.set_yticks(y_ticks)
            ax_inset.set_yticklabels(y_labels)
            ax_inset.tick_params(axis = 'x', which = 'major', labelsize = 6, length = 2, pad = 1, width = 0.6)
            ax_inset.tick_params(axis = 'x', which = 'minor', length = 1, width = 0.4)
            ax_inset.tick_params(axis = 'y', which = 'major', labelsize = 6, length = 2, labelrotation = 90, pad = 2, width = 0.6)
            ax_inset.tick_params(axis = 'y', which = 'minor', length = 1, width = 0.4)
            for label in ax_inset.get_yticklabels():
                label.set_verticalalignment('center')
            ax_inset.patch.set_alpha(0.5)

            import matplotlib.ticker as ticker
            ax_inset.xaxis.set_major_locator(ticker.LogLocator(numticks = 9))
            ax_inset.xaxis.set_minor_locator(ticker.LogLocator(subs = (0.2, 0.4, 0.6, 0.8), numticks = 9))
            ax_inset.yaxis.set_major_locator(ticker.LogLocator(numticks = 9))
            ax_inset.yaxis.set_minor_locator(ticker.LogLocator(subs = (0.2, 0.4, 0.6, 0.8), numticks = 9))

    # Add violin plot
    if apply_violin == True:
        for i in range(len(inlist)):

            zoi_res = zoi_res_list[i]
            zoi_cov = zoi_cov_list[i]
            doi_res = doi_res_list[i]
            doi_cov = doi_cov_list[i]

            # Add violin plots
            import seaborn as sns
            import pandas as pd
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes
            
            df_a = pd.DataFrame({
                'value': np.concatenate([zoi_res, doi_res]),
                'group': ['A'] * (len(zoi_res) + len(doi_res)),
                'category': np.concatenate([
                    np.repeat('X', len(zoi_res)),
                    np.repeat('Y', len(doi_res))
                ])
            })
            df_b = pd.DataFrame({
                'value': np.concatenate([zoi_cov, doi_cov]),
                'group': ['B'] * (len(zoi_cov) + len(doi_cov)),
                'category': np.concatenate([
                    np.repeat('X', len(zoi_cov)),
                    np.repeat('Y', len(doi_cov))
                ])
            }) 
            
            ax_inset = ax[i].inset_axes([0.8, 0, 0.2, 0.45])
            ax_inset.patch.set_alpha(0.5)
            ax_inset.set_ylim(scatplot_ymin, scatplot_ymax * 1.3) # * 1.3 to get some space at the upper end
            ax_inset.set_yscale('log')
            ax_inset.set_ylabel('Res. [$\Omega$m]', fontsize = 6, labelpad = 2, ha = 'center', va = 'bottom')
            ax_inset.tick_params(axis = 'y', which = 'major', labelsize = 6, length = 2, labelrotation = 90, pad = 1, width = 0.6)
            ax_inset.tick_params(axis = 'y', which = 'minor', length = 1, width = 0.4)
            ax_inset.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            ax_inset.yaxis.set_minor_locator(LogLocator(base = 10, subs = (0.2, 0.4, 0.6, 0.8)))
            for label in ax_inset.get_yticklabels():
                label.set_verticalalignment('center')

            import matplotlib.ticker as ticker
            ax_inset.yaxis.set_major_locator(ticker.LogLocator(numticks = 9))
            ax_inset.yaxis.set_minor_locator(ticker.LogLocator(subs = (0.2, 0.4, 0.6, 0.8), numticks = 9))
            ax_inset.set_yticks(y_ticks)
            ax_inset.set_yticklabels(y_labels)
            ax_inset_b = ax_inset.twinx() # second y axis
            ax_inset_b.set_ylim(scatplot_xmin, scatplot_xmax)
            ax_inset_b.set_yscale('log')
            ax_inset_b.set_ylabel('Sensitivity', fontsize = 6, labelpad = 2, ha = 'center', va = 'top')
            ax_inset_b.tick_params(axis = 'y', which = 'major', labelsize = 6, length = 2, labelrotation = 90, pad = 2, width = 0.6)
            ax_inset_b.tick_params(axis = 'y', which = 'minor', length = 1, width = 0.4)
            ax_inset_b.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            ax_inset_b.yaxis.set_minor_locator(LogLocator(base = 10, subs = (0.2, 0.4, 0.6, 0.8)))
            for label in ax_inset_b.get_yticklabels():
                label.set_verticalalignment('center')
            ax_inset_b.yaxis.set_major_locator(ticker.LogLocator(numticks = 9))
            ax_inset_b.yaxis.set_minor_locator(ticker.LogLocator(subs = (0.2, 0.4, 0.6, 0.8), numticks = 9))
            ax_inset_b.set_yticks(x_ticks)
            ax_inset_b.set_yticklabels(x_labels)

            sns.violinplot(data = df_a, x = "group", y = "value", hue = "category", split = True, ax = ax_inset, linewidth = 0.8, cut = 0, density_norm = 'width')
            sns.violinplot(data = df_b, x = "group", y = "value", hue = "category", split = True, ax = ax_inset_b, linewidth = 0.8, cut = 0, density_norm = 'width')
            
            ax_inset.legend_.remove()
            ax_inset_b.legend_.remove()

            ax_inset.text(-0.25, ax_inset.get_ylim()[1] * 1.5, 'ZOI', ha = 'center', va = 'bottom', fontsize = 6) # scale is logarithmic --> * 0.1
            ax_inset.text(0.25, ax_inset.get_ylim()[1] * 1.5, 'TES', ha = 'center', va = 'bottom', fontsize = 6)
            ax_inset.text(0.75, ax_inset.get_ylim()[1] * 1.5, 'ZOI', ha = 'center', va = 'bottom', fontsize = 6)
            ax_inset.text(1.25, ax_inset.get_ylim()[1] * 1.5, 'TES', ha = 'center', va = 'bottom', fontsize = 6)

    # Common colorbar
    fig_height = fig.get_size_inches()[1]
    cbar_height = 0.25 * cm
    cbar_ax = fig.add_axes([0.15, -0.0125, 0.7, cbar_height / fig_height]) # new space for plotting colorbar (axes): x-pos: at 15% of plot width, y-pos: 1.25% of plot height below plot, width: 0.7 -> 70% width, height: 2.5% of plot height
    cbar_params2 = cbar_params
    cbar_params2.pop('label') # so far labelsize of legend label can't be adjusted --> remove label from cbar_params and add as title of legend
    cbar = pg.viewer.mpl.createColorBarOnly(ax = cbar_ax, **cbar_params2)
    cbar.tick_params(labelsize = 8, width = 0.6)
    cbar.minorticks_on()
    cbar.set_title('Resistivity [$\Omega$m]', fontsize = 8)
    for spine in cbar.spines.values():
        spine.set_linewidth(0.6)

    from matplotlib.ticker import FuncFormatter
    def number_format(x, pos):
        s = '{:,.0f}'.format(x) # format the number with commas then replace commas with spaces
        return s.replace(',', ' ')
    cbar.xaxis.set_major_formatter(FuncFormatter(number_format))

    fig.savefig(infolder + '/' + inlist[i].split('_')[-1] + '_figure.png', format = 'png', bbox_inches = 'tight', dpi = 600)
    plt.close(fig)

import os
import multiprocess
from multiprocess import Pool
max_proc = os.cpu_count() - 1

with Pool(max_proc) as pool:

    pool.starmap(pyGIMLi_visu_2D, inparams)
