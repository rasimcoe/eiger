import numpy as np

def get_rectified_2d_spectrum(wfss_sci, wfss_var, x_trace, y_trace, w_grid, module='A'):

    if w_grid[0]>3.0 or w_grid[-1]<4.21:
        raise ValueError('Need w_grid[0]<=3.0 and w_grid[1]>=4.21.')    
    
    print(wfss_sci.shape)
    ### Mod A, xlim0 = blue, xlim1 = red
    ### Mod B, xlim0 = red, xlim1 = blue
    xlim0 = np.int32(1 + np.min(x_trace))
    xlim1 = np.int32(np.max(x_trace))

    ### Mod A, xmin = blue, xmax = red
    ### Mod B, xmin = red, xmin = blue
    xmin = np.max([0, xlim0])
    xmax = np.min([wfss_sci.shape[1]-1, xlim1])

    if xmax<0 or xmin>wfss_sci.shape[1]-1:
        print('xmin, xmax: ', xmin, xmax)
        raise ValueError('Spectrum is out of fiels.')
    
    ### X-pixel axis
    ### Mod A, blue --> red
    ### Mod B, red --> blue
    xaxis = xmin + np.arange(xmax - xmin + 1)

    ### "trace" at xaxis
    if module=='A':
        y_trace_xaxis = np.interp(xaxis, x_trace, y_trace)
        wav_xaxis = np.interp(xaxis, x_trace, w_grid)
    elif module=='B':
        y_trace_xaxis = np.interp(xaxis, np.flip(x_trace), np.flip(y_trace))
        wav_xaxis = np.interp(xaxis, np.flip(x_trace), np.flip(w_grid))
    else:
        raise ValueError('Module must be either A or B.')
    

    ###
    ylim0 = np.int32(np.min(y_trace_xaxis)) - 20
    ylim1 = np.int32(np.max(y_trace_xaxis)) + 20

    ymin = np.max([0, ylim0])
    ymax = np.min([wfss_sci.shape[0]-1, ylim1])

    if ymax<0 or ymin>wfss_sci.shape[0]-1:
        raise ValueError('Spectrum is out of fiels.')

    yaxis = ymin + np.arange(ymax - ymin + 1)

    ## Keep the dispersion direction in the module A and B
    spc_2d = np.zeros((31, xaxis.size))
    var_2d = np.zeros((31, xaxis.size))

    for ix in range(xaxis.size):

        flx_col = wfss_sci[ymin: ymax+1, xaxis[ix]]
        var_col = wfss_var[ymin: ymax+1, xaxis[ix]]
        bad_col = var_col*0.
        bad_col[np.where((var_col<=0.0)|(np.isnan(var_col))|(np.isnan(flx_col)))[0]]=1.0

        ### Cumulative sum
        flx_col_cum = np.nancumsum(flx_col)
        var_col_cum = np.nancumsum(var_col)
        bad_col_cum = np.nancumsum(bad_col)

        ### New axis
        yaxis_new = y_trace_xaxis[ix]-15+np.arange(31)

        flx_col_new0 = np.interp(yaxis_new - 0.5, yaxis+0.5, flx_col_cum, left=np.nan, right=np.nan)
        flx_col_new1 = np.interp(yaxis_new + 0.5, yaxis+0.5, flx_col_cum, left=np.nan, right=np.nan)
        var_col_new0 = np.interp(yaxis_new - 0.5, yaxis+0.5, var_col_cum, left=np.nan, right=np.nan)
        var_col_new1 = np.interp(yaxis_new + 0.5, yaxis+0.5, var_col_cum, left=np.nan, right=np.nan)
        bad_col_new0 = np.interp(yaxis_new - 0.5, yaxis+0.5, bad_col_cum, left=np.nan, right=np.nan)
        bad_col_new1 = np.interp(yaxis_new + 0.5, yaxis+0.5, bad_col_cum, left=np.nan, right=np.nan)


        bad_col_new = bad_col_new1-bad_col_new0
        flx_col_new = bad_col_new * 0.
        var_col_new = bad_col_new * 0.

        idx_bad  = np.where(bad_col_new>=0.5)[0]
        idx_good = np.where(bad_col_new<0.5)[0]

        bad_col_new[idx_bad]=1.0

        flx_col_new[idx_good] = (flx_col_new1-flx_col_new0)[idx_good]/(1-bad_col_new[idx_good])
        flx_col_new[idx_bad] = np.nan

        var_col_new[idx_good] = (var_col_new1-var_col_new0)[idx_good]/(1-bad_col_new[idx_good])
        var_col_new[idx_bad] = np.nan

        spc_2d[:,ix] = flx_col_new
        var_2d[:,ix] = var_col_new

    ### Linearly-spaced wavelength grid    
    dwav = 0.000975
    waxis_lin = 3.00 + dwav * np.arange(1240)

    ### Both mod A and mod B is in the same wavelength order
    spc_w2d = np.zeros((31, waxis_lin.size))
    var_w2d = np.zeros((31, waxis_lin.size))

    ### xaxis interpolated at the linearly-spaced wave grid
    ### Mod A, x_waxis is in ascending order (x_waxis_lin0 < x_waxis_lin1)
    ### Mod B, x_waxis is in descending order (x_waxis_lin0 > x_waxis_lin1)
    x_waxis_lin  = np.interp(waxis_lin,         w_grid, x_trace)
    x_waxis_lin0 = np.interp(waxis_lin-dwav/2., w_grid, x_trace)
    x_waxis_lin1 = np.interp(waxis_lin+dwav/2., w_grid, x_trace)

    for iy in range(31):
    #for iy in [15]:

        flx_row = spc_2d[iy, :]
        var_row = var_2d[iy, :]
        bad_row = var_row*0.
        bad_row[np.where((var_row<=0.0)|(np.isnan(var_row))|(np.isnan(flx_row)))[0]]=1.0

        if module=='B':
            flx_row = np.flip(flx_row)
            var_row = np.flip(var_row)
            bad_row = np.flip(bad_row)

        flx_row_cum = np.nancumsum(flx_row)
        var_row_cum = np.nancumsum(var_row)
        bad_row_cum = np.nancumsum(bad_row)

        ### Interpolation
        if module=='A':
            flx_row_new0 = np.interp(x_waxis_lin0, xaxis+0.5, flx_row_cum, left=np.nan, right=np.nan)
            flx_row_new1 = np.interp(x_waxis_lin1, xaxis+0.5, flx_row_cum, left=np.nan, right=np.nan)
            var_row_new0 = np.interp(x_waxis_lin0, xaxis+0.5, var_row_cum, left=np.nan, right=np.nan)
            var_row_new1 = np.interp(x_waxis_lin1, xaxis+0.5, var_row_cum, left=np.nan, right=np.nan)
            bad_row_new0 = np.interp(x_waxis_lin0, xaxis+0.5, bad_row_cum, left=np.nan, right=np.nan)
            bad_row_new1 = np.interp(x_waxis_lin1, xaxis+0.5, bad_row_cum, left=np.nan, right=np.nan)

        elif module=='B':
            flx_row_new0 = np.interp(-x_waxis_lin0, -np.flip(xaxis)+0.5, flx_row_cum, left=np.nan, right=np.nan)
            flx_row_new1 = np.interp(-x_waxis_lin1, -np.flip(xaxis)+0.5, flx_row_cum, left=np.nan, right=np.nan)
            var_row_new0 = np.interp(-x_waxis_lin0, -np.flip(xaxis)+0.5, var_row_cum, left=np.nan, right=np.nan)
            var_row_new1 = np.interp(-x_waxis_lin1, -np.flip(xaxis)+0.5, var_row_cum, left=np.nan, right=np.nan)
            bad_row_new0 = np.interp(-x_waxis_lin0, -np.flip(xaxis)+0.5, bad_row_cum, left=np.nan, right=np.nan)
            bad_row_new1 = np.interp(-x_waxis_lin1, -np.flip(xaxis)+0.5, bad_row_cum, left=np.nan, right=np.nan)

        bad_row_new = bad_row_new1-bad_row_new0
        flx_row_new = bad_row_new * 0.
        var_row_new = bad_row_new * 0.

        idx_bad  = np.where(bad_row_new>=0.5)[0]
        idx_good = np.where(bad_row_new<0.5)[0]

        bad_row_new[idx_bad]=1.0

        flx_row_new[idx_good] = (flx_row_new1-flx_row_new0)[idx_good]/(1-bad_row_new[idx_good])
        flx_row_new[idx_bad] = np.nan

        var_row_new[idx_good] = (var_row_new1-var_row_new0)[idx_good]/(1-bad_row_new[idx_good])
        var_row_new[idx_bad] = np.nan


        spc_w2d[iy,:] = flx_row_new
        var_w2d[iy,:] = var_row_new
    
    return waxis_lin, spc_w2d, var_w2d
    