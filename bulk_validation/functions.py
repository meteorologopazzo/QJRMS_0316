#### Authors: Alessandro Storer, Matteo Borgnino, Agostino Niyonkuru Meroni
#### date   : 2nd September 2025
#### this file contains two functions which act on lat-lon numpy.ndarrays


# slopes_r_p_mix()  
# computes the least-square regression line between y and x ;
# the degrees of freedom used to assess the significance of the estimated slope 
# are obtained through subsampling the data: 
# 1 point every nt points in time ; 1 point every nskip points in space
# nt and nskip were obtained from the computation of the autocorrelation function of the different fields
# for details check: Borgnino et al. 2025, https://doi.org/10.1029/2024GL112294

def slopes_r_p_mix(x, y, nt, nskip, ls=False):
    from scipy import stats
    import numpy as np
    
    xx = x[::nt,::nskip,::nskip]
    yy = y[::nt,::nskip,::nskip]
    
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    
    xx = xx[~np.isnan(xx)]
    yy = yy[~np.isnan(yy)]
    
    linreg = stats.linregress(x,y)
    corr_coeff, trash = stats.spearmanr(x,y)
    
    df = np.size(xx)-2
    mean_x = np.mean(x);  mean_x2 = np.mean(x**2); lever_arm = mean_x2-mean_x**2
    
    sigma_y = np.sqrt( np.sum(  (y-linreg.slope*x-linreg.intercept)**2 )/df  )
    sigma_slope = sigma_y/( np.sqrt(np.size(xx)*(lever_arm) ) )
    sigma_intercept = sigma_y*np.sqrt(mean_x2/( np.size(xx)*(lever_arm) ))
    
    sigmas = (sigma_slope, sigma_intercept)
    
    t_value_cannelli = linreg.slope/sigma_slope     # SOMETHING MISSING?
    p_value_cannelli = 2*(1 - stats.t.cdf(t_value_cannelli,df=df))
    
    # to ADD: scipy.stats.chisquare(f_obs, f_exp=None) --> not working as expected
    # chisq = np.sum(  (y-linreg.slope*x-linreg.intercept)**2 / std_y**2 )
    # in realtà mi servirebbero le dev std delle diverse osservazioni y --> calcolo X2 solo per i percentili e sto contento
    

    t_value = np.abs(corr_coeff)*np.sqrt((df)/(1-corr_coeff**2))
    p_value = 2*(1 - stats.t.cdf(t_value,df=df))
    
    if ls:
        return [linreg, corr_coeff, p_value, p_value_cannelli, sigmas]
    else:
        return linreg, corr_coeff, p_value, p_value_cannelli, sigmas

    
    


def div_sphere(field_a, field_b, llon, llat):
    """
    Function to calculate the divergence of a 2D vectorial field over a sphere, given the coordinates in degrees on 
    the same 2D grid. The derivatives are taken as second-order differences in the interior, and first-order 
    (forward or backward) on the edges.
    """
    import numpy as np
    R = 6371.0e3 # Earth radius in km.
    
    field_a = np.double(field_a)
    field_b = np.double(field_b)
    llon = np.double(llon)
    llat = np.double(llat)
    
    costheta = np.cos(llat*np.pi/180)

    div_a = field_a-field_a
    div_a[:,1:-1] = (field_a[:,2:]-field_a[:,:-2])/(R*costheta[:,1:-1]*(llon[:,2:]-llon[:,:-2])*np.pi/180)
    div_a[:,0] = (field_a[:,1]-field_a[:,0])/(R*costheta[:,0]*(llon[:,1]-llon[:,0])*np.pi/180)
    div_a[:,-1] = (field_a[:,-1]-field_a[:,-2])/(R*costheta[:,-1]*(llon[:,-1]-llon[:,-2])*np.pi/180)
    
    div_b = field_b-field_b
    div_b[1:-1,:] = (field_b[2:,:]*costheta[2:,:]-field_b[:-2,:]*costheta[:-2,:])/(R*costheta[1:-1,:]*(llat[2:,:]-llat[:-2,:])*np.pi/180)
    div_b[0,:] = (field_b[1,:]*costheta[1,:]-field_b[0,:]*costheta[0,:])/(R*costheta[0,:]*(llat[1,:]-llat[0,:])*np.pi/180)
    div_b[-1,:] = (field_b[-1,:]*costheta[-1,:]-field_b[-2,:]*costheta[-2,:])/(R*costheta[-1,:]*(llat[-1,:]-llat[-2,:])*np.pi/180)
        
    div = div_a + div_b
    return div




# nan_gaussian_filter()
# takes the instantaneous lat-lon field and applies an isotropic 2D Gaussian smoothing
# the width of the filter is set by its standard deviation sigma, which is expfessed in terms of number of points

def nan_gaussian_filter(field,sigma):
    """
    Function to smooth the field ignoring the NaNs.
    I follow the first answer here 
    https://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
    By default, the filter is truncated at 4 sigmas.
    If the sigma provided is zero, the function just returns the input field (by Ale 26.07.24)
    """
    import numpy as np
    from scipy.ndimage import gaussian_filter
    
    field = np.double(field)
    
    # Take the original field and replace the NaNs with zeros.
    field0 = field.copy()
    field0[np.isnan(field)] = 0
    
    if sigma == 'inf':
        return np.nanmean(field, axis=(0,1))
    
    elif sigma > 0:
        ff = gaussian_filter(field0, sigma=sigma)
    
        # Create the smoothed weight field.
        weight = 0*field.copy()+1
        weight[np.isnan(field)] = 0
        ww = gaussian_filter(weight, sigma=sigma)

        zz = ff/(ww*weight) # This rescale for the actual weights used in the filter and set to NaN where the field
                            # was originally NaN.
        #zz[zz == np.inf] = np.nan
        zz[np.isinf(zz)] = np.nan
        return zz
    
    elif sigma == 0:
        return field
    

    
    
# boxcar_filter()
# takes the instantaneous lat-lon field and applies a uniform 2D Gaussian smoothing
# the width of the filter (size) is expfessed in terms of number of points

def boxcar_filter(field,size):
    """
    Function to smooth the field ignoring the NaNs.
    I follow the first answer here 
    https://stackoverflow.com/questions/18697532/gaussian-filtering-a-image-with-nan-in-python
    By default, the filter is truncated at 4 sigmas.
    If the sigma provided is zero, the function just returns the input field (by Ale 26.07.24)
    """
    import numpy as np
    from scipy.ndimage import uniform_filter
    
    field = np.double(field)
    
    # Take the original field and replace the NaNs with zeros.
    field0 = field.copy()
    field0[np.isnan(field)] = 0
    
    if size == 'inf':
        return np.nanmean(field, axis=(0,1))
    
    elif size > 0:
        ff = uniform_filter(field0, size=size)
    
        # Create the smoothed weight field.
        weight = 0*field.copy()+1
        weight[np.isnan(field)] = 0
        ww = uniform_filter(weight, size=size)

        zz = ff/(ww*weight) # This rescale for the actual weights used in the filter and set to NaN where the field
                            # was originally NaN.
        #zz[zz == np.inf] = np.nan
        zz[np.isinf(zz)] = np.nan
        return zz
    
    elif size == 0:
        return field
    