def seascyc_var_TLL(var1: str,
                    path: str,
                    begi: str = 'beg',
                    endi: str = 'end',
                    var2: str = None, 
                    math: str = None, 
                    level: int = None, 
                    plev: np.ndarray = default_levels,
                    mult: float = 1.) -> xr.DataArray: 

  '''Reads in up to two variables from a netCDF file (that includes a dimension for time) and creates
  Optional Parameters:
  --------------------------------------------------------------------------------------------------------
  begi, endi: class 'string' or 'int', Determines the first and last time element to be included in 
                                       analysis. Defaults 'beg' and 'end' read in the entire length of
                                       the file. Specifying an integer will read in time from begi as int 
                                       to endi as int. If integers, (endi-begi)/12 must be multiple of 12.  
  _______________________________________________________________________________________________________
  '''

  # Read in dataset 
  data = xr.open_dataset(path)

  # Assign values for begi and endi
  begi = 0 if begi == 'beg' else begi
  endi = len(data.time) if endi == 'end' else endi

  # Read in variable(s)
  var1_data = data[var1][begi:endi,:,:]*mult
  if var2 != None:
   var2_data = data[var2][begi:endi,:,:]*mult
 
