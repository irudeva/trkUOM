from netCDF4 import Dataset
import numpy as np
# from windspharm.standard import VectorWind
# from windspharm.tools import prep_data, recover_data, order_latdim
# from scipy import interpolate
import datetime as datetime
from datetime import timedelta  # Python standard library datetime  module
# import time as ftime
# import scipy.ndimage.filters as filters
# import scipy.ndimage.morphology as morphology
# import scipy.ndimage as sp
import calendar


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def dist(plat,plon,lat,lon):

    rad  =6371032 #in m
    pi   = np.pi
    dtr  = pi/180.
    rtd  = 180./pi

    dist=0.
    lat1=plat*dtr
    lon1=plon*dtr

    lat2=lat*dtr
    lon2=lon*dtr

    if lon1==lon2 :
        dist=abs(lat1-lat2)
        dist=rad*dist
    elif lat1==lat2 :
        dist=np.arccos(np.cos(lat1)**2*(np.cos(lon2-lon1)-1.)+1.)
        dist=rad*dist
    else :
        dist=np.arccos(np.sin(lat1)*np.sin(lat2)+np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1))
        dist=rad*dist

    return dist


R  =6371032 #in m
pi   = np.pi
deg = R*pi/180

#create lon/lat grid
grdres = 1  # grid reolution
lons = np.arange(0,360,grdres)
lats = np.arange(90,-91,-grdres)

cycrad = 2.5*deg   # cyclone area
#Parameter ird
#   ird = 0 => cyc area = cycrad
#   ird = 1 => cyc area = rd from trk file

ird  = 1

dset = "erain"
# scheme = "UOM"

#IMPORTANT: for trk files from 16 Nov YR-1 to 31 Jan YR+1
# netcdf reference time
time0 = 'hours since 1900-01-01 00:00:0.0'
dt0 = datetime.datetime(1900, 1, 1, 0)

max_ntrk =  20000   # max number of tracks per year
max_trklength = 200      # max length of a track

# creat arrays for recording tracks
trklen  = np.empty(max_ntrk,np.int16); trklen.fill(0)
trktime = np.empty([max_ntrk,max_trklength],np.int); trktime.fill(0)
trklon  = np.empty([max_ntrk,max_trklength]); trklon.fill(0)
trklat  =  np.copy(trklon); trklat.fill(0)
trkslp  =  np.copy(trklon); trklat.fill(0)
trkrd   =  np.copy(trklon); trklat.fill(0)
trklpl  =  np.copy(trklon); trklat.fill(0)
trkdp   =  np.copy(trklon); trklat.fill(0)
trkiop  = np.copy(trktime); trklat.fill(0)



for yr in range(1980,2017):
    #create time variable
    dt_nc = [ datetime.datetime(yr-1, 11, 16, 0)+i*timedelta(hours=6) for i in range(1704)]
    if calendar.isleap(yr):
        dt_nc = [ datetime.datetime(yr-1, 11, 16, 0)+i*timedelta(hours=6) for i in range(1708)]

    dt_hours = np.empty(len(dt_nc)); dt_hours.fill(0)
    for it, dt in enumerate(dt_nc) :
        tdelta =  dt - datetime.datetime(1900, 1, 1, 0)
        dt_hours[it] = divmod(tdelta.total_seconds(),3600)[0]


    #grid for cyclone records
    cycgrd_area = np.empty([len(lons),len(lats),len(dt_nc)],np.int16); cycgrd_area.fill(0)
    # cycgrd_rad  = np.copy(cycgrd_area)    #for radius from cyc tracking
    cycgrd_cen  = np.copy(cycgrd_area)


    #read trk
    ftrk = "../ERAint/trk/ntrkdat.%s.%d"%(dset,yr)
    print "ftrk =",ftrk
    print "Reading TRK:"



    print 'ftrk reading:'

    npnt  = np.zeros(max_ntrk,dtype = np.int16)
    clon  = np.zeros([max_ntrk,max_trklength])
    clat  = np.zeros_like(clon)
    cslp = np.zeros_like(clon)
    crad = np.zeros_like(clon)
    cdp  = np.zeros_like(clon)
    clpl = np.zeros_like(clon)
    date = np.zeros_like(clon,dtype = np.int)
    ntime  = np.zeros_like(date)
    ctime  = np.empty_like(date); ctime.fill(-1)
    cyr    = np.empty_like(date); cyr.fill(-1)
    cmon   = np.empty_like(date); cmon.fill(-1)
    cdate  = np.empty_like(date); cdate.fill(-1)
    chour  = np.empty_like(date); chour.fill(-1)
    ciop   = np.zeros_like(date)  #open/close marker

    f = open(ftrk, 'r')

    for line in range(1,max_ntrk):
        emptyline = f.readline()
        header    = f.readline()
        emptyline = f.readline()
        smthline = f.readline()

        # fout.write(emptyline)
        # fout.write(header)
        # fout.write(emptyline)
        # fout.write(smthline)

        if header == '':
         print ' ftrk is at the eof'
         break
        else :
         print ' header: ', header.strip()

        columns = header.split()
        tmp=columns[1]
        # print header
        if(tmp[-1] == "(") :
            ntrk = int(tmp[:-1])
            nit=int(columns[14])
        elif (tmp[-1] == ":") :
            ntrk = int(tmp[:-8])
            nit=int(columns[13])

        if nit > max_trklength :
            print " ERROR!!! Track length is more than max_trklength: nit = %d, max_trklength = %d" %(nit,max_trklength)
            quit()

        npnt[ntrk-1]=nit
        for n in range(0,nit):
            l = f.readline()
            columns = l.split()
            ciop[ntrk-1,n]=float(columns[5])
            clon[ntrk-1,n]=float(columns[7])
            clat[ntrk-1,n]=float(columns[8])
            cslp[ntrk-1,n]=float(columns[9])
            clpl[ntrk-1,n]=float(columns[10])
            cdp [ntrk-1,n]=float(columns[11])
            crad[ntrk-1,n]=float(columns[12])
            date[ntrk-1,n]=columns[1]
            ntime[ntrk-1,n]=columns[2]

            # for ind,year in enumerate(np.arange(yr-1,yr+1)):
            #     if str(year) == str(date[ntrk-1,n])[0:4]:
            #         iyr[ntrk-1,n] = ind
            # if iyr[ntrk-1,n] == -1:
            #     print "!!!! Check years in trk file"
            #     quit()
            cyr[ntrk-1,n]    = int(str(date[ntrk-1,n])[0:4])
            # to process the current year only
            # if cyr[ntrk-1,n] != yr:
            #     continue
            cmon[ntrk-1,n]   = int(str(date[ntrk-1,n])[4:6])
            cdate[ntrk-1,n]  = int(str(date[ntrk-1,n])[6:8])
            chour[ntrk-1,n]  = int(str(ntime[ntrk-1,n]/100))
            for ind,t in enumerate(dt_nc) :
                if t == datetime.datetime(cyr[ntrk-1,n],cmon[ntrk-1,n],cdate[ntrk-1,n],chour[ntrk-1,n]):
                    ctime[ntrk-1,n] = ind
                    break
            if ctime[ntrk-1,n] == -1:
                print "!!!! Check the time in trk file"
                quit()

            ####################################################
            # cyclones onto grid

            #find the closest grid value to the cyclone center
            # clon[ntrk-1,n] = 1 # test
            # clat[ntrk-1,n] = 55.  #test!

            if ird == 1:
                cycrad = crad[ntrk-1,n]*deg

            if clon[ntrk-1,n] < 0:
                 clon[ntrk-1,n]=clon[ntrk-1,n]+360
            ilon = find_nearest(lons,clon[ntrk-1,n])
            ilat = find_nearest(lats,clat[ntrk-1,n])

            cycgrd_area[ilon,ilat,ctime[ntrk-1,n]] = ntrk*10
            # cycgrd_rad[ilon,ilat,ctime[ntrk-1,n]] = ntrk*10
            cycgrd_cen[ilon,ilat,ctime[ntrk-1,n]]  = ntrk


            #check if cyclone area is close to poles
            if clat[ntrk-1,n] >= 0 :
                cdist_pole = (90 - clat[ntrk-1,n])*deg
            else :
                cdist_pole = (90 + clat[ntrk-1,n])*deg

            # print 'cyclone at ', ilat,ilon
            # print 'dist_pole', cdist_pole

            for nlat in np.arange(ilat,ilat+10) :
                if nlat <= len(lats)-1 :
                    # cdist = dist(lats[nlat],lons[ilon],clat[ntrk-1,n],clon[ntrk-1,n])
                    cdist1 = np.absolute(lats[nlat] - clat[ntrk-1,n])*deg
                        # print 'nlat=',nlat
                    if cdist1 <= cycrad :
                        for nlon in np.arange(ilon,ilon+len(lons)-1) :
                            if nlat == ilat and nlon == ilon:
                                continue
                            nnlon = nlon
                            if nlon >= len(lons):
                                nnlon = nlon - len(lons)
                                # print lons[nnlon]
                            # print 'nnlon=',nnlon
                            # print " // ++ ", nlat,nnlon
                            cdist = dist(lats[nlat],lons[nnlon],clat[ntrk-1,n],clon[ntrk-1,n])
                            # print cdist,cycrad
                            if cdist <= cycrad :
                                # print " ++ ", nlat, nnlon
                                cycgrd_area[nnlon,nlat,ctime[ntrk-1,n]] = cycgrd_area[nnlon,nlat,ctime[ntrk-1,n]] + 1
                            else :
                                if cdist_pole > cycrad:
                                    break
                        if cdist_pole > cycrad :
                            for nlon in np.arange(ilon-1,ilon-len(lons)+1,-1) :
                                nnlon = nlon
                                # cdist = dist(lats[nlat],lons[nnlon],clat[ntrk-1,n],clon[ntrk-1,n])
                                cdist = dist(lats[nlat],lons[nnlon],lats[ilat],lons[ilon])
                                if cdist <= cycrad :
                                    # print " +- ", nlat, nnlon
                                    cycgrd_area[nnlon,nlat,ctime[ntrk-1,n]] = cycgrd_area[nnlon,nlat,ctime[ntrk-1,n]] + 1
                                else :
                                    break
                    else :
                        # print "dist along longitude > cycdist => no need to check 1 "
                        break #dist along longitude > cycdist => no need to check
                else :
                    # print "out of bounds"
                    break #out of bounds

            for nlat in np.arange(ilat-1,ilat-10,-1) :
                if nlat >= 0:
                    # cdist = dist(lats[nlat],lons[ilon],clat[ntrk-1,n],clon[ntrk-1,n])
                    # cdist1 = dist(lats[nlat],lons[ilon],lats[ilat],lons[ilon])
                    # cdist1 = np.absolute(lats[nlat] - clat[ntrk-1,n])*deg
                    cdist1 = np.absolute(lats[nlat] - lats[ilat])*deg
                    if cdist1 <= cycrad :
                        for nlon in np.arange(ilon,ilon+len(lons)-1) :
                            nnlon = nlon
                            if nlon >= len(lons):
                                nnlon = nlon - len(lons)
                            # print " // -+ ", nlat,nnlon
                            cdist = dist(lats[nlat],lons[nnlon],clat[ntrk-1,n],clon[ntrk-1,n])
                            if cdist <= cycrad :
                                # print " -+ ", nlat, nnlon
                                cycgrd_area[nnlon,nlat,ctime[ntrk-1,n]] = cycgrd_area[nnlon,nlat,ctime[ntrk-1,n]] + 1
                            else :
                                if cdist_pole > cycrad:
                                    break
                        if cdist_pole > cycrad :
                            for nlon in np.arange(ilon-1,ilon-len(lons)+1,-1) :
                                nnlon = nlon
                                # print " // -- ", nlat,nnlon
                                cdist = dist(lats[nlat],lons[nnlon],clat[ntrk-1,n],clon[ntrk-1,n])
                                if cdist <= cycrad :
                                    # print " -- ", nlat, nnlon
                                    cycgrd_area[nnlon,nlat,ctime[ntrk-1,n]] = cycgrd_area[nnlon,nlat,ctime[ntrk-1,n]] + 1
                                else :
                                    break
                    else :
                        # print "dist along longitude > cycdist => no need to check 2"
                        break #dist along longitude > cycdist => no need to check
                else :
                    # print "out of bounds"
                    break #out of bounds

            # break #test - first point of the first cyclone onto the grid

        ####################################################
        # record tracks

        trklen[ntrk-1] = nit; print "trklen =", trklen[ntrk-1]," ntrk = ", ntrk
        trklon [ntrk-1,:nit]  = clon[ntrk-1,:nit]
        trklat [ntrk-1,:nit]  = clat[ntrk-1,:nit]
        trktime[ntrk-1,:nit]  = dt_hours[ctime[ntrk-1,:nit]]
        trkslp [ntrk-1,:nit]  = cslp[ntrk-1,:nit]
        trkrd  [ntrk-1,:nit]  = crad[ntrk-1,:nit]
        trkdp  [ntrk-1,:nit]  = cdp [ntrk-1,:nit]
        trklpl [ntrk-1,:nit]  = clpl[ntrk-1,:nit]
        trkiop [ntrk-1,:nit]  = ciop[ntrk-1,:nit]


        # break  #test - first cyclone onto the grid
    # print cyr[ntrk-1,n],cmon[ntrk-1,n],cdate[ntrk-1,n],chour[ntrk-1,n], n, ntrk
    # quit()
    #
    # break #test - first cyclone onto the grid


    f.close()
    print 'ftrk closed'

    #---Cyclone number to NetCDF write---------------------------------------------------------------
    print("Start NetCDF writing")

    # varlist = np.zeros(16, dtype = {'names': ['name', 'outname', 'data', 'scale'],
    #                                 'formats': ['a5', 'a5', '(241,480)f4', 'f4']} )
    #
    #
    # varlist[0] = ("u","u",u,1)
    # varlist[1] = ("v","v",v,1)
    # for iv in range(varlist['name'].size) :


    # ncvar = "cycnum"
    # print 'ncvar=',ncvar
    fcyc = '../ERAint/trkgrid/cycloc_rd.%d.nc' % (yr)
    nccyc = Dataset(fcyc, 'w', format='NETCDF4')
    nccyc.description = "Cyclone centers from  %s. Cyclone center number corresponds \
                         to the tracking number in the source file." % (ftrk)

    dimnam=('lon','lat','time')
    varnam=['longitude','latitude','time','cycloc','cycarea']

    #dimensions
    nccyc.createDimension(dimnam[0], lons.size)
    nccyc.createDimension(dimnam[1], lats.size)
    nccyc.createDimension(dimnam[2], None)


    #variables
    # for nv in range(0, 3) :
    #     # nccyc_var = nccyc.createVariable(varnam[nv], nc.variables[varnam[nv]].dtype,dimnam[nv])
    #     nccyc_var = nccyc.createVariable(varnam[nv], 'f',dimnam[nv])
        # for ncattr in nc.variables[varnam[nv]].ncattrs():
            # nccyc_var.setncattr(ncattr, nc.variables[varnam[nv]].getncattr(ncattr))
    #print(nc.variables['latitude'].ncattrs())

    #create  variables
     #lons
    nv=0
    nccyc_var = nccyc.createVariable(varnam[nv], 'f',dimnam[nv])
    if varnam[nv] == 'longitude' :
        nccyc_var.long_name = varnam[nv]
        nccyc_var.units = 'degrees_east'
        nccyc.variables[varnam[nv]][:] = lons

     #lats
    nv=1
    nccyc_var = nccyc.createVariable(varnam[nv], 'f',dimnam[nv])
    if varnam[nv] == 'latitude' :
        nccyc_var.long_name = varnam[nv]
        nccyc_var.units = 'degrees_north'
        nccyc.variables[varnam[nv]][:] = lats

     #time
    nv=2
    nccyc_var = nccyc.createVariable(varnam[nv], 'f',dimnam[nv])
    nccyc_var.long_name = varnam[nv]
    if varnam[nv] == 'time' :
        nccyc_var.calendar = 'gregorian'
        # nccyc_var.units = 'hours since 1900-01-01 00:00:0.0'
        nccyc_var.units = time0
        #for one time step - test!
        # tdelta =  datetime.datetime(cyr[ntrk-1,n],cmon[ntrk-1,n],cdate[ntrk-1,n],chour[ntrk-1,n]) - datetime.datetime(1900, 1, 1, 0)
        # nccyc.variables[varnam[nv]][:] = divmod(tdelta.total_seconds(),3600)[0]
        nccyc.variables[varnam[nv]][:] = dt_hours


    #cyclone centers to netcdf
    nv = 3
    nccyc_var = nccyc.createVariable(varnam[nv], 'f',dimnam[::-1])
    # print nccyc_var.shape
    nccyc_var.long_name = 'location of cyclone centers'
    # nccyc_var.scale_factor = varlist["scale"][iv]
    # nccyc_var.add_offset   = 0.
    # nccyc_var.units        = 'scale   %s' % varlist["scale"][iv]

    #print qx.shape
    #print nccyc_var.shape
    # nccyc_var[:,:,:] = np.swapaxes(gridcyc,0,2)
    nccyc_var[:,:,:] = np.swapaxes(cycgrd_cen[:,:,:],0,2)

    #cyclone area to netcdf
    nv =4
    nccyc_var = nccyc.createVariable(varnam[nv], 'f',dimnam[::-1])
    # print nccyc_var.shape
    nccyc_var.long_name = 'cyclone area'
    if ird == 0 :
        nccyc_var.radius = cycrad
    elif ird == 1:
        nccyc_var.radius = "rd from original trk file"
    # nccyc_var.scale_factor = varlist["scale"][iv]
    # nccyc_var.add_offset   = 0.
    # nccyc_var.units        = 'scale   %s' % varlist["scale"][iv]

    #print qx.shape
    #print nccyc_var.shape
    # nccyc_var[:,:,:] = np.swapaxes(gridcyc,0,2)
    nccyc_var[:,:,:] = np.swapaxes(cycgrd_area[:,:,:],0,2)



    nccyc.close()

    ##---Tracks to NetCDF write---------------------------------------------------------------


    ftrk = '../ERAint/trkgrid/trk.%d.nc' % (yr)
    nctrk = Dataset(ftrk, 'w', format='NETCDF4')
    nctrk.description = "Tracks from  %s" % (ftrk)

    #dimensions
    dimnam=('n','ntrk')

    nctrk.createDimension(dimnam[0], max_trklength )
    nctrk.createDimension(dimnam[1], None)

    #variables
    #  1D variable
    ncvar = nctrk.createVariable('trklen', 'i2',dimnam[1])
    ncvar.long_name = 'length of track'
    ncvar[:]    = trklen[:ntrk]
    #for one time step - test!
    # tdelta =  datetime.datetime(cyr[ntrk-1,n],cmon[ntrk-1,n],cdate[ntrk-1,n],chour[ntrk-1,n]) - datetime.datetime(1900, 1, 1, 0)
    # nccyc.variables[varnam[nv]][:] = divmod(tdelta.total_seconds(),3600)[0]

    #  2D variables
    frmt = "(%s,%s)f4" % (ntrk,max_trklength)
    varlist = np.zeros( 8, dtype = {  'names': ['name', 'long_name', 'dtype',  'units', 'data'],
                                    'formats': [  'a7',       'a31',    'a5',    'a44',   frmt]})


    # varlist[0] = ( "trklen",                     "length of track", 'i2',                       '',  trklen[:ntrk])
    varlist[1] = ("trktime",     "time of the current track point",  'i',                    time0, trktime[:ntrk,:])
    varlist[0] = ( "trklon","longitude of the current track point",  'f',           'degrees_east',  trklon[:ntrk,:])
    varlist[2] = ( "trklat", "latitude of the current track point",  'f',          'degrees_north',  trklat[:ntrk,:])
    varlist[3] = ( "trkslp", "slp in cyclone center from trk file",  'f',                    'hPa',  trkslp[:ntrk,:])
    varlist[4] = ( "trklpl",             "Laplacian from trk file",  'f','hPa/degrees_latitude**2',  trklpl[:ntrk,:])
    varlist[5] = ( "trkdp",                  "depth from trk file",  'f',                    'hPa',   trkdp[:ntrk,:])
    varlist[6] = ( "trkrad",                "radius from trk file",  'f',       'degrees_latitude',   trkrd[:ntrk,:])
    varlist[7] = ( "trkiop",                      "cyclone status", 'i2','0->strong/1->weak/10->str_open/11->weak_open',  trkiop[:ntrk,:])


    # varnam=['trklen','trktime','trklon','trklat','trkslp','trkrad','trklpl','trkdp','trkiop']

    for iv in range(varlist['name'].size) :
        # print varlist["name"][iv]
        ncvar = nctrk.createVariable(varlist["name"][iv], varlist["dtype"][iv],dimnam[::-1])
        ncvar.long_name = varlist["long_name"][iv]
        ncvar.units     = varlist[    "units"][iv]
        if varlist["name"][iv] == "trktime": ncvar.calendar = 'gregorian'
            # print varlist[     "data"][iv]
            # print trktime[:ntrk,:]
        ncvar[:]        = varlist[     "data"][iv]


    # nv=0
    # nccyc_var = nctrk.createVariable(varnam[nv], 'i2',dimnam[1])
    # nccyc_var.long_name = varnam[nv]
    # if varnam[nv] == 'trklen' :
    #     nccyc_var.long_name = 'length of track'
    #     # nccyc_var.units = 'hours since 1900-01-01 00:00:0.0'
    #     # nccyc_var.units = time0
    #     #for one time step - test!
    #     # tdelta =  datetime.datetime(cyr[ntrk-1,n],cmon[ntrk-1,n],cdate[ntrk-1,n],chour[ntrk-1,n]) - datetime.datetime(1900, 1, 1, 0)
    #     # nccyc.variables[varnam[nv]][:] = divmod(tdelta.total_seconds(),3600)[0]
    #     nccyc_var[:] = trklen[:ntrk]
    #
    # nv=1
    # nccyc_var = nctrk.createVariable(varnam[nv], 'i',dimnam[::-1])
    # nccyc_var.long_name = varnam[nv]
    # if varnam[nv] == 'trktime':
    #     nccyc_var.long_name = 'time of the current track point'
    #     nccyc_var.calendar = 'gregorian'
    #     # nccyc_var.units = 'hours since 1900-01-01 00:00:0.0'
    #     nccyc_var.units = time0
    #     #for one time step - test!
    #     # tdelta =  datetime.datetime(cyr[ntrk-1,n],cmon[ntrk-1,n],cdate[ntrk-1,n],chour[ntrk-1,n]) - datetime.datetime(1900, 1, 1, 0)
    #     # nccyc.variables[varnam[nv]][:] = divmod(tdelta.total_seconds(),3600)[0]
    #     nccyc_var[:,:] = trktime[:ntrk,:]
    #
    #
    #
    # nv=2
    # nccyc_var = nctrk.createVariable(varnam[nv], 'f',dimnam[::-1])
    # if varnam[nv] == 'trklon' :
    #     nccyc_var.long_name = 'longitude of the current track point'
    #     nccyc_var.units = 'degrees_east'
    #     # nctrk.variables[varnam[nv]][:] = lats
    #     nccyc_var[:,:] = trklon[:ntrk,:]
    #
    # nv=3
    # nccyc_var = nctrk.createVariable(varnam[nv], 'f',dimnam[::-1])
    # if varnam[nv] == 'trklat' :
    #     nccyc_var.long_name = 'latitude of the current track point'
    #     nccyc_var.units = 'degrees_north'
    #     # nctrk.variables[varnam[nv]][:] = lats
    #     nccyc_var[:,:] = trklat[:ntrk,:]
    #
    # nv=4
    # nccyc_var = nctrk.createVariable(varnam[nv], 'f',dimnam[::-1])
    # if varnam[nv] == 'trkslp' :
    #     nccyc_var.long_name = 'slp in cyclone center from trk file'
    #     nccyc_var.units = 'hPa'
    #     # nctrk.variables[varnam[nv]][:] = lats
    #     nccyc_var[:,:] = trkslp[:ntrk,:]
    #
    # nv=5
    # nccyc_var = nctrk.createVariable(varnam[nv], 'f',dimnam[::-1])
    # if varnam[nv] == 'trkrad' :
    #     nccyc_var.long_name = 'radius from trk file'
    #     nccyc_var.units = 'degrees_of_latitude'
    #     # nctrk.variables[varnam[nv]][:] = lats
    #     nccyc_var[:,:] = trkrd[:ntrk,:]
    #
    nctrk.close()


    ##---End NetCDF write---------------------------------------------------------------
