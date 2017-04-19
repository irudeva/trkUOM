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



# def polerot(plat,plon,ilat,ilon):
# # the location of the new pole is (plat,plon)
# # the new coordinates are in a Cartesian c. system where:
# # real North pole is (y=90,x=0) in the new system if plat > 0 and the East is (y=0,x=90)
# # real South pole is (90,0) in the new system if plat < 0 and the East is (y=0,x=-90)
#
#
#     import numpy as np
#     rad  =6371032 #in m
#     pi   = np.pi
#     dtr  = pi/180.
#     rtd  = 180./pi
#
#     plat0 = np.copy(plat)
#     plon0 = np.copy(plon)
#
#     plat = (plat-90.)*dtr
#     plon = plon*dtr
#
#     ilon[:] = [x*dtr for x in ilon]
#     ilat[:] = [x*dtr for x in ilat]
#
#     nlat = np.zeros_like(ilat)
#     nlon = np.zeros_like(ilon)
#
#     for index, (lon,lat) in enumerate(zip(ilon, ilat)):
#
#         # convert to cartesian coordinates
#         # from bronstein p.217
#
#         cos_plon = np.cos(plon)
#         if plon == 90*dtr or plon == -90*dtr:
#             cos_plon = 0
#         cos_lon = np.cos(lon)
#         if lon == 90*dtr or lon == -90*dtr:
#             cos_lon = 0
#
#         x=rad*np.cos(lat)*cos_lon
#         y=rad*np.cos(lat)*np.sin(lon)
#         z=rad*np.sin(lat)
#
#
#         # turn XY
#         xn=x*cos_plon+y*np.sin(plon)
#         yn=x*(-np.sin(plon))+y*cos_plon
#         zn=z
#
#         # tilt XZ
#         xnn=xn*np.cos(plat)+zn*np.sin(plat)
#         ynn=yn
#         znn=xn*(-np.sin(plat))+zn*np.cos(plat)
#
#         # convert back to polar coordinates
#         rr=np.sqrt(xnn*xnn+ynn*ynn+znn*znn)
#         nlat[index]=np.arcsin(znn/rr)*rtd
#
#
#         nlon[index]=np.arctan2(ynn,xnn)*rtd
#
#         # print "res----", nlat[index],nlon[index]
#
#     return nlat,nlon
#
#
# def rotated_grid_transform(plat, plon, ilat, ilon, option):
#
#     import numpy as np
#
#     rad  =6371032 #in m
#     pi   = np.pi
#     dtr  = pi/180.
#     rtd  = 180./pi
#
#     ilon[:] = [x*dtr for x in ilon]
#     ilat[:] = [x*dtr for x in ilat]
#
#     nlat = np.zeros_like(ilat)
#     nlon = np.zeros_like(ilon)
#
#     theta = plat-90 # Rotation around y-axis
#     phi = plon #  Rotation around z-axis
#
#     phi = phi*dtr
#     theta = theta*dtr
#     if option == 2: # Rotated -> Regular
#      phi = -phi
#      theta = -theta
#
#     cos_phi = np.cos(phi)
#     if phi == 90*dtr or phi == -90*dtr:
#         cos_phi = 0
#     sin_phi = np.sin(phi)
#     if phi == 0*dtr or phi == 180*dtr or phi == -180*dtr:
#         sin_phi = 0
#     cos_theta = np.cos(theta)
#     if theta == 90*dtr or theta == -90*dtr:
#         cos_theta = 0
#     sin_theta = np.sin(theta)
#     if theta == 0*dtr :
#         sin_theta = 0
#
#     # print cos_theta,cos_phi, sin_phi,sin_theta
#
#     for index, (lon,lat) in enumerate(zip(ilon, ilat)):
#
#         cos_lon = np.cos(lon)
#         if lon == 90*dtr or lon == -90*dtr:
#             cos_lon = 0
#         sin_lon = np.sin(lon)
#         if lon == 0*dtr or lon == 180*dtr or lon == -180*dtr:
#             sin_lon = 0
#         cos_lat = np.cos(lat)
#         if lat == 90*dtr or lat == -90*dtr:
#             cos_lat = 0
#         sin_lat = np.sin(lat)
#         if lat == 0*dtr :
#             sin_lat = 0
#
#         x = cos_lat*cos_lon
#         y = cos_lat*sin_lon
#         z = sin_lat
#
#         # print cos_lat,cos_lon, sin_lon
#
#
#         if option == 1: # Regular -> Rotated
#
#             x_new = cos_theta*cos_phi*x + cos_theta*sin_phi*y + sin_theta*z
#             y_new = -sin_phi*x + cos_phi*y
#             z_new = -sin_theta*cos_phi*x - sin_theta*sin_phi*y + cos_theta*z
#
#
#         elif option == 2: # Rotated -> Regular
#
#             # phi = -phi
#             # theta = -theta
#
#             x_new = cos_theta*cos_phi*x + sin_phi*y + sin_theta*cos_phi*z
#             y_new = -cos_theta*sin_phi*x + cos_phi*y - sin_theta*sin_phi*z
#             z_new = -sin_theta*x + cos_theta*z
#
#
#         lon_new = np.arctan2(y_new,x_new) # % Convert cartesian back to spherical coordinates
#         lat_new = np.arcsin(z_new)
#
#         nlon[index] = lon_new*rtd  # % Convert radians back to degrees
#         nlat[index] = lat_new*rtd
#
#     # print 'end rotation'
#     return nlat,nlon
#
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

# def detect_local_minima(arr):
#     # http://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
#     """
#     Takes an array and detects the troughs using the local maximum filter.
#     Returns a boolean mask of the troughs (i.e. 1 when
#     the pixel's value is the neighborhood maximum, 0 otherwise)
#     """
#     # define an connected neighborhood
#     # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
#     neighborhood = morphology.generate_binary_structure(len(arr.shape),2)
#     # apply the local minimum filter; all locations of minimum value
#     # in their neighborhood are set to 1
#     # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#minimum_filter
#     local_min = (filters.minimum_filter(arr, footprint=neighborhood)==arr)
#     # local_min is a mask that contains the peaks we are
#     # looking for, but also the background.
#     # In order to isolate the peaks we must remove the background from the mask.
#     #
#     # we create the mask of the background
#     background = (arr==0)
#     #
#     # a little technicality: we must erode the background in order to
#     # successfully subtract it from local_min, otherwise a line will
#     # appear along the background border (artifact of the local minimum filter)
#     # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
#     eroded_background = morphology.binary_erosion(
#         background, structure=neighborhood, border_value=1)
#     #
#     # we obtain the final mask, containing only peaks,
#     # by removing the background from the local_min mask
#     detected_minima = local_min - eroded_background
#     return np.where(detected_minima)
#
# def local_minima(array2d):
#     return ((array2d <= np.roll(array2d,  1, 0)) &
#             (array2d <= np.roll(array2d, -1, 0)) &
#             (array2d <= np.roll(array2d,  1, 1)) &
            # (array2d <= np.roll(array2d, -1, 1)))

# test!!!!

# Rotate Pole
# nplat, nplon = rotated_grid_transform(11,165,[88],[50],1)
# print nplat, nplon
# Reverse rotation
# nplat, nplon = rotated_grid_transform(11,165,[nplat],[nplon],2)
# print nplat, nplon
# quit()

# end test!!!

R  =6371032 #in m
pi   = np.pi
deg = R*pi/180

#create lon/lat grid
lons = np.arange(0,360,1)
lats = np.arange(90,-91,-1)

cycrad = 2.*deg

dset = "erain"
# scheme = "UOM"


for yr in range(2005,2006):
    dt_nc = [ datetime.datetime(yr, 1, 1, 0)+i*timedelta(hours=6) for i in range(1460)]
    if calendar.isleap(yr):
        dt_nc = [ datetime.datetime(yr, 1, 1, 0)+i*timedelta(hours=6) for i in range(1464)]

    cycgrd = np.empty([len(lons),len(lats),len(dt_nc)]); cycgrd.fill(0)

    ftrk = "../ERAint/trk/ntrkdat.%s.%d"%(dset,yr)
    print "ftrk =",ftrk
    # fout = open(frad, 'wb')

    # for iyr,cyr in enumerate(np.arange(yr,yr+1)):
    # for iy,cyr in enumerate(np.arange(yr-1,yr+2)):
#     for iy,cyr in enumerate([yr-1,yr+1,yr]):
#     #  print "loop", iy, cyr
#      fnc  = "/Users/Irina/work/DATA/%st/erain.mslp.%d.nc"%(dset,cyr)
#      print "fnc  =",fnc
#
#     # read netcdfF
#      print 'fnc reading...'
#      dimnam=('longitude','latitude','time')
#      varnam=['longitude','latitude','time','msl']
#
#      nc = Dataset(fnc, 'r')
#      v=0
#      for var in varnam:
#         if nc.variables[varnam[v]].name != var:
#             print "Variables don't agree", var, nc.variables[varnam[v]].name, v
#             exit()
#         v += 1
#
#      lons = nc.variables[varnam[0]][:]
#      lats = nc.variables[varnam[1]][:]
#      time = nc.variables[varnam[2]][:]
#      mslp = nc.variables[varnam[3]][:]/100.
#      #mslp0=mslp/100.
#
#      #smoothing
#      mslp = sp.filters.gaussian_filter(mslp, sigma = 2, order = 0)
#
#      # fix for leap years !!!
#      if iy == 0:
#          slp=np.tile(mslp,(3,1,1,1))
#      elif iy == 1 :
#          gridcyc = np.zeros_like(mslp)
#      slp[iy,:time.size,:,:] = mslp
#
# # create an array with 3 time series
#      dt_nc[0:time.size] = [datetime.datetime(1900, 1, 1, 0) + datetime.timedelta(hours=int(t))\
#            for t in time]
#      if iy == 0:
#         dtall=np.tile(dt_nc,(3,1))
#
#      dtall[iy,:] = dt_nc
#     #  print len(dtall)
#     #  print dtall[:,:10]

    print "Reading TRK"


    #read trk

    print 'ftrk reading:'
    max_ntrk =  505 #20000
    max_trklength = 200
    npnt  = np.zeros(max_ntrk,dtype = np.int)
    clon  = np.zeros((max_ntrk,max_trklength))
    clat  = np.zeros_like(clon)
    cslp = np.zeros_like(clon)
    date = np.zeros_like(clon,dtype = np.int)
    trktime  = np.zeros_like(date)
    ctime  = np.empty_like(date); ctime.fill(-1)
    cyr    = np.empty_like(date); cyr.fill(-1)
    cmon   = np.empty_like(date); cmon.fill(-1)
    cdate  = np.empty_like(date); cdate.fill(-1)
    chour  = np.empty_like(date); chour.fill(-1)

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

        #print 'ntrk=',ntrk, 'nit=',nit

        npnt[ntrk-1]=nit
        for n in range(0,nit):
            l = f.readline()
            # print "l: ",l
            # fout.write(l)
            columns = l.split()
            clon[ntrk-1,n]=float(columns[7])
            clat[ntrk-1,n]=float(columns[8])
            cslp[ntrk-1,n]=float(columns[9])
            date[ntrk-1,n]=columns[1]
            trktime[ntrk-1,n]=columns[2]

            # for ind,year in enumerate(np.arange(yr-1,yr+1)):
            #     if str(year) == str(date[ntrk-1,n])[0:4]:
            #         iyr[ntrk-1,n] = ind
            # if iyr[ntrk-1,n] == -1:
            #     print "!!!! Check years in trk file"
            #     quit()
            cyr[ntrk-1,n] = int(str(date[ntrk-1,n])[0:4])
            if cyr[ntrk-1,n] != yr:
                continue
            cmon[ntrk-1,n]  = int(str(date[ntrk-1,n])[4:6])
            cdate[ntrk-1,n]  = int(str(date[ntrk-1,n])[6:8])
            chour[ntrk-1,n]   = int(str(trktime[ntrk-1,n]/100))
            print ctime[ntrk-1,n]
            for ind,t in enumerate(dt_nc) :
                if t == datetime.datetime(cyr[ntrk-1,n],cmon[ntrk-1,n],cdate[ntrk-1,n],chour[ntrk-1,n]):
                    ctime[ntrk-1,n] =ind
                    break
            if ctime[ntrk-1,n] == -1:
                print "!!!! Check the time in trk file"
                quit()

            # record netcdf file

            #find the closest grid value to the cyclone center
            clon[ntrk-1,n] = 1 # test
            clat[ntrk-1,n] = 89.  #test!

            if clon[ntrk-1,n] < 0:
                 clon[ntrk-1,n]=clon[ntrk-1,n]+360
            ilon = find_nearest(lons,clon[ntrk-1,n])
            ilat = find_nearest(lats,clat[ntrk-1,n])

            cycgrd[ilon,ilat,ctime[ntrk-1,n]] = 22

            #circle around cyc center
            # TO DO!!! add lons>360,lat>90...

            #check if cyclone area is close to poles
            # cdist1 = dist(90,0,clat[ntrk-1,n],clon[ntrk-1,n])
            # cdist2 = dist(-90,0,clat[ntrk-1,n],clon[ntrk-1,n])
            # print cdist1
            # cdist_pole = np.amin([cdist1,cdist2])
            if clat[ntrk-1,n] >= 0 :
                cdist_pole = (90 - clat[ntrk-1,n])*deg
            else :
                cdist_pole = (90 + clat[ntrk-1,n])*deg

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
                            print " // ++ ", nlat,nnlon
                            cdist = dist(lats[nlat],lons[nnlon],clat[ntrk-1,n],clon[ntrk-1,n])
                            print cdist,cycrad
                            if cdist <= cycrad :
                                print " ++ ", nlat, nnlon
                                cycgrd[nnlon,nlat,ctime[ntrk-1,n]] = cycgrd[nnlon,nlat,ctime[ntrk-1,n]] + 1
                            else :
                                if cdist_pole > cycrad:
                                    print 'break'
                                    break
                        if cdist_pole > cycrad :
                            for nlon in np.arange(ilon-1,ilon-len(lons)+1,-1) :
                                nnlon = nlon
                                # if nlon <0 :
                                #     nnlon = nlon + len(lons)
                                print " // +- ", nlat, nnlon
                                # cdist = dist(lats[nlat],lons[nnlon],clat[ntrk-1,n],clon[ntrk-1,n])
                                cdist = dist(lats[nlat],lons[nnlon],lats[ilat],lons[ilon])
                                print cdist,cycrad
                                if cdist <= cycrad :
                                    print " +- ", nlat, nnlon
                                    cycgrd[nnlon,nlat,ctime[ntrk-1,n]] = cycgrd[nnlon,nlat,ctime[ntrk-1,n]] + 1
                                else :
                                    print 'break'
                                    break
                    else :
                        print "dist along longitude > cycdist => no need to check 1 "
                        break #dist along longitude > cycdist => no need to check
                else :
                    print "out of bounds"
                    break #out of bounds

            for nlat in np.arange(ilat-1,ilat-10,-1) :
                if nlat >= 0:
                    # cdist = dist(lats[nlat],lons[ilon],clat[ntrk-1,n],clon[ntrk-1,n])
                    # cdist1 = dist(lats[nlat],lons[ilon],lats[ilat],lons[ilon])
                    # cdist1 = np.absolute(lats[nlat] - clat[ntrk-1,n])*deg
                    cdist1 = np.absolute(lats[nlat] - lats[ilat])*deg
                    if cdist1 <= cycrad :
                        for nlon in np.arange(ilon+1,ilon+len(lons)-1) :
                            nnlon = nlon
                            if nlon >= len(lons):
                                nnlon = nlon - len(lons)
                            print " // -+ ", nlat,nnlon
                            cdist = dist(lats[nlat],lons[nnlon],clat[ntrk-1,n],clon[ntrk-1,n])
                            print cdist,cycrad
                            if cdist <= cycrad :
                                print " -+ ", nlat, nnlon
                                cycgrd[nnlon,nlat,ctime[ntrk-1,n]] = cycgrd[nnlon,nlat,ctime[ntrk-1,n]] + 1
                            else :
                                if cdist_pole > cycrad:
                                    print 'break'
                                    break
                        if cdist_pole > cycrad :
                            for nlon in np.arange(ilon-1,ilon-len(lons)+1,-1) :
                                nnlon = nlon
                                # if nlon < 0 :
                                #     nnlon = nlon + len(lons)
                                print " // -- ", nlat,nnlon
                                cdist = dist(lats[nlat],lons[nnlon],clat[ntrk-1,n],clon[ntrk-1,n])
                                print cdist,cycrad
                                if cdist <= cycrad :
                                    print " -- ", nlat, nnlon
                                    cycgrd[nnlon,nlat,ctime[ntrk-1,n]] = cycgrd[nnlon,nlat,ctime[ntrk-1,n]] + 1
                                else :
                                    print 'break'
                                    break
                    else :
                        print "dist along longitude > cycdist => no need to check 2"
                        print cdist1, cycrad
                        if cdist1 == cycrad :
                            print " == "
                        else :
                            print " != "
                            print cdist1 - cycrad
                        print nlat,ilat
                        break #dist along longitude > cycdist => no need to check
                else :
                    print "out of bounds"
                    break #out of bounds

            print cyr[ntrk-1,n],cmon[ntrk-1,n],cdate[ntrk-1,n],chour[ntrk-1,n], n, ntrk

            break #test - first cyclone onto the grid
    # print cyr[ntrk-1,n],cmon[ntrk-1,n],cdate[ntrk-1,n],chour[ntrk-1,n], n, ntrk
    # quit()
    #
    # break #test - first cyclone onto the grid


            #
            #
            # # start radius estimation
            # plon = [lon[ntrk-1,n]]
            # plat = [lat[ntrk-1,n]]
            #
            # #find the closest grid value to the cyclone center
            # iplon = find_nearest(lons,plon[0])
            # iplat = find_nearest(lats,plat[0])
            #
            # iplon0 = iplon
            # iplat0 = iplat
            # #
            # # print "initial iplon,iplat = ", iplon,iplat
            #
            # #check if it is a local minimum
            #
            # # create an extended array
            # slpext = np.zeros((lats.size+4,lons.size+4))
            # slpext[2:-2,2:-2] = slp[iyr[ntrk-1,n],it[ntrk-1,n],:,:]
            #
            #
            # slpext[:,:2] = slpext[:,-4:-2]
            # slpext[:,-2:] = slpext[:,2:4]
            #
            #
            # lon180 = np.where(lons==180)[0]
            #
            # slpext[1,2:lon180+2] = [slpext[3,ii] for ii in range(lon180+2,lons.size+2)]
            # slpext[0,2:lon180+2] = [slpext[4,ii] for ii in range(lon180+2,lons.size+2)]
            #
            # slpext[1,lon180+2:-2] = [slpext[3,ii] for ii in range(2,lon180+2)]
            # slpext[0,lon180+2:-2] = [slpext[4,ii] for ii in range(2,lon180+2)]
            #
            # slpext[:2,:2] = slpext[:2,-4:-2]
            # slpext[:2,-2:] = slpext[:2,2:4]
            #
            # # check if cyclone center is a loc min
            #
            # # slp16 = slp[iyr[ntrk-1,n],it[ntrk-1,n],iplat-2:iplat+3,plon_g-2:plon_g+3]
            #
            # #fix iplat and iplon for the extended slp matrix
            # iplatext = iplat+2
            # iplonext = iplon+2
            # slp16 = slpext[iplatext-2:iplatext+3,iplonext-2:iplonext+3]
            #
            # print slp16
            #
            # locmin = local_minima(slp16)
            # print locmin
            #
            # #check if the central point is a loc min
            # if locmin[2,2]:
            #     print 'good to go - the center is in the right place'
            # # check if one of the 8 surrounding points is a loc min
            # else:
            #     nmin = np.count_nonzero(locmin[1:4,1:4])
            #     if nmin == 1 :
            #         locmin2 = np.where(locmin[1:4,1:4])
            #         iplonext = iplonext-1+locmin2[1][0]
            #         iplatext = iplatext-1+locmin2[0][0]
            #     #if there are 2 loc mins within the surrounding 8 points take the slp min
            #     elif nmin == 2 :
            #         locmin2 = np.where(locmin[1:4,1:4])
            #         iplonext1 = iplonext-1+locmin2[1][0]
            #         iplatext1 = iplatext-1+locmin2[0][0]
            #         iplonext2 = iplonext-1+locmin2[1][1]
            #         iplatext2 = iplatext-1+locmin2[0][1]
            #         # if slp[iyr[ntrk-1,n],it[ntrk-1,n],iplat,iplon]>slp[iyr[ntrk-1,n],it[ntrk-1,n],iplat1,iplon1] :
            #         if slpext[iplatext2,iplonext2]>slpext[iplatext1,iplonext1] :
            #             iplonext = iplonext1
            #             iplatext = iplatext1
            #         else:
            #             iplonext = iplonext2
            #             iplatext = iplatext2
            #     else :
            #         print '!!!! ERROR: No loc min around the cyclone center, nmin =', nmin, "time step =",n,"/",nit
            #         miss   = -99.9
            #         cglon  = miss
            #         cglat  = miss
            #         cgslp  = miss
            #         fslp   = miss
            #         effrad = miss
            #         fout.write(" {:>14}{:8.3f}{:>7} {:7.3f}{:>7} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format
            #                   ("cglon=",lons[iplon],"cglat=",lats[iplat],"cgslp=",cgslp,"fslp=",fslp,"effrad=",effrad))
            #         fout.write(" {:>14}{:7.3f}\n".format("effrad=",effrad))
            #         continue
            #
            # iplon = iplonext - 2
            # iplat = iplatext - 2
            #
            # if iplon >= lons.size :
            #     iplon = iplon - lons.size
            # if iplon < 0 :
            #     iplon = iplon + lons.size
            #
            # if iplat < 0 :
            #     iplat = -iplat
            #     iplon = iplon + 180
            #     if iplon >= lons.size:
            #          iplon = iplon - lons.size
            # if iplat >= lats.size :
            #     iplat = 2*lats.size-iplat
            #     iplon = iplon + 180
            #     if iplon >= lons.size:
            #          iplon = iplon - lons.size
            #
            # print "initial cyc center slp= ", slp[iyr[ntrk-1,n],it[ntrk-1,n],iplat0,iplon0]
            # print "grid    cyc center slp= ", slp[iyr[ntrk-1,n],it[ntrk-1,n],iplat,iplon]
            # print "initial iplat =",iplat0,"final iplon=",iplon0
            # print "grid iplat =",iplat,"final plon=",iplon


            # if np.abs(iplat-iplat0)>1 or (np.abs(iplon-iplon0)>1 and np.abs(iplat-iplat0)<359):
            #     print "!!!! ERROR: check cyclone center location"
            #     print "Tracking cyclone center (lon,lat) = (",plon,",",plat,")"
            #     print "Current cyclone center (lon,lat) = (",lons[iplon],",",lats[iplat],")"
            #     quit()


            # alternative  approach may be
            # local_minima_locations = detect_local_minima(slp16)
            # print slp16[local_minima_locations[0][0],local_minima_locations[1][0]]

            #continue rad estimation using the new center
            #  plon = [162]
            #  plat = [33]
            # nplat0 = plat[0]
            # nplon0 = plon[0]

            # grid around cyclone center
            # dlon = 10
            # dlat = 0.5
            # dlon = 20
        #     dlat = 1.
        #     lonrange = np.arange(0., 360., dlon)
        #     # latrange = np.arange(90., 64.5, -dlat)
        #     latrange = np.arange(90., 64., -dlat)
        #
        #     rslp     = np.zeros_like(latrange)
        #     dslp     = np.zeros_like(latrange[:-1])
        #     lslp     = np.zeros_like(lonrange)
        #     flslp    = np.copy(lslp)
        #     flat     = np.copy(lslp)
        #     flon     = np.copy(lslp)
        #     rad      = np.copy(lslp)
        #
        #     nlon     = np.zeros((latrange.size,lonrange.size))
        #     nlat     = np.zeros_like(nlon)
        #
        #     mr = 2  # min radius = mr*dlat (deg.lat)
        #
        #     for i,ilon in enumerate(lonrange) :
        #         #  print "t=",n,"/",nit,"   ilon=",ilon
        #
        #          gridlat = np.copy(latrange)
        #          gridlon = np.zeros_like(gridlat)+ilon
        #
        #         #  nlat[:,i], nlon[:,i] = rotated_grid_transform(plat[0],plon[0],gridlat,gridlon,2)
        #          nlat[:,i], nlon[:,i] = rotated_grid_transform(lats[iplat],lons[iplon],gridlat,gridlon,2)
        #
        #          for j in range(nlon[:,i].size):
        #              if nlon[j,i] < 0:
        #                  nlon[j,i] = nlon[j,i]+360.
        #
        #          slpint = interpolate.interp2d(lons, lats, slp[iyr[ntrk-1,n],it[ntrk-1,n],:,:], kind='cubic')
        #
        #          for j,jlat in enumerate(latrange[:-1]) :
        #              rslp[j] = slpint(nlon[j,i],nlat[j,i])[0]
        #              dslp[j] = slpint(nlon[j+1,i],nlat[j+1,i]) - slpint(nlon[j,i],nlat[j,i])
        #
        #         #  print rslp
        #
        #         #  if any(dslp[0:mr-1]) < 0. :
        #         #      print "!!!! slp bug"
        #         #      print dslp
        #         #     #  ftime.sleep(20)
        #         #      quit()
        #
        #          for j in range(mr,latrange.size-1):
        #              if dslp[j] < 0 or j == latrange.size-2:
        #                  lslp[i] = rslp[j]
        #                 #  fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:4.1f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format("angle=",ilon,"lon=",nlon[j+1,i],"lat=",nlat[j+1,i],"rad=",(j+1)*dlat,"cslp=",slp0[0],"fslp=",lslp[i]))
        #                 #  ftime.sleep(10)
        #                  break
        #
        #
        #     # find the last closed isobar (fslp)
        #     fslp = np.amin(lslp)
        #     cgslp = slp[iyr[ntrk-1,n],it[ntrk-1,n],iplat,iplon]
        #     cycdepth = fslp-cgslp
        #
        #     # set rad = 0 for weak cyclones
        #     if cycdepth >= 0.5 :
        #         # find where last closed isobar fslp cross each radius
        #         gridlat = np.copy(latrange)
        #         for i,ilon in enumerate(lonrange) :
        #             for j in range(mr,latrange.size-1):
        #              slp2 = slpint(nlon[j+1,i],nlat[j+1,i])[0]
        #              slp1 = slpint(nlon[j,i],nlat[j,i])[0]
        #
        #              if slp1 <= fslp and slp2 > fslp :
        #                 rad[i] = j*dlat + (fslp-slp1)/(slp2-slp1)*dlat
        #                 # flat[i], flon[i] = rotated_grid_transform(plat[0],plon[0],[90-rad[i]],[ilon],2)
        #                 flat[i], flon[i] = rotated_grid_transform(lats[iplat],lons[iplon],[90-rad[i]],[ilon],2)
        #                 if flon[i] < 0:
        #                     flon[i] = flon[i]+360
        #                 # print "flat=",flat, "flon=",flon
        #                 flslp[i] = slpint(flon[0],flat[0])[0]
        #                 #  fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format("angle=",ilon,"lon=",flon[0],"lat=",flat[0],"rad=",rad[i],"cslp=",slp0[0],"fslp=",flslp[i]))
        #                 break
        #
        #              if j == latrange.size-1:
        #                 if slpint(nlon[mr,i],nlat[mr,i])[0] > fslp :
        #                     rad[i] = mr*dlat
        #                     # flat[i], flon[i] = rotated_grid_transform(plat[0],plon[0],rad[i],[ilon],2)
        #                     flat[i], flon[i] = rotated_grid_transform(lat[iplat],lons[iplon],rad[i],[ilon],2)
        #                     if flon[i] < 0:
        #                         flon[i] = flon[i]+360
        #                     flslp[i] = slpint(flon[0],flat[0])[0]
        #                     #  fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format("angle=",ilon,"lon=",flon[0],"lat=",flat[0],"rad=",rad[i],"cslp=",slp0[0],"fslp=",flslp[i]))
        #                 elif slpint(nlon[j],nlat[j])[0] < fslp :
        #                     print "!!!! fslp bug"
        #                     print fslp, slpint(nlon[j],nlat[j])[0]
        #                     #  ftime.sleep(20)
        #                     quit()
        #
        #             # fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format
        #             #           ("angle=",ilon,"lon=",flon[i],"lat=",flat[i],"rad=",rad[i],"cslp=",slp0[0],"fslp=",flslp[i]))
        #
        #
        #     else :
        #         for i,ilon in enumerate(lonrange) :
        #             rad[i]   = 0
        #             flat[i]  =  plat[0]
        #             flon[i]  =  plon[0]
        #             flslp[i]    = fslp
        #             # fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format
        #             #           ("angle=",ilon,"lon=",flon[i],"lat=",flat[i],"rad=",rad[i],"cslp=",slp0[0],"fslp=",flslp[i]))
        #
        #     Area = 0
        #     for i,ilon in enumerate(lonrange) :
        #         Area = Area + rad[i]**2
        #     Area = Area/lonrange.size
        #
        #     if all(rad) > 0 :
        #         effrad = np.sqrt(Area)
        #         fout.write(" {:>14}{:8.3f}{:>6} {:7.3f}{:>7} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format
        #                   ("cglon=",lons[iplon],"cglat=",lats[iplat],"cgslp=",cgslp,"fslp=",fslp,"effrad=",effrad))
        #
        #         for i,ilon in enumerate(lonrange) :
        #             fout.write(" {:>14}{:4.0f}{:>5}{:8.3f}{:>5} {:7.3f}{:>5} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format
        #                   ("angle=",ilon,"lon=",flon[i],"lat=",flat[i],"rad=",rad[i],"cslp=",cgslp,"fslp=",flslp[i]))
        #
        #     else:
        #         effrad = 0
        #         fout.write(" {:>14}{:8.3f}{:>7} {:7.3f}{:>7} {:7.3f} {:>5} {:7.3f} {:>5} {:7.3f}\n".format
        #                   ("cglon=",lons[iplon],"cglat=",lats[iplat],"cgslp=",cgslp,"fslp=",fslp,"effrad=",effrad))
        #
        #    # cyclones onto the regular grid
        #
        #     print "rad = ",rad
        #     print "effrad = ", effrad
        #
        #     # if effrad != 0 and iyr[ntrk-1,n]==yr:
        #     if effrad != 0 :
        #
        #        maxrad = np.amax(rad)+1
        #
        #        minlat = lats[iplat]-maxrad
        #        if minlat < -90 :
        #            minlat = -90.
        #            j1 = np.where(lats==-90.)[0]
        #        else :
        #            j1 = find_nearest(lats,minlat)
        #
        #        maxlat = plat[0]+maxrad
        #        if maxlat > 90 :
        #            maxlat = 90.
        #            j2 = np.where(lats==90.)[0]
        #        else :
        #            j2 = find_nearest(lats,maxlat)
        #
        #        if j1 < j2:
        #            jrange = range(j1,j2+1)
        #        else:
        #            jrange = range(j2,j1+1)
        #
        #        for j in jrange:
        #             for i,ilon in enumerate(lons):
        #
        #                 # gdist = dist(plat[0],plon[0],lats[j],ilon)
        #                 gdist = dist(lats[iplat],lons[iplon],lats[j],ilon)
        #                 gdist = gdist/deg
        #
        #                 if gdist <= maxrad:
        #                     # glat, glon = rotated_grid_transform(plat[0],plon[0],[lats[j]],[ilon],1)
        #                     glat, glon = rotated_grid_transform(lats[iplat],lons[iplon],[lats[j]],[ilon],1)
        #
        #                     if glon[0] < 0 :
        #                         glon[0] = glon[0] + 360
        #
        #                     indlon = find_nearest(lonrange,glon[0])
        #
        #                     if gdist <= rad[indlon] :
        #                         gridcyc[it[ntrk-1,n],j,i] = 1
        #
        #                     # print "to grid:",it[ntrk-1,n],plat[0],plon[0],effrad
        #
        #     #else:
        #         # j = find_nearest(lats,plat[0])
        #         # i = find_nearest(lons,plon[0])
        #
        #         # gridcyc[it[ntrk-1,n],j,i] = 2
        #
        #         # gridcyc[it[ntrk-1,n],iplat,iplon] = 2
        #
        #     gridcyc[it[ntrk-1,n],iplat,iplon] = 2
        #
        #
        #     # end radius estimation
        #     # break  # for the first step of the track
        #
        # # break  # for testing - first tack only


    f.close()
    print 'ftrk closed'

    #---NetCDF write---------------------------------------------------------------
print("Start NetCDF writing")

# varlist = np.zeros(16, dtype = {'names': ['name', 'outname', 'data', 'scale'],
#                                 'formats': ['a5', 'a5', '(241,480)f4', 'f4']} )
#
#
# varlist[0] = ("u","u",u,1)
# varlist[1] = ("v","v",v,1)
# for iv in range(varlist['name'].size) :


ncvar = "cycnum"
print 'ncvar=',ncvar
fcyc = '../ERAint/trkgrid/trk.%d.nc' % (yr)
ncout = Dataset(fcyc, 'w', format='NETCDF4')
ncout.description = "Cyclone number from  %s" % (ftrk)

# Using our previous dimension info, we can create the new time dimension
# Even though we know the size, we are going to set the size to unknown

dimnam=('lon','lat,'time')
varnam=['longitude','latitude','time',ncvar]

ncout.createDimension(dimnam[0], lons.size)
ncout.createDimension(dimnam[1], lats.size)
ncout.createDimension(dimnam[2], None)



# for nv in range(0, 3) :
#     # ncout_var = ncout.createVariable(varnam[nv], nc.variables[varnam[nv]].dtype,dimnam[nv])
#     ncout_var = ncout.createVariable(varnam[nv], 'f',dimnam[nv])
    # for ncattr in nc.variables[varnam[nv]].ncattrs():
        # ncout_var.setncattr(ncattr, nc.variables[varnam[nv]].getncattr(ncattr))
#print(nc.variables['latitude'].ncattrs())

#create  variables
 #lons
nv=0
ncout_var = ncout.createVariable(varnam[nv], 'f',dimnam[nv])
if varnam[nv] == 'longitude' :
    ncout_var.long_name = varnam[nv]
    ncout_var.units = 'degrees_east'
    ncout.variables[dimnam[0]][:] = lons

 #lats
nv=1
ncout_var = ncout.createVariable(varnam[nv], 'f',dimnam[nv])
if varnam[nv] == 'latitude' :
    ncout_var.long_name = varnam[nv]
    ncout_var.units = 'degrees_north'
    ncout.variables[dimnam[1]][:] = lats

 #time
nv=2
ncout_var = ncout.createVariable(varnam[nv], 'f',dimnam[nv])
ncout_var.long_name = varnam[nv]
if varnam[nv] == 'time' :
    ncout_var.calendar = 'gregorian'
    ncout_var.units = 'hours since 1900-01-01 00:00:0.0'
    #for one time step - test!
    tdelta =  datetime.datetime(cyr[ntrk-1,n],cmon[ntrk-1,n],cdate[ntrk-1,n],chour[ntrk-1,n]) - datetime.datetime(1900, 1, 1, 0)
    ncout.variables[dimnam[2]][:] = divmod(tdelta.total_seconds(),3600)[0]


print cycgrd[:,:,:].shape
print cycgrd[:,:,0].shape
print lons.size
print lats.size
# print time.size

ncout_var = ncout.createVariable(ncvar, 'f',dimnam[::-1])
print ncout_var.shape
ncout_var.long_name = 'number of cylones'
# ncout_var.scale_factor = varlist["scale"][iv]
# ncout_var.add_offset   = 0.
# ncout_var.units        = 'scale   %s' % varlist["scale"][iv]

#print qx.shape
#print ncout_var.shape
# ncout_var[:,:,:] = np.swapaxes(gridcyc,0,2)
ncout_var[0,:,:] = cycgrd[:,:,0]

ncout.close()

##---End NetCDF write---------------------------------------------------------------
