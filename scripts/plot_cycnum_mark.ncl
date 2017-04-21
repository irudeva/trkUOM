begin
;************************************************
; read in netCDF file
;************************************************
yr = 2005
t = 0
  slpyr = yr
  a = addfile("/Users/irudeva/work/DATA/ERAint/erain.mslp."+slpyr+".nc","r")
  slp = short2flt(a->msl(t,:,:))
  slp1 =0.01*slp
  copy_VarMeta(slp, slp1)



  b = addfile("../ERAint/trkgrid/trk."+yr+".nc","r")
  var = b->cycnum(t,:,:)
  lat = b->latitude
  lon = b->longitude

  printVarSummary(var)
  printVarSummary(slp1)
  ;printVarSummary(lat)
  ;print (var(:,359))
  ;print(lat)
  ;print(var(2,:))
  copy_VarMeta(slp, var)


  wks = gsn_open_wks("png","../ERAint/trkgrid/cycnum_m.t"+t)                  ; send graphics to PNG file

res = True
res@gsnDraw              = False             ; do not draw the plot
res@gsnFrame             = False             ; do not advance the frame

res@cnFillOn            = True        ; Turn on contour fill
res@cnLinesOn           = False       ; Turn off contour lines
res@cnFillPalette       = "BlueDarkRed18"  ; Set color map

res@gsnPolar   = "NH"

;res@mpLimitMode="LatLon"
res@mpMinLatF            = 45      ; range to zoom in on
res@mpMaxLatF            =  90.
;res@mpMinLonF            =  -60.
;res@mpMaxLonF            = 120.


res@tiMainString        = "t = " + t

;plot1 = gsn_csm_contour_map(wks,slp1,res)
plot1 = gsn_csm_map(wks,res)


mres = True

mres@gsMarkerIndex = 16
mres@gsMarkerSizeF = .003

;mres@gsnDraw              = False             ; do not draw the plot
;mres@gsnFrame             = False             ; do not advance the frame

;mres@cnFillOn            = False        ; Turn on contour fill
;mres@cnLinesOn           = True       ; Turn off contour lines
;mres@cnMonoFillColor     = True         ; Use single pattern color;
;mres@cnLineColor         = "green"      ; Turn off contour lines


do j = 0,dimsizes(var(:,0))-1
;print ("j="+j)

 indices = ind(var(j,:).gt.0)
 indcyc = ind(var(j,:).eq.22)
 ;print (indices)

 ;print(ismissing(indices))
 if(.not.any(ismissing(indices)))then
  mlon = lon(indices)
  mlat = new(dimsizes(mlon),float)
  mlat = lat(j)

  print (j+" mlon="+mlon)
  print (j+" mlat="+mlat)

  mres@gsMarkerColor      = 20

  str = unique_string("poly")            ;-- result is poly0-poly15
  plot1@$str$ = gsn_add_polymarker(wks,plot1,mlon,mlat,mres)
  ;gsn_polymarker(wks,plot1,mlon,mlat,mres)
  ;draw(plot1)
  ;frame(wks)
  ;status_exit(1)

  if(.not.any(ismissing(indcyc)))then


   cyclon = lon(indcyc)
   cyclat = lat(j)
   print (j+" cyclon="+cyclon)
   print (j+" cyclat="+cyclat)

   mres@gsMarkerColor   = 200
   str = unique_string("poly")            ;-- result is poly0-poly15
   plot1@$str$ = gsn_add_polymarker(wks,plot1,cyclon,cyclat,mres)
   ;gsn_polymarker(wks,plot1,cyclon,cyclat,mres)

   delete(cyclon)
   delete(cyclat)

  end if

  delete(mlon)
  delete(mlat)
 end if
 delete(indices)
 delete(indcyc)

end do

;status_exit(1)




;res1=True

;res1@gsnDraw              = False             ; do not draw the plot
;res1@gsnFrame             = False             ; do not advance the frame

;res1@cnFillOn            = False        ; Turn on contour fill
;res1@cnLinesOn           = True       ; Turn off contour lines
;res1@cnMonoFillColor     = True         ; Use single pattern color
;res1@cnLineColor         = "green"      ; Turn off contour lines

;res1@cnLevelSelectionMode = "ManualLevels"
;res1@cnMinLevelValF = 1.
;res1@cnMaxLevelValF = 23.
;res1@cnLevelSpacingF = 2.

;plot2 = gsn_csm_contour(wks,var,res1)         ; create a default plot


;overlay(plot1,plot2)



draw(plot1)
frame(wks)


end
