begin
;************************************************
; Selection
;************************************************
month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)

; time
year = (/2016,2016/)
yr = year(0)
mm = (/1,12/)  ;(/9,9/)
dd = (/0,31/)
hr = (/0,18/)

;ssn   / use "YYY" for whole year
ssn="YYY"

;ntrk
trk1 = 0; 5518 ;0
trk2 = 20000; 5518 ;20000

;reg
hs = "SH"
;lon = (/100-2.5,100+2.5/)
;lat = (/-39.-2.5,-39+2.5/)
;reg = "indoc"
;lon = (/0,360/)
;lat = (/-90,0/)
;reg = "hs"
lon = (/139-2.5,139+2.5/)
lat = (/-39.-2.5,-35+2.5/)
reg = "adel"

;reg2
ngen = 3
lon2 = (/0,120/)
lat2 = (/-60,-35/)
reg2 = "IO"
;lon2 = (/100-2.5,100+2.5/)
;lat2 = (/-39.-2.5,-39+2.5/)
;reg2 = "indoc"


selyr = 2016
seltrk = 5517

;************************************************
; define some variables
;************************************************

ntrk_total = 0

slparr1 = new(1, float)
slparr1 = 0.
slparr = slparr1

;************************************************
; read in netCDF files
;************************************************
; 1. slp
;************************************************

;t = 0
  ;slpyr = yr
  ;a = addfile("/Users/irudeva/work/DATA/ERAint/erain.mslp."+slpyr+".nc","r")

  ;;slp = short2flt(a->msl(t,:,:))
  ;slp = short2flt(a->msl(:,:,:))
  ;slp1 =0.01*slp
  ;copy_VarMeta(slp, slp1)

  ;timeslp = a->time
  ;dateslp = cd_calendar(timeslp, 0)

  ;print(dateslp)
  ;status_exit(1)

  ;do nt=0, dimsizes(timeslp)-1
  ;  if (dateslp(nt,0).eq.yr.and.dateslp(nt,1).eq.mm.and.dateslp(nt,2).eq.dd.and.dateslp(nt,3).eq.hr)then
  ;   tslp = nt
  ;   break
  ;  end if
  ;end do

  ;print(n)
  ;print (cd_calendar(timeslp(n), 0))
  ;status_exit(1)



;************************************************
; Plotting - preparation
;************************************************

date_str1 = sprinti("%0.4i", year(0))+sprinti("%0.2i", mm(0))+ \
             sprinti("%0.2i", dd(0)) +"_"+sprinti("%0.2iZ", hr(0))
date_str2 = sprinti("%0.4i", year(1))+sprinti("%0.2i", mm(1))+ \
             sprinti("%0.2i", dd(1)) +"_"+sprinti("%0.2iZ", hr(1))

wks = gsn_open_wks("png","../ERAint/trkgrid/tracks."+ssn+"."+date_str1+"_"+date_str2+"."+reg+"_"+reg2+"."+hs )                  ; send graphics to PNG file

res = True
res@gsnDraw              = False             ; do not draw the plot
res@gsnFrame             = False             ; do not advance the frame

;res@cnFillOn            = True        ; Turn on contour fill
res@cnLinesOn           = True       ; Turn on/off contour lines
res@cnFillPalette       = "BlueDarkRed18"  ; Set color map

res@gsnPolar   = hs

;res@mpLimitMode="LatLon"
res@mpMinLatF            = -90      ; range to zoom in on
res@mpMaxLatF            =  -10.
;res@mpMinLonF            =  -60.
;res@mpMaxLonF            = 120.

;res@tiMainString        = "From "+date_str1+" to "+date_str2 + "~C~~Z90~ reg lon: ("+lon(0)+"E to "+lon(1)+"E)  lat: ("+lat(0)+"N to "+lat(1)+"N)" + "~C~~Z90~ reg2 lon: ("+lon2(0)+"E to "+lon2(1)+"E)  lat: ("+lat2(0)+"N to "+lat2(1)+"N)"
res@tiMainFontHeightF   = 0.015

;plot1 = gsn_csm_contour_map(wks,slp1,res)
;plot1 = gsn_csm_contour_map_polar(wks,slp1(tslp,:,:),res)
plot1 = gsn_csm_map(wks,res)


;draw(plot1)
;frame(wks)


;status_exit(1)

plres = True

 plres@gsLineColor      = "black"                  ; color of lines
 plres@gsLineThicknessF = 2.0                    ; thickness of lines

pl1res = True   ; for a selected track

  pl1res@gsLineColor      = "red"                  ; color of lines
  pl1res@gsLineThicknessF = 3.0                    ; thickness of lines

mres = True

 mres@gsMarkerIndex = 16
 mres@gsMarkerSizeF = .003

txres               = True
 txres@txFontHeightF = 0.01

cycres             = True
 cycres@gsLineColor      = "red"                  ; color of lines
 cycres@gsLineThicknessF = 2.0
 ;************************************************
 ; read in netCDF files
 ;************************************************
 ; 2. cyclones
 ;************************************************

 do iyr =year(0),year(1)

 ftrk = "../ERAint/trkgrid/trk."+iyr+".nc"
 print("ftrk= " + ftrk)
 b = addfile(ftrk,"r")

 trklen = b->trklen
 trktime= b->trktime
 trklon = b->trklon
 trklat = b->trklat
 trkslp = b->trkslp
 trkrad = b->trkrad


 ; what's in netcdf:
 ;varlist[1] = ("trktime",     "time of the current track point",  'i',                    time0, trktime[:ntrk,:])
 ;varlist[0] = ( "trklon","longitude of the current track point",  'f',           'degrees_east',  trklon[:ntrk,:])
 ;varlist[2] = ( "trklat", "latitude of the current track point",  'f',          'degrees_north',  trklat[:ntrk,:])
 ;varlist[3] = ( "trkslp", "slp in cyclone center from trk file",  'f',                    'hPa',  trkslp[:ntrk,:])
 ;varlist[4] = ( "trklpl",             "Laplacian from trk file",  'f','hPa/degrees_latitude**2',  trklpl[:ntrk,:])
 ;varlist[5] = ( "trkdp",                  "depth from trk file",  'f',                    'hPa',   trkdp[:ntrk,:])
 ;varlist[6] = ( "trkrad",                "radius from trk file",  'f',       'degrees_latitude',   trkrd[:ntrk,:])
 ;varlist[7] = ( "trkiop",                      "cyclone status", 'i2','0->strong/1->weak/10->str_open/11->weak_open',  trkiop[:ntrk,:])



;************************************************
; masking
;************************************************


;create mask
  trkmask = trktime
  trkmask = 0

do ntrk = 0,dimsizes(trklen)-1
  ;print ("ntrk="+ntrk+" n= "+ trklen(ntrk))
  if (ntrk+1.ge.trk1.and.ntrk+1.le.trk2)then

   nit = trklen(ntrk)
   trkdate = cd_calendar(trktime(ntrk,:nit-1), 0)
   do n = 0,nit-1
    ;print(trkdate(n,1)+" ")
    ;detect season
    if(any((/1,2,12/).eq.trkdate(n,1)))then
      trkssn = "DJF"
    end if
    if(any((/3,4,5/).eq.trkdate(n,1)))then
      trkssn = "MAM"
    end if
    if(any((/6,7,8/).eq.trkdate(n,1)))then
      trkssn = "JJA"
    end if
    if(any((/9,10,11/).eq.trkdate(n,1)))then
      trkssn = "SON"
    end if

    if (trkdate(n,0).ge.year(0).and.trkdate(n,0).le.year(1).and.trkdate(n,0).eq.iyr)then
     if (trkdate(n,1).ge.mm(0).and.trkdate(n,1).le.mm(1))then
      if (trkdate(n,2).ge.dd(0).and.trkdate(n,2).le.dd(1))then
       if (trkdate(n,3).ge.hr(0).and.trkdate(n,3).le.hr(1))then

        if (ssn.eq."YYY".or.trkssn.eq.ssn)then

          if (trklat(ntrk,n).ge.lat(0).and.trklat(ntrk,n).le.lat(1))then
           if (trklon(ntrk,n).ge.lon(0).and.trklon(ntrk,n).le.lon(1)) then
            trkmask(ntrk,n) = 1
            ;print (n+" for reg  "+ntrk)
            ;print (n+"  "+ntrk)
            ;print (trklon(ntrk,n))
            ;print (trklat(ntrk,n))
            ;status_exit(1)
           end if  ;lon
          end if  ;lat
          ; for generation (n<=ngen) in reg 2
          if (n.le.ngen)then
            if (trklat(ntrk,n).ge.lat2(0).and.trklat(ntrk,n).le.lat2(1))then
             if (trklon(ntrk,n).ge.lon2(0).and.trklon(ntrk,n).le.lon2(1)) then
              trkmask(ntrk,n) = 2
              ;print (n+" for reg2  "+ntrk)
              ;print (trklon(ntrk,n))
              ;print (trklat(ntrk,n))
             end if  ;lon
            end if  ;lat
          end if

        end if

       end if ;hr
      end if ;dd
     end if ;mm
    end if ;year
   end do  ;nit

   delete(trkdate)

  end if  ; trk1/trk2
end do  ; ntrk


;************************************************
; Plotting
;************************************************

do ntrk = 0,dimsizes(trklen)-1
 nit = trklen(ntrk)
 if(any(trkmask(ntrk,:).eq.1).and.any(trkmask(ntrk,:).eq.2))then
  print ("  "+ntrk)
  ntrk_total = ntrk_total+1

  ;************************************************
  ; add polyline
  ;************************************************

  ;print (ntrk)
  ;print (trklon(ntrk,:nit-1))
  ;print (trklat(ntrk,:nit-1))
  str = unique_string("poly")            ;-- result is poly0-poly15
  if(selyr.eq.2016.and.ntrk.eq.seltrk)then
   plot1@$str$ = gsn_add_polyline(wks,plot1,trklon(ntrk,:nit-1),trklat(ntrk,:nit-1),pl1res)
  else
   plot1@$str$ = gsn_add_polyline(wks,plot1,trklon(ntrk,:nit-1),trklat(ntrk,:nit-1),plres)
  end if

  ;************************************************
  ; add polymarker
  ;************************************************
  mres@gsMarkerColor      = 20
  str = unique_string("marker")            ;-- result is poly0-poly15
  ;plot1@$str$ = gsn_add_polymarker(wks,plot1,trklon(ntrk,:nit-1),trklat(ntrk,:nit-1),mres)

  ;************************************************
  ; add date
  ;************************************************

  nt = 0  ;which date to show
  trkdate = cd_calendar(trktime(ntrk,:nit-1), 0)

  yy     = tointeger(trkdate(nt,0))    ; Convert to integer for
  month  = tointeger(trkdate(nt,1))    ; use sprinti
  day    = tointeger(trkdate(nt,2))
  hour   = tointeger(trkdate(nt,3))
  date2plot = sprinti("%0.4i", yy)+month_abbr(month)+ \
             sprinti("%0.2i", day) +"_"+sprinti("%0.2iZ", hour)

  txt = unique_string("txt")
  ;plot1@$txt$ = gsn_add_text(wks,plot1,date2plot,trklon(ntrk,nt)+3.,trklat(ntrk,nt)-3.,txres)

  ;************************************************
  ; pdf stat
  ;************************************************
   trkdr = new(nit-4,float)
   do in =0, nit-5
    trkdr(in) = (trkslp(ntrk,in+4)-trkslp(ntrk,in))*sin(60*3.14/180)/sin(3.14*abs(avg((/trklat(ntrk,in+4),trklat(ntrk,in)/)))/180)
    tmp = trkslp(ntrk,in+4)-trkslp(ntrk,in)
    print("dr: "+in+"  "+trkdr(in)+ " "+tmp)
    print("slp: "+in+"  "+trkslp(ntrk,in)+"  "+trkslp(ntrk,in+4))
   end do


   if(any(slparr1.ne.0))then
    delete(slparr)
    delete(drarr)
    slparr = array_append_record(slparr1, min(trkslp(ntrk,ind(trkmask(ntrk,:nit-1).eq.1))), 0)
    drarr = array_append_record(drarr1,min(trkdr), 0)
    if(selyr.eq.2016.and.ntrk.eq.seltrk)then
     seldr  = min(trkdr)
     selslp = min(trkslp(ntrk,ind(trkmask(ntrk,:nit-1).eq.1)))
    end if
    delete(slparr1)
    delete(drarr1)
    slparr1 = slparr
    drarr1  = drarr
   else
    delete(slparr1)
    slparr1 = trkslp(ntrk,ind(trkmask(ntrk,:nit-1).eq.1))
    drarr1 = min(trkdr)
    drarr  = drarr1

   end if
  ;************************************************
  print ("  "+ntrk+"  dr="+min(trkdr))
  delete(trkdate)
  delete(trkdr)

  ;status_exit(1)
  ;break
 end if

end do  ;ntrk

delete(trklen)
delete(trktime)
delete(trklon)
delete(trklat)
delete(trkslp)
delete(trkrad)
delete(trkmask)

end do   ; yr

;draw(plot1)
;frame(wks)

;************************************************
; pdf stat
;************************************************

opt         = True
  opt@bin_min = 976.
  opt@bin_max = 1024.
  nbins = (opt@bin_max - opt@bin_min)/2.

slppdf = pdfx(slparr, toint(nbins), opt)

print("PDF slp:")
print(slppdf@bin_center+" "+slppdf)

; deepening rate
opt@bin_min = -38.
opt@bin_max = 0.
nbins = (opt@bin_max - opt@bin_min)/2

drpdf = pdfx(drarr(ind(drarr.lt.0)), toint(nbins), opt)

print("PDF dr:")
print(drpdf@bin_center+" "+drpdf)


;------------------------------------------------
; Plot each PDF
;------------------------------------------------
  plot    = new ( 2, "graphic")

  ;wks1  = gsn_open_wks ("png","../ERAint/trkgrid/pdf."+ssn+"."+date_str1+"_"+date_str2+"."+reg+"_"+reg2+"."+hs )             ; send graphics to PNG file

  zres  = True
  zres@gsnDraw                = False
  zres@gsnFrame               = False

  zres@vpHeightF = 0.4

  zres@gsnXYBarChart          = True            ; Create bar plot
  zres@gsnXYBarChartBarWidth  = 1.5

  zres@xyLineThicknessF       = 2

  zres@gsnLeftStringFontHeightF     = 0.02
  zres@gsnRightStringFontHeightF    = 0.02
  zres@tiYAxisFontHeightF           = 0.015


  zres@gsnLeftString          = "SLP"
  zres@tiYAxisString          = "PDF (%)"
  zres@tiXAxisString          = "hPa"
  zres@gsnRightString         = "N of trkacks "+ntrk_total
  plot(0) = gsn_csm_xy (wks, slppdf@bin_center, slppdf, zres)


  ; adding a polyline

    xval = selslp

    getvalues plot(0)
       "trYMinF": ymin
       "trYMaxF": ymax
    end getvalues

   xbars = (/xval, xval/)
   ybars = (/ymin, ymax/)

   slpline = gsn_add_polyline(wks,plot(0),xbars,ybars,cycres)

   ; end add polyline

   ; plotting DR
  zres@gsnXYBarChartBarWidth  = 1.5
  zres@gsnYRefLine  = 0.
  zres@gsnLeftString          = "Maximum deepening rate"
  zres@tiYAxisString          = ""
  zres@gsnRightString         = "N trkacks "+num(drarr.lt.0)
  zres@tiXAxisString          = "hPa / 24 hours"
  plot(1) = gsn_csm_xy (wks, drpdf@bin_center, drpdf, zres)

  ; adding a polyline

    xval = seldr

    getvalues plot(1)
       "trYMinF": ymin
       "trYMaxF": ymax
    end getvalues

   xbars = (/xval, xval/)
   ybars = (/ymin, ymax/)

   drline = gsn_add_polyline(wks,plot(1),xbars,ybars,cycres)

   ; end add polyline

; panel plot
  resP    = True
  resP@gsnMaximize      = True                   ; use full page
  resP@gsnPanelBottom   = 0.06                   ; save space at bottom
  resP@gsnPanelLeft     = 0.06                   ; save space at bottom
  resP@gsnPanelYWhiteSpacePercent  = 5.


  ;resP@txString        = res@tiMainString
  resP@gsnPanelRowSpec = True
  ;gsn_panel(wks,plot(:1),(/2,0/),resP)
  ;gsn_panel(wks1,plot(:1),2,resP)
  gsn_panel(wks,(/plot1,plot(0),plot(1)/),(/1,2/),resP)


end
