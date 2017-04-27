begin
;************************************************
; Selection
;************************************************
month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)

; time
year = (/1979,2017/)
yr = 2016 ; year(0)
mm = (/0,12/)  ;(/9,9/)
dd = (/0,31/)
hr = (/0,18/)

;ssn   / use "YYY" for whole year
ssn="SON"

;ntrk
trk1 = 0; 5518 ;0
trk2 = 20000; 5518 ;20000

;reg
hs = "SH"
lon = (/100-2.5,100+2.5/)
lat = (/-39.-2.5,-39+2.5/)
reg = "indoc"


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
  ; 2. cyclones
  ;************************************************

  b = addfile("../ERAint/trkgrid/trk."+yr+".nc","r")

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
; Plotting - preparation
;************************************************

;create mask
  trkmask = trktime
  trkmask = 0



date_str1 = sprinti("%0.4i", year(0))+sprinti("%0.2i", mm(0))+ \
             sprinti("%0.2i", dd(0)) +"_"+sprinti("%0.2iZ", hr(0))
date_str2 = sprinti("%0.4i", year(1))+sprinti("%0.2i", mm(1))+ \
             sprinti("%0.2i", dd(1)) +"_"+sprinti("%0.2iZ", hr(1))

wks = gsn_open_wks("png","../ERAint/trkgrid/tracks."+ssn+"."+date_str1+"_"+date_str2+"."+reg+"."+hs )                  ; send graphics to PNG file

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

res@tiMainString        = "From "+date_str1+" to "+date_str2 + "~C~~Z90~ reg lon: ("+lon(0)+"E to "+lon(1)+"E)  lat: ("+lat(0)+"N to "+lat(1)+"N)"
res@tiMainFontHeightF   = 0.015

;plot1 = gsn_csm_contour_map(wks,slp1,res)
;plot1 = gsn_csm_contour_map_polar(wks,slp1(tslp,:,:),res)
plot1 = gsn_csm_map(wks,res)


;draw(plot1)
;frame(wks)


;status_exit(1)

plres = True

 plres@gsLineColor      = "red"                  ; color of lines
 plres@gsLineThicknessF = 2.0                    ; thickness of lines

mres = True

 mres@gsMarkerIndex = 16
 mres@gsMarkerSizeF = .003

txres               = True
 txres@txFontHeightF = 0.01

;************************************************
; masking
;************************************************

do ntrk = 0,dimsizes(trklen)-1
  ;print ("ntrk="+ntrk+" n= "+ trklen(ntrk))
  if (ntrk+1.ge.trk1.and.ntrk+1.le.trk2)then

   nit = trklen(ntrk)
   trkdate = cd_calendar(trktime(ntrk,:nit-1), 0)
   do n = 0,nit-1
    print(trkdate(n,1)+" ")
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

    if (trkdate(n,0).ge.year(0).and.trkdate(n,0).le.year(1))then
     if (trkdate(n,1).ge.mm(0).and.trkdate(n,1).le.mm(1))then
      if (trkdate(n,2).ge.dd(0).and.trkdate(n,2).le.dd(1))then
       if (trkdate(n,3).ge.hr(0).and.trkdate(n,3).le.hr(1))then

        if (ssn.eq."YYY".or.trkssn.eq.ssn)then

          if (trklat(ntrk,n).ge.lat(0).and.trklat(ntrk,n).le.lat(1))then
           if (trklon(ntrk,n).ge.lon(0).and.trklon(ntrk,n).le.lon(1)) then
            trkmask(ntrk,n) = 1
            ;print (n+"  "+ntrk)
            ;print (trklon(ntrk,n))
            ;print (trklat(ntrk,n))
            ;status_exit(1)
           end if  ;lon
          end if  ;lat

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
 if(any(trkmask(ntrk,:).eq.1))then
  print ("  "+ntrk)

  ;************************************************
  ; add polyline
  ;************************************************

  ;print (ntrk)
  ;print (trklon(ntrk,:nit-1))
  ;print (trklat(ntrk,:nit-1))
  str = unique_string("poly")            ;-- result is poly0-poly15
  plot1@$str$ = gsn_add_polyline(wks,plot1,trklon(ntrk,:nit-1),trklat(ntrk,:nit-1),plres)

  ;************************************************
  ; add polymarker
  ;************************************************
  mres@gsMarkerColor      = 20
  str = unique_string("marker")            ;-- result is poly0-poly15
  plot1@$str$ = gsn_add_polymarker(wks,plot1,trklon(ntrk,:nit-1),trklat(ntrk,:nit-1),mres)

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

  delete(trkdate)

  ;status_exit(1)
  ;break
 end if

end do  ;ntrk

draw(plot1)
frame(wks)

end
