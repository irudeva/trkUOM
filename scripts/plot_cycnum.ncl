begin
;************************************************
; Selection
;************************************************
month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                    "Oct","Nov","Dec"/)

; time

 ;ssn   / use "YYY" for whole year
  ssn="DJF"

 year = (/2016,2016/)
 yr = year(0)
 mm = (/1,12/)  ;(/9,9/)
 dd = (/0,31/)
 hr = (/0,18/)

 if(mm(0).ne.1.and.ssn.eq."YYY")then
  print("WARNING: Check the season and time settings!")
  print("  season = "+ssn)
  print("  month = "+month_abbr(mm(0))+" to "+month_abbr(mm(1)))
  print("  Check 'cyr' var!!!")
  sleep(1)
 end if


;ntrk
trk1 = 0; 5518 ;0
trk2 = 20000; 5518 ;20000

;reg

hs ="SH"

rlon = (/-181,361/)
rlat = (/-90,90/)
reg = ""


;rlon = (/100,150/)
;rlat = (/-60,-30/)
;reg = "AU"+"."

;************************************************
; Plotting - preparation
;************************************************

date_str1 = sprinti("%0.4i", year(0))+sprinti("%0.2i", mm(0))+ \
             sprinti("%0.2i", dd(0)) +"_"+sprinti("%0.2iZ", hr(0))
date_str2 = sprinti("%0.4i", year(1))+sprinti("%0.2i", mm(1))+ \
             sprinti("%0.2i", dd(1)) +"_"+sprinti("%0.2iZ", hr(1))

regstr = reg+hs
wks = gsn_open_wks("png","../ERAint/trkgrid/cycnum."+ssn+"."+date_str1+"_"+date_str2+"."+regstr )                  ; send graphics to PNG file

res = True
res@gsnDraw              = False             ; do not draw the plot
res@gsnFrame             = False             ; do not advance the frame


res@cnFillOn            = True        ; Turn on contour fill
res@cnLinesOn           = True       ; Turn on/off contour lines
res@cnFillPalette       = "BlueDarkRed18"  ; Set color map

res@gsnPolar   = hs

;res@mpLimitMode="LatLon"
res@mpMinLatF            = -90      ; range to zoom in on
res@mpMaxLatF            = -10.
;res@mpMinLonF           = -60.
;res@mpMaxLonF           = 120.


res@cnLevelSelectionMode = "ManualLevels"
res@cnMinLevelValF   = 0
res@cnMaxLevelValF   = .7
res@cnLevelSpacingF  = .05

res@tiMainString        = date_str1 +" to "+ date_str2


;************************************************
; read in netCDF files
;************************************************

  do iyr =year(0),year(1)
  fin = "../ERAint/trkgrid/cycloc_rd."+iyr+".nc"
  print("fin = "+fin)

  b = addfile(fin,"r")
  cycloc = b->cycloc(:,:,:)
  cycarea = b->cycarea(:,:,:)
  lat = b->latitude
  lon = b->longitude

  timecyc = b->time
  datecyc = cd_calendar(timecyc, 0)

  ;print(datecyc)
  ;status_exit(1)

  if(iyr.eq.year(0))then

   cycgrd = new((/dimsizes(lat),dimsizes(lon)/),float)
   cycgrd = 0
   cycgrd!1 = "lon"
   cycgrd!0 = "lat"

   cycgrd&lon = lon
   cycgrd&lat = lat

   ntstp = 0

  end if

  do nt=0, dimsizes(timecyc)-1
    cyr = datecyc(nt,0)
    if(any((/1,2,12/).eq.datecyc(nt,1)))then
      if(datecyc(nt,1).eq.12.and.ssn.ne."YYY") then
       cyr = cyr+1
      end if
      cssn = "DJF"
    end if
    if(any((/3,4,5/).eq.datecyc(nt,1)))then
      cssn = "MAM"
    end if
    if(any((/6,7,8/).eq.datecyc(nt,1)))then
      cssn = "JJA"
    end if
    if(any((/9,10,11/).eq.datecyc(nt,1)))then
      cssn = "SON"
    end if
  ;print(datecyc(nt,:))

    if (cyr.ge.year(0).and.cyr.le.year(1).and.cyr.eq.iyr)then
     if (datecyc(nt,1).ge.mm(0).and.datecyc(nt,1).le.mm(1))then
      if (datecyc(nt,2).ge.dd(0).and.datecyc(nt,2).le.dd(1))then
       if (datecyc(nt,3).ge.hr(0).and.datecyc(nt,3).le.hr(1))then

        if (ssn.eq."YYY".or.cssn.eq.ssn)then

          cycgrd = cycgrd+where(cycarea(nt,:,:).gt.0, 1, 0)
          ntstp  = ntstp + 1

         end if

        end if ;hr
       end if ;dd
      end if ;mm
     end if ;year
    ;end do  ;nit
  ;  end if ;
  ;  print ("ku-ku")
  end do

  delete(cycloc)
  delete(cycarea)
  delete(timecyc)
  delete(datecyc)

  end do  ; year



  ;************************************************
  ; 3. Plotting
  ;************************************************


 ;tstr = cd_calendar(timecyc(tcyc), 3)

 ;month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
 ;                     "Oct","Nov","Dec"/)
 ;year   = tointeger(datecyc(tcyc,0))    ; Convert to integer for
 ;month  = tointeger(datecyc(tcyc,1))    ; use sprinti
 ;day    = tointeger(datecyc(tcyc,2))
 ;hour   = tointeger(datecyc(tcyc,3))
 ;date_strm = sprinti("%0.4i", year)+month_abbr(month)+ \
 ;            sprinti("%0.2i", day) +"_"+sprinti("%0.2iZ", hour)
 ;date_str = sprinti("%0.4i", year)+sprinti("%0.2i", month)+ \
;            sprinti("%0.2i", day) +"_"+sprinti("%0.2iZ", hour)

 ;print(date_str + "")
 ;status_exit(1)

cycgrd = cycgrd/ntstp

regcycgrd = cycgrd
regcycgrd = 0
regcycgrd(ind(lat.gt.rlat(0).and.lat.le.rlat(1)),ind(lon.gt.rlon(0).and.lon.le.rlon(1))) = cycgrd(ind(lat.gt.rlat(0).and.lat.le.rlat(1)),ind(lon.gt.rlon(0).and.lon.le.rlon(1)))

plot1 = gsn_csm_contour_map_polar(wks,regcycgrd,res)

draw(plot1)
frame(wks)

status_exit(1)

end
