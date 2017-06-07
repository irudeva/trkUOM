begin

ncyc = (/3,2,3,6,24,29,30,36,28,14,5,3/)
mon = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

; print (sum(ncyc))

;********************************************************
; create the plot
;********************************************************
  wks = gsn_open_wks("png","../ERAint/output/AdelStorm/Ncyc_mon.Adel")                  ; send graphics to PNG file
  res                       = True                 ; plot mods desired
  res@gsnFrame              = False                ; don't advance frame yet
  res@gsnXYBarChart         = True                 ; turn on bar chart
  res@gsnXYBarChartBarWidth = 0.75                 ; change bar widths
  ; res@gsnXYBarChartColors = (/"firebrick","red","orange","green", \
  ; "navy","blue","SkyBlue","SlateBlue"/)            ; choose colors


  ; res@tmXBOn                = False                ; turn off tickmarks at bot
  res@trYMinF               = 0                    ; bring bars down to zero
  ; res@trXMinF               = 0                    ; adds space on either end
  ; res@trXMaxF               = 9                    ; of the 1st and last bars

  ; res@tiMainString          = " "

  res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
  res@tmXBValues    = ispan(0,11,1)
  res@tmXBLabels    =mon

  res@tmXBLabelAngleF = 45


  plot = gsn_csm_xy (wks,ispan(0,11,1),ncyc,res)                  ; create plot

frame(wks)



end
