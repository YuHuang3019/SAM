; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "open_wks.ncl"

;************************************************
begin
;************************************************
; read in data
;************************************************
 filename = "./../../../OUT_STAT/GATE_IDEAL_S_2048x2048x256_100m_2s.nc"

;wks = open_wks("NCGM")
;wks = open_wks("X11")
wks = open_wks("PS")

time_max = 24.

 f     = addfile (filename , "r")             ; add file

;************************************************
; plotting parameters
;************************************************

 res                   = True                       ; plot mods desired

; res@tiMainString      = fld@long_name             ; add title
 res@trYReverse        = False                       ; reverse Y-axis
; res@xyDashPatterns    = 1                         ; choose dash patterns
 res@pmLegendDisplayMode    = "Always"
 res@pmLegendSide           = "Top"                 ; Change location of 
 res@pmLegendParallelPosF   = .60                   ; move units right
 res@pmLegendOrthogonalPosF = -0.6                  ; more neg = down
 
 res@pmLegendWidthF         = 0.12                  ; Change width and
 res@pmLegendHeightF        = 0.15                  ; height of legend.
 res@lgLabelFontHeightF     = .016                   ; change font height
 res@lgPerimOn              = False                 ; no box around
 res@xyLineThicknessF = 3.0 

 field_names = (/"PREC","LHF","SHF"/)
 res@xyExplicitLegendLabels = (/"Rainfall","LHF","SHF"/)

 data = new((/dimsizes(field_names),dimsizes(f->time({:time_max/24.}))/),float)

 do i=0,dimsizes(field_names)-1
   data(i,:) = f->$field_names(i)$({:time_max/24.})
 end do
 res@xyLineColors    =   (/1,2,4,3,5,6,7/) 
 plot  = gsn_csm_xy (wks,data&time*24.,data,res) ; create plot



end

