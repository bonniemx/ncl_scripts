; 原来的脚本来自pigliu。感谢师兄，
; 2020年10月21日 马欣
; 内容：xy线图，缺省值的检测和替换，年平均和区域平均的计算。
; ____________________


begin
dir = "/home/max/cmip6/fracLUT/all/"
cesm = addfile(dir+"fracLut_Emon_CESM2_land-hist_r1i1p1f1_gn_185001-201512.nc", "r")    ;NOTE 时间范围
ipsl = addfile(dir+"fracLut_Emon_IPSL-CM6A-LR_land-hist_r1i1p1f1_gr_185001-201412.nc", "r")
cmcc = addfile(dir+"fracLut_Emon_CMCC-ESM2_land-hist_r1i1p1f1_gn_185001-201412.nc", "r")
cnrm = addfile(dir+"fracLut_Emon_CNRM-ESM2-1_land-hist_r1i1p1f2_gr_185001-201412.nc", "r")
gfdl = addfile(dir+"fracLut_Emon_GFDL-ESM4_land-hist_r1i1p1f1_gr1_180001-201412.nc", "r")     ;NOTE 时间范围
mpi = addfile(dir+"fracLut_Emon_MPI-ESM1-2-LR_land-hist_r1i1p1f1_gn_185001-201412.nc", "r")      


cesm_fra = cesm ->fracLut
;print(cesm_fra&landuse)
cesm_fracLut = cesm_fra(0:1979,:,:,:) ;cesm:1992,but ipsl:1980
ipsl_fracLut = ipsl ->fracLut
cmcc_fracLut = cmcc ->fracLut
cnrm_fracLut = cnrm ->fracLut
gfdl_fra = gfdl ->fracLut
gfdl_fracLut = gfdl_fra(600:,:,:,:)   ;600:就是从600到结尾的意思。
mpi_fracLut = mpi ->fracLut  

delete(cesm_fra)
delete(gfdl_fra)
delete(cesm)
delete(ipsl)
delete(cmcc)
delete(cnrm)
delete(gfdl)
delete(mpi)

;------------------- 将nan转换成missingvalue -------------------
if (any(isnan_ieee(cesm_fracLut))) then
      value = 1.e20
      replace_ieeenan (cesm_fracLut, value, 0)
      cesm_fracLut@_FillValue = value
end if
 ;print(num(ismissing(fracLut(0,2,:,:))))

if (any(isnan_ieee(ipsl_fracLut))) then
      value = 1.e20
      replace_ieeenan (ipsl_fracLut, value, 0)
      ipsl_fracLut@_FillValue = value
end if

if (any(isnan_ieee(cmcc_fracLut))) then
      value = 1.e20
      replace_ieeenan (cmcc_fracLut, value, 0)
      cmcc_fracLut@_FillValue = value
end if

if (any(isnan_ieee(cnrm_fracLut))) then
      value = 1.e20
      replace_ieeenan (cnrm_fracLut, value, 0)
      cnrm_fracLut@_FillValue = value
end if

;--------gfdl这货跟别人不一样,缺省是0但是写的是1e20-----------
;print(num(ismissing(gfdl_fracLut(0,0,:,:))))
gfdl_fracLut@_FillValue = 0.0 
replace_ieeenan(gfdl_fracLut,1.e20,0)
gfdl_fracLut@_FillValue = 1.e20

if (any(isnan_ieee(mpi_fracLut))) then
      value = 1.e20
      replace_ieeenan (mpi_fracLut, value, 0)
      mpi_fracLut@_FillValue = value
end if

;------------------- 计算年平均 -------------------
yyyymm = yyyymm_time(1850, 2014, "integer")             ;nyear=2015-1850+1=166,而不是我想象里的165
cesm_fracLut_ann = month_to_annual_weighted(yyyymm, cesm_fracLut, 1)      ;1 for 加权的结果除以12得到年平均的结果
ipsl_fracLut_ann = month_to_annual_weighted(yyyymm, ipsl_fracLut, 1)      
cmcc_fracLut_ann = month_to_annual_weighted(yyyymm, cmcc_fracLut, 1)      
cnrm_fracLut_ann = month_to_annual_weighted(yyyymm, cnrm_fracLut, 1)      
gfdl_fracLut_ann = month_to_annual_weighted(yyyymm, gfdl_fracLut, 1)      
mpi_fracLut_ann = month_to_annual_weighted(yyyymm, mpi_fracLut, 1)      

delete(cesm_fracLut)
delete(ipsl_fracLut)
delete(cmcc_fracLut)
delete(cnrm_fracLut)
delete(gfdl_fracLut)
delete(mpi_fracLut)

landuse=ispan(0,3,1)    ;用来给landuse赋值

;------------------- 计算cesm区域平均 -------------------
lat=tofloat(cesm_fracLut_ann&lat)
lon=tofloat(cesm_fracLut_ann&lon)
ny=dimsizes(lat)
area0 = new(ny,float)
R=6371.22
rad=4.0*atan(1.0)/180.0
dy=abs(lat(1)-lat(2))
dx=abs(lon(1)-lon(2))
do i=0,ny-1
      area0(i)=rad*(R^2)*abs(sin(rad*(lat(i)+dy/2))-sin(rad*(lat(i)-dy/2)))*dx
      area0(i)=area0(i)/1000000    ;*10^6 km2         
end do

cesm_fracLut_mean = wgt_areaave_Wrap(cesm_fracLut_ann, area0, 1.0, 1) 
printVarSummary(cesm_fracLut_mean)
cesm_fracLut_mean1 = transpose(cesm_fracLut_mean)     ;4*165

delete(cesm_fracLut_ann)
delete(cesm_fracLut_mean)

;------------------- 计算ipsl区域平均 -------------------
lat:=tofloat(ipsl_fracLut_ann&lat)
lon:=tofloat(ipsl_fracLut_ann&lon)
ny:=dimsizes(lat)
area0 := new(ny,float)
R=6371.22
rad:=4.0*atan(1.0)/180.0
dy:=abs(lat(1)-lat(2))
dx:=abs(lon(1)-lon(2))
do i=0,ny-1
      area0(i)=rad*(R^2)*abs(sin(rad*(lat(i)+dy/2))-sin(rad*(lat(i)-dy/2)))*dx
      area0(i)=area0(i)/1000000    ;*10^6 km2         
end do

ipsl_fracLut_mean = wgt_areaave_Wrap(ipsl_fracLut_ann, area0, 1.0, 1) 
ipsl_fracLut_mean&landuse = landuse      ;因为print出来没有landuse的坐标信息，赋值
printVarSummary(ipsl_fracLut_mean)
ipsl_fracLut_mean1 = transpose(ipsl_fracLut_mean)    ;4*165

delete(ipsl_fracLut_ann)
delete(area0)
delete(ipsl_fracLut_mean)
landuse=ispan(0,3,1)

;------------------- 计算cmcc区域平均 -------------------
lat:=tofloat(cmcc_fracLut_ann&lat)
lon:=tofloat(cmcc_fracLut_ann&lon)
ny:=dimsizes(lat)
area0 := new(ny,float)
R=6371.22
rad:=4.0*atan(1.0)/180.0
dy:=abs(lat(1)-lat(2))
dx:=abs(lon(1)-lon(2))
do i=0,ny-1
      area0(i)=rad*(R^2)*abs(sin(rad*(lat(i)+dy/2))-sin(rad*(lat(i)-dy/2)))*dx
      area0(i)=area0(i)/1000000    ;*10^6 km2         
end do

cmcc_fracLut_mean = wgt_areaave_Wrap(cmcc_fracLut_ann, area0, 1.0, 1) 
cmcc_fracLut_mean&landuse = landuse      ;因为print出来没有landuse的坐标信息，赋值
printVarSummary(cmcc_fracLut_mean)
cmcc_fracLut_mean1 = transpose(cmcc_fracLut_mean)    ;4*165

delete(cmcc_fracLut_ann)
delete(area0)
delete(cmcc_fracLut_mean)

;------------------- 计算cnrm区域平均 -------------------
lat:=tofloat(cnrm_fracLut_ann&lat)
lon:=tofloat(cnrm_fracLut_ann&lon)
ny:=dimsizes(lat)
area0 := new(ny,float)
R=6371.22
rad:=4.0*atan(1.0)/180.0
dy:=abs(lat(1)-lat(2))
dx:=abs(lon(1)-lon(2))
do i=0,ny-1
      area0(i)=rad*(R^2)*abs(sin(rad*(lat(i)+dy/2))-sin(rad*(lat(i)-dy/2)))*dx
      area0(i)=area0(i)/1000000    ;*10^6 km2         
end do

cnrm_fracLut_mean = wgt_areaave_Wrap(cnrm_fracLut_ann, area0, 1.0, 1) 
cnrm_fracLut_mean&landUse = landuse      ;因为print，landUse的坐标信息是1-4，不是0-3，重新赋值
printVarSummary(cnrm_fracLut_mean)
cnrm_fracLut_mean1 = transpose(cnrm_fracLut_mean)    ;4*165

delete(cnrm_fracLut_ann)
delete(area0)
delete(cnrm_fracLut_mean)

;------------------- 计算gfdl区域平均 -------------------
lat:=tofloat(gfdl_fracLut_ann&lat)
lon:=tofloat(gfdl_fracLut_ann&lon)
ny:=dimsizes(lat)
area0 := new(ny,float)
R=6371.22
rad:=4.0*atan(1.0)/180.0
dy:=abs(lat(1)-lat(2))
dx:=abs(lon(1)-lon(2))
do i=0,ny-1
      area0(i)=rad*(R^2)*abs(sin(rad*(lat(i)+dy/2))-sin(rad*(lat(i)-dy/2)))*dx
      area0(i)=area0(i)/1000000    ;*10^6 km2         
end do

gfdl_fracLut_mean = wgt_areaave_Wrap(gfdl_fracLut_ann, area0, 1.0, 1) 
printVarSummary(gfdl_fracLut_mean)
gfdl_fracLut_mean1 = transpose(gfdl_fracLut_mean)    ;4*165

delete(gfdl_fracLut_ann)
delete(area0)


;------------------- 计算mpi区域平均 -------------------
lat:=tofloat(mpi_fracLut_ann&lat)
lon:=tofloat(mpi_fracLut_ann&lon)
ny:=dimsizes(lat)
area0 := new(ny,float)
R=6371.22
rad:=4.0*atan(1.0)/180.0
dy:=abs(lat(1)-lat(2))
dx:=abs(lon(1)-lon(2))
do i=0,ny-1
      area0(i)=rad*(R^2)*abs(sin(rad*(lat(i)+dy/2))-sin(rad*(lat(i)-dy/2)))*dx
      area0(i)=area0(i)/1000000    ;*10^6 km2         
end do

mpi_fracLut_mean = wgt_areaave_Wrap(mpi_fracLut_ann, area0, 1.0, 1) 
copy_VarCoords(gfdl_fracLut_mean, mpi_fracLut_mean)
printVarSummary(mpi_fracLut_mean)
mpi_fracLut_mean1 = transpose(mpi_fracLut_mean)    ;4*165

delete(gfdl_fracLut_mean)
delete(mpi_fracLut_ann)
delete(area0)
delete(mpi_fracLut_mean)

;--------------------- 放到一个数组里 ------------------------
frac=new((/6,4,165/), float)
frac(0,:,:) = cesm_fracLut_mean1
frac(1,:,:) = cmcc_fracLut_mean1
frac(2,:,:) = cnrm_fracLut_mean1
frac(3,:,:) = gfdl_fracLut_mean1
frac(4,:,:) = ipsl_fracLut_mean1
frac(5,:,:) = mpi_fracLut_mean1

delete(cesm_fracLut_mean1)
delete(cmcc_fracLut_mean1)
delete(cnrm_fracLut_mean1)
delete(gfdl_fracLut_mean1)
delete(ipsl_fracLut_mean1)
delete(mpi_fracLut_mean1)

;----------------------------------------------------------
wks = gsn_open_wks("eps", "fracLUT_psl_1850-2014_all_v10")

 res                    = True                      ; plot mods desired
 res@gsnDraw            = False
 res@gsnFrame           = False
 res@tiMainString       = "psl:1850-2014"       ; add title
 res@tiMainFont         = 25       ; {25,    "times-roman"}



 ;------------------- x和y轴的范围 ------------------- 
  ;res@trYMinF                  = min(frac(:,0,:))    :大概是48
  res@trYMinF                  = 40     ;留点地方给label
  res@trYMaxF                  = 100        ;某些地表会有总和超过1的情况  
  res@trXMinF                  = 1850
  res@trXMaxF                  = 2015

  res@tmLabelAutoStride        = True     ;当此资源设置为真时，勾选标记标签将在绘制之前检查是否有重叠
  res@tmXBLabelFont = "times-roman"
  res@tmYLLabelFont = "times-roman"
  
  res@lbTitleFont = "times-roman"
  res@lbLabelFont = "times-roman"
 ;------------------- 线条设置 ------------------- 
  res@xyLineColors       = (/"DodgerBlue","tomato","dark slate grey","SpringGreen3","sienna3","DarkGoldenrod1"/)          ; change line color  ,
  res@xyLineThicknesses  = (/3.,3.,3.,3.,3.,3./)     ;改变线条粗细  ,3.,3.
  res@xyMonoDashPattern      = True     ; force all solid lines


 ;----------------- legend 设置 --------------------
  ;res@pmLegendDisplayMode    = "Always"              ; 打开legend
  ;res@pmLegendWidthF         = 0.12                  ; Change width and
  ;res@pmLegendHeightF        = 0.15                  ; height of legend.
  ;res@pmLegendOrthogonalPosF = -0.08                 ; move up slightly 
  ;res@lgLabelFontHeightF     = .011                  ; change font height
  ;res@lgPerimOn              = False                 ; no box around决定是否在图例的周长周围画线。

  ;res@xyExplicitLegendLabels = " " + (/"CESM2","CMCC-ESM2","CNRM-ESM2-1","GFDL-ESM4","IPSL-CM6A-LR","MPI-ESM1-2-LR"/)         ;

  plot  = gsn_csm_xy (wks,frac&year,frac(:,0,:),res) ; create plot

; Manually create and attach legend--来自xy18
;***************************************************** 
   res_text                    = True                  ; text mods desired
   res_text@txFontHeightF      = 0.015                 ; change text size
   res_text@txJust             = "CenterLeft"          ; 垂直对齐，水平对齐、
   
   res_lines                   = True                  ; polyline mods desired
   res_lines@gsLineDashPattern = 0.                    ; solid line
   res_lines@gsLineThicknessF  = 3.                    ; line thicker
   res_lines@gsLineLabelFont   = 125                   ; 字体
   ;res_lines@gsLineThicknessF  = 0.5

   res_lines@gsLineColor       = "DodgerBlue"                ; line color
   xx = (/1860,1875/)
   yy1 = (/58,58/)
   dum1 = gsn_add_polyline(wks,plot,xx,yy1,res_lines)              ; add polyline
   dum2 = gsn_add_text(wks,plot,"CESM2",1880,58,res_text)        ; add text
  
   yy2 = (/55,55/)
   res_lines@gsLineColor       = "tomato"                                 ; change to blue
   dum3 = gsn_add_polyline(wks,plot,xx,yy2,res_lines)                ; add polyline
   dum4 = gsn_add_text(wks,plot,"CMCC-ESM2",1880,55,res_text)       ; add text
  
   yy3 = (/52,52/)
   res_lines@gsLineColor       = "dark slate grey"                                ; change to black
   dum5 = gsn_add_polyline(wks,plot,xx,yy3,res_lines)                ; add polyline
   dum6 = gsn_add_text(wks,plot,"CNRM-ESM2-1",1880,52,res_text) ; add text
  
   yy4 = (/49,49/)
   res_lines@gsLineColor       = "SpringGreen3"                                 ; change to blue
   dum7 = gsn_add_polyline(wks,plot,xx,yy4,res_lines)                ; add polyline
   dum8 = gsn_add_text(wks,plot,"GFDL-ESM4",1880,49,res_text)       ; add text
  
   yy5 = (/46,46/)
   res_lines@gsLineColor       = "sienna3"                                 ; change to blue
   dum9 = gsn_add_polyline(wks,plot,xx,yy5,res_lines)                ; add polyline
   dum10 = gsn_add_text(wks,plot,"IPSL-CM6A-LR",1880,46,res_text)       ; add text
   
   yy6 = (/43,43/)
   res_lines@gsLineColor       = "DarkGoldenrod1"                                 ; change to blue
   dum11 = gsn_add_polyline(wks,plot,xx,yy6,res_lines)                ; add polyline
   dum12 = gsn_add_text(wks,plot,"MPI-ESM1-2-LR",1880,43,res_text)       ; add text


  draw(plot)
  frame(wks)
end 


