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

;------------------- 处理missingvalue -------------------
 ;print(num(ismissing(cesm_fracLut(0,0,:,:))))   --34283
if (any(isnan_ieee(cesm_fracLut))) then
      value = 1.e20
      replace_ieeenan (cesm_fracLut, value, 0)
      cesm_fracLut@_FillValue = value
end if
 ;print(num(ismissing(cesm_fracLut(0,0,:,:))))   --34283

;--------------------------------------------
if (any(isnan_ieee(ipsl_fracLut))) then
      value = 1.e20
      replace_ieeenan (ipsl_fracLut, value, 0)
      ipsl_fracLut@_FillValue = value
end if
;--------------------------------------------
if (any(isnan_ieee(cmcc_fracLut))) then
      value = 1.e20
      replace_ieeenan (cmcc_fracLut, value, 0)
      cmcc_fracLut@_FillValue = value
end if
;--------------------------------------------
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

;--------------------------统计结果,165年的最大，最小，均值，标准差，和最大最小值的和（试图看一下是不是等于1虽然知道这样算不太对）-----------------------------------
; stas = new((/6,4,5/), float)
; do i = 0, 5 ;不同模式
;   do j = 0, 3 ;不同土地类型
;     stas(i,j,0) = max(frac(i,j,:))
;     stas(i,j,1) = min(frac(i,j,:))
;     stas(i,j,2) = avg(frac(i,j,:))
;     stas(i,j,3) = stddev(frac(i,j,:))
;     stas(i,j,4) = max(frac(i,j,:))+min(frac(i,j,:))

;     print_table([/"mode" + i + j/], "%s")
;     print(stas(i,j,:))  
;   end do 
; end do

ver = "v2"
;***************** 0层，primary and secondary land ***********************
wks1 = gsn_open_wks("eps", "frac_pst_1850-2014_" + ver)

 res1                    = True                      ; plot mods desired
 res1@gsnDraw            = False
 res1@gsnFrame           = False
 res1@tiMainString       = "Primary and Secondary Land:1850-2014"       ; add title
 res1@tiMainFont         = 26       ; {26,    "times-bold"}  :https://www.ncl.ucar.edu/Document/HLUs/Classes/Base.shtml#NhlTFont

 ;------------------- x和y轴的范围 ------------------- 
  res1@trYMinF                  = min(frac(:,0,:)) - 2   
  ;res1@trYMinF                  = min(value)     ;留点地方给label
  res1@trYMaxF                  = max(frac(:,0,:)) + 2       ;某些地表会有总和超过1的情况  
  res1@trXMinF                  = 1850
  res1@trXMaxF                  = 2015

  res1@tmLabelAutoStride        = True     ;当此资源设置为真时，勾选标记标签将在绘制之前检查是否有重叠
  res1@tmXBLabelFont = 26
  res1@tmYLLabelFont = 26
  
  res1@lbTitleFont = 26
  res1@lbLabelFont = 26
 ;------------------- 线条设置 ------------------- 
  res1@xyLineColors       = (/"DodgerBlue","RosyBrown1","dark slate grey","SpringGreen3","sienna3","DarkGoldenrod1"/)          ; change line color  ,
  res1@xyLineThicknesses  = (/3.,3.,3.,3.,3.,3./)     ;改变线条粗细  ,3.,3.
  res1@xyMonoDashPattern      = True     ; force all solid lines

  plot1  = gsn_csm_xy (wks1,frac&year,frac(:,0,:),res1) ; create plot

;-----------------------Label,放到图的下面,参考黄平书p161-----------------------------
  ndcres1                  =   True
  ndcres1@txFontHeightF    =   0.015
  ndcres1@txFont           =   26    ;times-bold
  ndcres1@txJust           =   "CenterLeft"

  string1                  =   "____"
  ndcres1@txFontColor      =   "DodgerBlue" 
  gsn_text_ndc(wks1, string1, 0.14, 0.105, ndcres1) ;调整一下，往左一点点，往上一点点

  string2                  =   "CESM2"
  ndcres1@txFontColor      =   "black" 
  gsn_text_ndc(wks1, string2, 0.20, 0.10, ndcres1)

  string3                  =   "____"
  ndcres1@txFontColor      =   "RosyBrown1" 
  gsn_text_ndc(wks1, string3, 0.39, 0.105, ndcres1)

  string4                  =   "CMCC-ESM2"
  ndcres1@txFontColor      =   "black" 
  gsn_text_ndc(wks1, string4, 0.45, 0.10, ndcres1)

  string5                  =   "____"
  ndcres1@txFontColor      =   "dark slate grey"
  gsn_text_ndc(wks1, string5, 0.64, 0.105, ndcres1)

  string6                  =   "CNRM-ESM2-1"
  ndcres1@txFontColor      =   "black"
  gsn_text_ndc(wks1, string6, 0.70, 0.10, ndcres1)

  string7                  =   "____"
  ndcres1@txFontColor      =   "SpringGreen3"
  gsn_text_ndc(wks1, string7, 0.14, 0.055, ndcres1)

  string8                  =   "GFDL-ESM4"
  ndcres1@txFontColor      =   "black"
  gsn_text_ndc(wks1, string8, 0.20, 0.05, ndcres1)

  string9                  =   "____"
  ndcres1@txFontColor      =   "sienna3"
  gsn_text_ndc(wks1, string9, 0.39, 0.055, ndcres1)

  string10                 =   "IPSL-CM6A-LR"
  ndcres1@txFontColor      =   "black"
  gsn_text_ndc(wks1, string10, 0.45, 0.05, ndcres1)

  string11                 =   "____"
  ndcres1@txFontColor      =   "DarkGoldenrod1"
  gsn_text_ndc(wks1, string11, 0.64, 0.055, ndcres1)

  string12                 =   "MPI-ESM1-2-LR"
  ndcres1@txFontColor      =   "black"
  gsn_text_ndc(wks1, string12, 0.70, 0.05, ndcres1)

  draw(plot1)
  frame(wks1)

;********************* 1 层，crop land *******************************
wks2 = gsn_open_wks("eps", "frac_crp_1850-2014_" + ver)

 res2                    = True                      ; plot mods desired
 res2@gsnDraw            = False
 res2@gsnFrame           = False
 res2@tiMainString       = "Cropland:1850-2014"       ; add title
 res2@tiMainFont         = 26       ; {26,    "times-bold"}  :https://www.ncl.ucar.edu/Document/HLUs/Classes/Base.shtml#NhlTFont

 ;------------------- x和y轴的范围 ------------------- 
  res2@trYMinF                  = min(frac(:,1,:)) - 2   
  ;res2@trYMinF                  = min(value)     ;留点地方给label
  res2@trYMaxF                  = max(frac(:,1,:)) + 2      ;某些地表会有总和超过1的情况  
  res2@trXMinF                  = 1850
  res2@trXMaxF                  = 2015

  res2@tmLabelAutoStride        = True     ;当此资源设置为真时，勾选标记标签将在绘制之前检查是否有重叠
  res2@tmXBLabelFont = 26
  res2@tmYLLabelFont = 26
  
  res2@lbTitleFont = 26
  res2@lbLabelFont = 26
 ;------------------- 线条设置 ------------------- 
  res2@xyLineColors       = (/"DodgerBlue","RosyBrown1","dark slate grey","SpringGreen3","sienna3","DarkGoldenrod1"/)          ; change line color  ,
  res2@xyLineThicknesses  = (/3.,3.,3.,3.,3.,3./)     ;改变线条粗细  ,3.,3.
  res2@xyMonoDashPattern      = True     ; force all solid lines

  plot2  = gsn_csm_xy (wks2,frac&year,frac(:,1,:),res2) ; create plot

;-----------------------Label,放到图的下面,参考黄平书p161-----------------------------
  ndcres2                  =   True
  ndcres2@txFontHeightF    =   0.015
  ndcres2@txFont           =   26    ;times-bold
  ndcres2@txJust           =   "CenterLeft"

  string13                 =   "____"
  ndcres2@txFontColor      =   "DodgerBlue" 
  gsn_text_ndc(wks2, string13, 0.14, 0.105, ndcres2) ;调整一下，往左一点点，往上一点点

  string14                =   "CESM2"
  ndcres2@txFontColor      =   "black" 
  gsn_text_ndc(wks2, string14, 0.20, 0.10, ndcres2)

  string15                 =   "____"
  ndcres2@txFontColor      =   "RosyBrown1" 
  gsn_text_ndc(wks2, string15, 0.39, 0.105, ndcres2)

  string16                =   "CMCC-ESM2"
  ndcres2@txFontColor      =   "black" 
  gsn_text_ndc(wks2, string16, 0.45, 0.10, ndcres2)

  string17                 =   "____"
  ndcres2@txFontColor      =   "dark slate grey"
  gsn_text_ndc(wks2, string17, 0.64, 0.105, ndcres2)

  string18                =   "CNRM-ESM2-1"
  ndcres2@txFontColor      =   "black"
  gsn_text_ndc(wks2, string18, 0.70, 0.10, ndcres2)

  string19                 =   "____"
  ndcres2@txFontColor      =   "SpringGreen3"
  gsn_text_ndc(wks2, string19, 0.14, 0.055, ndcres2)

  string20                =   "GFDL-ESM4"
  ndcres2@txFontColor      =   "black"
  gsn_text_ndc(wks2, string20, 0.20, 0.05, ndcres2)

  string21                 =   "____"
  ndcres2@txFontColor      =   "sienna3"
  gsn_text_ndc(wks2, string21, 0.39, 0.055, ndcres2)

  string22                =   "IPSL-CM6A-LR"
  ndcres2@txFontColor      =   "black"
  gsn_text_ndc(wks2, string22, 0.45, 0.05, ndcres2)

  string23                 =   "____"
  ndcres2@txFontColor      =   "DarkGoldenrod1"
  gsn_text_ndc(wks2, string23, 0.64, 0.055, ndcres2)

  string24                =   "MPI-ESM1-2-LR"
  ndcres2@txFontColor      =   "black"
  gsn_text_ndc(wks2, string24, 0.70, 0.05, ndcres2)

;------------------------------------------------------------
  draw(plot2)
  frame(wks2)

;********************* 2 层，pastureland *******************************
wks3 = gsn_open_wks("eps", "frac_pas_1850-2014_" + ver)

 res3                    = True                      ; plot mods desired
 res3@gsnDraw            = False
 res3@gsnFrame           = False
 res3@tiMainString       = "Pastureland:1850-2014"       ; add title
 res3@tiMainFont         = 26       ; {26,    "times-bold"}  :https://www.ncl.ucar.edu/Document/HLUs/Classes/Base.shtml#NhlTFont

 ;------------------- x和y轴的范围 ------------------- 
  res3@trYMinF                  = min(frac(:,2,:)) - 2    
  ;res3@trYMinF                  = min(value)     ;留点地方给label
  res3@trYMaxF                  = max(frac(:,2,:)) + 2        ;某些地表会有总和超过1的情况  
  res3@trXMinF                  = 1850
  res3@trXMaxF                  = 2015

  res3@tmLabelAutoStride        = True     ;当此资源设置为真时，勾选标记标签将在绘制之前检查是否有重叠
  res3@tmXBLabelFont = 26
  res3@tmYLLabelFont = 26
  
  res3@lbTitleFont = 26
  res3@lbLabelFont = 26
 ;------------------- 线条设置 ------------------- 
  res3@xyLineColors       = (/"DodgerBlue","RosyBrown1","dark slate grey","SpringGreen3","sienna3","DarkGoldenrod1"/)          ; change line color  ,
  res3@xyLineThicknesses  = (/3.,3.,3.,3.,3.,3./)     ;改变线条粗细  ,3.,3.
  res3@xyMonoDashPattern      = True     ; force all solid lines

  plot3  = gsn_csm_xy (wks3,frac&year,frac(:,2,:),res3) ; create plot

;-----------------------Label,放到图的下面,参考黄平书p161-----------------------------
  ndcres3                  =   True
  ndcres3@txFontHeightF    =   0.015
  ndcres3@txFont           =   26    ;times-bold
  ndcres3@txJust           =   "CenterLeft"

  string25                 =   "____"
  ndcres3@txFontColor      =   "DodgerBlue" 
  gsn_text_ndc(wks3, string25, 0.14, 0.105, ndcres3) ;调整一下，往左一点点，往上一点点

  string26                 =   "CESM2"
  ndcres3@txFontColor      =   "black" 
  gsn_text_ndc(wks3, string26, 0.20, 0.10, ndcres3)

  string27                 =   "____"
  ndcres3@txFontColor      =   "RosyBrown1" 
  gsn_text_ndc(wks3, string27, 0.39, 0.105, ndcres3)

  string28                 =   "CMCC-ESM2"
  ndcres3@txFontColor      =   "black" 
  gsn_text_ndc(wks3, string28, 0.45, 0.10, ndcres3)

  string29                 =   "____"
  ndcres3@txFontColor      =   "dark slate grey"
  gsn_text_ndc(wks3, string29, 0.64, 0.105, ndcres3)

  string30                 =   "CNRM-ESM2-1"
  ndcres3@txFontColor      =   "black"
  gsn_text_ndc(wks3, string30, 0.70, 0.10, ndcres3)

  string31                 =   "____"
  ndcres3@txFontColor      =   "SpringGreen3"
  gsn_text_ndc(wks3, string31, 0.14, 0.055, ndcres3)

  string32                 =   "GFDL-ESM4"
  ndcres3@txFontColor      =   "black"
  gsn_text_ndc(wks3, string32, 0.20, 0.05, ndcres3)

  string33                 =   "____"
  ndcres3@txFontColor      =   "sienna3"
  gsn_text_ndc(wks3, string33, 0.39, 0.055, ndcres3)

  string34                 =   "IPSL-CM6A-LR"
  ndcres3@txFontColor      =   "black"
  gsn_text_ndc(wks3, string34, 0.45, 0.05, ndcres3)

  string35                 =   "____"
  ndcres3@txFontColor      =   "DarkGoldenrod1"
  gsn_text_ndc(wks3, string35, 0.64, 0.055, ndcres3)

  string36                 =   "MPI-ESM1-2-LR"
  ndcres3@txFontColor      =   "black"
  gsn_text_ndc(wks3, string36, 0.70, 0.05, ndcres3)

;------------------------------------------------------------
  draw(plot3)
  frame(wks3)

;********************* 3 层，urban *******************************
wks4 = gsn_open_wks("eps", "frac_urb_1850-2014_" + ver)

 res4                    = True                      ; plot mods desired
 res4@gsnDraw            = False
 res4@gsnFrame           = False
 res4@tiMainString       = "Pastureland:1850-2014"       ; add title
 res4@tiMainFont         = 26       ; {26,    "times-bold"}  :https://www.ncl.ucar.edu/Document/HLUs/Classes/Base.shtml#NhlTFont

 ;------------------- x和y轴的范围 ------------------- 
  res4@trYMinF                  = min(frac(:,3,:)) - 2    
  ;res4@trYMinF                  = min(value)     ;留点地方给label
  res4@trYMaxF                  = max(frac(:,3,:)) + 2         ;某些地表会有总和超过1的情况  
  res4@trXMinF                  = 1850
  res4@trXMaxF                  = 2015

  res4@tmLabelAutoStride        = True     ;当此资源设置为真时，勾选标记标签将在绘制之前检查是否有重叠
  res4@tmXBLabelFont = 26
  res4@tmYLLabelFont = 26
  
  res4@lbTitleFont = 26
  res4@lbLabelFont = 26
 ;------------------- 线条设置 ------------------- 
  res4@xyLineColors       = (/"DodgerBlue","RosyBrown1","dark slate grey","SpringGreen3","sienna3","DarkGoldenrod1"/)          ; change line color  ,
  res4@xyLineThicknesses  = (/3.,3.,3.,3.,3.,3./)     ;改变线条粗细  ,3.,3.
  res4@xyMonoDashPattern      = True     ; force all solid lines

  plot4  = gsn_csm_xy (wks4,frac&year,frac(:,3,:),res4) ; create plot

;-----------------------Label,放到图的下面,参考黄平书p161-----------------------------
  ndcres4                  =   True
  ndcres4@txFontHeightF    =   0.015
  ndcres4@txFont           =   26    ;times-bold
  ndcres4@txJust           =   "CenterLeft"

  string37                 =   "____"
  ndcres4@txFontColor      =   "DodgerBlue" 
  gsn_text_ndc(wks4, string37, 0.14, 0.105, ndcres4) ;调整一下，往左一点点，往上一点点

  string38                 =   "CESM2"
  ndcres4@txFontColor      =   "black" 
  gsn_text_ndc(wks4, string38, 0.20, 0.10, ndcres4)

  string39                 =   "____"
  ndcres4@txFontColor      =   "RosyBrown1" 
  gsn_text_ndc(wks4, string39, 0.39, 0.105, ndcres4)

  string40                 =   "CMCC-ESM2"
  ndcres4@txFontColor      =   "black" 
  gsn_text_ndc(wks4, string40, 0.45, 0.10, ndcres4)

  string41                 =   "____"
  ndcres4@txFontColor      =   "dark slate grey"
  gsn_text_ndc(wks4, string41, 0.64, 0.105, ndcres4)

  string42                 =   "CNRM-ESM2-1"
  ndcres4@txFontColor      =   "black"
  gsn_text_ndc(wks4, string42, 0.70, 0.10, ndcres4)

  string43                 =   "____"
  ndcres4@txFontColor      =   "SpringGreen3"
  gsn_text_ndc(wks4, string43, 0.14, 0.055, ndcres4)

  string44                 =   "GFDL-ESM4"
  ndcres4@txFontColor      =   "black"
  gsn_text_ndc(wks4, string44, 0.20, 0.05, ndcres4)

  string45                 =   "____"
  ndcres4@txFontColor      =   "sienna3"
  gsn_text_ndc(wks4, string45, 0.39, 0.055, ndcres4)

  string46                 =   "IPSL-CM6A-LR"
  ndcres4@txFontColor      =   "black"
  gsn_text_ndc(wks4, string46, 0.45, 0.05, ndcres4)

  string47                 =   "____"
  ndcres4@txFontColor      =   "DarkGoldenrod1"
  gsn_text_ndc(wks4, string47, 0.64, 0.055, ndcres4)

  string48                 =   "MPI-ESM1-2-LR"
  ndcres4@txFontColor      =   "black"
  gsn_text_ndc(wks4, string48, 0.70, 0.05, ndcres4)

;------------------------------------------------------------
  draw(plot4)
  frame(wks4)

end 


