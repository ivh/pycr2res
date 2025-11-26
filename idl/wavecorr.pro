function wavecorr_ref,wl,sp,xref,wref,WINDOW=win,POWER=pow,COEFFICIENTS=coe $
                        ,PLOT=plot,MODEL=mode,FILTER=filter

  if keyword_set(win) then window=win else window=21
  if keyword_set(pow) then power=pow else power=1
  if keyword_set(mode) then model=strlowcase(mode) else model='gauss+lorentz'
  if keyword_set(filter) then filter_width=filter else filter_width=60.d0
  nlines=n_elements(xref)
  nx=n_elements(wl)

  if nlines le 1 then begin
    o=sp
    oo=middle(sp,filter_width)
    o[where(o gt oo-stddev(o-sp))]=max(o)

    o=max(o)-o
    i=where(o gt 0, n0)

    i1=[i[0],i[where(i[1:*]-i[0:*] gt 1,n1)+1]]
    n1=n1+1
    i2=[i[where(i[1:*]-i[0:*] gt 1,n2)],i[n0-1]]
    n2=n2+1

    ii=(i1+i2)/2

    i=where(i2-i1 ge 8 and i1 gt window/2 and $
            i2 lt nx-1-window/2 and sp[ii] lt 0.9*oo[ii], nn)

    i1=i1[i]
    i2=i2[i]

;    i1=ii-window/2
;    i2=ii+window/2

    ii=(i1+i2)/2

    nn=n_elements(ii)

    xref=dblarr(nn)
    wref=dblarr(nn)
    o=max(sp)-sp

    dev=dblarr(nn)
    for ilin=0,nn-1 do begin
      xx=indgen(i2[ilin]-i1[ilin]+1)+i1[ilin]
      case model of
      'gauss':    begin
                    dummy1=disp_gaussfit(   xx,o[xx],a,xx,/NOSLOP)
                    dummy=dummy1
                    xref[ilin]=a[1]
                    dev[ilin]=total((dummy-o[xx])^2)
                  end
      'lorentz':  begin
                    dummy2=disp_lorentzfit( xx,o[xx],b,xx,/NOSLOP)
                    dummy=dummy2
                    xref[ilin]=b[1]
                    dev[ilin]=total((dummy-o[xx])^2)
                  end
      'gaussbox': begin
                    dummy3=disp_gaussboxfit(xx,o[xx],c,xx,/NOSLOP)
                    dummy=dummy3
                    xref[ilin]=c[1]
                    dev[ilin]=total((dummy-o[xx])^2)
                  end
            else: begin
                    dummy1=disp_gaussfit(   xx,o[xx],a,xx,/NOSLOP)
                    dummy2=disp_lorentzfit(xx,o[xx],b,xx,/NOSLOP)
                    dummy3=disp_gaussboxfit(xx,o[xx],c,xx,/NOSLOP)
                    std1=total((dummy1-o[xx])^2)
                    std2=total((dummy2-o[xx])^2)
                    std3=total((dummy3-o[xx])^2)
                    std4=total(((dummy1+dummy2)*0.5-o[xx])^2)
                    std=[std1,std2,std3,std4]
                    if min(std) eq std1 then begin
                      dummy=dummy1
                      xref[ilin]=a[1]
                      dev[ilin]=std1
                    endif else if min(std) eq std2 then begin
                      dummy=dummy2
                      xref[ilin]=b[1]
                      dev[ilin]=std2
                    endif else if min(std) eq std3 then begin
                      dummy=dummy3
                      xref[ilin]=c[1]
                      dev[ilin]=std3
                    endif else if min(std) eq std4 then begin
                      dummy=(dummy1+dummy2)*0.5d0
                      xref[ilin]=(a[1]+b[1])*0.5d0
                      dev[ilin]=std4
                    endif
                  end
      endcase

      wref[ilin]=spline(xx,wl[xx],xref[ilin])
      if keyword_set(plot) then begin
        plot,wl[xx],sp[xx],xs=3,ys=3,tit='Line #'+strtrim(ilin,2)
        oplot,wl[xx],max(sp)-dummy,col=c24(2)
        oplot,wref[ilin]+[0,0],!y.crange,col=c24(4),line=2
        ssss=get_kbrd(1)
;stop
      endif
    endfor
    jgood=where(dev lt middle(dev,10.d0)+stddev(dev-middle(dev,10.d0)),nn)

    xref=xref[jgood]
    wref=wref[jgood]
    i1=i1[jgood]
    i2=i2[jgood]

    if keyword_set(plot) then begin
      plot,wl,sp,xs=1,ys=3
      for j=0,nn-1 do begin
        oplot,wl[i1[j]]+[0,0],!y.crange,line=2
        oplot,wl[i2[j]]+[0,0],!y.crange,line=3
        oplot,wref[j]+[0,0],!y.crange,col=c24(4)
      endfor
    endif

    return,wl
  endif

; Finding match to the selected reference lines in the current spectrum.


  o=sp
  oo=middle(sp,filter_width)
  o[where(o gt oo-stddev(o-sp))]=max(o)

  o=max(o)-o
  i=where(o gt 0, n0)

  o=max(sp)-sp

  nn=n_elements(xref)
  i1=intarr(nn)
  i2=intarr(nn)

  for jj=0,nn-1 do begin
    j=round(xref[jj])
    j1=j-window/2
    j2=j+window/2
    dummy=max(o[j1:j2],jmax)
    jmax=jmax+j1
    if abs(jmax-j) le window/4 then begin
      i1[jj]=j-window/2
      i2[jj]=j+window/2
    endif else begin
      print,'Reference line #'+strtrim(jj,2)+' not found in the current spectrum'
      stop
    endelse
  endfor

  n1=nn
  n2=nn

  ii=(i1+i2)/2

  nlines=nn

  o=max(sp)-sp
  x0=dblarr(nn)
  dev=x0
  w0=x0

  for ilin=0,nlines-1 do begin
    xx=indgen(i2[ilin]-i1[ilin]+1)+i1[ilin]
    case model of
    'gauss':    begin
                  dummy1=disp_gaussfit(   xx,o[xx],a,xx,/NOSLOP)
                  dummy=dummy1
                  x0[ilin]=a[1]
                end
    'lorentz':  begin
                  dummy2=disp_lorentzfit( xx,o[xx],b,xx,/NOSLOP)
                  dummy=dummy2
                  x0[ilin]=b[1]
                end
    'gaussbox': begin
                  dummy3=disp_gaussboxfit(xx,o[xx],c,xx,/NOSLOP)
                  dummy=dummy3
                  x0[ilin]=c[1]
                end
          else: begin
                  dummy1=disp_gaussfit(   xx,o[xx],a,xx,/NOSLOP)
                  dummy2=disp_lorentzfit(xx,o[xx],b,xx,/NOSLOP)
                  dummy3=disp_gaussboxfit(xx,o[xx],c,xx,/NOSLOP)
                  std1=total((dummy1-o[xx])^2)
                  std2=total((dummy2-o[xx])^2)
                  std3=total((dummy3-o[xx])^2)
                  std4=total(((dummy1+dummy2)*0.5-o[xx])^2)
                  std=[std1,std2,std3,std4]
                  if min(std) eq std1 then begin
                    dummy=dummy1
                    x0[ilin]=a[1]
                  endif else if min(std) eq std2 then begin
                    dummy=dummy2
                    x0[ilin]=b[1]
                  endif else if min(std) eq std3 then begin
                    dummy=dummy3
                    x0[ilin]=c[1]
                  endif else if min(std) eq std4 then begin
                    dummy=(dummy1+dummy2)*0.5d0
                    x0[ilin]=(a[1]+b[1])*0.5d0
                  endif
                end
    endcase

    w0[ilin]=spline(xx,wl[xx],x0[ilin])
    if keyword_set(plot) then begin
      plot,wl[xx],sp[xx],xs=3,tit='Line #'+strtrim(ilin,2)
      oplot,wl[xx],max(sp)-dummy,col=c24(2)
      oplot,w0[ilin]+[0,0],!y.crange,col=c24(2),line=2
      oplot,wref[ilin]+[0,0],!y.crange,col=c24(4),line=2
;      stop
;      ssss=get_kbrd(1)
    endif
  endfor

  if keyword_set(plot) then begin
    plot,wl,sp,xs=1;,xr=[2072,2074.5]
    for j=0,nlines-1 do begin
      oplot,wl[i1[j]]+[0,0],!y.crange,line=2
      oplot,wl[i2[j]]+[0,0],!y.crange,line=3
      oplot,w0  [j]+[0,0],!y.crange,col=c24(3),line=2
      oplot,wref[j]+[0,0],!y.crange,col=c24(4),line=2
    endfor
  endif

  wmean=mean(w0)
;  coe=poly_fit(w0-wmean,wref/w0,power,/DOUBLE)
;  wl_new=poly(wl-wmean,coe)*wl

  x=dindgen(nx)
  coe=poly_fit(x0,wref/w0,power,/DOUBLE)
  wl_new=poly(x,coe)*wl

  return,wl_new
end

pro wavecorr,wave1,obs1,unc1,wave2,obs2,unc2,wave3,obs3,unc3 $
            ,REF_ORDER=ref_ord,REF_PHASE=ref_phase,PLOT=plot

  nphase=n_elements(obs1[0,0,*])
  nord=n_elements(wave1[0,*])
  nwl1=n_elements(wave1[*,0])
  x1=dindgen(nwl1)
  nwl2=n_elements(wave2[*,0])
  x2=dindgen(nwl2)+nwl1
  nwl3=n_elements(wave3[*,0])
  x3=dindgen(nwl3)+nwl1+nwl2

  if keyword_set(ref_ord) then iord=ref_ord else iord=1
  if keyword_set(ref_phase) then ref_ph=ref_phase else ref_ph=1

  ww=wave1[*,iord]
  oo=obs1[*,iord,ref_ph]
  wl=wavecorr_ref(ww,oo,xref1,wref1)
  if keyword_set(plot) then begin
    plot,ww,oo,xs=1
    oplot,ww,middle(oo,30.d0),col=c24(2)
    for i=0,n_elements(xref1)-1 do oplot,wref1[i]+[0,0],!y.crange
    ssss=get_kbrd(1)
  endif

  ww=wave2[*,iord]
  oo=obs2[*,iord,ref_ph]
  wl=wavecorr_ref(ww,oo,xref2,wref2)
  if keyword_set(plot) then begin
    plot,ww,oo,xs=1
    oplot,ww,middle(oo,30.d0),col=c24(2)
    for i=0,n_elements(xref2)-1 do oplot,wref2[i]+[0,0],!y.crange
    ssss=get_kbrd(1)
  endif

  ww=wave3[*,iord]
  oo=obs3[*,iord,ref_ph]
  wl=wavecorr_ref(ww,oo,xref3,wref3)
  if keyword_set(plot) then begin
    plot,ww,oo,xs=1
    oplot,ww,middle(oo,30.d0),col=c24(2)
    for i=0,n_elements(xref3)-1 do oplot,wref3[i]+[0,0],!y.crange
    ssss=get_kbrd(1)
  endif

;---------------------------------
  xref=[xref1,xref2+nwl1,xref3+nwl1+nwl2]
  wref=[wref1,wref2,wref3]
  oo=[obs1[*,iord,ref_ph],obs2[*,iord,ref_ph],obs3[*,iord,ref_ph]]
  w1=wave1[*,iord]
  w2=wave2[*,iord]
  w3=wave3[*,iord]
  w =[w1,w2,w3]

  oo1=oo[        0:          nwl1-1] 
  oo2=oo[     nwl1:     nwl2+nwl1-1]
  oo3=oo[nwl2+nwl1:nwl3+nwl2+nwl1-1]

  for jphase=0,nphase-1 do begin
    o=[obs1[*,iord,jphase],obs2[*,iord,jphase],obs3[*,iord,jphase]]
    wl=wavecorr_ref(w,o,xref,wref,COE=coe,POWER=2)
    onew=spline(wl,o,w,/DOUBLE)

    if keyword_set(plot) then begin
      plot,w,oo,xs=1,ys=3,tit='Phase '+strtrim(jphase,2)
;      plot,w,oo,xs=1,xr=[2051,2054],ys=3,tit='Det 1, Phase '+strtrim(jphase,2),/NODATA
      plot,w,oo,xs=1,xr=[2068,2071],ys=3,tit='Det 2, Phase '+strtrim(jphase,2);,/NODATA
;      plot,w,oo,xs=1,xr=[2084,2087],ys=3,tit='Det 3, Phase '+strtrim(jphase,2),/NODATA

      o1=o[        0:          nwl1-1] & onew1=onew[        0:          nwl1-1]
      o2=o[     nwl1:     nwl2+nwl1-1] & onew2=onew[     nwl1:     nwl2+nwl1-1]
      o3=o[nwl2+nwl1:nwl3+nwl2+nwl1-1] & onew3=onew[nwl2+nwl1:nwl3+nwl2+nwl1-1]

      a1=total(o1*(oo1-mean(oo1)))/(total(o1*o1)-mean(o1)*total(o1))
      b1=mean(oo1)-a1*mean(o1)
      a2=total(o2*(oo2-mean(oo2)))/(total(o2*o2)-mean(o2)*total(o2))
      b2=mean(oo2)-a2*mean(o2)
      a3=total(o3*(oo3-mean(oo3)))/(total(o3*o3)-mean(o3)*total(o3))
      b3=mean(oo3)-a3*mean(o3)
      oplot,w,[oo1*a1+b1,oo2*a2+b2,oo3*a3+b3],col=c24(2)

      a1=total(onew1*(oo1-mean(oo1)))/(total(onew1*onew1)-mean(onew1)*total(onew1))
      b1=mean(oo1)-a1*mean(onew1)
      a2=total(onew2*(oo2-mean(oo2)))/(total(onew2*onew2)-mean(onew2)*total(onew2))
      b2=mean(oo2)-a2*mean(onew2)
      a3=total(onew3*(oo3-mean(oo3)))/(total(onew3*onew3)-mean(onew3)*total(onew3))
      b3=mean(oo3)-a3*mean(onew3)
      oonew=[onew1*a1+b1,onew2*a2+b2,onew3*a3+b3]
      oplot,w,oonew,col=c24(4)
      wait,0.02
;      ssss=get_kbrd(1)
    endif
    dev2=stddev(oo-oonew)

    brd1=0.d0
    brd2=0.d0
    brd3=2.d-4

    onew1=gaussbroad(w1,onew[        0:          nwl1-1],brd3*mean(w1)/mean(w2))
    onew2=gaussbroad(w2,onew[     nwl1:     nwl2+nwl1-1],brd3)
    onew3=gaussbroad(w3,onew[nwl2+nwl1:nwl3+nwl2+nwl1-1],brd3*mean(w3)/mean(w2))

    a1=total(onew1*(oo1-mean(oo1)))/(total(onew1*onew1)-mean(onew1)*total(onew1))
    b1=mean(oo1)-a1*mean(onew1)
    a2=total(onew2*(oo2-mean(oo2)))/(total(onew2*onew2)-mean(onew2)*total(onew2))
    b2=mean(oo2)-a2*mean(onew2)
    a3=total(onew3*(oo3-mean(oo3)))/(total(onew3*onew3)-mean(onew3)*total(onew3))
    b3=mean(oo3)-a3*mean(onew3)
    oonew=[onew1*a1+b1,onew2*a2+b2,onew3*a3+b3]
    dev3=stddev(oo-oonew)

    for jord=0,nord-1 do begin
      wl_new=poly(x1,coe)*wave1[*,jord]
      onew=spline(wl_new,obs1[*,jord,jphase],wave1[*,jord],/DOUBLE)
      obs1[*,jord,jphase]=onew
      unew=spline(wl_new,unc1[*,jord,jphase],wave1[*,jord],/DOUBLE)
      unc1[*,jord,jphase]=unew

      wl_new=poly(x2,coe)*wave2[*,jord]
      onew=spline(wl_new,obs2[*,jord,jphase],wave2[*,jord],/DOUBLE)
      obs2[*,jord,jphase]=onew
      unew=spline(wl_new,unc2[*,jord,jphase],wave2[*,jord],/DOUBLE)
      unc2[*,jord,jphase]=unew

      wl_new=poly(x3,coe)*wave3[*,jord]
      onew=spline(wl_new,obs3[*,jord,jphase],wave3[*,jord],/DOUBLE)
      obs3[*,jord,jphase]=onew
      unew=spline(wl_new,unc3[*,jord,jphase],wave3[*,jord],/DOUBLE)
      unc3[*,jord,jphase]=unew
    endfor
  endfor

end
