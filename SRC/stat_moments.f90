! Module to compute to compute statistics over user specified area
! domains (i.e. goal is to compute moments for each 1 km x 1 km area for the
! projected super run).  The computation of these moments/statistics can
! easily be turned off by changing a parameter in the setparm.f90 file
! Original implementation by Peter Bogenschutz, U. Uta, March 2008

module stat_moments

use vars
use params
use rad, only: qrad

implicit none

integer, parameter :: nxm = max(1,(nx_gl/nsubdomains_x)/navgmom_x)
integer, parameter :: nym = max(1,(ny_gl/nsubdomains_y)/navgmom_y)

real moment1(nzm,5,nxm,nym), moment2(nzm,5,nxm,nym), moment3(nzm,9,nxm,nym)
real mom1cld(nzm,8,nxm,nym), moment5(nzm,8,nxm,nym), moment4(nzm,3,nxm,nym)

! added Gentine
real qrad_avg(nxm,nym,nzm)
real wt_avg(nxm,nym),wq_avg(nxm,nym),wu_avg(nxm,nym),wv_avg(nxm,nym)
real ddx_ut(nxm,nym,nzm),ddx_uq(nxm,nym,nzm),ddx_uqn(nxm,nym,nzm)
real ddy_vt(nxm,nym,nzm),ddy_vq(nxm,nym,nzm),ddy_vqn(nxm,nym,nzm)
real ddz_wt(nxm,nym,nzm),ddz_wq(nxm,nym,nzm),ddz_wqn(nxm,nym,nzm)
real ddx_uu(nxm,nym,nzm),ddx_uv(nxm,nym,nzm),ddx_uw(nxm,nym,nzm)
real ddy_vu(nxm,nym,nzm),ddy_vv(nxm,nym,nzm),ddy_vw(nxm,nym,nzm)
real ddz_wu(nxm,nym,nzm),ddz_wv(nxm,nym,nzm),ddz_ww(nxm,nym,nzm)
real udtdx(nxm,nym,nzm),udqdx(nxm,nym,nzm),udqndx(nxm,nym,nzm)
real vdtdy(nxm,nym,nzm),vdqdy(nxm,nym,nzm),vdqndy(nxm,nym,nzm)
real wdtdz(nxm,nym,nzm),wdqdz(nxm,nym,nzm),wdqndz(nxm,nym,nzm)
real ududx(nxm,nym,nzm),udvdx(nxm,nym,nzm),udwdx(nxm,nym,nzm)
real vdudy(nxm,nym,nzm),vdvdy(nxm,nym,nzm),vdwdy(nxm,nym,nzm)
real wdudz(nxm,nym,nzm),wdvdz(nxm,nym,nzm),wdwdz(nxm,nym,nzm)
real qcl_avg(nxm,nym,nzm),qci_avg(nxm,nym,nzm),qpl_avg(nxm,nym,nzm)
real qpi_avg(nxm,nym,nzm),qsatw_avg(nxm,nym,nzm),Mc_avg(nxm,nym,nzm)
real hydro3D(nxm,nym,nzm),mcup(nxm,nym,nzm),mcdns(nxm,nym,nzm)
real mcdnu(nxm,nym,nzm),mc(nxm,nym,nzm),cldd_avg(nxm,nym,nzm)
real mcrdns(nxm,nym,nzm),mcrdnu(nxm,nym,nzm),mcr(nxm,nym,nzm)
real mcrup(nxm,nym,nzm),condavg_mask_avg(nxm,nym,nzm,5)
real wdse_avg(nxm,nym,nzm),wmse_avg(nxm,nym,nzm),wfmse_avg(nxm,nym,nzm)
real wfmsecla_avg(nxm,nym,nzm),dse_avg(nxm,nym,nzm),mse_avg(nxm,nym,nzm)
real fmse_avg(nxm,nym,nzm),fmsecla_avg(nxm,nym,nzm),qv_avg(nxm,nym,nzm)
real wqv_avg(nxm,nym,nzm),wqcl_avg(nxm,nym,nzm),qt_avg(nxm,nym,nzm)
real wqci_avg(nxm,nym,nzm),wqpl_avg(nxm,nym,nzm),wqpi_avg(nxm,nym,nzm)
real prec_xy_avg(nxm,nym),shf_xy_avg(nxm,nym)
real lhf_xy_avg(nxm,nym),tau_xy_avg(nxm,nym)


CONTAINS

subroutine compmoments()

! SUBROUTINE TO COMPUTE STATISTICAL PROFILES 
! Profiles computed for given domain area (controlled in domain.f90 file)
! Output can be used to construct moments and covariances for selected 
! variables.  

logical condition_cl, condition
integer i, j, k, ii,jj, starti, startj,nn
real(8) thlfac, qwfac, qtfac, wfac, divfac,coef
real(8) wfac2, wfac3, wfac4
real(8) thetal(nzm,4), totalqw(nzm,4), vertw(nzm,12), qclavg(nzm,5), cldcnt(nzm)
real(8)  uavg(nzm,3), vavg(nzm,3), micro(nzm,6)
real auto, accre, rflux, evapa, evapb
real(8) temp, temp2, Dz, tvz_k, tvirt
real(8) dsefac,msefac,ssefac,fmsefac,fmseclafac

! added Gentine

call t_startf ('moments')

do jj=1,nym
do ii=1,nxm
wt_avg(i,j)= 0.
wq_avg(i,j)= 0.
wu_avg(i,j)= 0.
wv_avg(i,j)= 0.
prec_xy_avg(i,j)= 0.
shf_xy_avg(i,j)= 0.
lhf_xy_avg(i,j)= 0.
tau_xy_avg(i,j)= 0.
do k=1,nzm
qrad_avg(i,j,k)= 0.
ddx_ut(i,j,k)= 0.
ddx_uq(i,j,k)= 0.
ddx_uqn(i,j,k)= 0.
ddy_vt(i,j,k)= 0.
ddy_vq(i,j,k)= 0.
ddy_vqn(i,j,k)= 0.
ddz_wt(i,j,k)= 0.
ddz_wq(i,j,k)= 0.
ddz_wqn(i,j,k)= 0.
ddx_uu(i,j,k)= 0.
ddx_uv(i,j,k)= 0.
ddx_uw(i,j,k)= 0.
ddy_vu(i,j,k)= 0.
ddy_vv(i,j,k)= 0.
ddy_vw(i,j,k)= 0.
ddz_wu(i,j,k)= 0.
ddz_wv(i,j,k)= 0.
ddz_ww(i,j,k)= 0.
udtdx(i,j,k)= 0.
udqdx(i,j,k)= 0.
udqndx(i,j,k)= 0.
vdtdy(i,j,k)= 0.
vdqdy(i,j,k)= 0.
vdqndy(i,j,k)= 0.
wdtdz(i,j,k)= 0.
wdqdz(i,j,k)= 0.
wdqndz(i,j,k)= 0.
ududx(i,j,k)= 0.
udvdx(i,j,k)= 0.
udwdx(i,j,k)= 0.
vdudy(i,j,k)= 0.
vdvdy(i,j,k)= 0.
vdwdy(i,j,k)= 0.
wdudz(i,j,k)= 0.
wdvdz(i,j,k)= 0.
wdwdz(i,j,k)= 0.
qcl_avg(i,j,k)= 0.
qci_avg(i,j,k)= 0.
qpl_avg(i,j,k)= 0.
qpi_avg(i,j,k)= 0.
qsatw_avg(i,j,k)= 0.
Mc_avg(i,j,k)= 0.
hydro3D(i,j,k)= 0.
mcup(i,j,k)= 0.
mcdns(i,j,k)= 0.
mcdnu(i,j,k)= 0.
mc(i,j,k)= 0.
cldd_avg(i,j,k)= 0.
mcrdns(i,j,k)= 0.
mcrdnu(i,j,k)= 0.
mcr(i,j,k)= 0.
mcrup(i,j,k)= 0.
wdse_avg(i,j,k)= 0.
wmse_avg(i,j,k)= 0.
wfmse_avg(i,j,k)= 0.
wfmsecla_avg(i,j,k)= 0.
dse_avg(i,j,k)= 0.
mse_avg(i,j,k)= 0.
fmse_avg(i,j,k)= 0.
fmsecla_avg(i,j,k)= 0.
qv_avg(i,j,k)= 0.
wqv_avg(i,j,k)= 0.
wqcl_avg(i,j,k)= 0.
qt_avg(i,j,k)= 0.
wqci_avg(i,j,k)= 0.
wqpl_avg(i,j,k)= 0.
wqpi_avg(i,j,k)= 0.
do nn=1,5
condavg_mask_avg(i,j,k,nn)= 0.
enddo
enddo
enddo
enddo


startj=1
do jj=1,nym
  starti=1.
  do ii=1,nxm

    thetal(:,:)        = 0.
    totalqw(:,:)       = 0.
    vertw(:,:)         = 0.
    qclavg(:,:)        = 0.
    cldcnt(:)          = 0.
    uavg(:,:)          = 0.
    vavg(:,:)          = 0.
    micro(:,:)         = 0.
    wt_avg(:,:)        = 0.
    wq_avg(:,:)        = 0.
    wu_avg(:,:)        = 0.
    wv_avg(:,:)        = 0.
    prec_xy_avg(:,:)   = 0.
    shf_xy_avg(:,:)    = 0.
    lhf_xy_avg(:,:)    = 0.
    tau_xy_avg(:,:)    = 0.

    divfac=float(navgmom_x)*float(navgmom_y)
    !print*,'divfac = ',divfac 

    do k=1,nzm
      if(LES) then
          coef=0.
      else
          coef=min(1.e-5,0.01*qsatw(tabs0(k),pres(k)))
      endif
      tvz_k = 0. 
      do j=startj,startj+(navgmom_y-1)
        do i=starti,starti+(navgmom_x-1)
          tvz_k=tvz_k+tvirt/(divfac)
        enddo
      enddo
      
      do j=startj,startj+(navgmom_y-1)
        do i=starti,starti+(navgmom_x-1)


        if(k .eq. 1) then
            ! added Gentine
            wt_avg(ii,jj)      = wt_avg(ii,jj) + fluxbt(i,j)/(divfac - float(navgmom_y))
            wq_avg(ii,jj)      = wq_avg(ii,jj) + fluxbq(i,j)/(divfac)
            wu_avg(ii,jj)      = wu_avg(ii,jj) + fluxbu(i,j)/(divfac)
            wv_avg(ii,jj)      = wv_avg(ii,jj) + fluxbv(i,j)/(divfac)
            prec_xy_avg(ii,jj) = prec_xy_avg(ii,jj)+ prec_xy(i,j)/(divfac)
            shf_xy_avg(ii,jj)  = shf_xy_avg(ii,jj) + shf_xy(i,j)/(divfac)
            lhf_xy_avg(ii,jj)  = lhf_xy_avg(ii,jj) + lhf_xy(i,j)/(divfac)
            tau_xy_avg(ii,jj)  = tau_xy_avg(ii,jj) + tau_xy(i,j)/(divfac)
        endif
        ! Liquid water potential temperature
          thlfac=(tabs(i,j,k)*prespot(k))*(1-fac_cond*(qcl(i,j,k)+qci(i,j,k))/tabs(i,j,k))
        !  water vapor mixing ratio
          qwfac=qv(i,j,k)
          ! Total water mixing ratio
          qtfac=qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)
        !if(i.eq.2 .and. j.eq.2 .and. k.eq. 1) then
        !        print*,'qv(i,j,k)=',qv(i,j,k)
        !        print*,'qcl(i,j,k)=',qcl(i,j,k)
        !        print*,'qci(i,j,k)=',qci(i,j,k)
        !        print*,'qwfac = ',qwfac
        !        print*,'qtfac = ',qtfac
        !end if
        ! Vertical Velocity
          wfac=(w(i,j,k)+w(i,j,k+1))/2.
        ! Virtual temperature
          tvirt=tabs(i,j,k)*prespot(k)*(1.+epsv*qv(i,j,k)-(qcl(i,j,k)+qci(i,j,k))-(qpl(i,j,k)+qpi(i,j,k)))   
           
        ! added Pierre Gentine
        dsefac=tabs(i,j,k)+gamaz(k)
        msefac=dsefac+fac_cond*qv(i,j,k)
        ssefac=tabs(i,j,k)+gamaz(k)+fac_cond*qsatw(tabs(i,j,k),pres(k))
        fmsefac=t(i,j,k)+fac_cond*(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)+qpl(i,j,k)+qpi(i,j,k))
        fmseclafac=t(i,j,k)-t0(k) &
            +fac_cond*(qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)+qpl(i,j,k)+qpi(i,j,k)-q0(k)-qp0(k))
        qv_avg(ii,jj,k)   = qv_avg(ii,jj,k)  + qwfac / (divfac)
        qt_avg(ii,jj,k)   = qt_avg(ii,jj,k)  + qtfac / (divfac)
        qcl_avg(ii,jj,k)  = qcl_avg(ii,jj,k)  + qcl(i,j,k) / (divfac)
        qci_avg(ii,jj,k)  = qci_avg(ii,jj,k)  + qci(i,j,k) / (divfac)
        qpl_avg(ii,jj,k)  = qpl_avg(ii,jj,k)  + qpl(i,j,k) / (divfac)
        qpi_avg(ii,jj,k)  = qpi_avg(ii,jj,k)  + qpi(i,j,k) / (divfac)
        qsatw_avg(ii,jj,k)= qsatw_avg(ii,jj,k)+ qsatw(tabs(i,j,k),pres(k)) / (divfac)
          !if(i.eq.2 .and. j.eq.2 .and. k.eq. 1) then
          !      print*,'qpl = ',qpl(i,j,k)
          !      print*,'qpi = ',qpi(i,j,k)
          !      print*,'pres(k) = ',pres(k)
          !endif
        
        !-------------------------------------------------------------
        !	Mass flux, hydrometeor fraction statistics
        !-------------------------------------------------------------
        if(qcl(i,j,k)+qci(i,j,k).gt.coef) then
            hydro3D(ii,jj,k) = hydro3D(ii,jj,k) + 1 / (divfac)
        if(wfac.gt.0.) then
            mcup(ii,jj,k)=mcup(ii,jj,k)+rho(k)*wfac / (divfac)
        else
            mcdns(ii,jj,k)=mcdns(ii,jj,k)+rho(k)*wfac / (divfac)
        endif
        elseif(qpl(i,j,k)+qpi(i,j,k).gt.1.e-4) then
            hydro3D(ii,jj,k) = hydro3D(ii,jj,k) + 1. / (divfac)
        if(w(i,j,k)+w(i,j,k+1).lt.0.) &
            mcdnu(ii,jj,k)=mcdnu(ii,jj,k)+rho(k)*wfac / (divfac)
        endif

        mc(ii,jj,k)=mcup(ii,jj,k)+mcdns(ii,jj,k)+mcdnu(ii,jj,k)

        !-------------------------------------------------------------
        !	Updraft Core statistics:
        !-------------------------------------------------------------
        condition_cl = qcl(i,j,k)+qci(i,j,k).gt.coef
        condition = tvirt.gt.tvz_k
        if(CEM) condition=condition.and.w(i,j,k)+w(i,j,k+1).gt.2.
        if(LES) condition=condition_cl.and.condition
        if(condition) then
            mcrup(ii,jj,k)=mcrup(ii,jj,k)+rho(k)*wfac / (divfac)
            if(condition_cl) then
               cldd_avg(ii,jj,k)=cldd_avg(ii,jj,k)+1 / (divfac)
            end if
        endif

        !-------------------------------------------------------------
        !	Cloud Downdraft Core statistics:
        !-------------------------------------------------------------
        condition_cl = qcl(i,j,k)+qci(i,j,k).gt.coef .or. qpl(i,j,k)+qpi(i,j,k).gt.1.e-4
        condition = tvirt.lt.tvz_k
        if(CEM) condition=condition.and.w(i,j,k)+w(i,j,k+1).lt.-2.
        if(LES) condition=condition_cl.and.condition
        if(condition) then
        if(condition_cl) then
            mcrdns(ii,jj,k)=mcrdns(ii,jj,k)+rho(k)*wfac / (divfac)
            cldd_avg(ii,jj,k)=cldd_avg(ii,jj,k) + 1 / (divfac)
        else
            mcrdnu(ii,jj,k)=mcrdnu(ii,jj,k)+rho(k)*wfac / (divfac)
        end if
        endif
        mcr(ii,jj,k)=mcrup(ii,jj,k)+mcrdns(ii,jj,k)+mcrdnu(ii,jj,k)
	Mc_avg(ii,jj,k)   = Mc_avg(ii,jj,k)   + rho(k)*wfac / (divfac)

        wdse_avg(ii,jj,k) = wdse_avg(ii,jj,k) + wfac*dsefac / (divfac)
        wmse_avg(ii,jj,k) = wmse_avg(ii,jj,k) + wfac*msefac / (divfac)
        wfmse_avg(ii,jj,k)= wfmse_avg(ii,jj,k)+ wfac*fmsefac / (divfac)
        wfmsecla_avg(ii,jj,k)= wfmsecla_avg(ii,jj,k) + wfac*fmseclafac / (divfac)
        dse_avg(ii,jj,k)  = dse_avg(ii,jj,k) + dsefac / (divfac)
        mse_avg(ii,jj,k)  = mse_avg(ii,jj,k) + msefac / (divfac)
        fmse_avg(ii,jj,k) = fmse_avg(ii,jj,k) + fmsefac / (divfac)
        fmsecla_avg(ii,jj,k)= fmsecla_avg(ii,jj,k) + fmseclafac / (divfac)
        wqv_avg(ii,jj,k)  = wqv_avg(ii,jj,k)  + wfac*qv(i,j,k) / (divfac)
        wqcl_avg(ii,jj,k) = wqcl_avg(ii,jj,k) + wfac*qcl(i,j,k) / (divfac)
        wqci_avg(ii,jj,k) = wqci_avg(ii,jj,k) + wfac*qci(i,j,k) / (divfac)
        wqpl_avg(ii,jj,k) = wqpl_avg(ii,jj,k) + wfac*qpl(i,j,k) / (divfac)
        wqpi_avg(ii,jj,k) = wqpi_avg(ii,jj,k) + wfac*qpi(i,j,k) / (divfac)

         ! advection and convergence, Pierre Gentine
        if (i .le. (nx-1)) then
            temp = tabs(i+1,j,k)+gamaz(k)
            temp2 = (temp*u(i+1,j,k) - dsefac*u(i,j,k))/dx
            ddx_ut(ii,jj,k) = ddx_ut(ii,jj,k) + temp2 / (divfac - float(navgmom_y))
            temp2 = (u(i,j,k)+u(i+1,j,k))/2 * (temp - dsefac)/dx
            udtdx(ii,jj,k)  = udtdx(ii,jj,k)  + temp2 / (divfac - float(navgmom_y))

            !print*,'udT/dx = ',temp2
          
            temp = qv(i+1,j,k)
            temp2 = (temp*u(i+1,j,k) - qwfac*u(i,j,k))/dx
            ddx_uq(ii,jj,k) = ddx_uq(ii,jj,k) + temp2 / (divfac - float(navgmom_y))
            temp2 = (u(i,j,k)+u(i+1,j,k))/2 * (temp - qwfac)/dx
            udqdx(ii,jj,k)  = udqdx(ii,jj,k)  + temp2 / (divfac - float(navgmom_y))
            if(qv(i+1,j,k) .lt. 0 .or. qv(i,j,k) .lt. 0) then
                print*,'qv(i+1,j,k) = ',qv(i+1,j,k)
                print*,'qv(i,j,k) = ',qv(i,j,k)
            endif

            !print*,'udq/dx = ',temp2
                
            temp = qcl(i+1,j,k)+qci(i+1,j,k)
            temp2 = (temp*u(i+1,j,k) - (qcl(i,j,k)+qci(i,j,k))*u(i,j,k))/dx
            ddx_uqn(ii,jj,k) = ddx_uqn(ii,jj,k) + temp2 / (divfac - float(navgmom_y))
            temp2 = (u(i,j,k)+u(i+1,j,k))/2 * (temp - (qcl(i,j,k)+qci(i,j,k)))/dx
            udqndx(ii,jj,k)  = udqndx(ii,jj,k)  + temp2 / (divfac - float(navgmom_y))

            temp = u(i+1,j,k)
            temp2 = (temp*u(i+1,j,k) - (u(i,j,k))*u(i,j,k))/dx
            ddx_uu(ii,jj,k) = ddx_uu(ii,jj,k) + temp2 / (divfac - float(navgmom_y))
            temp2 = (u(i,j,k)+u(i+1,j,k))/2 * (temp - u(i,j,k))/dx
            ududx(ii,jj,k)  = ududx(ii,jj,k)  + temp2 / (divfac - float(navgmom_y))

            temp = v(i+1,j,k)
            temp2 = (temp*u(i+1,j,k) - v(i,j,k)*u(i,j,k))/dx
            ddx_uu(ii,jj,k) = ddx_uu(ii,jj,k) + temp2 / (divfac - float(navgmom_y))
            temp2 = (u(i,j,k)+u(i+1,j,k))/2 * (temp - v(i,j,k))/dx
            udvdx(ii,jj,k)  = udvdx(ii,jj,k)  + temp2 / (divfac - float(navgmom_y))

            temp = w(i+1,j,k)
            temp2 = (temp*u(i+1,j,k) - w(i,j,k)*u(i,j,k))/dx
            ddx_uu(ii,jj,k) = ddx_uu(ii,jj,k) + temp2 / (divfac - float(navgmom_y))
            temp2 = (u(i,j,k)+u(i+1,j,k))/2 * (temp - w(i,j,k))/dx
            udwdx(ii,jj,k)  = udwdx(ii,jj,k)  + temp2 / (divfac - float(navgmom_y))
        end if

        if (j .le. (ny-1)) then
            temp = tabs(i,j+1,k)+gamaz(k)
            temp2 = (temp*v(i,j+1,k) - dsefac*v(i,j,k))/dy
            ddy_vt(ii,jj,k) = ddy_vt(ii,jj,k) + temp2 / (divfac - float(navgmom_x))
            temp2 = (v(i,j,k)+v(i,j+1,k))/2 * (temp - dsefac)/dy
            vdtdy(ii,jj,k)  = vdtdy(ii,jj,k)  + temp2 / (divfac - float(navgmom_x))

            temp = qv(i,j+1,k)
            temp2 = (temp*v(i,j+1,k) - qwfac*v(i,j,k))/dy
            ddy_vq(ii,jj,k) = ddy_vq(ii,jj,k) + temp2 / (divfac - float(navgmom_x))
            temp2 = (v(i,j,k)+v(i,j+1,k))/2 * (temp - qwfac)/dy
            vdqdy(ii,jj,k)  = vdqdy(ii,jj,k)  + temp2 / (divfac - float(navgmom_x))

            temp = qcl(i,j+1,k)+qci(i,j+1,k)
            temp2 = (temp*v(i,j+1,k) - (qcl(i,j,k)+qci(i,j,k))*v(i,j,k))/dy
            ddy_vqn(ii,jj,k) = ddy_vqn(ii,jj,k) + temp2 / (divfac - float(navgmom_x))
            temp2 = (v(i,j,k)+v(i,j+1,k))/2 * (temp - (qcl(i,j,k)+qci(i,j,k)))/dy
            vdqndy(ii,jj,k)  = vdqndy(ii,jj,k)  + temp2 / (divfac - float(navgmom_x))

            temp = u(i,j+1,k)
            temp2 = (temp*v(i,j+1,k) - u(i,j,k)*v(i,j,k))/dy
            ddy_vu(ii,jj,k) = ddy_vu(ii,jj,k) + temp2 / (divfac - float(navgmom_x))
            temp2 = (v(i,j,k)+v(i,j+1,k))/2 * (temp - u(i,j,k))/dy
            vdudy(ii,jj,k)  = vdudy(ii,jj,k)  + temp2 / (divfac - float(navgmom_x))

            temp = v(i,j+1,k)
            temp2 = (temp*v(i,j+1,k) - v(i,j,k)*v(i,j,k))/dy
            ddy_vv(ii,jj,k) = ddy_vv(ii,jj,k) + temp2 / (divfac - float(navgmom_x))
            temp2 = (v(i,j,k)+v(i,j+1,k))/2 * (temp - v(i,j,k))/dy
            vdvdy(ii,jj,k)  = vdvdy(ii,jj,k)  + temp2 / (divfac - float(navgmom_x))

            temp = w(i,j+1,k)
            temp2 = (temp*v(i,j+1,k) - w(i,j,k)*v(i,j,k))/dy
            ddy_vw(ii,jj,k) = ddy_vw(ii,jj,k) + temp2 / (divfac - float(navgmom_x))
            temp2 = (v(i,j,k)+v(i,j+1,k))/2 * (temp - w(i,j,k))/dy
            vdwdy(ii,jj,k)  = vdwdy(ii,jj,k)  + temp2 / (divfac - float(navgmom_x))
        end if

        if (k .le. (nzm-1)) then
            Dz   = z(k+1)-z(k)

            temp = tabs(i,j,k+1)+gamaz(k+1)
            temp2 = (temp*w(i,j,k+1) - dsefac*w(i,j,k))/Dz
            ddz_wt(ii,jj,k) = ddz_wt(ii,jj,k) + temp2 / (divfac)
            temp2 = (w(i,j,k)+w(i,j,k+1))/2 * (temp - dsefac)/Dz
            wdtdz(ii,jj,k)  = wdtdz(ii,jj,k)  + temp2 / (divfac)

            temp = qv(i,j,k+1)
            temp2 = (temp*w(i,j,k+1) - qwfac*w(i,j,k))/Dz
            ddz_wq(ii,jj,k) = ddz_wq(ii,jj,k) + temp2 / (divfac)
            temp2 = (w(i,j,k)+w(i,j,k+1))/2 * (temp - qwfac)/Dz
            wdqdz(ii,jj,k)  = wdqdz(ii,jj,k)  + temp2 / (divfac)

            temp = qcl(i,j,k+1)+qci(i,j,k+1)
            temp2 = (temp*w(i,j,k+1) - (qcl(i,j,k)+qci(i,j,k))*w(i,j,k))/Dz
            ddz_wqn(ii,jj,k) = ddz_wqn(ii,jj,k) + temp2 / (divfac)
            temp2 = wfac * (temp - (qcl(i,j,k)+qci(i,j,k)))/Dz
            wdqndz(ii,jj,k)  = wdqndz(ii,jj,k)  + temp2 / (divfac)

            temp = u(i,j,k+1)
            temp2 = (temp*w(i,j,k+1) - u(i,j,k)*w(i,j,k))/Dz
            ddz_wu(ii,jj,k) = ddz_wq(ii,jj,k) + temp2 / (divfac)
            temp2 = (w(i,j,k)+w(i,j,k+1))/2 * (temp - u(i,j,k))/Dz
            wdudz(ii,jj,k)  = wdudz(ii,jj,k)  + temp2 / (divfac)

            temp = v(i,j,k+1)
            temp2 = (temp*w(i,j,k+1) - v(i,j,k+1)*w(i,j,k))/Dz
            ddz_wv(ii,jj,k) = ddz_wv(ii,jj,k) + temp2 / (divfac)
            temp2 = (w(i,j,k)+w(i,j,k+1))/2 * (temp - v(i,j,k+1))/Dz
            wdvdz(ii,jj,k)  = wdvdz(ii,jj,k)  + temp2 / (divfac)

            temp = w(i,j,k+1)
            temp2 = (temp*w(i,j,k+1) - w(i,j,k)*w(i,j,k))/Dz
            ddz_ww(ii,jj,k) = ddz_ww(ii,jj,k) + temp2 / (divfac)
            temp2 = (w(i,j,k)+w(i,j,k+1))/2 * (temp - w(i,j,k))/Dz
            wdwdz(ii,jj,k)  = wdwdz(ii,jj,k)  + temp2 / (divfac)

        end if

	call autoconversion(qcl(i,j,k),auto)
        call accretion(qcl(i,j,k),qpl(i,j,k),accre)
        call evaporation(qpl(i,j,k),tabs(i,j,k),pres(k),qv(i,j,k),evapa,evapb)
        call rainflux(rho(k),qpl(i,j,k),rflux)

          qclavg(k,1)=qclavg(k,1)+qcl(i,j,k)	! First mom. cloud water
          uavg(k,1)=uavg(k,1)+u(i,j,k)		! First mom. u-wind
          vavg(k,1)=vavg(k,1)+v(i,j,k)		! First mom. v-wind

          if (qcl(i,j,k)+qci(i,j,k) .gt. 0) then
            cldcnt(k)=cldcnt(k)+1.
          endif  

          thetal(k,1)=thetal(k,1)+thlfac	! First mom. theta_l
          totalqw(k,1)=totalqw(k,1)+qwfac	! First mom. total water
          vertw(k,1)=vertw(k,1)+wfac 		! First mom. vert. vel.

          thetal(k,2)=thetal(k,2)+thlfac**2.	! Sum squares theta_l
          totalqw(k,2)=totalqw(k,2)+qwfac**2.	! Sum squares total water
          vertw(k,2)=vertw(k,2)+wfac**2.	! Sum squares vert. vel.
          uavg(k,2)=uavg(k,2)+u(i,j,k)**2.	! Sum squares u-wind
          vavg(k,2)=vavg(k,2)+v(i,j,k)**2.	! Sum squares v-wind

          thetal(k,3)=thetal(k,3)+(thlfac*wfac)	! Flux theta_l and vert. vel.	
          totalqw(k,3)=totalqw(k,3)+(thlfac*qwfac) ! Flux theta_l and tot wat
          vertw(k,3)=vertw(k,3)+(wfac*qwfac)	! Flux tot wat and vert vel
          uavg(k,3)=uavg(k,3)+(wfac*u(i,j,k))	! Flux vert vel and u-wind
          vavg(k,3)=vavg(k,3)+(wfac*v(i,j,k))   ! Flux vert vel and v-wind

          qclavg(k,2)=qclavg(k,2)+(qcl(i,j,k)*wfac)	! Flux qc and vert vel
          qclavg(k,3)=qclavg(k,3)+(wfac*wfac*qcl(i,j,k)) ! Flux w^2 and qc
          qclavg(k,4)=qclavg(k,4)+(thlfac*qcl(i,j,k))	! Flux theta_l and qc
          qclavg(k,5)=qclavg(k,5)+(qwfac*qcl(i,j,k))	! flux qw and qc

          thetal(k,4)=thetal(k,4)+thlfac**3.	! Sum cubes theta_l
          totalqw(k,4)=totalqw(k,4)+qwfac**3.	! Sum cubes tot water
          vertw(k,4)=vertw(k,4)+wfac**3.	! Sum cubes vert vel

          vertw(k,5)=vertw(k,5)+wfac**4.	! Sum quad vert vel
          
          vertw(k,6)=vertw(k,6)+(wfac*wfac*thlfac) ! Flux w^2 and theta_l
          vertw(k,7)=vertw(k,7)+(wfac*wfac*qwfac) ! Flux w^2 and tot wat
          vertw(k,8)=vertw(k,8)+(wfac*thlfac*thlfac) ! Flux theta_l^2 and w
          vertw(k,9)=vertw(k,9)+(wfac*qwfac*qwfac) ! Flux qw^2 and vert vel
          vertw(k,10)=vertw(k,10)+(wfac*qwfac*thlfac) ! Flux w, theta_l, qw
          vertw(k,11)=vertw(k,11)+(wfac*u(i,j,k)*u(i,j,k)) ! Flux u^2, w
          vertw(k,12)=vertw(k,12)+(wfac*v(i,j,k)*v(i,j,k)) ! Flux v^2, w

          micro(k,1)=micro(k,1)+qpl(i,j,k)  ! Precipitation (liquid)
          micro(k,2)=micro(k,2)+auto
          micro(k,3)=micro(k,3)+accre
          micro(k,4)=micro(k,4)+rflux
          micro(k,5)=micro(k,5)+evapa
          micro(k,6)=micro(k,6)+evapb




        enddo
      enddo
    enddo

    starti=starti+navgmom_x

  enddo
  !print*,'qv_avg = ',qv_avg(ii,jj,k)
          
  startj=startj+navgmom_y
enddo

! Call function to output the moments
call write_moments()

call t_stopf ('moments')

return

end subroutine compmoments

!----------------------------------------------------------------------
subroutine autoconversion(qc1,auto)

real qc1, auto
real, parameter :: alpha = 0.001  ! autoconversion rate
real, parameter :: q_co =  0.001

auto=max(0.,alpha*(qc1-q_co))

return

end subroutine autoconversion

!----------------------------------------------------------------------
subroutine accretion(qc1,qr1,accre)

real qc1, qr1, accre
real, parameter :: br = 0.8

accre = qc1*qr1**((3.+br)/4.)

return

end subroutine accretion

!----------------------------------------------------------------------
subroutine evaporation(qr1,tabs1,pres1,qv1,evapa,evapb)

real qr1, tabs1, pres1, qv1, evapa, evapb
real S
real, parameter :: br = 0.8

S = qv1/qsatw(tabs1,pres1)

evapa=qr1**(1./2.)*(S-1.)
evapb=qr1**((5.+br)/8.)*(S-1.)

return

end subroutine evaporation

!----------------------------------------------------------------------
subroutine rainflux(rho1,qr1,rflux)

real rho1, qr1, rflux
real, parameter :: br = 0.8

rflux = (rho1*qr1)**((1.+br)/4.)

return

end subroutine rainflux


subroutine write_moments()

implicit none

character *120 filename
character *80 long_name
character *8 name
character *10 timechar
character *4 rankchar
character *5 sepchar
character *6 filetype
character *10 units
character *10 myString
character *12 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
integer i,j,k,nfields,nfields1,nn
real(4) tmp(nxm,nym,nzm)

nfields=120

if(masterproc.or.output_sep) then

  if(output_sep) then
     write(rankchar,'(i4)') rank
     sepchar="_"//rankchar(5-lenstr(rankchar):4)
  else
     sepchar=""
  end if
  write(rankchar,'(i4)') nsubdomains
  write(timechar,'(i10)') nstep
  do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
  end do

  if(RUN3D) then
    if(savemombin) then
      filetype = '.bin3D'
    else
      filetype = '.com3D'
    end if
    filename='./OUT_MOMENTS/'//trim(case)//'_moments_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
    open(46,file=filename,status='unknown',form='unformatted')

  else
    if(savemombin) then
     if(savemomsep) then
       filetype = '.bin3D'
     else
       filetype = '.bin2D'
     end if
    else
     if(savemomsep) then
       filetype = '.com3D'
     else
       filetype = '.com2D'
     end if
    end if
    if(savemomsep) then
      filename='./OUT_MOMENTS/'//trim(case)//'_moments_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
      open(46,file=filename,status='unknown',form='unformatted')	
    else
      filename='./OUT_MOMENTS/'//trim(case)//'_moments_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//filetype//sepchar
      if(nrestart.eq.0.and.notopenedmom) then
         open(46,file=filename,status='unknown',form='unformatted')	
      else
         open(46,file=filename,status='unknown', &
                              form='unformatted', position='append')
      end if
      notopenedmom=.false.
    end if  

  end if

  if(masterproc) then

   if(savemombin) then

     write(46) nxm,nym,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields
     print*,nxm,nym,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields
     do k=1,nzm
       write(46) z(k) 
     end do
     do k=1,nzm
       write(46) pres(k)
     end do
     write(46) dx*float(navgmom_x)
     write(46) dy*float(navgmom_x)
     write(46) nstep*dt/(3600.*24.)+day0

   else

     write(long_name,'(8i4)') nxm,nym,nzm,nsubdomains, &
                                    nsubdomains_x,nsubdomains_y,nfields
     do k=1,nzm
        write(c_z(k),'(f12.3)') z(k)
     end do
     do k=1,nzm
        write(c_p(k),'(f12.3)') pres(k)
     end do
     write(c_dx,'(f12.5)') dx*float(navgmom_x)
     write(c_dy,'(f12.5)') dy*float(navgmom_y)
     write(c_time,'(f12.5)') nstep*dt/(3600.*24.)+day0
	
     write(46) long_name(1:32)
     write(46) c_time,c_dx,c_dy, (c_z(k),k=1,nzm),(c_p(k),k=1,nzm)

    end if ! savemombin

   end if ! msterproc 

end if ! masterproc.or.output_sep

nfields1=0




nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=qrad_avg(i,j,k)
        enddo
    enddo
enddo
name='QRAD'
long_name='Radiative heating rate'
units='K/day'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=ddx_ut(i,j,k)
        enddo
    enddo
enddo
name='D_UDSE_DX'
long_name='x-horizontal DSE divergence'
units='K m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=ddy_vt(i,j,k)
        enddo
    enddo
enddo
name='D_VDSE_DY'
long_name='y-horizontal DSE divergence'
units='K m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=ddz_wt(i,j,k)
        enddo
    enddo
enddo
name='D_WDSE_DX'
long_name='w DSE divergence'
units='K m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=ddx_uq(i,j,k)
        enddo
    enddo
enddo
name='D_UQ_DX'
long_name='x-horizontal q divergence'
units='kg/kg m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=ddy_vq(i,j,k)
        enddo
    enddo
enddo
name='D_VQ_DY'
long_name='y-horizontal q divergence'
units='kg/kg m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=ddz_wq(i,j,k)
        enddo
    enddo
enddo
name='D_WQ_Dz'
long_name='w q divergence'
units='kg/kg m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=ddx_uqn(i,j,k)
        enddo
    enddo
enddo
name='D_UQN_DX'
long_name='x-horizontal qn divergence'
units='kg/kg m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=ddy_vqn(i,j,k)
        enddo
    enddo
enddo
name='D_VQN_DY'
long_name='y-horizontal qn divergence'
units='kg/kg m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=ddz_wqn(i,j,k)
        enddo
    enddo
enddo
name='D_WQN_Dz'
long_name='w qn divergence'
units='kg/kg m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

!10

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=udtdx(i,j,k)
        enddo
    enddo
enddo
name='UD_DSE_DX'
long_name='x-horizontal DSE advection'
units='K m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=vdtdy(i,j,k)
        enddo
    enddo
enddo
name='VD_DSE_DY'
long_name='y-horizontal DSE advection'
units='K m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=wdtdz(i,j,k)
        enddo
    enddo
enddo
name='WD_DSE_DX'
long_name='w DSE advection'
units='K m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=udqdx(i,j,k)
        enddo
    enddo
enddo
name='UDQ_DX'
long_name='x-horizontal q advection'
units='kg/kg m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=vdqdy(i,j,k)
        enddo
    enddo
enddo
name='VDQ_DY'
long_name='y-horizontal q advection'
units='kg/kg m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=wdqdz(i,j,k)
        enddo
    enddo
enddo
name='WDQ_Dz'
long_name='w q advection'
units='kg/kg m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=udqndx(i,j,k)
        enddo
    enddo
enddo
name='UD_QN_DX'
long_name='x-horizontal qn advection'
units='kg/kg m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=vdqndy(i,j,k)
        enddo
    enddo
enddo
name='VD_QN_DY'
long_name='y-horizontal qn advection'
units='kg/kg m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=wdqndz(i,j,k)
        enddo
    enddo
enddo
name='WD_QN_Dz'
long_name='w qn advection'
units='kg/kg m / s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=qcl_avg(i,j,k)
        enddo
    enddo
enddo
name='QCL'
long_name='liquid condensate mixing ratio'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=qci_avg(i,j,k)
        enddo
    enddo
enddo
name='QCI'
long_name='ice condensate mixing ratio'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=qpl_avg(i,j,k)
        enddo
    enddo
enddo
name='QPL'
long_name='precipitate liquid  mixing ratio'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=qpl_avg(i,j,k)
        enddo
    enddo
enddo
name='QPLBIS'
long_name='precipitate liquid  mixing ratio'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=qpi_avg(i,j,k)
        enddo
    enddo
enddo
name='QPI'
long_name='precipitate ice  mixing ratio'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=qsatw_avg(i,j,k)
        enddo
    enddo
enddo
name='QSAT'
long_name='Saturation mixing ratio'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=Mc_avg(i,j,k)
        enddo
    enddo
enddo
name='Mc'
long_name='Mass flux'
units='kg m-2 s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=condavg_mask_avg(i,j,k,1)
        enddo
    enddo
enddo
name='FAVG1'
long_name='Conditional averaging 1 - fraction'
units=''
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=condavg_mask_avg(i,j,k,2)
        enddo
    enddo
enddo
name='FAVG2'
long_name='Conditional averaging 2 - fraction'
units=''
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=condavg_mask_avg(i,j,k,3)
        enddo
    enddo
enddo
name='FAVG3'
long_name='Conditional averaging 3 - fraction'
units=''
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=condavg_mask_avg(i,j,k,4)
        enddo
    enddo
enddo
name='FAVG4'
long_name='Conditional averaging 4 - fraction'
units=''
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)
!30
nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=condavg_mask_avg(i,j,k,5)
        enddo
    enddo
enddo
name='FAVG5'
long_name='Conditional averaging 5 - fraction'
units=''
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)



! end Gentine

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment1(k,2,i,j)
    enddo
  enddo
enddo
name='THL1'
long_name='First Mom Liquid Wat Pot Tmp'
units='K'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment2(k,2,i,j)
    enddo
  enddo
enddo
name='THL2'
long_name='Second Mom Liquid Wat Pot Tmp'
units='K^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment3(k,2,i,j)
    enddo
  enddo
enddo
name='THLW'
long_name='Flux THL and W'
units='K m/s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment4(k,2,i,j)
    enddo
  enddo
enddo
name='THL3'
long_name='Third Mom Liquid Wat Pot Tmp'
units='K^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment1(k,3,i,j)
    enddo
  enddo
enddo
name='QW1'
long_name='First Mom Total Water'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment2(k,3,i,j)
    enddo
  enddo
enddo
name='QW2'
long_name='Second Mom Total Water'
units='kg/kg ^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment3(k,3,i,j)
    enddo
  enddo
enddo
name='THLQW'
long_name='Flux THL and QW'
units='kg/kg K'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment4(k,3,i,j)
    enddo
  enddo
enddo
name='QW3'
long_name='Third Mom QW'
units='kg/kg ^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment1(k,1,i,j)
    enddo
  enddo
enddo
name='W1'
long_name='First Mom Vert Vel'
units='m/s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)
!40
nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment2(k,1,i,j)
    enddo
  enddo
enddo
name='W2'
long_name='Second Mom Vert Vel'
units='m/s ^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment3(k,1,i,j)
    enddo
  enddo
enddo
name='WQW'
long_name='Flux W and QW'
units='kg/kg m/s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1 
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment4(k,1,i,j)
    enddo
  enddo
enddo
name='W3'
long_name='Third Mom Vert Vel'
units='m/s ^3'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment1(k,4,i,j)
    enddo
  enddo
enddo
name='U1'
long_name='First Mom U-Wind'
units='m/s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment2(k,4,i,j)
    enddo
  enddo
enddo
name='U2'
long_name='Second Mom U-Wind'
units='m/s ^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment3(k,4,i,j)
    enddo
  enddo
enddo
name='UW'
long_name='Flux U and W'
units='m/s ^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment1(k,5,i,j)
    enddo
  enddo
enddo
name='V1'
long_name='First Mom V-Wind'
units='m/s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment2(k,5,i,j)
    enddo
  enddo
enddo
name='V2'
long_name='Second Mom V-Wind'
units='m/s ^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment3(k,5,i,j)
    enddo
  enddo
enddo
name='VW'
long_name='Flux V and W'
units='m/s ^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1 
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=mom1cld(k,1,i,j)
    enddo
  enddo
enddo
name='QL'
long_name='QL AVG'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains) 

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=mom1cld(k,2,i,j)
    enddo
  enddo
enddo
name='CDFRC'
long_name='Cloud Fraction'
units='-'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment3(k,6,i,j)
    enddo
  enddo
enddo
name='WQL'
long_name='Flux W and QL'
units='m/s kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment3(k,7,i,j)
    enddo
  enddo
enddo
name='W2QL'
long_name='Flux W^2 and QL'
units='m/s ^2 kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment3(k,8,i,j)
    enddo
  enddo
enddo
name='THLQL'
long_name='Flux THL and QL'
units='K kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment3(k,9,i,j)
    enddo
  enddo
enddo
name='QTQL'
long_name='Flux QT and QL'
units='kg/kg ^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment5(k,1,i,j)
    enddo
  enddo
enddo
name='W4'
long_name='Fourth Moment W'
units='m/s ^4'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment5(k,2,i,j)
    enddo
  enddo
enddo
name='W2THL'
long_name='Flux W2 and THL'
units='m/s ^2 K'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment5(k,3,i,j)
    enddo
  enddo
enddo
name='W2QT'
long_name='Flux W2 and QT'
units='m/s ^2 kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment5(k,4,i,j)
    enddo
  enddo
enddo
name='WTHL2'
long_name='Flux W and THL2'
units='m/s K^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment5(k,5,i,j)
    enddo
  enddo
enddo
name='WQT2'
long_name='Flux W and QT2'
units='m/s kg/kg ^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)
!60
nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment5(k,6,i,j)
    enddo
  enddo
enddo
name='WQTTHL'
long_name='Flux W and QT and THL'
units='m/s kg/kg K'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment5(k,7,i,j)
    enddo
  enddo
enddo
name='WU2'
long_name='Flux W and U2'
units='m/s ^3'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

!write(*,*) 'here I am', nfields1

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=moment5(k,8,i,j)
    enddo
  enddo
enddo
name='WV2'
long_name='Flux W and V2'
units='m/s ^3'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=mom1cld(k,3,i,j)
    enddo
  enddo
enddo
name='QR1'
long_name='Rain Water Mixing Ratio'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=mom1cld(k,4,i,j)
    enddo
  enddo
enddo
name='AUTO1'
long_name='Autoconversion '
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=mom1cld(k,5,i,j)
    enddo
  enddo
enddo
name='ACRE1'
long_name='Accretion '
units='kg/kg ^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=mom1cld(k,6,i,j)
    enddo
  enddo
enddo
name='RFLUX'
long_name='Rain Flux'
units='kg m-3'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=mom1cld(k,7,i,j)
    enddo
  enddo
enddo
name='EVAPa'
long_name='Rain Evap Prof 1'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=mom1cld(k,8,i,j)
    enddo
  enddo
enddo
name='EVAPb'
long_name='Rain Evap Prof 2'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

! added Gentine
! mass fluxes

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=hydro3D(i,j,k)
        enddo
    enddo
enddo
name='HYDRO'
long_name='hydro fraction'
units='-'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=mcup(i,j,k)
        enddo
    enddo
enddo
name='MCUP'
long_name='Updraft mass flux'
units='kg m-2 s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=mcdns(i,j,k)
        enddo
    enddo
enddo
name='MCDNS'
long_name='Downdraft saturated mass flux'
units='kg m-2 s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=mcdnu(i,j,k)
        enddo
    enddo
enddo
name='MCDNU'
long_name='Downdraft unsaturated mass flux'
units='kg m-2 s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=cldd_avg(i,j,k)
        enddo
    enddo
enddo
name='CLDD'
long_name='CLD fraction'
units='-'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=mcrup(i,j,k)
        enddo
    enddo
enddo
name='MCRUP'
long_name='Core Updraft mass flux'
units='kg m-2 s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=mcrdns(i,j,k)
        enddo
    enddo
enddo
name='MCRDNS'
long_name='Core Downdraft saturated mass flux'
units='kg m-2 s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=mcrdnu(i,j,k)
        enddo
    enddo
enddo
name='MCRDNU'
long_name='Core Downdraft unsaturated mass flux'
units='kg m-2 s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=mcr(i,j,k)
        enddo
    enddo
enddo
name='MCR'
long_name='Core mass flux'
units='kg m-2 s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
    do j=1,nym
        do i=1,nxm
            tmp(i,j,k)=mc(i,j,k)
        enddo
    enddo
enddo
name='MC_avg'
long_name='Mass flux'
units='kg m-2 s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
   do j=1,nym
      do i=1,nxm
         tmp(i,j,k)=wdse_avg(i,j,k)
      enddo
   enddo
enddo
name='WDSE'
long_name='W Dry Static Energy'
units='K m s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)
!80
nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=wmse_avg(i,j,k)
enddo
enddo
enddo
name='WMSE'
long_name='W Moist Static Energy'
units='K m s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=wfmse_avg(i,j,k)
enddo
enddo
enddo
name='WFMSE'
long_name='W Ice Frozen Static Energy '
units='K m s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=wfmsecla_avg(i,j,k)
enddo
enddo
enddo
name='WFMSEcla'
long_name='W Ice Frozen Static Energy Cla'
units='K m s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=dse_avg(i,j,k)
enddo
enddo
enddo
name='DSE'
long_name='Dry static energy'
units='K'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=mse_avg(i,j,k)
enddo
enddo
enddo
name='MSE'
long_name='Moist static energy'
units='K'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=fmse_avg(i,j,k)
enddo
enddo
enddo
name='FMSE'
long_name='Ice Frozen Moist static energy'
units='K'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=fmsecla_avg(i,j,k)
enddo
enddo
enddo
name='FMSE_cla'
long_name='Ice Frozen Moist static energy Cla'
units='K'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=qv_avg(i,j,k)
enddo
enddo
enddo
name='QV'
long_name='Vapor mixing ratio'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=qt_avg(i,j,k)
enddo
enddo
enddo
name='QT'
long_name='Total mixing ratio'
units='kg/kg'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=wqv_avg(i,j,k)
enddo
enddo
enddo
name='WQV'
long_name='W QV flux'
units='kg/kg m s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=wqcl_avg(i,j,k)
enddo
enddo
enddo
name='WQCL'
long_name='W QCL flux'
units='kg/kg m s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=wqci_avg(i,j,k)
enddo
enddo
enddo
name='WQCI'
long_name='W QCI flux'
units='kg/kg m s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=wqpl_avg(i,j,k)
enddo
enddo
enddo
name='WQPL'
long_name='W QPL flux'
units='kg/kg m s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=wqpi_avg(i,j,k)
enddo
enddo
enddo
name='WQPI'
long_name='W QPI flux'
units='kg/kg m s-1'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

!write 2 D fields in 3 D outputs so that we can read them again
nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=wt_avg(i,j)
    enddo
  enddo
enddo
name='WTHETAS'
long_name='Surface wtheta flux'
units='K m /s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=wq_avg(i,j)
    enddo
  enddo
enddo
name='WQS'
long_name='Surface wq flux'
units='kg/kg m /s'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=wu_avg(i,j)
    enddo
  enddo
enddo
name='WU'
long_name='Surface wu flux'
units='m^2 / s^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=wv_avg(i,j)
    enddo
  enddo
enddo
name='WV'
long_name='Surface wv flux'
units='m^2 / s^2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=prec_xy_avg(i,j)/dt*86400
    enddo
  enddo
enddo
name='PRECIPS'
long_name='Surface Precip. Rate'
units='mm/day'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=shf_xy_avg(i,j)*rhow(1)*cp
    enddo
  enddo
enddo
name='SHFS'
long_name='Sensible Heat Flux'
units='W/m2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=lhf_xy_avg(i,j)*rhow(1)*lcond
    enddo
  enddo
enddo
name='LHFS'
long_name='Surface Latent Heat Flux'
units='W/m2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)
!100
nfields1=nfields1+1
do k=1,nz-1
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=tau_xy_avg(i,j)
    enddo
  enddo
enddo
name='TAUS'
long_name='Surface shear'
units='m2/s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)
!101











nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=ddx_uu(i,j,k)
enddo
enddo
enddo
name='D_UU_DX'
long_name='x-horizontal U divergence'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=ddy_vu(i,j,k)
enddo
enddo
enddo
name='D_VU_DY'
long_name='y-horizontal U divergence'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=ddz_wu(i,j,k)
enddo
enddo
enddo
name='D_U_DX'
long_name='w U divergence'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=ddx_uv(i,j,k)
enddo
enddo
enddo
name='D_UV_DX'
long_name='x-horizontal v divergence'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=ddy_vv(i,j,k)
enddo
enddo
enddo
name='D_VV_DY'
long_name='y-horizontal v divergence'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=ddz_wv(i,j,k)
enddo
enddo
enddo
name='D_WV_Dz'
long_name='w v divergence'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=ddx_uw(i,j,k)
enddo
enddo
enddo
name='D_UW_DX'
long_name='x-horizontal w divergence'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=ddy_vw(i,j,k)
enddo
enddo
enddo
name='D_VW_DY'
long_name='y-horizontal w divergence'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=ddz_wqn(i,j,k)
enddo
enddo
enddo
name='D_WW_Dz'
long_name='w qn divergence'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=ududx(i,j,k)
enddo
enddo
enddo
name='UD_U_DX'
long_name='x-horizontal U advection'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=vdudy(i,j,k)
enddo
enddo
enddo
name='VD_U_DY'
long_name='y-horizontal U advection'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=wdudz(i,j,k)
enddo
enddo
enddo
name='WD_U_DX'
long_name='w U advection'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)


nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=udvdx(i,j,k)
enddo
enddo
enddo
name='UDV_DX'
long_name='x-horizontal v advection'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=vdvdy(i,j,k)
enddo
enddo
enddo
name='VDV_DY'
long_name='y-horizontal v advection'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=wdvdz(i,j,k)
enddo
enddo
enddo
name='WDV_Dz'
long_name='w v advection'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=udwdx(i,j,k)
enddo
enddo
enddo
name='UD_W_DX'
long_name='x-horizontal w advection'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=vdwdy(i,j,k)
enddo
enddo
enddo
name='VD_W_DY'
long_name='y-horizontal w advection'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

nfields1=nfields1+1
do k=1,nz-1
do j=1,nym
do i=1,nxm
tmp(i,j,k)=wdwdz(i,j,k)
enddo
enddo
enddo
name='WD_W_Dz'
long_name='w w advection'
units='m2 /s2'
call compress3D(tmp,nxm,nym,nzm,name,long_name,units,savemombin,dompi,rank,nsubdomains)

!119






call task_barrier()

if (nfields .ne. nfields1) then
  if(masterproc) then
        print*,'write_moments error: nfields'
        write(myString,'(i10)') nfields        
        print*,myString
        write(myString,'(i10)') nfields1
        print*,myString
  endif
  call task_abort()
endif
if (masterproc) then
  close(46)
endif

if(nfields.ne.nfields1) then
    if(masterproc) print*,'write_fields3D error: nfields'
    call task_abort()
end if
if(masterproc) then
    if(RUN3D.or.savemomsep) then
       if(dogzip3D) call systemf('gzip -f '//filename)
       print*, 'Writting statistical-moments data. file:'//filename
    else
       print*, 'Appending statistical-moments data. file:'//filename
    end if
  endif

end subroutine write_moments

end module stat_moments
