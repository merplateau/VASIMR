subroutine fdfd_ini
    use the_whole_varibles
    implicit none

    integer :: nmlid

    namelist /fdfdConfig/ &
        & mTruncation, &
        & iswitch_specific_mode, &
        & m_specific, &
        & iswitch_vacum_dielectric_type
    
    namelist /boundaryConfig/ &
        & switch_plasma_vacuum_boundary, &
        & switch_z_boundary, &
        & switch_boundary_half_variables, &
        & switch_r_0_boundary, &
        & switch_source_type, &
        & switch_divE

    open(newunit=nmlid, file=trim(nmlfile), status="old")
    read(nmlid, nml=fdfdConfig)
    close(nmlid)

    if(iswitch_antenna_type==0)then
        m_start=0
        m_end=-m_start
        m_delta=+2
    else
        m_start=-mTruncation
        m_end=-m_start
        m_delta=+1

        if(iswitch_specific_mode == 1)then
            m_start = m_specific
            m_end = m_specific
            m_delta=+2
        end if
    endif

    open(newunit=nmlid, file=trim(nmlfile), status="old")
    read(nmlid, nml=boundaryConfig)
    close(nmlid)

    !switch_plasma_vacuum_boundary=43 !switch boundary of plasma-vacuum rp; 0:nothing more; 1:tangent continuous ;2:tangent continuous + Div(E)=0(vacuum);
    !3 :tangent continuous + Div(E)=0 + Jpn=0
    !4 :tangent continuous + Div(E)=0 + pdv{Dn}{r}continuous
    !5 :tangent continuous + Div(E)=0 + Dn condinuous
    !7 :tangent continuous + Div(E)=0[outside] + Div(D)=0[inside]
    !8 :tangent continuous + B_th continuous
    !1X : only Ez continuous + X[don't use]
    !2X : X + Ez continuous
    !switch_z_boundary=1 ! switch boundary condition of z=0/l; 1:conductor Boundary  ;  ; 3: z=0:conductor;z=zl:absorbing(For Vacuum Benchmark)
    !switch_boundary_half_variables=1 !Switch how to handle variables on half Grid in boundaries.
    !0 : by eqs(inside) or extrapolation (outside);
    !1: Interpolation/Extrapolation to find it's quantity or derivative at boundary
    !switch_r_0_boundary=2   !Switch how to handle the boundary condition of r=0. 1:old 2:new;
    !switch_source_type=1 !0: z=0 Source(For Benchmark) 1:Antenna
    !switch_divE=2!  0: ignore div(E) ; 1:not ignore div(E)(div(E) neq ===0); 2:ignore div(E) in Vacuum and nod ignore in Plasma Area
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!MUST NOT SET IT TO 0 (Except TESTing)
end subroutine fdfd_ini

subroutine maxwell_FDFD
    use the_whole_varibles
    IMPLICIT NONE
    call find_func_cputime_1_of_2  !-------------------1/2


    if(iswitch_RF2==  1  .and. mod(maxwellcount,2)==1)then
        om=2*pi*frequency2;
    else
        om=2*pi*frequency;
    endif
    omc2=(om/c)**2
    iom=i*om*mu0

    if (dev_Only4CalE == 1) then
        call dev_FDFD_pre
    end if

    call find_Dielectric_tensor
    call set_0


    do m= m_start,m_end,m_delta
        m_real=dble(m)
        im=cmplx(0.0d0, 1.0d0)*m
        call antenna
        call set_eq_b
        call set_A
        call do_reallocate
        call complex_coo_pardiso_solve(A_helicon,eq_xlast(m,:),eq_b)
        call totalsolve
    enddo
    !call find_irf_now
    call update_and_calculate
    maxwellcount=maxwellcount+1
    call record_loading_history(2)

    if (dev_Only4CalE == 1) then
        call dev_FDFD_post
    end if

    call find_func_cputime_2_of_2(func_time(9))  !-----2/2
end subroutine maxwell_FDFD

subroutine dev_FDFD_pre
    use the_whole_varibles
    implicit none

    

end subroutine dev_FDFD_pre

subroutine dev_FDFD_post
    use the_whole_varibles
    implicit none

    call record_dev_FDFD
    stop 'stop after dev_FDFD_post for only calculating E and dielectric tensor'

end subroutine dev_FDFD_post

subroutine do_reallocate
    use the_whole_varibles
    IMPLICIT NONE
    integer i_reallocate

    if(iswitch_RF2 ==  1 )then
        i_reallocate=1
    else
        i_reallocate=2
    endif

    if(maxwellcount>i_reallocate) return   !only for the first calculation of FDFD
    if(maxwellcount==0)then
        nzmax_pt=merge(A_helicon.nzmax,nzmax_pt,A_helicon.nzmax>nzmax_pt)
        call reallocate(A_helicon,2*A_helicon.nzmax)
        ! *2 to make sure the A_helicon of other m modes would not overflow!!!
    endif
    if(maxwellcount==i_reallocate) call reallocate(A_helicon,nzmax_pt)

end subroutine do_reallocate

subroutine find_Dielectric_tensor
    use the_whole_varibles
    implicit none
    integer*8::n_err,n_sec
    real*8::ikap,vl,vl_pic,fcr,err_kap
    real*8::tpp(1:2,1:2)! parallel or perp; i or e ; eV
    real*8::ti_ll,ti_perp,te_ll,te_perp
    complex*16::err_now,kap1,kap2,kapt,scat
    complex*16::kap_ll,func,fc1,fc2,kap_rec(1:nr,1:nz)
    complex*16::Ak(1:3,1:2),Bk(1:2) !(N=-1,0,+1;L=i,e),n_B=0
    Complex*16 :: omcr,omcz
    Complex*16 :: sigma_0,sigma_1,ve_iom2,ve_iom
    Complex*16 :: omce,omci,qim,qem !Defined and just used in find_Dielectric_tensor,mei
    Complex*16 :: ompi2,ompe2,omci2,omce2 !Defined and just used in find_Dielectric_tensor
    Complex*16 K_perp,K_phi,K_ll,coeff_ve_om,coeff_vi_om
    real*8     b0r,b0z,b0_abs,b02,om2,ve,vi,coeff_ompi2,coeff_ompe2
    !--------------------------coefficient--------------------------!

    ! change
    if (iswitch_analytic_profile == 1 .or. iswitch_analytic_profile == 2) then
        ni_midVal=ni ! 把真 ni 的值先存起来，解析的算完之后再赋值回去，最小改动源代码
        ni = ni_read
    end if
    ! change

    coeff_ompe2= qe_abs*qe_abs/(epson0*me)
    coeff_ompi2= qi*qi/(epson0*mi)
    !mei=me/mp
    qem=qe_abs/me
    qim=q_mass_i
    om2=om*om;

    ep=0.0D0*i
    si=0.0D0*i
    ep(:,:,1)=1.0D0+0.0D0*i
    ep(:,:,5)=1.0D0+0.0D0*i
    ep(:,:,9)=1.0D0+0.0D0*i


    do ir=1,nr
        do iz=1,nz

            if (ir > nrp) then
                if (iswitch_vacum_dielectric_type == 1 .or. iswitch_vacum_dielectric_type == 11) cycle
            else
                if (iswitch_vacum_dielectric_type == 11) cycle
            end if

                

            call find_ve_vi(ve,vi)
            
            if (iswitch_dielectric==3 .or. iswitch_dielectric==4) then
                ve=0.0d0;vi=0.0d0
            end if

            b0r=b0_DC(ir,iz,1)
            b0z=b0_DC(ir,iz,3)
            b0_abs=b0_DC(ir,iz,4)
            b02=b0_abs*b0_abs
            coeff_ve_om=1.0D0+i*(ve/om)
            ompe2=coeff_ompe2*ni(ir,iz)/coeff_ve_om
            omce=qem*b0_abs/coeff_ve_om
            omce2=omce*omce
            coeff_vi_om=1.0D0+i*(vi/om)
            ompi2=coeff_ompi2*ni(ir,iz)/coeff_vi_om
            if(iswitch_dielectric==2)ompi2=0.0d0
            omci=qim*b0_abs/coeff_vi_om
            omci2=omci*omci

            !-------------start kinetics Kpp---------------!

            !vl = 5e4 ! m/s !test
            vl = 0.0d0
            if (iswitch_dielectric==3 .or. iswitch_dielectric==4) then
                if(iswitch_v_closure_type==1)then
                    vl = vl_in_FDFD
                elseif(iswitch_v_closure_type==2)then
                    vl_pic = vl_closure_prev(ir,iz)
                    if(statistic_ready==1 .and. n_depo(ir,iz)>=1.0d-10) vl_pic = nu_depo(ir,iz,6)
                    vl_closure_prev(ir,iz) = (1.0d0-alpha_v)*vl_closure_prev(ir,iz) + alpha_v*vl_pic
                    vl = vl_closure_prev(ir,iz)
                endif
            endif

            ti_perp = Ek_ion_2D_r(ir,iz)/1.5;
            ti_ll = (Ek_ion_2D(ir,iz)-Ek_ion_2D_r(ir,iz))/1.5
            if(Ek_ion_2D(ir,iz)<1.5 .or. Ek_ion_2D_r(ir,iz)<1.) then
                ti_perp=1.5/1.5; ti_ll=1.5/3.0
            endif

            te_ll = te_in_fdfd(ir,iz)/3.0; te_perp = te_in_fdfd(ir,iz)/1.5
            !if(te(ir,iz)<fv_te_limit) then
            !    te_ll=fv_te_limit/3.0; te_perp=fv_te_limit/1.5
            !endif

            tpp(1,1) = ti_ll; tpp(2,1) = ti_perp
            tpp(1,2) = te_ll; tpp(2,2) = te_perp

            if (iswitch_T_closure_type == 1) then
                tpp(1,1) = till_uni; tpp(2,1) = tipp_uni
                tpp(1,2) = tell_uni; tpp(2,2) = t2pp_uni
            elseif (iswitch_T_closure_type == 2) then
                tpp(1,1) = TminLim
                tpp(2,1) = TminLim
                if (statistic_ready==1 .and. n_depo(ir,iz)>=1.0d-10) then
                    if (t_depo(ir,iz,1)>=TminLim) tpp(1,1) = t_depo(ir,iz,1)
                    if (t_depo(ir,iz,2)>=TminLim) tpp(2,1) = t_depo(ir,iz,2)
                endif
            end if

            if (iswitch_dielectric==3 .or. iswitch_dielectric==4) then
                if (k_closure_type==0) then
                    kap1=200.0d0
                    kap2=kap1*(1.0d0+1.0d-5)
                    kapt=kap2
                    err_kap=1.0d-3
                    err_now=err_kap+1.0d0
                    n_err=0
                    n_sec=0

                    do while (n_err<3)
                        n_sec=n_sec+1
                        call find_kfc(kap1,tpp,vl,Ak,Bk,omci,omce,ompi2,ompe2,fc1)
                        call find_kfc(kap2,tpp,vl,Ak,Bk,omci,omce,ompi2,ompe2,fc2)
                        fc1=real(fc1)
                        fc2=real(fc2)
                        scat=(fc1-fc2)/(kap1-kap2)
                        kapt=kapt-fc2/scat
                        err_now=(kap2-kapt)/(kap2+kapt)
                        kap1=kap2
                        kap2=kapt
                        if (abs(err_now)==0.0d0) kap1=kapt*(1.0d0+1.0d-5)
                        if (abs(err_now)<err_kap) then
                            n_err=n_err+1
                        else
                            n_err=0
                        endif
                    enddo

                    kap_ll=kap1
                elseif (k_closure_type==1) then
                    kap_ll=dcmplx(k_constant,0.0d0)
                elseif (k_closure_type==2) then
                    kap_ll = kap_ll_read(ir,iz)
                endif

                kap_rec(ir,iz)=kap_ll
                call find_kfc(kap_ll,tpp,vl,Ak,Bk,omci,omce,ompi2,ompe2,func)
            endif


            if(iswitch_dielectric==3)then
                ! ion kinetics and electron cold
                K_perp=1.0d0+ompi2/(2.0d0*om)*(Ak(1,1)+Ak(3,1))&
                    &-ompe2/om2*(om-kap_ll*vl)**2/((om-kap_ll*vl)**2-omce2)
                K_phi=ompi2/(2.0d0*om)*(Ak(1,1)-Ak(3,1))&
                    &-ompe2*omce/om2*(om-kap_ll*vl)/((om-kap_ll*vl)**2-omce2)
                 K_ll=1.0d0+2.0d0*ompi2/(kap_ll*2.0d0*qi*tpp(2,1)/mi)*(vl/om+Bk(1))&
                   &-ompe2/(om*(om-kap_ll*vl))
                !K_ll=1.0d0-ompi2/(om*om)&
                !    &-ompe2/(om*om)

            elseif(iswitch_dielectric==4)then
                ! ion and electron kinetics
                K_perp=1.0d0+ompi2/(2.0d0*om)*(Ak(1,1)+Ak(3,1))&
                    &+ompe2/(2.0d0*om)*(Ak(1,2)+Ak(3,2))
                K_phi=ompi2/(2.0d0*om)*(Ak(1,1)-Ak(3,1))&
                    &+ompe2/(2.0d0*om)*(Ak(1,2)-Ak(3,2))
                K_ll=1.0d0+2.0d0*ompi2/(kap_ll*2.0d0*qi*tpp(2,1)/mi)*(vl/om+Bk(1))&
                    &+2.0d0*ompe2/(kap_ll*2.0d0*qe_abs*tpp(2,2)/me)*(vl/om+Bk(2))

            else
                ! cold plasma
                K_perp=0.;K_phi=0.;K_ll=0.;

                coeff_ompi2= qi*qi/(epson0*mi)
                coeff_vi_om=1.0D0+i*(vi/om)
                ompi2=coeff_ompi2*ni(ir,iz)/coeff_vi_om
                if(iswitch_dielectric==2)ompi2=0.0d0
                omci=q_mass_i*b0_abs/coeff_vi_om
                omci2=omci*omci

                K_perp=K_perp-ompi2/(om2-omci2)
                K_phi=K_phi+ompi2*omci/(om*(om2-omci2))
                K_ll=K_ll-ompi2/om2

                K_perp=K_perp+1-ompe2/(om2-omce2);
                K_phi=K_phi-ompe2*omce/(om*(om2-omce2))
                K_ll=K_ll+1-ompe2/om2;
            endif


            si(ir,iz,1,1)=b0z**2/b02*K_perp+b0r**2/b02*K_ll-1.0D0
            si(ir,iz,1,2)=-b0z/b0_abs*i*K_phi
            si(ir,iz,1,3)=b0r*b0z/b02*(K_ll-K_perp)
            si(ir,iz,2,2)=K_perp-1.0D0
            si(ir,iz,2,3)=-b0r/b0_abs*i*K_phi
            si(ir,iz,3,3)=b0z**2/b02*K_ll+b0r**2/b02*K_perp-1.0D0
            si(ir,iz,2,1)=-1.0D0*si(ir,iz,1,2)
            si(ir,iz,3,1)=si(ir,iz,1,3)
            si(ir,iz,3,2)=-1.0D0*si(ir,iz,2,3)
        enddo
    enddo
    si=-si*i*epson0*om


    si(nrp,:,1,:)=0.0D0*i

    if (iswitch_vacum_dielectric_type == 1) si(nrp+1:nr,:,:,:)=0.0D0*i
    if (iswitch_vacum_dielectric_type == 11) si(:,:,:,:)=0.0D0*i

    ep(:,:,1)=1.0D0+i/(epson0*om)*si(:,:,1,1)
    ep(:,:,5)=1.0D0+i/(epson0*om)*si(:,:,2,2)
    ep(:,:,9)=1.0D0+i/(epson0*om)*si(:,:,3,3)
    ep(:,:,4)=i/(epson0*om)*si(:,:,2,1)
    ep(:,:,7)=i/(epson0*om)*si(:,:,3,1)
    ep(:,:,8)=i/(epson0*om)*si(:,:,3,2)
    ep(:,:,2)=-ep(:,:,4)
    ep(:,:,3)=ep(:,:,7)
    ep(:,:,6)=-ep(:,:,8)

    if (iswitch_analytic_profile == 1 .or. iswitch_analytic_profile == 2) then
        ni=ni_midVal
    end if
    continue
end subroutine find_Dielectric_tensor

subroutine set_0
    use the_whole_varibles
    implicit none
    ptotal=0.0D0
    ptotal_boundary=0.0D0
    ptotal_inside=0.0D0
    !ptotal_BI=0.0D0
    power_depo_2D=0.0D0*i
    power_depo_rthz=0.
    ptotm=0.0D0
    A_helicon.col=0
    A_helicon.row=0
    A_helicon.val=0.0D0*i
    e_m=0.0D0; e_int=0.0D0; ja_m_nor=0.0D0;  ja_record=0.0D0;
    e_AC=0.0D0;e_record=0.0D0;ja_AC=0.0D0;!ptotm=0.0D0
    b_m=0.0D0;b_AC=0.0D0;b_record=0.0D0
    jp_m=0.0D0;jp_record=0.0D0;jp_AC=0.0D0

end subroutine set_0

subroutine set_eq_b
    use the_whole_varibles
    IMPLICIT NONE
    eq_b(:)=0.0D0*i
    if(switch_source_type==0)then
        ir=1;i3=1;
        if(m==1.or.m==-1)then
            eq_b(3*(nz*ir+iz-nz)-3+i3)=1.0D0/(2.0D0*DBesselJZero(m))+0.0D0*i!Er(ir,1)
        else
            eq_b(3*(nz*ir+iz-nz)-3+i3)=0.0D0*i!Er(ir,1)
        endif
        i3=2
        eq_b(3*(nz*ir+iz-nz)-3+i3)=i/(2.0D0*DBesselJZero(m))*(Bessel_Jn(m-1,DBesselJZero(m)*r(ir)/rl)&
            &-Bessel_Jn(m+1,DBesselJZero(m)*r(ir)/rl))+0.0D0*i!Eth(ir,1)
        Do ir=2,nr
            i3=1
            eq_b(3*(nz*ir+iz-nz)-3+i3)=m_real*rl/(DBesselJZero(m)**2*r(ir))*Bessel_Jn(m,DBesselJZero(m)*r(ir)/rl)+0.0D0*i!Er(ir,1)
            i3=2
            eq_b(3*(nz*ir+iz-nz)-3+i3)=i/(2.0D0*DBesselJZero(m))*(Bessel_Jn(m-1,DBesselJZero(m)*r(ir)/rl)&
                &-Bessel_Jn(m+1,DBesselJZero(m)*r(ir)/rl))+0.0D0*i!Eth(ir,1)
        EndDo
    elseif(switch_source_type==1)then
        Do ir=1,nr
            Do iz=1,nz
                Do i3=1,n3
                    eq_b(3*(nz*ir+iz-nz)-3+i3)=iom*ja_m_nor(ir,iz,i3)
                EndDo
            EndDo
        EndDo
    endif
    continue
    !pause
end subroutine set_eq_b

function DBesselJZero(nbessel)
    implicit none
    integer*4 :: nbessel
    real*8 :: DBesselJZero
    real*8 :: listpt(6)=(/ 1.841D0, 3.054D0, 4.201D0, &
        &5.317D0, 6.415D0, 7.501D0 /)
    if(nbessel==0)then
        DBesselJZero=3.831705970207512D0
    elseif(nbessel>=1 .and. nbessel<=6)then
        DBesselJZero=listpt(nbessel)
    else
        write(*,*)'DBesselJZero(x_plasma) xmax=6!!!!!! Please Check'
        pause
        DBesselJZero=0
    endif
    return
end function DBesselJZero

subroutine set_A
    use the_whole_varibles
    implicit none
    integer*4 :: lcount,scount
    A_helicon.n=3*nr*nz
    A_helicon.nzmax=0
    A_helicon.row=0
    A_helicon.col=0
    A_helicon.val=0.0D0*i
    iswitch_output_Region=1


    if(iswitch_output_Region==1)open(unit=545,file=trim(outputDir)//trim('Region.txt'),status='REPLACE',iostat=ierror)
    Call Find_Region(nr,nz,n_vac,nr_vac,nr_met,nz_vac,FindRegion) !,nrz_diploe,iswtich_inner_dipole
    do ir=1,nr
        do iz=1,nz
            !!!plasma region(2 <= ir <= nrp-1)
            Region=FindRegion(ir,iz)
            if(iswitch_output_Region==1)then
                write(545,"(I5,I5,I5)")ir,iz,Region!TEST
            endif
            if(Region==1)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(ir,iz)
                r(:)=r(:)+1.0D0/2.0D0 * dr
                call set(lcount,0,m_real**2/r(ir)**2+2.0e0/dz22-omc2*(ep(ir+1,iz,1)+ep(ir,iz,1))/2.0e0+0.0e0*i,A_helicon)!Er(ir,iz)
                call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz+1)
                call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz-1)
                call set(lcount,1,im/r(ir)**2 *(1.0D0/2.0D0-r(ir)/dr)-omc2*ep(ir,iz,2)/2.0D0+0.0D0*i,A_helicon)!Eth(ir,iz)
                call set(lcount,1+3*nz,im/r(ir)**2 *(1.0D0/2.0D0+r(ir)/dr)-omc2*ep(ir+1,iz,2)/2.0D0+0.0D0*i,A_helicon)!Eth(ir+1,iz)
                call set(lcount,2+3*nz,1.0D0/(dr*dz)-omc2*(ep(ir+1,iz,3)+ep(ir+1,iz+1,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                call set(lcount,2-3,1.0D0/(dr*dz)-omc2*(ep(ir,iz-1,3)+ep(ir+1,iz,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir,iz-1)
                call set(lcount,2,-1.0D0/(dr*dz)-omc2*(ep(ir,iz,3)+ep(ir,iz+1,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir,iz)
                call set(lcount,2+3*nz-3,-1.0D0/(dr*dz)-omc2*(ep(ir+1,iz-1,3)+ep(ir+1,iz,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir+1,iz-1)

                r(:)=r(:)-1.0D0/2.0D0 * dr


                lcount=3*(nz*ir+iz-nz)-1!!!Eth(ir,iz)
                call set(lcount,-1,im/r(ir)**2 *(-1.0D0/2.0D0+r(ir)/dr)-omc2*(ep(ir+1,iz,4)+ep(ir,iz,4))/4.0D0+0.0D0*i,A_helicon)!Er(ir,iz)
                call set(lcount,-1-3*nz,im/r(ir)**2 * (-1.0D0/2.0D0-r(ir)/dr)-omc2*(ep(ir,iz,4)+ep(ir-1,iz,4))/4.0D0+0.0D0*i,A_helicon)!Er(ir-1,iz)
                call set(lcount,0,1.0D0/r(ir)**2 + 2.0D0/dr22+2.0D0/dz22-omc2*ep(ir,iz,5)+0.0D0*i,A_helicon)!Eth(ir,iz)
                call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Eth(ir+1,iz)
                call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Eth(ir-1,iz)
                call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Eth(ir,iz+1)
                call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Eth(ir,iz-1)
                call set(lcount,1,im/(r(ir)*dz)-omc2*(ep(ir,iz,6)+ep(ir,iz+1,6))/4.0D0+0.0D0*i,A_helicon)!Ez(ir,iz)
                call set(lcount,1-3,-im/(r(ir)*dz)-omc2*(ep(ir,iz-1,6)+ep(ir,iz,6))/4.0D0+0.0D0*i,A_helicon)!Ez(ir,iz-1)


                lcount=3*(nz*ir+iz-nz)!!!Ez(ir,iz)
                call set(lcount,-2,-1.0D0/(2.0D0*r(ir)*dz)-1.0D0/(dr*dz)-omc2*(ep(ir,iz,7)+ep(ir+1,iz,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir,iz)
                call set(lcount,-2-3*nz,-1.0D0/(2.0D0*r(ir)*dz)+1.0D0/(dr*dz)-omc2*(ep(ir-1,iz,7)+ep(ir,iz,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir-1,iz)
                call set(lcount,-2+3,1.0D0/(2.0D0*r(ir)*dz)+1.0D0/(dr*dz)-omc2*(ep(ir,iz+1,7)+ep(ir+1,iz+1,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir,iz+1)
                call set(lcount,-2-3*nz+3,+1.0D0/(2.0D0*r(ir)*dz)-1.0D0/(dr*dz)-omc2*(ep(ir-1,iz+1,7)+ep(ir,iz+1,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir-1,iz+1)
                call set(lcount,-1,-im/(r(ir)*dz)-omc2*ep(ir,iz,8)/2.0D0+0.0D0*i,A_helicon)!Eth(ir,iz)
                call set(lcount,-1+3,im/(r(ir)*dz)-omc2*ep(ir,iz+1,8)/2.0D0+0.0D0*i,A_helicon)!Eth(ir,iz+1)
                call set(lcount,0,m_real**2/(r(ir)**2)+2.0D0/dr22-omc2*(ep(ir,iz,9)+ep(ir,iz+1,9))/2.0D0+0.0D0*i,A_helicon)!Ez(ir,iz)
                call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir-1,iz)

                !!!Vacuum Region(nrp <= ir <= nr-1)
            elseif(Region==2 .or. Region==9 .or. Region==14)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(ir,iz)
                r(:)=r(:)+1.0D0/2.0D0 * dr

                if(ir==nr_vac(1))then
                    call set(lcount,0+3*nz,1.0D0/(2.0D0*(rp+dr))+1.0D0/dr+0.0D0*i,A_helicon)!Er(nrp+1,iz)
                    call set(lcount,0,1.0D0/(2.0D0*(rp+dr))-1.0D0/dr+0.0D0*i,A_helicon)!Er(nrp,iz)
                    call set(lcount,1,im/rp+0.0D0*i,A_helicon)!Eth(nrp+1,iz)
                    call set(lcount,2,1.0D0/dz+0.0D0*i,A_helicon)!Ez(nrp+1,iz)
                    call set(lcount,2-3,-1.0D0/dz+0.0D0*i,A_helicon)!Ez(nrp+1,iz-1)
                else
                    call set(lcount,0,(1.0D0+m_real**2)/r(ir)**2+2.0D0/dz22+2.0D0/dr22-0.5D0*omc2*(ep(ir,iz,1)+ep(ir+1,iz,1))+0.0D0*i,A_helicon)!Er(ir,iz)
                    call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz+1)
                    call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz-1)
                    call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Er(ir+1,iz)
                    call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Er(ir-1,iz)
                    call set(lcount,1,im/r(ir)**2-0.25D0*omc2*(ep(ir,iz,2)+ep(ir+1,iz,2))+0.0D0*i,A_helicon)!Eth(ir,iz)
                    call set(lcount,1+3*nz,im/r(ir)**2-0.25D0*omc2*(ep(ir,iz,2)+ep(ir+1,iz,2))+0.0D0*i,A_helicon)!Eth(ir+1,iz)
                endif
                r(:)=r(:)-1.0D0/2.0D0 * dr

                lcount=3*(nz*ir+iz-nz)-1!!!Eth(ir,iz)
                call set(lcount,-1,-im/r(ir)**2-0.5D0*omc2*ep(ir,iz,4)+0.0D0*i,A_helicon)!Er(ir,iz)
                call set(lcount,-1-3*nz,-im/r(ir)**2-0.5D0*omc2*ep(ir,iz,4)+0.0D0*i,A_helicon)!Er(ir-1,iz)
                call set(lcount,0,(1.0D0+m_real**2)/r(ir)**2 + 2.0D0/dr22+2.0D0/dz22-omc2*ep(ir,iz,5)+0.0D0*i,A_helicon)!Eth(ir,iz)
                call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Eth(ir+1,iz)
                call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Eth(ir-1,iz)
                call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Eth(ir,iz+1)
                call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Eth(ir,iz-1)

                lcount=3*(nz*ir+iz-nz)!!!Ez(ir,iz)
                if(Region==2 .or. Region==9)then
                    call set(lcount,0,m_real**2/(r(ir)**2)+2.0D0/dr22+2.0D0/dz22-0.5D0*omc2*(ep(ir,iz,9)+ep(ir,iz+1,9))+0.0D0*i,A_helicon)!Ez(ir,iz)
                    call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                    call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir-1,iz)
                    call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Ez(ir,iz+1)
                    call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Ez(ir,iz-1)
                elseif(Region==14)then
                    call set(lcount,0,m_real**2/(r(ir)**2)+2.0D0/dr22+1.0D0/(3.0D0*dz22)-0.5D0*omc2*(ep(ir,iz,9)+ep(ir,iz+1,9))+0.0D0*i,A_helicon)!Ez(ir,iz)
                    call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                    call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir-1,iz)
                    call set(lcount,0-3,-2.0D0/(3.0D0*dz22)+0.0D0*i,A_helicon)!Ez(ir,iz-1)
                    call set(lcount,0-2*3,1.0D0/(3.0D0*dz22)+0.0D0*i,A_helicon)!Ez(ir,iz-2)
                endif
                !!!!!!!!!!!!!!!!!!!! Boundary Condition !!!!!!!!!!!!!!!!!!!!!!!!!!
                !----------------------   z=0 [z=0] (iz=1)   ----------------------!
            elseif(Region==4.or.Region==12.or.Region==22)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(ir,1)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(ir,1)
                lcount=3*(nz*ir+iz-nz)-1!!!Eth(ir,1)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(ir,1)
                lcount=3*(nz*ir+iz-nz)!!!Ez(ir,1)
                if(switch_boundary_half_variables==0)then
                    call set(lcount,-2,-1.0D0/(2.0D0*r(ir)*dz)-1.0D0/(dr*dz)-omc2*(ep(ir,iz,7)+ep(ir+1,iz,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir,iz)
                    call set(lcount,-2-3*nz,-1.0D0/(2.0D0*r(ir)*dz)+1.0D0/(dr*dz)-omc2*(ep(ir-1,iz,7)+ep(ir,iz,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir-1,iz)
                    call set(lcount,-2+3,1.0D0/(2.0D0*r(ir)*dz)+1.0D0/(dr*dz)-omc2*(ep(ir,iz+1,7)+ep(ir+1,iz+1,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir,iz+1)
                    call set(lcount,-2-3*nz+3,+1.0D0/(2.0D0*r(ir)*dz)-1.0D0/(dr*dz)-omc2*(ep(ir-1,iz+1,7)+ep(ir,iz+1,7))/8.0D0+0.0D0*i,A_helicon)!Er(ir-1,iz+1)
                    call set(lcount,-1,-im/(r(ir)*dz)-omc2*ep(ir,iz,8)/2.0D0+0.0D0*i,A_helicon)!Eth(ir,iz)
                    call set(lcount,-1+3,im/(r(ir)*dz)-omc2*ep(ir,iz+1,8)/2.0D0+0.0D0*i,A_helicon)!Eth(ir,iz+1)
                    call set(lcount,0,m_real**2/(r(ir)**2)+2.0D0/dr22-omc2*(ep(ir,iz,9)+ep(ir,iz+1,9))/2.0D0+0.0D0*i,A_helicon)!Ez(ir,iz)
                    call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                    call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir-1,iz)
                elseif(switch_boundary_half_variables==1)then
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(ir,1)
                    if(jieshu==1)then
                        call set(lcount,0+3,-1.0D0+0.0D0*i,A_helicon)!Ez(ir,2)
                    else
                        call set(lcount,0+3,-26765.0d0/9129.0d0+0.0D0*i,A_helicon)!Ez(ir,2)
                        call set(lcount,0+3*2,11630.0d0/3043.0d0+0.0D0*i,A_helicon)!Ez(ir,3)
                        call set(lcount,0+3*3,-8590.0d0/3043.0D0+0.0D0*i,A_helicon)!Ez(ir,4)
                        call set(lcount,0+3*4,10205.0d0/9129.0d0+0.0D0*i,A_helicon)!Ez(ir,5)
                        call set(lcount,0+3*5,-563.0d0/3043.0d0+0.0D0*i,A_helicon)!Ez(ir,6)
                    endif
                else
                    write(*,*)'switch_boundary_half_variables ERROR PLEASECHECK'
                    pause
                endif
            elseif(Region==7)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(ir,1)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(ir,1)
                lcount=3*(nz*ir+iz-nz)-1!!!Eth(ir,1)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(ir,1)
                lcount=3*(nz*ir+iz-nz)!!!Ez(ir,1)
                if(switch_boundary_half_variables==0)then
                    call set(lcount,0,m_real**2/(r(ir)**2)+2.0D0/dr22+2.0D0/dz22-0.5D0*omc2*(ep(ir,iz,9)+ep(ir,iz+1,9))+0.0D0*i,A_helicon)!Ez(ir,iz)
                    call set(lcount,0+3*nz,-1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                    call set(lcount,0-3*nz,1.0D0/(2.0D0*r(ir)*dr)-1.0D0/dr22+0.0D0*i,A_helicon)!Ez(ir-1,iz)
                    call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Ez(ir,iz+1)
                    call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Ez(ir,iz-1)
                elseif(switch_boundary_half_variables==1)then
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(ir,1)
                    if(jieshu==1)then
                        call set(lcount,0+3,-1.0D0+0.0D0*i,A_helicon)!Ez(ir,2)
                    else
                        call set(lcount,0+3,-26765.0d0/9129.0d0+0.0D0*i,A_helicon)!Ez(ir,2)
                        call set(lcount,0+3*2,11630.0d0/3043.0d0+0.0D0*i,A_helicon)!Ez(ir,3)
                        call set(lcount,0+3*3,-8590.0d0/3043.0D0+0.0D0*i,A_helicon)!Ez(ir,4)
                        call set(lcount,0+3*4,10205.0d0/9129.0d0+0.0D0*i,A_helicon)!Ez(ir,5)
                        call set(lcount,0+3*5,-563.0d0/3043.0d0+0.0D0*i,A_helicon)!Ez(ir,6)
                    endif
                else
                    write(*,*)'switch_boundary_half_variables ERROR PLEASECHECK'
                    pause
                endif
                !----------------------   z=0 [z=0] (iz=1)   ----------------------!



                !----------------------   z=zl [z=l] (iz=nz)   ----------------------!
            elseif(Region==5 .or. Region==8 .or. Region==13.or.Region==24)then
                if(switch_z_boundary==3)then
                    call set(3*nz*ir-2,0,23.0D0/24.0D0-i*dz/2.0D0*sqrt((om/c)**2-(DBesselJZero(m)/rl)**2)+0.0D0*i,A_helicon)!!!Er(ir,nz)
                    call set(3*nz*ir-2,-3,-7.0D0/8.0D0-i*dz/2.0D0*sqrt((om/c)**2-(DBesselJZero(m)/rl)**2)+0.0D0*i,A_helicon)!Er(ir,nz-1)
                    call set(3*nz*ir-2,-2*3,-1.0D0/8.0D0+0.0D0*i,A_helicon)!Er(ir,nz-2)
                    call set(3*nz*ir-2,-3*3,1.0D0/24.0D0+0.0D0*i,A_helicon)!Er(ir,nz-3)

                    call set(3*nz*ir-1,0,11.0D0/6.0D0-i*dz*sqrt((om/c)**2-(DBesselJZero(m)/rl)**2)+0.0D0*i,A_helicon)!!!Eth(ir,nz)
                    call set(3*nz*ir-1,-3,-3.0D0+0.0D0*i,A_helicon)!Eth(ir,nz-1)
                    call set(3*nz*ir-1,-2*3,3.0D0/2.0D0+0.0D0*i,A_helicon)!Eth(ir,nz-2)
                    call set(3*nz*ir-1,-3*3,-1.0D0/3.0D0+0.0D0*i,A_helicon)!Eth(ir,nz-3)
                else
                    lcount=3*(nz*ir+iz-nz)-2!!!Er(ir,nz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(ir,nz)
                    lcount=3*(nz*ir+iz-nz)-1!!!Eth(ir,nz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(ir,nz)
                endif

                lcount=3*(nz*ir+iz-nz)!!!Ez(ir,nz)
                if(switch_boundary_half_variables==0)then
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(ir,nz)
                    call set(lcount,0-3,-2.0D0+0.0D0*i,A_helicon)!Ez(ir,nz-1)
                    call set(lcount,0-2*3,1.0D0+0.0D0*i,A_helicon)!Ez(ir,nz-2)
                elseif(switch_boundary_half_variables==1)then
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(ir,nz)
                    if(jieshu==1)then
                        call set(lcount,0-3,-1.0D0+0.0D0*i,A_helicon)!Ez(ir,nz-1)
                    else
                        call set(lcount,0-3,-335.0d0/563.0d0+0.0D0*i,A_helicon)!Ez(ir,nz-1)
                        call set(lcount,0-3*2,-1430.0d0/1689.0d0+0.0D0*i,A_helicon)!Ez(ir,nz-2)
                        call set(lcount,0-3*3,370.0d0/563.0d0+0.0D0*i,A_helicon)!Ez(ir,nz-3)
                        call set(lcount,0-3*4,-145.0d0/563.0d0+0.0D0*i,A_helicon)!Ez(ir,nz-4)
                        call set(lcount,0-3*5,71.0d0/1689.0d0+0.0D0*i,A_helicon)!Ez(ir,nz-5)
                    endif
                endif
                !----------------------   z=zl [z=l] (iz=nz)   ----------------------!



                !----------------------   r=rl [r=R] (ir=nr)   ----------------------!
            elseif(Region==3.or.Region==6.or.Region==21.or.Region==25.or.Region==26)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(nr,iz)
                if(switch_boundary_half_variables==0)then
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(nr,iz)
                    call set(lcount,0-3*nz,-2.0D0+0.0D0*i,A_helicon)!Er(nr-1,iz)
                    call set(lcount,0-2*3*nz,1.0D0+0.0D0*i,A_helicon)!Er(nr-2,iz)
                elseif(switch_boundary_half_variables==1)then
                    call set(lcount,0,1.0D0/(2.0D0*rl) + 1.0D0/dr+0.0D0*i,A_helicon) !Er(nr,iz)
                    call set(lcount,0-3*nz,1.0D0/(2.0D0*rl) - 1.0D0/dr+0.0D0*i,A_helicon) !Er(nr-1,iz)
                endif

                lcount=3*(nz*ir+iz-nz)-1!!!Eth(nr,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(nr,iz)
                lcount=3*(nz*ir+iz-nz)!!!Ez(nr,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(nr,iz)
                !----------------------   r=rl [r=R] (ir=nr)   ----------------------!



                !----------------------   r=0 [r=0] (ir=1)   ----------------------!
            elseif(Region==10)then
                if(m==0)then
                    lcount=3*(nz*ir+iz-nz)-2!!!Er(1,iz)
                    if(switch_boundary_half_variables==1)then
                        call set(lcount,0,3.0D0/2.0D0+0.0D0*i,A_helicon)!Er(1,iz)
                        call set(lcount,0+3*nz,-1.0D0/2.0D0+0.0D0*i,A_helicon)!Er(2,iz)
                    endif
                    lcount=3*(nz*ir+iz-nz)-1!!!Eth(1,iz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(1,iz)
                    lcount=3*(nz*ir+iz-nz)!!!Ez(1,iz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(1,iz)
                    call set(lcount,0+3*nz,-1.0D0+0.0D0*i,A_helicon)!Ez(2,iz)
                elseif(m==-1.or.m==1)then
                    lcount=3*(nz*ir+iz-nz)-2!!!Er(1,iz)
                    if(switch_r_0_boundary==2)then
                        if(iz==1 .or. iz==nz)then
                            call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(1,1)
                        else
                            r(:)=r(:)+1.0D0/2.0D0 * dr
                            call set(lcount,0,m_real**2/r(ir)**2+2.0D0/dz22-0.5D0*omc2*(ep(ir,iz,1)+ep(ir+1,iz,1))+0.0D0*i,A_helicon)!Er(ir,iz)
                            call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz+1)
                            call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz-1)
                            call set(lcount,1,im/r(ir)**2 *(1.0D0/2.0D0-r(ir)/dr)-0.25D0*omc2*(ep(ir,iz,2)+ep(ir+1,iz,2))+0.0D0*i,A_helicon)!Eth(ir,iz)
                            call set(lcount,1+3*nz,im/r(ir)**2 *(1.0D0/2.0D0+r(ir)/dr)-0.25D0*omc2*(ep(ir,iz,2)+ep(ir+1,iz,2))+0.0D0*i,A_helicon)!Eth(ir+1,iz)
                            call set(lcount,2+3*nz,1.0D0/(dr*dz)+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                            call set(lcount,2-3,1.0D0/(dr*dz)+0.0D0*i,A_helicon)!Ez(ir,iz-1)
                            call set(lcount,2,-1.0D0/(dr*dz)+0.0D0*i,A_helicon)!Ez(ir,iz)
                            call set(lcount,2+3*nz-3,-1.0D0/(dr*dz)+0.0D0*i,A_helicon)!Ez(ir+1,iz-1)
                            r(:)=r(:)-1.0D0/2.0D0 * dr
                        endif
                    elseif(switch_r_0_boundary==1)then
                        if(switch_boundary_half_variables==1)then
                            call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(1,iz)
                            call set(lcount,0+3*nz,-1.0D0+0.0D0*i,A_helicon)!Er(2,iz)
                        endif
                    else
                        write(*,*)'switch_r_0_boundary ERROR please check(=1 or 2)'
                    endif
                    lcount=3*(nz*ir+iz-nz)-1!!!Eth(1,iz)
                    if(switch_r_0_boundary==2)then
                        call set(lcount,0,m_real*i*(1.0D0+0.0D0*i),A_helicon)!Eth(1,iz)
                        call set(lcount,-1,3.0D0/2.0D0+0.0D0*i,A_helicon)!Er(1,iz)
                        call set(lcount,-1+3*nz,-1.0D0/2.0D0+0.0D0*i,A_helicon)!Er(2,iz)
                    elseif(switch_r_0_boundary==1)then
                        call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(1,iz)
                        call set(lcount,0+3*nz,-1.0D0+0.0D0*i,A_helicon)!Eth(2,iz)
                    endif
                    lcount=3*(nz*ir+iz-nz)!!!Ez(1,iz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(1,iz)
                else
                    lcount=3*(nz*ir+iz-nz)-2!!!Er(1,iz)
                    if(switch_boundary_half_variables==1)then
                        call set(lcount,0,3.0D0/2.0D0+0.0D0*i,A_helicon)!Er(1,iz)
                        call set(lcount,0+3*nz,-1.0D0/2.0D0+0.0D0*i,A_helicon)!Er(2,iz)
                    endif
                    lcount=3*(nz*ir+iz-nz)-1!!!Eth(1,iz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(1,iz)
                    lcount=3*(nz*ir+iz-nz)!!!Ez(1,iz)
                    call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(1,iz)
                endif
            elseif(Region==11)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(ir,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Er(ir,iz)
                lcount=3*(nz*ir+iz-nz)-1!!!Eth(ir,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Erh(ir,iz)
                lcount=3*(nz*ir+iz-nz)!!!Ez(ir,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(ir,iz)

                !----------------------   dipole upper   ----------------------!
            elseif(Region==23)then
                lcount=3*(nz*ir+iz-nz)-2!!!Er(nr,iz)
                r(:)=r(:)+1.0D0/2.0D0 * dr
                call set(lcount,0,m_real**2/r(ir)**2+2.0D0/dz22-omc2*(ep(ir+1,iz,1)+ep(ir,iz,1))/2.0D0+0.0D0*i,A_helicon)!Er(ir,iz)
                call set(lcount,0+3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz+1)
                call set(lcount,0-3,-1.0D0/dz22+0.0D0*i,A_helicon)!Er(ir,iz-1)
                call set(lcount,1,im/r(ir)**2 *(1.0D0/2.0D0-r(ir)/dr)-omc2*ep(ir,iz,2)/2.0D0+0.0D0*i,A_helicon)!Eth(ir,iz)
                call set(lcount,1+3*nz,im/r(ir)**2 *(1.0D0/2.0D0+r(ir)/dr)-omc2*ep(ir+1,iz,2)/2.0D0+0.0D0*i,A_helicon)!Eth(ir+1,iz)
                call set(lcount,2+3*nz,1.0D0/(dr*dz)-omc2*(ep(ir+1,iz,3)+ep(ir+1,iz+1,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir+1,iz)
                call set(lcount,2-3,1.0D0/(dr*dz)-omc2*(ep(ir,iz-1,3)+ep(ir+1,iz,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir,iz-1)
                call set(lcount,2,-1.0D0/(dr*dz)-omc2*(ep(ir,iz,3)+ep(ir,iz+1,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir,iz)
                call set(lcount,2+3*nz-3,-1.0D0/(dr*dz)-omc2*(ep(ir+1,iz-1,3)+ep(ir+1,iz,3))/8.0D0+0.0D0*i,A_helicon)!Ez(ir+1,iz-1)
                r(:)=r(:)-1.0D0/2.0D0 * dr
                lcount=3*(nz*ir+iz-nz)-1!!!Eth(nr,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Eth(nr,iz)
                lcount=3*(nz*ir+iz-nz)!!!Ez(nr,iz)
                call set(lcount,0,1.0D0+0.0D0*i,A_helicon)!Ez(nr,iz)
                !----------------------   dipole upper   ----------------------!

            endif
        enddo
    enddo
    !pause

    if(iswitch_output_Region==1) close(545)!pause;!TEST
    call findnzmax(A_helicon)
    !pause
end subroutine set_A

subroutine set(l1,sml1,Als,A1)
    use Types
    implicit none
    integer*4 :: l1,sml1
    Complex*16 :: Als
    type(sparse_complex) :: A1
    call setTriplet(l1,l1+sml1,Als,A1)
end subroutine set

subroutine setTriplet(icount,jcount,value,A)
    use Types
    implicit none
    type(sparse_complex) :: A
    integer*4 :: icount,jcount
    Complex*16 :: value
    !type(sparse_complex) :: A(:)

    A.nzmax=A.nzmax+1
    A.row(A.nzmax)=icount
    A.col(A.nzmax)=jcount
    A.val(A.nzmax)=value

end subroutine setTriplet

subroutine findnzmax(A)
    use types
    implicit none
    type(sparse_complex) :: A
    integer*4 :: icount2
    Do icount2=1,A.n**2
        if(icount2>size(A.col))then
            A.nzmax=icount2-1
            EXIT
        elseif(A.col(icount2)==0)then
            A.nzmax=icount2-1
            EXIT
        end if
    End Do
    continue
    !write(*,*)A.nzmax;pause!TEST!!!FOR TEST ONLY
end subroutine findnzmax


recursive subroutine quicksort(A1, A2, A3, first, last)
    implicit none
    integer*4 ::  A1(*),A2(*), x1, x2, t1, t2
    Complex*16 :: A3(*),t3
    integer*4 :: first, last
    integer*4 :: i00, j00
    ! quicksort.f -*-f90-*-
    ! Author: t-nissie
    ! License: GPLv3
    ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
    !!
    x1 = A1( (first+last) / 2 )
    x2 = A2( (first+last) / 2 )
    i00 = first
    j00 = last
    do
        do while ((A1(i00) < x1) .or. ((A1(i00)==x1).and.(A2(i00)<x2)))
            i00=i00+1
        end do
        do while ((A1(j00) > x1) .or. ((A1(j00)==x1).and.(A2(j00)>x2)))
            j00=j00-1
        end do
        if (i00 >= j00) exit
        t1 = A1(i00);  A1(i00) = A1(j00);  A1(j00) = t1
        t2 = A2(i00);  A2(i00) = A2(j00);  A2(j00) = t2
        t3 = A3(i00);  A3(i00) = A3(j00);  A3(j00) = t3
        i00=i00+1
        j00=j00-1
    end do
    if (first < i00-1) call quicksort(A1,A2,A3, first, i00-1)
    if (j00+1 < last)  call quicksort(A1,A2,A3, j00+1, last)
end subroutine quicksort

subroutine complex_coo_pardiso_solve(A,x_plasma,B_rf)
    use types
    IMPLICIT NONE
    type(sparse_complex) :: A
    Complex*16 :: x_plasma(A.n),B_rf(A.n)
    type(sparse_complex) :: A_csr3
    integer*8 :: pt(64)
    integer*4 :: maxfct,mnum,mtype,phase,msglvl,nrhs,iparm(64),error
    integer*4,allocatable :: perm(:)
    allocate (A_csr3.row(1:A.n+1),A_csr3.col(1:A.nzmax),A_csr3.val(1:A.nzmax))
    allocate (perm(1:A.n))
    call findnzmax(A)
    A_csr3.n=A.n;A_csr3.nzmax=A.nzmax
    call coo_to_csr(A,A_csr3)
    !maxfct=1;mnum=1;mtype=13;phase=13;nrhs=1;msglvl=0
    !call pardisoinit(pt, mtype, iparm)
    maxfct=1;mnum=1;mtype=13;phase=13;nrhs=1;iparm(:)=0;msglvl=0;pt(:)=0
    iparm(3)=4 !OMP_NUM_THREADS
    call pardiso(pt,maxfct,mnum,mtype,phase,A_csr3.n,A_csr3.val,A_csr3.row,A_csr3.col,perm,nrhs,iparm,msglvl,B_rf,x_plasma,error)
    ! 释放 pardiso 内部内存(否则每次求解都泄漏因子内存，累积耗尽内存崩溃)
    phase=-1
    call pardiso(pt,maxfct,mnum,mtype,phase,A_csr3.n,A_csr3.val,A_csr3.row,A_csr3.col,perm,nrhs,iparm,msglvl,B_rf,x_plasma,error)
    deallocate(perm)
    deallocate(A_csr3.row,A_csr3.col,A_csr3.val)
end subroutine complex_coo_pardiso_solve

subroutine coo_to_csr(A,A_csr3)
    use types
    implicit none
    type(sparse_complex) :: A,A_csr3
    integer*4 :: col1,row1,rowcount,ncount,count1
    Complex*16 :: val1
    integer*4 :: rownum(1:A.n+1)
    rownum=0
    call quicksort(A.row,A.col,A.val,1,A.nzmax)
    DO rowcount=1,A.nzmax
        rownum(A.row(rowcount))=rownum(A.row(rowcount))+1
    End DO
    A_csr3.row(1)=1
    do rowcount=1,A.n
        A_csr3.row(rowcount+1)=A_csr3.row(rowcount)+rownum(rowcount)
    end do
    A_csr3.col(1:A.nzmax)=A.col(1:A.nzmax)
    A_csr3.val(1:A.nzmax)=A.val(1:A.nzmax)
end subroutine coo_to_csr

subroutine reallocate(A,nzmax_input)
    use Types
    implicit none
    type(sparse_complex) :: A
    integer :: nzmax_input
    integer*4,allocatable :: Arow(:),Acol(:)
    Complex*16,allocatable :: Aval(:)
    call findnzmax(A)
    allocate(Arow(1:A.nzmax),Acol(1:A.nzmax),Aval(1:A.nzmax))
    Arow(1:A.nzmax)=A.row(1:A.nzmax)
    Acol(1:A.nzmax)=A.col(1:A.nzmax)
    Aval(1:A.nzmax)=A.val(1:A.nzmax)
    Deallocate(A.row,A.col,A.val)
    Allocate(A.row(1:nzmax_input),A.col(1:nzmax_input),A.val(1:nzmax_input))
    A.row(1:A.nzmax)=Arow(1:A.nzmax)
    A.col(1:A.nzmax)=Acol(1:A.nzmax)
    A.val(1:A.nzmax)=Aval(1:A.nzmax)
    A.row(A.nzmax+1:nzmax_input)=0
    A.col(A.nzmax+1:nzmax_input)=0
    A.val(A.nzmax+1:nzmax_input)=0.0D0+0.0*(0.0D0,1.0D0)
    Deallocate(Acol,Arow,Aval)
    call findnzmax(A)
end subroutine reallocate

! subroutine find_fvm
    ! use the_whole_varibles
    ! implicit none
    ! real*8 :: tev,vte,ane,vei,ve,ven,ug
    ! ug=10.9D0; !ug Coulomb logarithm;vte is the electron thermal velocity
    ! fvm=0.0D0
    ! do iz=1,nz
    ! do ir=1,nrp !test!!!!!!!!!!!!!
    ! tev=te_in_FDFD(ir,iz)*qe_abs;
    ! if(te_in_FDFD(ir,iz)<fv_te_limit)tev=fv_te_limit*qe_abs;
    ! vte=sqrt(8.0D0*tev/(pi*me));ane=ni(ir,iz);
    ! ven=nn*sig_en*vte;
    ! vei=(tev/qe_abs)**(-1.5D0)*ug*(2.91D-12)*ane;
    ! ve=ven+vei;
    ! fvm(ir,iz)=ve;
    ! end do
    ! end do
    ! continue
! end subroutine find_fvm

subroutine find_ve_vi(ve,vi)
    use the_whole_varibles
    implicit none
    real*8 tev,vte,ni_tp,vei,ve,ven,ug
    real*8 ti_tp,vti,vie,vin,vi
    ug=10.9D0; !ug Coulomb logarithm;vte is the electron thermal velocity

    ni_tp=ni(ir,iz);
    tev=te_in_FDFD(ir,iz)*qe_abs;
    !if(te_in_FDFD(ir,iz)<fv_te_limit)tev=fv_te_limit*qe_abs;
    vte=sqrt(8.0D0*tev/(pi*me));
    ven=nn*sig_en*vte;
    vei=(tev/qe_abs)**(-1.5D0)*ug*(2.91D-12)*ni_tp;
    ve=ven+vei;

    ti_tp=ek_ion_2D(ir,iz)/1.5d0 !eV
    if(ti_tp<1.)ti_tp=1.0d0
    vti=sqrt(8.d0*qe_abs*ti_tp/(pi*mi));

    !!Huba, NRL PLASMA FORMULARY,2002, Page 28
    !vie=4.80d-14*mass_number**(-0.5)*ti_tp**(-1.5)*ug*ni_tp;

    !Ma Tengcai, et,al, principle of plasma physics, 2012 revision,Page 127,eq (3.2.106)
    vie=2.*me/mi*(te_in_FDFD(ir,iz)/ti_tp)**(1.5)*vei

    vin=nn*sig_in*vti;  !vi=500 m/s, 500K
    vi=vie+vin;
end subroutine find_ve_vi

subroutine totalsolve
    use the_whole_varibles
    implicit none
    integer*4 :: kcount3,ndirec1,ndirec2
    real*8 :: ave_p
    ave_p=0.5D0

    DO kcount3=1, nr*nz
        e_m(((kcount3-1)/nz+1),(MOD(kcount3-1,nz)+1),1)=eq_xlast(m,3*kcount3-2)
        e_m((kcount3-1)/nz+1,MOD(kcount3-1,nz)+1,2)=eq_xlast(m,3*kcount3-1)
        e_m((kcount3-1)/nz+1,MOD(kcount3-1,nz)+1,3)=eq_xlast(m,3*kcount3)
    End DO

    e_int(2:nr,1:nz,1)=(0.5D0+0.0D0*i)*(e_m(2:nr,1:nz,1)+e_m(1:nr-1,1:nz,1))
    e_int(1,1:nz,1)=e_m(1,1:nz,1)

    e_int(1:nr,1:nz,2)=e_m(1:nr,1:nz,2)

    e_int(1:nr,2:nz,3)=(0.5D0+0.0D0*i)*(e_m(1:nr,2:nz,3)+e_m(1:nr,1:nz-1,3))
    e_int(1:nr,1,3)=e_m(1:nr,1,3)

    call findb
    call findjp
    call ptotalsolve
    do ith=1,nth
        e_record(:,ith,:,:)=e_int(:,:,:)*exp(im*th(ith))
        b_record(:,ith,:,:)=b_m(:,:,:)*exp(im*th(ith))
        ja_record(:,ith,:,:)=ja_m_nor(:,:,:)*exp(im*th(ith))
        jp_record(:,ith,:,:)=jp_m(:,:,:)*exp(im*th(ith))
    enddo

    e_AC=e_AC+e_record
    b_AC=b_AC+b_record
    ja_AC=ja_AC+ja_record
    jp_AC=jp_AC+jp_record


    e_output(m,:,:,:)=e_int(:,:,:)
    max_eth(m)=maxval(abs(e_output(m,1:nrp,:,2)))

end subroutine totalsolve

subroutine findb
    use the_whole_varibles
    implicit none

    do ir=1,nr-1
        do iz=1,nz-1
            if(ir/=1)then
                b_m(ir,iz,1)=im/r(ir)*e_m(ir,iz,3)-1.0D0/dz*e_m(ir,iz+1,2)+1.0D0/dz*e_m(ir,iz,2)
            endif
            b_m(ir,iz,2)=1.0D0/dz*e_m(ir,iz+1,1)-1.0D0/dz*e_m(ir,iz,1)-1.0D0/dr*e_m(ir+1,iz,3)+1.0D0/dr*e_m(ir,iz,3)
            b_m(ir,iz,3)=1.0D0/(dr*(0.5D0*dr+r(ir)))*(r(ir+1)*e_m(ir+1,iz,2)-r(ir)*e_m(ir,iz,2))-im*e_m(ir,iz,1)/(r(ir)+0.5D0*dr)
        enddo
        if(m/=0)then
            b_m(1,:,1)=b_m(2,:,1)
        else
            b_m(1,:,1)=0.0D0+0.0D0*i
        endif
    enddo
    iz=nz
    do ir=1,nr-1
        b_m(ir,iz,1)=2.0D0*b_m(ir,iz-1,1)-b_m(ir,iz-2,1)
        b_m(ir,iz,2)=2.0D0*b_m(ir,iz-1,2)-b_m(ir,iz-2,2)
        b_m(ir,iz,3)=1.0D0/(dr*(0.5D0*dr+r(ir)))*(r(ir+1)*e_m(ir+1,iz,2)-r(ir)*e_m(ir,iz,2))-im*e_m(ir,iz,1)/(r(ir)+0.5D0*dr)
    enddo
    ir=nr
    do iz=1,nz-1
        b_m(ir,iz,1)=im/r(ir)*e_m(ir,iz,3)-1.0D0/dz*e_m(ir,iz+1,2)+1.0D0/dz*e_m(ir,iz,2)
        b_m(ir,iz,2)=2.0D0*b_m(ir-1,iz,2)-b_m(ir-2,iz,2)
        b_m(ir,iz,3)=2.0D0*b_m(ir-1,iz,3)-b_m(ir-2,iz,3)
    enddo
    ir=nr;iz=nz
    b_m(ir,iz,:)=-b_m(ir-1,iz-1,:)+b_m(ir,iz-1,:)+b_m(ir-1,iz,:)

    b_m(:,:,:)=b_m(:,:,:)/(i*om)
end subroutine findb

subroutine findjp
    use the_whole_varibles
    implicit none
    integer*4 :: ndirec1,ndirec2
    jp_m=0.0D0
    Do ndirec1=1,3
        Do ndirec2=1,3
            jp_m(:,:,ndirec1)=jp_m(:,:,ndirec1)+si(:,:,ndirec1,ndirec2)*e_int(:,:,ndirec2)
        EndDo
    EndDo
end subroutine findjp

subroutine ptotalsolve
    use the_whole_varibles
    implicit none
    integer*4 :: ndirec1,ndirec2,n_boundary
    real*8 :: ave_p,ptotal_pt
    complex*16 :: power_mode_2D(nr,nz)     ! 本 m 这一份的 2D 沉积，只用于 ptotal 积分
    n_boundary=nint(0.8*nrp)
    ave_p=0.5D0
    power_depo_2D(:,:)=0.0D0
    power_mode_2D(:,:)=0.0D0
    do ndirec1=1,3
        do ndirec2=1,3
            power_depo_rthz(:,:,ndirec2)=power_depo_rthz(:,:,ndirec2)+ave_p*2.0D0*pi*(si(:,:,ndirec2,ndirec1)*e_int(:,:,ndirec1)*DCONJG(e_int(:,:,ndirec2)))
            power_mode_2D(:,:)=power_mode_2D(:,:)+ave_p*2.0D0*pi*(si(:,:,ndirec2,ndirec1)*e_int(:,:,ndirec1)*DCONJG(e_int(:,:,ndirec2)))
        end do
    end do

    power_depo_2D(:,:)=power_depo_rthz(:,:,1)+power_depo_rthz(:,:,2)+power_depo_rthz(:,:,3)


    ptotal_pt=ptotal
    do ir=1,nrp-1
        do iz=1,nz-1
            ptotal=ptotal+ds/4.0D0*real(power_mode_2D(ir,iz)*r(ir)+power_mode_2D(ir+1,iz)*r(ir+1)&
                &+power_mode_2D(ir,iz+1)*r(ir)+power_mode_2D(ir+1,iz+1)*r(ir+1))
        end do
    end do
    do ir=1,n_boundary-1
        do iz=1,nz-1
            ptotal_inside=ptotal_inside+ds/4.0D0*real(power_mode_2D(ir,iz)*r(ir)+power_mode_2D(ir+1,iz)*r(ir+1)&
                &+power_mode_2D(ir,iz+1)*r(ir)+power_mode_2D(ir+1,iz+1)*r(ir+1))
        end do
    end do
    do ir=n_boundary,nrp-1
        do iz=1,nz-1
            ptotal_boundary=ptotal_boundary+ds/4.0D0*real(power_mode_2D(ir,iz)*r(ir)+power_mode_2D(ir+1,iz)*r(ir+1)&
                &+power_mode_2D(ir,iz+1)*r(ir)+power_mode_2D(ir+1,iz+1)*r(ir+1))
        end do
    end do
    ptotal_m(m)=ptotal-ptotal_pt

    !do ir=1,nr
    !ptotal_BI=ptotal_BI+alpha_BI*(abs(B_m(ir,1,1))**2+abs(B_m(ir,1,2))**2&
    !    &+abs(B_m(ir,nz,1))**2+abs(B_m(ir,nz,2))**2)*r(ir)
    !enddo
end subroutine ptotalsolve


subroutine Find_Region(nr,nz,n_vac,nr_vac,nr_met,nz_vac,FindRegion)
    implicit none
    integer*4 :: nr,nz,ir,iz,nn,n_vac
    integer*4 :: nr_vac(n_vac),nz_vac(n_vac),nr_met(n_vac)!,nrz_diploe(n_vac),iswtich_inner_dipole
    integer*4 :: FindRegion(nr,nz),F1(nr,0:nz+1)
    integer*4 :: nz_vac1(0:n_vac)
    nz_vac1(0)=1;nz_vac1(1:n_vac)=nz_vac(:);

    FindRegion(:,:)=0;  F1(:,:)=0
    F1(:,0)=11; F1(:,nz+1)=11


    do nn=1,n_vac
        do iz=nz_vac1(nn-1),nz_vac1(nn)
            F1(1,iz)=10
            F1(2:nr_vac(nn)-1,iz)=1
            F1(nr_vac(nn),iz)=9
            F1(nr_vac(nn)+1:nr_met(nn)-1,iz)=2
            F1(nr_met(nn),iz)=6
            F1(nr_met(nn)+1:nr,iz)=11
            if(nr_vac(nn)==1)F1(1,iz)=10
            if(nr_met(nn)==nr)F1(nr,iz)=6
            if(nr_vac(nn)==nr_met(nn))F1(nr_vac(nn),iz)=3
        enddo
    enddo

    do nn=0,n_vac
        iz=nz_vac1(nn)
        do ir=2,nr

            if(F1(ir,iz-1)==1)then
                if(F1(ir,iz+1)==1)F1(ir,iz)=1
                if(F1(ir,iz+1)==9)F1(ir,iz)=7
                if(F1(ir,iz+1)==2)F1(ir,iz)=7
                if(F1(ir,iz+1)==6)F1(ir,iz)=6
                if(F1(ir,iz+1)==11)F1(ir,iz)=5
                if(F1(ir,iz+1)==3)F1(ir,iz)=13
            elseif(F1(ir,iz-1)==9)then
                if(F1(ir,iz+1)/=9)F1(ir,iz-1)=14!!!!!!TEST: Can it be deleted?
                if(F1(ir,iz+1)==1)F1(ir,iz)=4
                if(F1(ir,iz+1)==9)F1(ir,iz)=9
                if(F1(ir,iz+1)==2)F1(ir,iz)=7
                if(F1(ir,iz+1)==6)F1(ir,iz)=6
                if(F1(ir,iz+1)==11)F1(ir,iz)=11
                if(F1(ir,iz+1)==3)F1(ir,iz)=13
            elseif(F1(ir,iz-1)==2)then
                if(F1(ir,iz+1)/=2)F1(ir,iz-1)=14!!!!!!TEST: Can it be deleted?
                if(F1(ir,iz+1)==1)F1(ir,iz)=4
                if(F1(ir,iz+1)==9)F1(ir,iz)=7
                if(F1(ir,iz+1)==2)F1(ir,iz)=2
                if(F1(ir,iz+1)==6)F1(ir,iz)=6
                if(F1(ir,iz+1)==11)F1(ir,iz)=11
                if(F1(ir,iz+1)==3)F1(ir,iz)=13
            elseif(F1(ir,iz-1)==6)then
                if(F1(ir,iz+1)==1)F1(ir,iz)=4
                if(F1(ir,iz+1)==9)F1(ir,iz)=4
                if(F1(ir,iz+1)==2)F1(ir,iz)=7
                if(F1(ir,iz+1)==6)F1(ir,iz)=6
                if(F1(ir,iz+1)==11)F1(ir,iz)=11
                if(F1(ir,iz+1)==3)F1(ir,iz)=3
            elseif(F1(ir,iz-1)==11)then
                if(F1(ir,iz+1)==1)F1(ir,iz)=4
                if(F1(ir,iz+1)==9)F1(ir,iz)=4
                if(F1(ir,iz+1)==2)F1(ir,iz)=7
                if(F1(ir,iz+1)==6)F1(ir,iz)=6
                if(F1(ir,iz+1)==11)F1(ir,iz)=11
                if(F1(ir,iz+1)==3)F1(ir,iz)=3
            elseif(F1(ir,iz-1)==3)then
                if(F1(ir,iz+1)==1)F1(ir,iz)=12
                if(F1(ir,iz+1)==9)F1(ir,iz)=4
                if(F1(ir,iz+1)==2)F1(ir,iz)=12
                if(F1(ir,iz+1)==6)F1(ir,iz)=6
                if(F1(ir,iz+1)==11)F1(ir,iz)=3
                if(F1(ir,iz+1)==3)F1(ir,iz)=3
            endif
        enddo
    enddo

    FindRegion(:,:)=F1(1:nr,1:nz)
    if(minval(FindRegion(:,:))==0)then
        write(*,*)'FindRegion ERROR, Please Check'
        pause
    endif
end subroutine Find_Region

subroutine find_kfc(kap_ll,tpp,vl,Ak,Bk,omci,omce,ompi2,ompe2,func)
    use the_whole_varibles
    implicit none
    integer*4::ia,ib
    real*8::vl,zu,zv,zup,zvp
    real*8::n_l(1:3,1:2),vpp(1:2,1:2),tpp(1:2,1:2)! parallel or perp
    complex*16::xomci,xomce,xompi2,xompe2
    complex*16::omci,omce,omc,ompi2,ompe2,omp2
    complex*16::kap_ll,Z_tmp,Zp_tmp,func
    complex*16::zeta(1:3,1:2),Ak(1:3,1:2),Bk(1:2) !(n,l),n_B=0


    ! n_l(a,b)  a=1,2,3 for \pm; b=1,2 for i,e
    n_l(1,1)=-1.; n_l(2,1)=0.; n_l(3,1)=1.;
    n_l(1,2)=1.; n_l(2,2)=0.; n_l(3,2)=-1.;

    ! vpp(a,b)  a=1,2 for ll,perp; b=1,2 for i,e
    vpp(1:2,1)=sqrt(2*qe_abs*tpp(1:2,1)/mi);
    vpp(1:2,2)=sqrt(2*qe_abs*tpp(1:2,2)/me);

    ! zeta(a,b) a=1,2,3 for n=-1,0,1; b=1,2 for i,e
    zeta(1:3,1)=(om-kap_ll*vl-n_l(1:3,1)*omci)/(kap_ll*vpp(1,1))
    zeta(1:3,2)=(om-kap_ll*vl-n_l(1:3,2)*omce)/(kap_ll*vpp(1,2))

    ! A(a,b)    a=1,2,3 for n=-1,0,1; b=1,2 for i,e
    ! B(b)      n=0; b=1,2 for i,e
    do ib=1,2
        do ia=1,3
            if (REAL(kap_ll)>0) then
                call zetaJpole(REAL(zeta(ia,ib)),IMAG(zeta(ia,ib)),zu,zv,zup,zvp)
                Z_tmp=CMPLX(zu,zv)
            else
                call zetaJpole(-REAL(zeta(ia,ib)),-IMAG(zeta(ia,ib)),zu,zv,zup,zvp)
                Z_tmp=-CMPLX(zu,zv)
            endif

            if (ib==1) then
                omc = omci
            else
                omc = omce
            endif
            Ak(ia,ib) = (tpp(2,ib)-tpp(1,ib))/(om*tpp(1,ib))+( zeta(ia,ib)*tpp(2,ib)/(om*tpp(1,ib))+n_l(ia,ib)*omc/(om*kap_ll*vpp(1,ib)) )*Z_tmp
            if (ia==2) then
                Bk(ib) = (om*tpp(2,ib)-kap_ll*vl*tpp(1,ib))/(om*kap_ll*tpp(1,ib))+zeta(ia,ib)*tpp(2,ib)/(kap_ll*tpp(1,ib))*Z_tmp
            endif
        enddo
    enddo

    if (iswitch_dielectric== 3) then
        func = om*om - kap_ll**2 * c**2 + ompi2*om*Ak(3,1) - ompe2*om/(om-kap_ll*vl-n_l(3,2)*omce)
    elseif (iswitch_dielectric== 4) then
        func = om*om - kap_ll**2 * c**2 + ompi2*om*Ak(3,1) + ompe2*om*Ak(3,2)
    endif
    func = func/(om*om)

end subroutine find_kfc

subroutine zetaJpole(xi,yi,zu,zv,zup,zvp)
    implicit none
    integer :: zj
    real*8,parameter :: sqpi = 1.77245D0 !=sqrt(\pi)
    real*8::xi,yi,zu,zv,zup,zvp
    complex*16 :: bzj(8),czj(8),ztmp,zi

    bzj(1) = (-0.017340112270401, -0.046306439626294)
    bzj(2) = (-0.739917811220052, 0.839518284620274)
    bzj(3) = (5.840632105105495, 0.953602751322040)
    bzj(4) = (-5.583374181615043, -11.208550459628098)
    czj(1) = (2.237687725134293, -1.625941024120362)
    czj(2) = (1.465234091939142, -1.789620299603315)
    czj(3) = (0.839253966367922, -1.891995211531426)
    czj(4) = (0.273936218055381, -1.941787037576095)

    bzj(5:8) = conjg(bzj(1:4))
    czj(5:8) = -conjg(czj(1:4))
    ztmp = (0,0); zi = CMPLX(xi,yi)

    if (yi>=0) then
        do zj = 1, 8
            ztmp = ztmp + bzj(zj) / (zi - czj(zj))
        end do
    else
        do zj = 1, 8
            ztmp = ztmp + CONJG(bzj(zj) / (CONJG(zi) - czj(zj)))
        end do
        ztmp = ztmp + 2.0*sqpi * CMPLX(sin(2.0*xi*yi), cos(2.0*xi*yi)) * exp(-xi**2+yi**2)
    endif

    zu = real(ztmp)
    zv = imag(ztmp)
    zup = -2*(1+xi*zu-yi*zv)
    zvp = -2*(xi*zv+yi*zu)

end subroutine zetaJpole

subroutine find_te
    USE the_whole_varibles
    implicit none
    real*8::dt_te0,x_k1(1:nrp,1:nz),x_k2(1:nrp,1:nz),x_k3(1:nrp,1:nz),x_k4(1:nrp,1:nz)
    real*8:: x_k(1:nrp,1:nz),x0(1:nrp,1:nz),x_k_out(1:nrp,1:nz),pe(1:nrp,1:nz),dt_te,ni_bd
    integer irun
    call find_func_cputime_1_of_2

    ni(1:nrp,1:nz)=n_macro*(density_2D(1:nrp,1:nz)+1e-2*maxval(density_2D))
    call find_dTe(dt_te0)

    nrun=nint(dt/dt_te0)
    dt_te=dt/nrun
    do irun=1,nrun

        x0=ni*te_mhd;
        x_k=x0;             call te_RK4(x_k,x_k_out);  x_k1=x_k_out;
        pe(1:nrp,1:nz)=x0+dt_te*x_k1

        te_mhd(1:nrp,1:nz)=pe(1:nrp,1:nz)/ni(1:nrp,1:nz)
        te_mhd(1,:)=te_mhd(2,:)
        do iz=1,nz
            do ir=1,nrp
                if(te_mhd(ir,iz)<1.)te_mhd(ir,iz)=1.
                if(te_mhd(ir,iz)>1e4)te_mhd(ir,iz)=1e4
            end do
        enddo
    enddo


    te_2D=te_mhd
    call find_func_cputime_2_of_2(func_time(8))
    continue
end subroutine find_te

subroutine te_RK4(x_k,x_k_out)
    use the_whole_varibles
    implicit none
    real*8::x_k(1:nrp,1:nz),x_k_out(1:nrp,1:nz),te_k(1:nrp,1:nz),divqe(1:nrp,1:nz)
    real*8 :: dter,dtez,term1,term2,dner,dnez
    integer::nb_type
    real*8::ne_bd,te_bd
    integer::ib_type,ir1,ir2,iz1,iz2
    real*8::div_val,r_tp,dr_tp,dz_tp
    real*8:: nu_0(1:nrp,1:nz),ur_tp,uz_tp
    real*8::para_in(1:5),flux_qe(0:nrp,0:nz,1:3),flux_je(0:nrp,0:nz,1:3)

    divqe=0.;
    flux_je=0.;
    flux_qe=0.;
    te_k=abs(x_k/ni)

    do iz=1,nz
        do ir=1,nrp
            if(te_k(ir,iz)<1.)te_k(ir,iz)=1.
            if(te_k(ir,iz)>1e4)te_k(ir,iz)=1e4
        end do
    enddo

    !-----------------------------jr-----------------------------------------------!
    do iz=2,nz-1
        do ir=1,nrp-1
            para_in(1)=0.5*(b0_DC(ir,iz,1)+b0_DC(ir+1,iz,1))
            para_in(2)=0.5*(b0_DC(ir,iz,3)+b0_DC(ir+1,iz,3))
            para_in(3)=0.5*(te_k(ir,iz)+te_k(ir+1,iz))
            para_in(4)=0.5*(ni(ir,iz)+ni(ir+1,iz))
            para_in(5)=0.5*(ek_ion_2D(ir,iz)+ek_ion_2D(ir+1,iz))/1.5
            call MHD_coeff(para_in)        !calculate the coefficient
            !density flux and energy flux
            dner =(ni(ir+1,iz)-ni(ir,iz))/dr
            dter =(te_k(ir+1,iz)-te_k(ir,iz))/dr
            tp1 = (ni(ir,iz)+ni(ir+1,iz))/2.  !the corresponding ni on half grids
            tp2=ni(ir,iz-1)
            tp3=ni(ir,iz+1)
            dnez=(tp3-tp2)/d2z

            term1=-dner*Da_perp*coeff_1_br2        !--u by MHD----------!
            term2=-dnez*Da_perp*mue_i_brz
            flux_je(ir,iz,1)= (term1+term2)  !
            u_mhd(ir,iz,1)= (term1+term2)/tp1
            ur_tp=u_mhd(ir,iz,1)                   !--u by MHD----------!
            !ur_tp=0.5*(u_pic(ir,iz,1)+u_pic(ir+1,iz,1)) !--u by PIC----------!

            term1=-1.5*de*coeff_1_mue2_br2/coeff_1_mue2_b2*tp1*dter
            term2=2.5*ur_tp*tp1*(.5*(te_k(ir+1, iz)+te_k( ir,iz)))
            flux_qe(ir,iz,1)= term1+term2
        end do
    end do

    !-----------------------------jz----------------------------------------!
    do iz=1,nz-1
        do ir=2,nrp-1
            para_in(1)=0.5*(b0_DC(ir,iz,1)+b0_DC(ir,iz+1,1))
            para_in(2)=0.5*(b0_DC(ir,iz,3)+b0_DC(ir,iz+1,3))
            para_in(3)=0.5*(te_k(ir,iz)+te_k(ir,iz+1))
            para_in(4)=0.5*(ni(ir,iz)+ni(ir,iz+1))
            para_in(5)=0.5*(ek_ion_2D(ir,iz)+ek_ion_2D(ir,iz+1))/1.5
            call MHD_coeff(para_in)
            dnez=(ni(ir,iz+1)-ni(ir,iz))/dz
            dtez=(te_k(ir,iz+1)-te_k(ir,iz))/dz
            tp1 = (ni(ir,iz)+ni(ir,iz+1))/2.  !the corresponding ni on half grids
            tp2=ni(ir-1,iz)
            tp3=ni(ir+1,iz)
            dner=(tp3-tp2)/d2r
            !
            term1=-dner*Da_perp*mue_i_brz        !--u by MHD----------!
            term2=-dnez*Da_perp*coeff_1_bz2
            flux_je(ir,iz,3)= (term1+term2)
            u_mhd(ir,iz,2)= (term1+term2)/tp1
            uz_tp=u_mhd(ir,iz,2)                 !--u by MHD----------!
            !uz_tp=0.5*(u_pic(ir,iz,2)+u_pic(ir,iz+1,2)) !--u by PIC----------!

            term1=-1.5*de*coeff_1_mue2_bz2/coeff_1_mue2_b2*tp1*dtez
            term2=2.5*uz_tp*tp1*(.5*(te_k(ir,iz+1)+te_k(ir,iz)))
            flux_qe(ir,iz,3)= term1+term2
        end do
    end do

    !----------------jr jz boundary------------------!
    do ir=1,nrp
        do iz=1,nz
            ib_type=i_plasma_region(ir,iz)!-10
            if((ib_type>1 .and. ib_type<=4) .or. ib_type>=6 )call MHD_flux_boundary(ib_type,te_k,flux_qe) !,x_k,flux_je
        enddo
    enddo
    !!!--------------------find divJe , divJi , divqe start------------------------------!!
    do ir=1,nrp
        do iz=1,nz
            ib_type=i_plasma_region(ir,iz) !-10
            if(ib_type>=1)then
                call MHD_div(ib_type,r_tp,dr_tp,dz_tp)! .and. ib_type<=5
                tp1=   (flux_qe(ir,iz,1)-flux_qe(ir-1,iz,1))/dr_tp
                tp2=.5*(flux_qe(ir,iz,1)+flux_qe(ir-1,iz,1))/r_tp
                tp3=   (flux_qe(ir,iz,3)-flux_qe(ir,iz-1,3))/dz_tp
                divqe(ir,iz)=tp1+tp2+tp3
            endif
        enddo
    enddo


    nu_0(1:nrp,1:nz)=ni(1:nrp,1:nz)*coeff_Te_1/(qe_abs*te_k(1:nrp,1:nz))**1.5;
    Q_ie(1:nrp,1:nz)=(nu_0(1:nrp,1:nz)*coeff_Te_3*(1+coeff_Te_2*Ek_ion_2D(1:nrp,1:nz)/te_k(1:nrp,1:nz))**(-1.5))* &
        & (Ek_ion_2D(1:nrp,1:nz)/1.5-te_k(1:nrp,1:nz))
    x_k_out(1:nrp,1:nz)=2/3.*(-divqe(1:nrp,1:nz)+ni(1:nrp,1:nz)*Q_ie(1:nrp,1:nz)+real(power_depo_2D(1:nrp,1:nz))/qe_abs);
end subroutine te_RK4

subroutine MHD_flux_boundary(nb_type,te_k,flux_qe)
    use the_whole_varibles
    implicit none
    integer::irb,izb,nb_type
    real*8::te_k(1:nrp,1:nz),ur_tp,uz_tp,vte
    real*8::flux_qe(0:nrp,0:nz,1:3)
    !-----------nb_type----------------!
    !11,23: left boundary   |<-plasma
    !12,24: right boundary     plasma->|
    !13: up  boundary     -^plasma
    !14: bottom boundary  _vplasma
    !15: r=0 boundary
    !16: r and z shared boundary

    vte=coeff_Te_4*sqrt(te_k(ir,iz))

    ur_tp=vte
    uz_tp=vte

    !-----------nb_type=11----------------!
    if(nb_type==11 .or. nb_type==23 .or. nb_type==27  )then
        u_mhd(ir,iz-1,2)= -uz_tp
        flux_qe(ir,iz-1,3)=   -uz_tp*ni(ir,iz)*te_k(ir,iz)  !t2.0  ->2.5  or 5

        !-----------nb_type=2----------------!
    elseif(nb_type==12 .or. nb_type==24 .or. nb_type==28   )then
        u_mhd(ir,iz,2)= uz_tp
        flux_qe(ir,iz,3)=   uz_tp*ni(ir,iz)*te_k(ir,iz)  ! 2.0  ->2.5  or 5

        !-----------nb_type=3----------------!
    elseif(nb_type==13)then
        u_mhd(ir,iz,1)= ur_tp
        flux_qe(ir,iz,1)=   ur_tp*ni(ir,iz)*te_k(ir,iz)  !2.0  ->2.5  or 5

        !-----------nb_type=4----------------!
    elseif(nb_type==14)then
        u_mhd(ir-1,iz,1)= -ur_tp
        flux_qe(ir-1,iz,1)=   -ur_tp*ni(ir,iz)*te_k(ir,iz)  !2.0  ->2.5  or 5
    endif
end subroutine MHD_flux_boundary

subroutine MHD_div(nb_type,r_tp,dr_tp,dz_tp)
    use the_whole_varibles
    implicit none
    integer::nb_type !irb,izb,
    real*8::r_tp,dr_tp,dz_tp !,i1,i2,i4div_val,
    !-----------nb_type----------------!
    !1: left boundary   |<-plasma
    !2: right boundary     plasma->|
    !3: up  boundary     -^plasma
    !4: bottom boundary  _vplasma
    !5: r=0 boundary
    !6: r and z shared boundary

    dr_tp=dr;
    dz_tp=dz;
    r_tp=r(ir)

    !-----------nb_type=1 and nb_type=2----------------!
    if(nb_type==11 .or. nb_type==12 )then ! .or. nb_type==23 .or. nb_type==24
        dz_tp=dz_2

        !-----------nb_type=5----------------!
    elseif(nb_type==15 )then
        dr_tp=dr_2;
        r_tp=dr*0.25

        !-----------nb_type=3----------------!
    elseif(nb_type==13)then !.or. nb_type==23 .or. nb_type==24
        dr_tp=dr_2;
        r_tp=r(ir-1)+dr*0.75!r2(ir)

    elseif(nb_type==14)then
        dr_tp=dr_2;
        r_tp=r(ir)+dr*0.25!r2(ir)

    elseif(nb_type==27 .or.nb_type==28  )then
        dz_tp=dz_2
        dr_tp=dr_2;
        r_tp=r(ir)+dr*0.25!r2(ir)

    elseif( nb_type==25 .or. nb_type==26)then
        dr_tp=dr_2;
        dz_tp=dz_2
        r_tp=dr*0.25

    elseif( nb_type==23 .or. nb_type==24)then
        dr_tp=dr_2;
        dz_tp=dz_2
        r_tp=r(ir-1)+dr*0.75
    endif
end subroutine MHD_div

subroutine MHD_coeff(para_in)
    use the_whole_varibles
    implicit none
    real*8::para_in(1:5)
    real*8::br00,bz00,te00,ne00,bt2,vti,vi,vie,vin,vei,ven,ve,vte
    real*8 :: coeff_1_b2,mue_mui,mue2,di,mui,mue
    real*8:: te_tp,ni_tp,ti_tp
    br00=para_in(1)
    bz00=para_in(2)
    te_tp=para_in(3)
    ni_tp=para_in(4)
    ti_tp=para_in(5)

    bt2=br00**2+bz00**2+1e-10  !B_total^2
    vte=sqrt(8.*te_tp*qe_abs/(pi*me));
    ven=nn*sig_en*vte;
    vei=(te_tp)**(-1.5)*ln_A0*(2.91e-12)*ni_tp;
    ve=ven+vei+1e8; !to avoid De too low
    mue=qe_abs/(me*ve);
    de=te_tp*mue;

    ti_tp=ek_ion_2D(ir,iz)/1.5 !eV
    if(ti_tp<1.)ti_tp=1.
    vti=sqrt(8.*qe_abs*ti_tp/(pi*mi));
    vie=4.80e-14*mass_number**(-0.5)*ti_tp**(-1.5)*ln_A0*ni_tp;
    vin=nn*sig_in*vti;
    vi=vie+vin;
    mui=qe_abs/(mi*vi); !mui
    di=mui*ti_tp    !Di=qe_abs*Ti/(mi*vi);

    mue2=mue*mue
    coeff_1_mue2_br2=1+mue2*br00**2
    coeff_1_mue2_bz2=1+mue2*bz00**2
    coeff_1_mue2_b2=1+mue2*bt2

    mue_mui=mue*mui
    mue_i_brz=mue_mui*br00*bz00
    coeff_1_br2=1+mue_mui*br00**2
    coeff_1_bz2=1+mue_mui*bz00**2
    coeff_1_b2=1+mue_mui*bt2
    Da=(di*mue+de*mui)/(mui+mue) !Da ll
    Da_perp=Da/coeff_1_b2        !Da perp
end subroutine MHD_coeff


subroutine find_dTe(dt_te0)
    use the_whole_varibles
    implicit none
    real*8::para_in(1:5),dTe_2D(nrp,nz),dt_te0
    integer Nte_max
    Nte_max=300
    dTe_2D=sqrt(-1.)
    do iz=1,nz
        do ir=1,nrp
            if(i_plasma_region(ir,iz)>0.5)then
                para_in(1)=b0_DC(ir,iz,1)
                para_in(2)=b0_DC(ir,iz,3)
                para_in(3)=te_mhd(ir,iz)
                para_in(4)=ni(ir,iz)
                para_in(5)=ek_ion_2D(ir,iz)
                call MHD_coeff(para_in)
                dTe_2D(ir,iz)=3.*dr**2/(5.*De)
            endif
        enddo
    enddo

    dt_te0=0.5*minval(dTe_2D)    
    if(dt_te0>dt)dt_te0=dt
end subroutine find_dTe
