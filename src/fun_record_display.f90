subroutine display_main
    USE the_whole_varibles
    implicit none
    integer,parameter::n6=6
    character(len=*),parameter::esc=achar(27)
    logical,save::first_log_display=.true.
    real*8::r_max(1:n6),tau_eq(1:nr,1:nz),nu_0,veq,loss_power
    call find_func_cputime_1_of_2  !-------------------1/2

    loss_power=power_loss_ave(3)
    !every 10 s display
    call end_ftime
    run_time_1=time_run*60.  !s
    if( (run_time_1-run_time_0)>10  .and. (iswitch_display/=0))then
        run_time_0=run_time_1

        !if((mod(t,t_show)<dt .or. it==1) .and. (iswitch_display/=0))then

        Ek_i_ave=0.
        do ip=1,n_active
            Ek_i_ave=Ek_i_ave+mass_q_i_05*sum(v(ip,1:3)**2); !ev
        enddo
        Ek_e_ave=mass_q_e_05*sum(v_e(1:n_active,1:3)**2); !ev

        Ek_i_ave=Ek_i_ave/real(n_active) !ev
        Ek_e_ave=Ek_e_ave/real(n_active) !ev
        ti_ave=Ek_i_ave/1.5
        te_ave=Ek_e_ave/1.5

        tau_eq=sqrt(-1.)
        do ir=1,nr
            do iz=1,nz
                if(i_plasma_region(ir,iz)>0.5)then
                    nu_0=ni(ir,iz)*coeff_Te_1/(qe_abs*te_mhd(ir,iz))**1.5;
                    veq=nu_0*coeff_Te_3*(1+coeff_Te_2*Ek_ion_2D(ir,iz)/te_mhd(ir,iz))**(-1.5);
                    tau_eq(ir,iz)=1/veq
                endif
            enddo
        enddo

        tp1=t/td;

        if(iswitch_log == 0) then
                    
            write(*,*)'------------------------------------------------------------------------'
            write(*, "('--------- The start time is ', i4,'.',i2,'.',i2,' ',i2.2,':',i2.2,'-------')") date_time(1:3,1),date_time(5:6,1)
            write(*, "('--------- The time now is  ', i4,'.',i2,'.',i2,' ',i2.2,':',i2.2,'-------')") date_time(1:3,2),date_time(5:6,2)
            write(*, "('--- The calculation time is ', f8.2,'  min (',f8.2,' , hour)----')") time_run,time_run/60.
            write(*, "(' Program is running at', f12.2,3x,'%   (estimating',f12.2,' hours left )')") tp1*100, time_run/60./tp1*(1-tp1)
            write(*, "(' total records, already records', 3i5)") n_pic,ifig-1
            write(*,*)
            write(*, "(' t=', f7.1,'  trf   (',es9.2,' s)')") t/trf,t
            write(*, "(' nr, nz=', 2i5)") nr,nz
            if(iswitch_antenna_type==0) write(*,"('antenna type (',i2, ' ) : single loop ')")iswitch_antenna_type
            if(iswitch_antenna_type==-1)write(*,"('antenna type : left helical   (',i3, ' )')")iswitch_antenna_type
            if(iswitch_antenna_type==1) write(*,"('antenna type : right helical   (',i3, ' )')")iswitch_antenna_type
            if(iswitch_antenna_type==2) write(*,"('antenna type : half loop  (',i3, ' )')")iswitch_antenna_type
            if(iswitch_antenna_type==3) write(*,"('antenna type : Nagoya type-III  (',i3, ' )')")iswitch_antenna_type
            write(*, "('azimuthal mode m: [ ',i3,', ',i3,' ]')") m_start,m_end
            !write(*, "('mass_number=', f7.1,', f=',f7.2,' MHz')") mass_number ,frequency/1e6
            write(*, "('center B0 =', f7.0,', B_resonance=',f7.1,' G')") 1e4*bz_max, 1e4*B_resonance
            write(*, "('RF: ,f=',f7.2,' MHz')")frequency/1e6
            write(*, "('irf, resistance = ',f7.1,' A, ',es10.2,' ohm')")i_now,resistance
            if(iswitch_RF2 ==  1)then
                write(*,*)'The RF2 antenna is on:'
                write(*, "('RF2: ,f=',f7.2,' MHz')")frequency2/1e6
                write(*, "('irf, resistance = ',f7.1,' A, ',es10.2,' ohm')")irf2_set,resistance
            endif
            write(*, "('total power, power_loss, power_Joule= ', 3es10.2,' kW')") &
                &       (loss_power+power_Joule)/1e3,loss_power/1e3,power_Joule/1e3
            !统计能量吸收效率
            write(*, "('total absorbed power =', es10.2,' kW')") absorbed_power/1e3
            write(*, "('plasma absorped power efficiency =', f7.1,' %')") 100*absorbed_power/(absorbed_power+power_Joule+1e-10)

            write(*,*)
            write(*, "('triple product(n*T_i*Tau_E) =', es10.2,' m^-3*keV*s')")rec_t(40)
            write(*, "('Tau_E =', es10.2,' s,  ','Tau_p =', es10.2,' s')") tau_E_ave,tau_p_ave
            write(*,*) '   trf,     dt,    dt_mcc,  t_power_on, td'
            write(*, "(5es9.1, ' s')")   trf,   dt,   minval(dt_mcc),   t_power_on,   td
            write(*,*) '  i-i,     e-e,     i-e,     e-i   dt_mcc'
            write(*, "(4es9.1, ' s')")   dt_mcc
            write(*,*) 'tau_ii,   tau_ee,  tau_ie,   tau_ei, tau_eq (e-i thermal equilibrium time)'
            write(*, "(5es9.1,'  ~',es9.1, ' s')")tau_ii,tau_ee,tau_ie,tau_ei,minval(tau_eq),maxval(tau_eq)
            write(*, "('ln_A =', f10.2)") ln_A
            write(*,*)

            write(*, "('T_i =', es10.2,' eV,  ','T_e =', es10.2,' eV')") ti_ave,te_ave
            write(*, "('<Ek_ion> =', es10.2,' eV,  ','<Ek_ele> =', es10.2,' eV')") Ek_i_ave, Ek_e_ave
            write(*, "('max te_in_FDFD =', es10.2,' eV,  ','min T_e =', es10.2,' eV')") maxval(te_mhd),minval(te_mhd)
            write(*, "('average kinetic energy of emitted ions =', es10.2,' eV')") Ek_loss_ave(4)
            !write(*, "('propulsive efficiency =', f7.1,' %')") 100*(power_loss_ave(2)+power_loss_ave(4))/(sum(power_loss_ave(1:4))+1e-10)
            write(*,*)

            write(*, "(' n_active=',i7,', n_capacity=',i7,', deltaN=',i7,', n_macro=',es10.2)") &
                & n_active,n_capacity,deltaN,n_macro
            write(*, "(' ni_max=',es9.2,', ni_min=',es9.2,'  m-3')")maxval(density_2D)*n_macro,minval(density_2D)*n_macro
            write(*, "(' inject, loss number', 2i5)") num_inject,N_lost_particles
            !write(*, "('function of mcc  has been called  =', f10.1,' times')")run_mcc_times
            write(*, "('function of FDFD has been called  =', f10.1,' times')")real(maxwellcount)
            write(*, "('dt/dTe  =', i8)")nrun
            write(*,*)

            tp1=sum(func_time(1:9))/100.;
            !write(*, "('cputime of fuction 1 (particles_inject)      =', f7.2,' %')")func_time(1)/tp1
            write(*, "('cputime of fuction 2 (escape_and_inject and inject)  =', f7.2,' %')")func_time(2)/tp1
            write(*, "('cputime of fuction 3 (mover->density_and_Es)         =', f7.2,' %')")func_time(3)/tp1
            write(*, "('cputime of fuction 4 (mover->push_RK4)               =', f7.2,' %')")func_time(4)/tp1
            if(iswitch_mcc==1)then
                write(*, "('cputime of fuction 5 (mcc)                           =', f7.2,' %')")func_time(5)/tp1
            else
                write(*, "('cputime of fuction 5 (mcc off)                       =', f7.2,' %')")func_time(5)/tp1
            endif
            write(*, "('cputime of fuction 6 (record_profiles)               =', f7.2,' %')")func_time(6)/tp1
            write(*, "('cputime of fuction 7 (display)                       =', f7.2,' %')")func_time(7)/tp1
            if(i_switch_MHD_Te==1)then
                write(*, "('cputime of fuction 8 (MHD on)                       =', f7.2,' %')")func_time(8)/tp1
            else
                write(*, "('cputime of fuction 8 (MHD off)                   =', f7.2,' %')")func_time(8)/tp1
            endif
            write(*, "('cputime of fuction 9 (FDFD)                          =', f7.2,' %')")func_time(9)/tp1
            write(*,*)
            write(*,*)
            write(*,*)
            write(*,*)
        
        else if (iswitch_log == 1) then
            if (.not. first_log_display) then
                write(*,'(a)',advance='no') esc//'[5A'
            endif
            
            write(*, "(a,'------------------------------------------------------------------------')") esc//'[2K'
            write(*, "(a,'running at ', f0.2,' % (estimating ',f0.2,' hours left)')") &
                esc//'[2K', tp1*100, time_run/60./tp1*(1-tp1)
            write(*, "(a,'(absorbed)  plasma loading = ', I0, ' mOhm')") &
                esc//'[2K', nint(2*absorbed_power/(irf_set**2)*1e3)
            write(*, "(a,'(deposited) plasma loading = ', I0, ' mOhm')") &
                esc//'[2K', nint(ptotal/(irf_set**2)*1e3)
            write(*, "(a,'recorded ', I0, ' of ', I0)") &
                esc//'[2K', ifig-1,n_pic

            first_log_display=.false.
            
        end if

    endif
    call find_func_cputime_2_of_2(func_time(7))  !-----2/2
end subroutine display_main



subroutine rec_time_ave
    USE the_whole_varibles
    implicit none
    integer,parameter::nz_len=6
    integer::k1,k2,i1,i2
    INTEGER*4::ip2 ,ip_acceleration
    Real*8::ek_x,ek_y,ek_z,ek_acceleration,zb3,z_count,dzn,zk1,zk2,Ek_ion_tp
    real*8::ek_ave(1:nz_len),ek_tol(1:nz_len),density(1:nz_len),z_c(1:nz_len)
    real*8::ek_tp1=0, ek_tp2=0, ek_tp3=0, ek_tp4=0,ek_ave_
    real*8::triple_product
302 format(<nrec_t>(e14.7,' '))

    if(mod(t,t_rec_evolution)<dt .or. it==1)then
        z_count=zl/nz_len
        dzn=zl/(nz_len-1)
        do iz=1,nz_len
            z_c(iz)=(iz-1)*dzn
        enddo
        density=0;
        ek_ave=0.;
        ek_tol=0.;
        ek_x=0.;ek_y=0.;ek_z=0.;
        do ip=1,n_active
            k1=int((x(ip,3)-z_count)/dzn)+1
            k2=int((x(ip,3)+z_count)/dzn)
            if(k1<0.5)k1=1
            if(k2>nz_len)k2=nz_len
            density(k1:k2)=density(k1:k2)+1

            tp1=mass_q_i_05*v(ip,1)**2
            tp2=mass_q_i_05*v(ip,2)**2
            tp3=mass_q_i_05*v(ip,3)**2
            ek_tol(k1:k2)=ek_tol(k1:k2)+tp1+tp2+tp3

            ek_x=ek_x+tp1
            ek_y=ek_y+tp2
            ek_z=ek_z+tp3
        enddo
        ek_ave(:)=ek_tol(:)/(1e-5+density)
        ek_x=ek_x/real(n_active)
        ek_y=ek_y/real(n_active)
        ek_z=ek_z/real(n_active)

        Ek_e_ave=mass_q_e_05*sum(v_e(1:n_active,1:3)**2)/real(n_active); !ev

        ek_acceleration=0
        ip_acceleration=0
        ek_tp1=0
        ek_tp2=0
        ek_tp3=0
        ek_tp4=0
        do ip=1,n_active
            Ek_ion_tp=0.5/q_mass_i*(v(ip,1)**2+v(ip,2)**2+v(ip,3)**2) !eV
            if (Ek_ion_tp>100.)then
                ip_acceleration=ip_acceleration+1;
                ek_acceleration=ek_acceleration+Ek_ion_tp;
            endif
            if (Ek_ion_tp>100. .and. Ek_ion_tp<1e3 )ek_tp1=ek_tp1+Ek_ion_tp
            if (Ek_ion_tp>1e3  .and. Ek_ion_tp<5e3 )ek_tp2=ek_tp2+Ek_ion_tp
            if (Ek_ion_tp>5e3  .and. Ek_ion_tp<10e3)ek_tp3=ek_tp3+Ek_ion_tp
            if (Ek_ion_tp>10e3 )                    ek_tp4=ek_tp4+Ek_ion_tp
        enddo
        !ek_acceleration=ek_acceleration/(ip_acceleration+1e-5);

        ek_ave_=ek_x+ek_y+ek_z

        if(sum(power_loss_ave(1:4))>1e-3)then
            tau_E_ave=(ek_x+ek_y+ek_z)/1.5*n_active*n_macro*qe_abs/sum(power_loss_ave(1:4))
        else
            tau_E_ave=1e-6
        endif
        
        triple_product=maxval(density_2D*n_macro*Ek_ion_2D/(1.5*1e3)*tau_E_ave)


        rec_t(1)=t
        rec_t(2)=real(it)
        rec_t(3)=n_active
        rec_t(4)=ek_x
        rec_t(5)=ek_y
        rec_t(6)=ek_z
        rec_t(7)=ip_acceleration
        rec_t(8)=ek_acceleration
        rec_t(9:14)=ek_ave(:)
        rec_t(15:20)=density(:)
        rec_t(21:24)=power_loss_ave(1:4)
        rec_t(25:28)=Ek_loss_ave(1:4)
        rec_t(29)=Ek_e_ave
        rec_t(34)=i_now
        rec_t(35)=ek_tp1
        rec_t(36)=ek_tp2
        rec_t(37)=ek_tp3
        rec_t(38)=ek_tp4

        rec_t(39:39)=ek_ave_
        rec_t(40)=triple_product
        rec_t(41)=tau_E_ave
        rec_t(42)=tau_p_ave
        write (42,302)rec_t
    endif
end subroutine rec_time_ave

subroutine rec_time_trajectory
    USE the_whole_varibles
    implicit none
    real*8:: B_mover(1:3),E_mover(1:3)
    real*8::sinth,costh,xtp,ytp,rtp
303 format(<nrec_p>(e14.7,' '))

    if(mod(t,t_rec_trajectory)<dt )then
        ip=1
        call interpolation_E_B(B_mover,E_mover)
        xtp=x(ip,1)
        ytp=x(ip,2)
        rtp=sqrt(xtp**2+ytp**2)+1e-20
        sinth=ytp/rtp
        costh=xtp/rtp

        rec_p(1)=t/trf
        rec_p(2:4)=x(ip,1:3)
        rec_p(5:7)=v(ip,1:3)
        rec_p(8:10)=B_mover(1:3)
        !rec_p(8)=v(ip,1)*costh+v(ip,2)*sinth
        !rec_p(9)=-v(ip,1)*sinth+v(ip,2)*costh
        !rec_p(10)=E_mover(1)*costh+E_mover(2)*sinth  !Er
        !rec_p(11)=-E_mover(1)*sinth+E_mover(2)*costh !Eth
        rec_p(12:14)=E_mover(1:3)

        rec_p(15:20)=0.0
        if(n_active>=2)then
            ip=2
            rec_p(15:17)=x(ip,1:3)
            rec_p(18:20)=v(ip,1:3)
        endif
        write (43,303)rec_p
    endif
end subroutine rec_time_trajectory

subroutine record_profiles
    USE the_whole_varibles
    implicit none
    real*8::Es_z(1:nz) !v_abs(1:np_max),z_p(1:np_max),ion_r(1:np_max),
    real*8, allocatable, save :: density_in_z(:),ek_ave(:,:)
    integer*4, save :: trigger=0
    integer*4, save :: acc_number=0
    character*30 fname
    character*256 fullpath
    call find_func_cputime_1_of_2  !-------------------1/2

    if(.not. allocated(density_in_z))then
        allocate(density_in_z(1:nz))
        allocate(ek_ave(1:nz,1:3))
        density_in_z=0.
        ek_ave=0.
        acc_number=0
        trigger=0
    endif

    if((mod((t+trf),(t_rec_profiles))<dt) .and. it>1 .and. trigger==0)then
        trigger = 1
        acc_number=0
        ek_ave=0.
        density_in_z=0.
    endif
    call find_Ek_1D_new(density_in_z,ek_ave,acc_number,trigger)
    if(mod(t,t_rec_profiles)<dt .or. it==1)then
        trigger = 0
        !call escape_periodicity !energy conservation
300     format(<nr>(e12.5,' '))
301     format(<np_loss_rec>(e12.5,' '))    
        !---------------------------record the whole data at specific time---------------------------!
        write (fname,120)ifig
120     format('plasma_',i0,'.dat')
        fullpath=trim(outputDir)//trim(fname)
        open (unit=20,file=fullpath,status='unknown',iostat=ierror)
        write (20,300)real(i_plasma_region(1:nr,1:nz))
        write (20,300)ni(1:nr,1:nz)
        write (20,300)te_in_FDFD(1:nr,1:nz)
        write (20,300)power_depo_rthz(1:nr,1:nz,1)+power_depo_rthz(1:nr,1:nz,2)+power_depo_rthz(1:nr,1:nz,3)
        close (20)

        write (fname,121)ifig
121     format('rf_field_',i0,'.dat')
        fullpath=trim(outputDir)//trim(fname)
        open (unit=21,file=fullpath,status='unknown',iostat=ierror)
        do itp=1,6
            write (21,300)abs(Erf6(1:nr,1:nz,itp))
        enddo
        close (21)


        do iz=1,nz
            Es_z(iz)=sum(Es_2D(:,iz,2))/nr
        enddo
        !call find_Ek_1D(density_in_z,ek_ave)

        if(acc_number>0)then
            ek_ave=ek_ave/(real(acc_number))
            density_in_z=density_in_z/(real(acc_number))
        else
            call find_Ek_1D(density_in_z,ek_ave)
        endif
        trigger = 0


322     format(<nz>(e12.5,' '))
122     format('Ek_ave_z_',i0,'.dat')
        write (fname,122)ifig
        fullpath=trim(outputDir)//trim(fname)
        open (unit=22,file=fullpath,status='unknown',iostat=ierror)
        write (22,322)z
        write (22,322)density_in_z
        write (22,322)ek_ave
        write (22,322)Es_z
        close(22)
323     format(<n_active>(e12.5,' '))
123     format('ion_rz_',i0,'.dat')
        write (fname,123)ifig
        fullpath=trim(outputDir)//trim(fname)
        open (unit=23,file=fullpath,status='unknown',iostat=ierror)
        write (23,323)sqrt(x(1:n_active,1)**2+x(1:n_active,2)**2 )!ion_r(1:np_max)
        write (23,323)x(1:n_active,3)
        write (23,323)v(1:n_active,1)
        write (23,323)v(1:n_active,2)
        write (23,323)v(1:n_active,3)
        write (23,323)t-life_and_ek(1:n_active,1)
        write (23,323)life_and_ek(1:n_active,2)
        close(23)
        

124     format('density_Es_2D_',i0,'.dat')
        write (fname,124)ifig
        fullpath=trim(outputDir)//trim(fname)
        open (unit=24,file=fullpath,status='unknown',iostat=ierror)
        write (24,300)density_2D
        write (24,300)Es_2D(:,:,1)
        write (24,300)Es_2D(:,:,2)
        write (24,300)Ek_ion_2D
        write (24,300)te_2D
        write (24,300)u_mhd(1:nr,1:nz,1)
        write (24,300)u_mhd(1:nr,1:nz,2)
        write (24,300)ni
        write (24,300)te_mhd
        write (24,300)Q_ie(1:nr,1:nz)
        write (24,300)u_pic(1:nr,1:nz,1)
        write (24,300)u_pic(1:nr,1:nz,2)
        write (24,300)density_2D*n_macro*Ek_ion_2D/(1.5*1e3)*tau_E_ave !triple_product
        close(24)


125     format('Erf_all_m_',i0,'.dat')
        write (fname,125)ifig
        fullpath=trim(outputDir)//trim(fname)
        open (unit=25,file=fullpath,status='unknown',iostat=ierror)
        do m=m_start,m_end
            do itp=1,3
                write (25,300)real(e_output(m,1:nr,1:nz,itp))
                write (25,300)imag(e_output(m,1:nr,1:nz,itp))
            enddo
        enddo
        close(25)
        
126     format('xv_loss_',i0,'.dat')
        write (fname,126)ifig
        fullpath=trim(outputDir)//trim(fname)
        open (unit=26,file=fullpath,status='unknown',iostat=ierror)
        write (26,301)xv_loss
        close(26)
        
        

        call rec_para
        ifig=ifig+1
    endif
    call find_func_cputime_2_of_2(func_time(6))  !-----2/2
end subroutine record_profiles


subroutine find_Ek_1D(density_x,ek_x) !,x_p,v_innz,np_max,z,
    USE the_whole_varibles
    implicit none
    !integer::np_max,nzip,
    integer::ix1,ix2,I_x,i_comp
    real*8::ek_x(1:nz,1:3),density_x(1:nz) !,x_p(np_max),v_in(1:np_max,1:3),z(1:nz)dz,
    real*8::xtp,s1,s2,min_x,dx_tp1,dx_tp2
    ek_x=0.
    density_x=1e-5

    do ip=1,n_active
        xtp=x(ip,3)
        I_x=minloc(abs(xtp-z),1);
        min_x=xtp-z(I_x);
        if (min_x<0)then
            ix1=I_x-1;
            ix2=I_x;
        else
            ix1=I_x;
            ix2=I_x+1;
        endif
        if(ix1<1)ix1=1
        if(ix2>nz)ix2=nz
        dx_tp1=abs(xtp-z(ix1)); !note abs()
        s1=dx_tp1/dz;
        s2=1.-s1;

        ek_x(ix1,1:3)=ek_x(ix1,1:3)+s2*mass_q_i_05*v(ip,1:3)**2
        ek_x(ix2,1:3)=ek_x(ix2,1:3)+s1*mass_q_i_05*v(ip,1:3)**2
        density_x(ix1)=density_x(ix1)+s2
        density_x(ix2)=density_x(ix2)+s1
    enddo
    do ix1=1,3
        ek_x(:,ix1)=ek_x(:,ix1)/density_x
    enddo
    density_x=density_x/(pi*rp**2*dz)
end subroutine find_Ek_1D

subroutine find_Ek_1D_new(density_x,ek_x,acc_number,trigger) !,x_p,v_innz,np_max,z,
    USE the_whole_varibles
    implicit none
    !integer::np_max,nzip,
    integer::ix1,ix2,I_x,i_comp
    real*8::ek_x(1:nz,1:3),density_x(1:nz) !,x_p(np_max),v_in(1:np_max,1:3),z(1:nz)dz,
    real*8::ek_x_inter(1:nz,1:3),density_x_inter(1:nz)
    real*8::xtp,s1,s2,min_x,dx_tp1,dx_tp2
    integer*4 :: acc_number,trigger
    
    if(trigger == 1) then
    
    acc_number=acc_number+1
    
    ek_x_inter=0.
    density_x_inter=1e-5

    do ip=1,n_active
        xtp=x(ip,3)
        I_x=minloc(abs(xtp-z),1);
        min_x=xtp-z(I_x);
        if (min_x<0)then
            ix1=I_x-1;
            ix2=I_x;
        else
            ix1=I_x;
            ix2=I_x+1;
        endif
        if(ix1<1)ix1=1
        if(ix2>nz)ix2=nz
        dx_tp1=abs(xtp-z(ix1)); !note abs()
        s1=dx_tp1/dz;
        s2=1.-s1;

        ek_x_inter(ix1,1:3)=ek_x_inter(ix1,1:3)+s2*mass_q_i_05*v(ip,1:3)**2
        ek_x_inter(ix2,1:3)=ek_x_inter(ix2,1:3)+s1*mass_q_i_05*v(ip,1:3)**2
        density_x_inter(ix1)=density_x_inter(ix1)+s2
        density_x_inter(ix2)=density_x_inter(ix2)+s1
    enddo
    do ix1=1,3
        ek_x_inter(:,ix1)=ek_x_inter(:,ix1)/density_x_inter
    enddo
    density_x_inter=density_x_inter/(pi*rp**2*dz)
    
    ek_x=ek_x+ek_x_inter
    density_x=density_x+density_x_inter
    
    endif
end subroutine find_Ek_1D_new

subroutine start_ftime
    USE the_whole_varibles
    implicit none
    call date_and_time(timefunc_char(1,1),timefunc_char(2,1), timefunc_char(3,1), timefunc_int(:,1))
end subroutine start_ftime

subroutine end_ftime
    USE the_whole_varibles
    implicit none
    call date_and_time(timefunc_char(1,2),timefunc_char(2,2), timefunc_char(3,2), timefunc_int(:,2))
    timefunc_int_real=real(timefunc_int)
    tp1=(timefunc_int_real(2,1)-1)*30.*24.*60.*60.+(timefunc_int_real(3,1)-1)*24.*60.*60.+timefunc_int_real(5,1)*60.*60.+ &
        &                timefunc_int_real(6,1)*60.+timefunc_int_real(7,1)+timefunc_int_real(8,1)*1e-3
    tp2=(timefunc_int_real(2,2)-1)*30.*24.*60.*60.+(timefunc_int_real(3,2)-1)*24.*60.*60.+timefunc_int_real(5,2)*60.*60.+ &
        &                timefunc_int_real(6,2)*60.+timefunc_int_real(7,2)+timefunc_int_real(8,2)*1e-3
    time_run=(tp2-tp1)/60.
    date_time=timefunc_int
end subroutine end_ftime



subroutine find_func_cputime_1_of_2
    USE the_whole_varibles
    implicit none
    call CPU_TIME(time_diff(1))
end subroutine find_func_cputime_1_of_2

subroutine find_func_cputime_2_of_2(time_used)
    USE the_whole_varibles
    implicit none
    real*8 time_used
    call CPU_TIME(time_diff(2))
    time_used=time_used+time_diff(2)-time_diff(1)
end subroutine find_func_cputime_2_of_2
