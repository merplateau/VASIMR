subroutine interp_position(n_len,pdf_x,x_grid,x_out)
    implicit none
    integer::n_bef,n_len
    real*8::rand_0_1,x_out,pdf_x(1:n_len),x_grid(1:n_len)
    integer::I_x,ix_bef1,ix_bef2
    real*8::s1,s2,min_x,dx_tp1,dx_tp2

    call random_number(rand_0_1)
    I_x=minloc(abs(rand_0_1-pdf_x),1);
    min_x=rand_0_1-pdf_x(I_x);
    if (min_x<0)then
        ix_bef1=I_x-1;
        ix_bef2=I_x;
    else
        ix_bef1=I_x;
        ix_bef2=I_x+1;
    endif

    if(ix_bef1<1)then
        ix_bef1=1
        ix_bef2=2
    elseif(ix_bef2>n_len)then
        ix_bef2=n_len
        ix_bef1=ix_bef2-1
    endif

    dx_tp1=rand_0_1-pdf_x(ix_bef1);
    dx_tp2=pdf_x(ix_bef2)-rand_0_1;
    s1=dx_tp1/(dx_tp1+dx_tp2);
    s2=1.-s1;
    x_out=s2*x_grid(ix_bef1)+s1*x_grid(ix_bef2);
    continue
end subroutine interp_position


subroutine particles_initialization
    USE the_whole_varibles
    implicit none
    integer*4 :: nseed(1)=123456
    real*8::rand_0_1,rand_1_1,xtp,ytp,rtp,th_tp
    real*8::vr_tp,vx_tp,vy_tp,vz_tp,br_tp,bz_tp,b0_tp,sinth_b0,costh_b0
    real*8::ztp,vz_ll,coeff_vx

    call RANDOM_SEED(put=nseed)

    !change the seed for the pseudorandom number generator.
    !Uncomment and simulation won't repeat absolutly run to run.
    !call RANDOM_SEED()

    coeff_vx=1;
    ip=1
    do while (ip<=n_active)
        !use PDF and random number to generate specific particle distribution

        call interp_position(nr,pdf_ne_r,r,rtp)
        call interp_position(nz,pdf_ne_z,z,ztp)

        if(rtp<0)rtp=abs(rtp)
        if(rtp<rp )then ! .or. (rtp<rl .and. ztp>z_vac(1) .and. ztp<z_vac(2))
            call random_number(rand_0_1)
            th_tp=2*pi*rand_0_1
            xtp=rtp*cos(th_tp);
            ytp=rtp*sin(th_tp);

            if(ztp<zs)ztp=2*zs-ztp
            if(ztp>zd)ztp=2*zd-ztp

            x(ip,1)=xtp
            x(ip,2)=ytp
            x(ip,3)=ztp

            call random_number(tp1)
            call random_number(tp2)
            tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
            v(ip,1)=coeff_vx*vi_ex*tp3
            v_e(ip,1)=coeff_vx*ve_ex*tp3

            call random_number(tp1)
            call random_number(tp2)
            tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
            v(ip,2)=vi_ex*tp3
            v_e(ip,2)=ve_ex*tp3

            call random_number(tp1)
            call random_number(tp2)
            tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
            v(ip,3)=vi_ex*tp3+vz0
            v_e(ip,3)=ve_ex*tp3
            t_np(ip)=0.0d0
            life_and_ek(ip,1)=0.0d0
            life_and_ek(ip,2)=mass_q_i_05*sum(v(ip,1:3)**2)

            !np_max=np_max+1
            ip=ip+1;
        endif
    end do
    !x_e=x
    np_e=n_active

    continue
end subroutine particles_initialization

subroutine set_particle_velocity(ip_set)
    USE the_whole_varibles
    implicit none
    integer*4 :: ip_set

    call random_number(tp1)
    call random_number(tp2)
    tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
    v(ip_set,1)=vi_ex*tp3
    v_e(ip_set,1)=ve_ex*tp3

    call random_number(tp1)
    call random_number(tp2)
    tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
    v(ip_set,2)=vi_ex*tp3
    v_e(ip_set,2)=ve_ex*tp3

    call random_number(tp1)
    call random_number(tp2)
    tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
    v(ip_set,3)=vi_ex*tp3+vz0
    v_e(ip_set,3)=ve_ex*tp3
end subroutine set_particle_velocity


!subroutine particles_inject
    !USE the_whole_varibles
    !implicit none
    !integer::ip_inj
    !real*8::rand_0_1,rtp
    !!real*8::vr_tp,vx_tp,vy_tp,vz_tp,br_tp,bz_tp,b0_tp,sinth_b0,costh_b0
    !!real*8::ztp,vz_ll
    !
    !call find_func_cputime_1_of_2  !-------------------1/2
    !
    !!change the seed for the pseudorandom number generator.
    !!Uncomment and simulation won't repeat absolutly run to run.
    !!call RANDOM_SEED()
    !
    !!inject ion
    !ip_inj=0
    !num_inject=np_max-np_max
    !do while (ip_inj<num_inject)
    !    !call random_number(rand_0_1)
    !    call interp_position(nrp,pdf_ne_source_r,r_particle_inj,rtp)
    !    if(rtp<0)rtp=abs(rtp)
    !    if(rtp<rp)then
    !        ip_inj=ip_inj+1
    !        np_max=np_max+1
    !        ip=np_max;
    !        if( np_max>np_max) exit
    !        call injection_sub(rtp,x(ip,1:3),v(ip,1:3),vi_ex)
    !    endif
    !end do
    !N_inject_particles=num_inject
    !
    !call find_func_cputime_2_of_2(func_time(1))  !-----2/2
!end subroutine particles_inject

subroutine injection_sub(rtp,x3,v3,v_ex)
    USE the_whole_varibles
    implicit none
    real*8:: x3(1:3),v3(1:3),v_ex
    real*8::rand_0_1,xtp,ytp,ztp,rtp,th_tp


    call random_number(rand_0_1)
    th_tp=2*pi*rand_0_1
    xtp=rtp*cos(th_tp);
    ytp=rtp*sin(th_tp);

    !call random_number(rand_0_1)
    call interp_position(nz,pdf_ne_source_z,z,ztp)
    if(ztp<zs)ztp=2*zs-ztp
    if(ztp>zd)ztp=2*zd-ztp

    x3(1)=xtp
    x3(2)=ytp
    x3(3)=ztp

    call random_number(tp1)
    call random_number(tp2)
    tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
    v3(1)=v_ex*tp3

    call random_number(tp1)
    call random_number(tp2)
    tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
    v3(2)=v_ex*tp3

    call random_number(tp1)
    call random_number(tp2)
    tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
    v3(3)=v_ex*tp3+vz0
end subroutine injection_sub

subroutine injection_inlet_sub(ip_set)
    USE the_whole_varibles
    implicit none
    integer*4 :: ip_set
    real*8 :: rand_0_1, rtp, th_tp, z_layer

    call interp_position(nrp,pdf_ne_source_r,r_particle_inj,rtp)
    if(rtp<0)rtp=abs(rtp)

    call random_number(rand_0_1)
    th_tp=2*pi*rand_0_1
    x(ip_set,1)=rtp*cos(th_tp)
    x(ip_set,2)=rtp*sin(th_tp)

    call random_number(rand_0_1)
    z_layer=vz0*dt*rand_0_1
    x(ip_set,3)=zs+z_layer

    call set_particle_velocity(ip_set)
    t_np(ip_set)=t
    life_and_ek(ip_set,1)=t
    life_and_ek(ip_set,2)=mass_q_i_05*sum(v(ip_set,1:3)**2)
end subroutine injection_inlet_sub

subroutine mover
    USE the_whole_varibles
    implicit none
    real*8:: B_mover(1:3),E_mover(1:3) !,xv(1:3)

    call find_func_cputime_1_of_2  !-------------------1/2
    call find_x_to_grid
    call find_density_and_Es_2D
    call find_func_cputime_2_of_2(func_time(3))  !-----2/2

    call find_func_cputime_1_of_2  !-------------------1/2
    expt=exp(-i*2*pi*frequency*t)
    expt2=exp(-i*2*pi*frequency2*t)
    do ip=1,n_active
        call interpolation_E_B(B_mover,E_mover)
        call push_RK4(B_mover,E_mover)
    enddo
    call find_func_cputime_2_of_2(func_time(4))  !-----2/2
end subroutine mover

subroutine find_x_to_grid !here x means (x,y,z)
    USE the_whole_varibles
    implicit none
    integer::ir_tp,iz_tp
    real*8::r_p(1:n_active)

    !r: s1= x_to_grid(1:n_active,1),s2=1-s1; z: s3= x_to_grid(1:n_active,2),s4=1-s3;
    x_to_grid=0.
    r_p(1:n_active)=sqrt(x(1:n_active,1)**2+x(1:n_active,2)**2)
    ir1_iz1_grid(1:n_active,1)=int((r_p(1:n_active)-r(1))/dr)+1
    ir1_iz1_grid(1:n_active,2)=int((x(1:n_active,3)-z(1))/dz)+1
    do ip=1,n_active
        ir_tp=ir1_iz1_grid(ip,1)
        iz_tp=ir1_iz1_grid(ip,2)
        x_to_grid(ip,1)=(r_p(ip)-r(ir_tp))/dr
        x_to_grid(ip,2)=(x(ip,3)-z(iz_tp))/dz
    enddo
end subroutine find_x_to_grid

subroutine interpolation_E_B(B_mover,E_mover)
    USE the_whole_varibles
    implicit none
    INTEGER*2:: ir1,ir2,iz1,iz2
    INTEGER  :: ixtp,iytp
    real*8::sinth,costh,erthz(1:3),xtp,ytp,rtp,th_tp
    complex*16::c1_6(1:6),Em_tp(1:3)
    complex*16::phase_factor !相位因子变量
    real*8:: B_mover(1:3),E_mover(1:3)
    Real*8::s1,s2,s3,s4
    Real*8::br,bz
    Real*8::Es_dc(1:3)=0
    
    !B_mover,E_mover - Bx By Bz Ex Ey Ez

    rtp=sqrt(x(ip,1)**2+x(ip,2)**2)+1e-20;
    xtp=x(ip,1)
    ytp=x(ip,2)
    sinth=ytp/rtp
    costh=xtp/rtp

    s1=x_to_grid(ip,1)*x_to_grid(ip,2)
    s2=x_to_grid(ip,1)*(1.-x_to_grid(ip,2))
    s3=(1.-x_to_grid(ip,1))*(1.-x_to_grid(ip,2))
    s4=(1.-x_to_grid(ip,1))*x_to_grid(ip,2)

    ir1=ir1_iz1_grid(ip,1)
    iz1=ir1_iz1_grid(ip,2)
    ir2=ir1+1
    iz2=iz1+1

    br=s3*b0_DC(ir1,iz1,1)+s1*b0_DC(ir2,iz2,1)+s2*b0_DC(ir2,iz1,1)+s4*b0_DC(ir1,iz2,1);
    bz=s3*b0_DC(ir1,iz1,3)+s1*b0_DC(ir2,iz2,3)+s2*b0_DC(ir2,iz1,3)+s4*b0_DC(ir1,iz2,3);
    B_mover(1)=br*costh
    B_mover(2)=br*sinth;
    B_mover(3)=bz

    Es_dc(1)=s3*Es_2D(ir1,iz1,1)+s1*Es_2D(ir2,iz2,1)+s2*Es_2D(ir2,iz1,1)+s4*Es_2D(ir1,iz2,1);
    Es_dc(3)=s3*Es_2D(ir1,iz1,2)+s1*Es_2D(ir2,iz2,2)+s2*Es_2D(ir2,iz1,2)+s4*Es_2D(ir1,iz2,2);
    
    c1_6 = (0.0d0, 0.0d0)

    th_tp = atan2(ytp, xtp)

    do m = m_start, m_end
        Em_tp = s3*e_output(m,ir1,iz1,:) + s1*e_output(m,ir2,iz2,:) + &
                s2*e_output(m,ir2,iz1,:) + s4*e_output(m,ir1,iz2,:)
        
        phase_factor = exp(cmplx(0.0d0, 1.0d0) * m * th_tp) !计算相位因子 e^(im*theta)
        
        c1_6(1:3) = c1_6(1:3) + Em_tp * phase_factor
    enddo

    if(iswitch_RF2 == 1) then
        c1_6(4:6) = s3*Erf6(ir1,iz1,4:6)+s1*Erf6(ir2,iz2,4:6)+ &
                    s2*Erf6(ir2,iz1,4:6)+s4*Erf6(ir1,iz2,4:6)
    endif
    
    !state_power_on_off  on->1;  0->off
    erthz(1:3)=Es_dc(1:3)+state_power_on_off*(real(expt*c1_6(1:3))+real(expt2*c1_6(4:6)));

    ! (r,th,z) ->(x,y,z)
    E_mover(1)=erthz(1)*costh-erthz(2)*sinth
    E_mover(2)=erthz(1)*sinth+erthz(2)*costh
    E_mover(3)=erthz(3)
end subroutine interpolation_E_B

subroutine push_RK4(B_mover,E_mover)
    USE the_whole_varibles
    implicit none
    real*8:: B_mover(1:3),E_mover(1:3) 
    real*8:: k1(1:6),k2(1:6),k3(1:6),k4(1:6)
    real*8:: xv(1:6)


    !find k1
    xv(1:3)=x(ip,1:3)
    xv(4:6)=v(ip,1:3)
    k1(1:3)=xv(4:6)
    k1(4)=q_mass_i*(E_mover(1)+xv(5)*B_mover(3)-xv(6)*B_mover(2));
    k1(5)=q_mass_i*(E_mover(2)+xv(6)*B_mover(1)-xv(4)*B_mover(3))
    k1(6)=q_mass_i*(E_mover(3)+xv(4)*B_mover(2)-xv(5)*B_mover(1))

    !find k2
    xv(1:3)=x(ip,1:3)+dt*0.5*k1(1:3)
    xv(4:6)=v(ip,1:3)+dt*0.5*k1(4:6)
    k2(1:3)=xv(4:6)
    k2(4)=q_mass_i*(E_mover(1)+xv(5)*B_mover(3)-xv(6)*B_mover(2));
    k2(5)=q_mass_i*(E_mover(2)+xv(6)*B_mover(1)-xv(4)*B_mover(3))
    k2(6)=q_mass_i*(E_mover(3)+xv(4)*B_mover(2)-xv(5)*B_mover(1))

    !find k3
    xv(1:3)=x(ip,1:3)+dt*0.5*k2(1:3)
    xv(4:6)=v(ip,1:3)+dt*0.5*k2(4:6)
    k3(1:3)=xv(4:6)
    k3(4)=q_mass_i*(E_mover(1)+xv(5)*B_mover(3)-xv(6)*B_mover(2));
    k3(5)=q_mass_i*(E_mover(2)+xv(6)*B_mover(1)-xv(4)*B_mover(3))
    k3(6)=q_mass_i*(E_mover(3)+xv(4)*B_mover(2)-xv(5)*B_mover(1))

    !find k4
    xv(1:3)=x(ip,1:3)+dt*k3(1:3)
    xv(4:6)=v(ip,1:3)+dt*k3(4:6)
    k4(1:3)=xv(4:6)
    k4(4)=q_mass_i*(E_mover(1)+xv(5)*B_mover(3)-xv(6)*B_mover(2));
    k4(5)=q_mass_i*(E_mover(2)+xv(6)*B_mover(1)-xv(4)*B_mover(3))
    k4(6)=q_mass_i*(E_mover(3)+xv(4)*B_mover(2)-xv(5)*B_mover(1))

    x(ip,1:3)=x(ip,1:3)+dt/6.*(k1(1:3)+2.*k2(1:3)+2.*k3(1:3)+k4(1:3))
    v(ip,1:3)=v(ip,1:3)+dt/6.*(k1(4:6)+2.*k2(4:6)+2.*k3(4:6)+k4(4:6))
end subroutine push_RK4

subroutine find_density_and_Es_2D
    USE the_whole_varibles
    implicit none
    real*8::dni_ni,ti_tmp,te_min_set,density_cri,N_te,Ek_ion_center,ratio_te !te_max_set,
    !real*8::ni_int
    integer n_sm,ifilter,ir1,ir2,iz1,iz2,i0,i1,i2,iEs !,itp

    te_min_set=2. !eV
    call density_Ek_2D_sub
    !itp=2
    do itp=1,2
        if(itp==1)n_sm=3   !1st smooth  3
        if(itp==2)n_sm=2   !2st smooth 2

        ir1=1;ir2=nr_vac(1);iz1=1;iz2=nz;
        call smooth_2d(n_sm,2*n_sm,ir1,ir2,iz1,iz2,density_2D(ir1:ir2,iz1:iz2))
        call smooth_2d(n_sm,2*n_sm,ir1,ir2,iz1,iz2, Ek_ion_2D(ir1:ir2,iz1:iz2))
        call smooth_2d(n_sm,2*n_sm,ir1,ir2,iz1,iz2, u_pic(ir1:ir2,iz1:iz2,1))
        call smooth_2d(n_sm,2*n_sm,ir1,ir2,iz1,iz2, u_pic(ir1:ir2,iz1:iz2,2))

        !ir1=1;ir2=nr_vac(2);iz1=nz_vac(1);iz2=nz;
        !call smooth_2d(n_sm,2*n_sm,ir1,ir2,iz1,iz2,density_2D(ir1:ir2,iz1:iz2))
        !call smooth_2d(n_sm,2*n_sm,ir1,ir2,iz1,iz2, Ek_ion_2D(ir1:ir2,iz1:iz2))
        !call smooth_2d(n_sm,2*n_sm,ir1,ir2,iz1,iz2, u_pic(ir1:ir2,iz1:iz2,1))
        !call smooth_2d(n_sm,2*n_sm,ir1,ir2,iz1,iz2, u_pic(ir1:ir2,iz1:iz2,2))
    enddo

    
    if(i_switch_MHD_Te==0)then
        te_2D=te_ave
    elseif(i_switch_MHD_Te==2)then
        te_2D=Te_set
    else
        te_2D=te_mhd
    endif

    !Es=-te_in_FDFD*grad(ni)/ni with ni=ni
    Es_2D=0.

    ir1=2;ir2=nr_vac(1)-1;iz1=1;iz2=nz;
    call find_Es_center(ir1,ir2,iz1,iz2,1) !find Er
    ir1=1;ir2=nr_vac(1);iz1=2;iz2=nz-1;
    call find_Es_center(ir1,ir2,iz1,iz2,2) !find Ez

    !ir1=nr_vac(1);ir2=nr_vac(2)-1;iz1=nz_vac(1);iz2=nz_vac(2);
    !call find_Es_center(ir1,ir2,iz1,iz2,1) !find Er
    !ir1=nr_vac(1);ir2=nr_vac(2);iz1=nz_vac(1)+1;iz2=nz_vac(2)-1;
    !call find_Es_center(ir1,ir2,iz1,iz2,2) !find Ez


    i0=1;i1=1;i2=nz; !r=0
    call find_Es_bdry(i0,i1,i2,1)

    i0=nr_vac(1);i1=1;i2=nz_vac(1); !r=r_vac(1)
    call find_Es_bdry(i0,i1,i2,2)
    !i0=nr_vac(2);i1=nz_vac(2);i2=nz_vac(2); !r=r_vac(2)
    !call find_Es_bdry(i0,i1,i2,2)

    i0=1;i1=1;i2=nr_vac(1); !z=0
    call find_Es_bdry(i0,i1,i2,3)
    !i0=nz_vac(1);i1=nr_vac(1);i2=nr_vac(2); !z=z_vac(1)
    !call find_Es_bdry(i0,i1,i2,3)

    i0=nz;i1=1;i2=nr_vac(1); !z=nz
    call find_Es_bdry(i0,i1,i2,4)
    !i0=nz;i1=1;i2=nr_vac(2); !z=nz
    !call find_Es_bdry(i0,i1,i2,4)
end subroutine find_density_and_Es_2D

subroutine find_Es_center(ir1,ir2,iz1,iz2,iEs)
    USE the_whole_varibles
    implicit none
    integer::ir1,ir2,iz1,iz2,iEs
    !iEs=1 ->find Er
    !iEs=2 ->find Ez

    if(iEs==1)then
        Es_2D(ir1:ir2,iz1:iz2,1)=-te_2D(ir1:ir2,iz1:iz2)*(density_2D(ir1+1:ir2+1,iz1:iz2)- &
            & density_2D(ir1-1:ir2-1,iz1:iz2))/(2*dr*density_2D(ir1:ir2,iz1:iz2))
    endif

    if(iEs==2)then
        Es_2D(ir1:ir2,iz1:iz2,2)=-te_2D(ir1:ir2,iz1:iz2)*(density_2D(ir1:ir2,iz1+1:iz2+1)- &
            & density_2D(ir1:ir2,iz1-1:iz2-1))/(2*dz*density_2D(ir1:ir2,iz1:iz2))
    endif
end subroutine find_Es_center

subroutine find_Es_bdry(i0,i1,i2,iEs)
    USE the_whole_varibles
    implicit none
    integer::i1,i2,i0,iEs
    !iEs=1 ->find low BD for Er
    !iEs=2 ->find up BD for Er
    !iEs=3 ->find left BD for Ez
    !iEs=4 ->find right BD for Ez

    if(iEs==1)then
        Es_2D(i0,i1:i2,1)=-te_2D(i0,i1:i2)*(density_2D(i0+1,i1:i2)-density_2D(i0,i1:i2))/(dr*density_2D(i0,i1:i2))
    endif

    if(iEs==2)then
        Es_2D(i0,i1:i2,1)=-te_2D(i0,i1:i2)*(density_2D(i0,i1:i2)-density_2D(i0-1,i1:i2))/(dr*density_2D(i0,i1:i2))
    endif

    if(iEs==3)then
        Es_2D(i1:i2,i0,2)=-te_2D(i1:i2,i0)*(density_2D(i1:i2,i0+1)-density_2D(i1:i2,i0))/(dr*density_2D(i1:i2,i0))
    endif

    if(iEs==4)then
        Es_2D(i1:i2,i0,2)=-te_2D(i1:i2,i0)*(density_2D(i1:i2,i0)-density_2D(i1:i2,i0-1))/(dr*density_2D(i1:i2,i0))
    endif
end subroutine find_Es_bdry

subroutine density_Ek_2D_sub
    USE the_whole_varibles
    implicit none
    real*8::sinth,costh,rtp,ur_tp
    real*8::s1,s2,s3,s4,vtp2
    integer::ir2,ir1,iz2,iz1

    density_2D=0.2 !set a small value to avoid zero or very small density
    Ek_ion_2D=0.
    Ek_ion_2D_r=0.
    u_pic=0.

    !s1-s4, squre of bottom left, bottom right, upper right, upper left
    ! ^
    ! |  Radial coordinates r, to up
    ! -> Axial coordinates  z, to right
    !  |s4| |s3|
    !   x(r,z)
    !  |s1| |s2|

    do ip=1,n_active
        s1=x_to_grid(ip,1)*x_to_grid(ip,2)
        s2=x_to_grid(ip,1)*(1.-x_to_grid(ip,2))
        s3=(1.-x_to_grid(ip,1))*(1.-x_to_grid(ip,2))
        s4=(1.-x_to_grid(ip,1))*x_to_grid(ip,2)

        ir1=ir1_iz1_grid(ip,1)
        iz1=ir1_iz1_grid(ip,2)
        ir2=ir1+1
        iz2=iz1+1

        density_2D(ir1,iz1)=density_2D(ir1,iz1)+s3
        density_2D(ir2,iz2)=density_2D(ir2,iz2)+s1
        density_2D(ir2,iz1)=density_2D(ir2,iz1)+s2
        density_2D(ir1,iz2)=density_2D(ir1,iz2)+s4

        vtp2=mass_q_i_05*(v(ip,1)**2+v(ip,2)**2+v(ip,3)**2)
        Ek_ion_2D(ir1,iz1)=Ek_ion_2D(ir1,iz1)+s3*vtp2
        Ek_ion_2D(ir2,iz2)=Ek_ion_2D(ir2,iz2)+s1*vtp2
        Ek_ion_2D(ir2,iz1)=Ek_ion_2D(ir2,iz1)+s2*vtp2
        Ek_ion_2D(ir1,iz2)=Ek_ion_2D(ir1,iz2)+s4*vtp2

        vtp2=mass_q_i_05*(v(ip,1)**2+v(ip,2)**2)
        Ek_ion_2D_r(ir1,iz1)=Ek_ion_2D_r(ir1,iz1)+s3*vtp2
        Ek_ion_2D_r(ir2,iz2)=Ek_ion_2D_r(ir2,iz2)+s1*vtp2
        Ek_ion_2D_r(ir2,iz1)=Ek_ion_2D_r(ir2,iz1)+s2*vtp2
        Ek_ion_2D_r(ir1,iz2)=Ek_ion_2D_r(ir1,iz2)+s4*vtp2

        rtp=sqrt(x(ip,1)**2+x(ip,2)**2)+1e-20;
        sinth=x(ip,2)/rtp
        costh=x(ip,1)/rtp
        ur_tp=v(ip,1)*costh+v(ip,2)*sinth
        u_pic(ir1,iz1,1)=u_pic(ir1,iz1,1)+s3*ur_tp
        u_pic(ir2,iz2,1)=u_pic(ir2,iz2,1)+s1*ur_tp
        u_pic(ir2,iz1,1)=u_pic(ir2,iz1,1)+s2*ur_tp
        u_pic(ir1,iz2,1)=u_pic(ir1,iz2,1)+s4*ur_tp

        u_pic(ir1,iz1,2)=u_pic(ir1,iz1,2)+s3*v(ip,3)
        u_pic(ir2,iz2,2)=u_pic(ir2,iz2,2)+s1*v(ip,3)
        u_pic(ir2,iz1,2)=u_pic(ir2,iz1,2)+s2*v(ip,3)
        u_pic(ir1,iz2,2)=u_pic(ir1,iz2,2)+s4*v(ip,3)

    enddo
    !Ek_ion_2D must be divided by the density without column coordinate coefficients (2pi*rdr)
    Ek_ion_2D=Ek_ion_2D/density_2D; 
    Ek_ion_2D_r=Ek_ion_2D_r/density_2D;

    !Note that densities in column coordinates need to be divided by the (2pi*rdr)
    do ir=2,nr!nrp
        density_2D(ir,:)=density_2D(ir,:)/(pi*(r2(ir)**2-r2(ir-1)**2)*dz)
    enddo
    ir=1;
    density_2D(ir,:)=density_2D(ir+1,:)
end subroutine density_Ek_2D_sub

!subroutine smooth_1d(n_sm,n_data,data_in)
    !implicit none
    !integer::n_sm,n_data,n_total,k1,k2,k3,ix_sm
    !real*8::data_in(1:n_data),data_out(1:n_data)
    !
    !n_total=2*n_sm+1
    !do ix_sm=n_sm+1,n_data-n_sm
    !    k1=ix_sm-n_sm
    !    k2=ix_sm+n_sm
    !    data_out(ix_sm)=sum(data_in(k1:k2))/n_total
    !enddo
    !do ix_sm=2,n_sm
    !    k1=1
    !    k2=(ix_sm-1)*2+1
    !    k3=k2-k1+1
    !    data_out(ix_sm)=sum(data_in(k1:k2))/k3
    !enddo
    !do ix_sm=n_data-n_sm+1,n_data-1
    !    k1=ix_sm-(n_data-ix_sm)
    !    k2=n_data
    !    k3=k2-k1+1
    !    data_out(ix_sm)=sum(data_in(k1:k2))/k3
    !enddo
    !data_in(2:n_data-1)=data_out(2:n_data-1)
!end subroutine smooth_1d

subroutine smooth_2d(n_smr,n_smz,nr_sta,nr_end,nz_sta,nz_end,data_in)
    implicit none
    integer  n_smr,n_smz,nr_sta,nr_end,nz_sta,nz_end,n_total,k1,k2,k3,k4,ir_sm,iz_sm
    integer  ir_reg1,ir_reg2,iz_reg1,iz_reg2
    real*8   data_in(nr_sta:nr_end,nz_sta:nz_end),data_out(nr_sta:nr_end,nz_sta:nz_end)

    ir_reg1=n_smr+nr_sta-1
    ir_reg2=nr_end-n_smr
    iz_reg1=n_smz+nz_sta-1
    iz_reg2=nz_end-n_smz
    do ir_sm=nr_sta,nr_end
        do iz_sm=nz_sta,nz_end
            if( ir_sm<=ir_reg1)then
                k1=1+nr_sta-1;k2=ir_sm+(ir_sm-k1)
            elseif(ir_sm>=ir_reg2)then
                k2=nr_end;k1=ir_sm-(nr_end-ir_sm);
            else
                k1=ir_sm-n_smr
                k2=ir_sm+n_smr
            endif

            if( iz_sm<=iz_reg1)then
                k3=1+nz_sta-1;k4=iz_sm+(iz_sm-k3)
            elseif(iz_sm>=iz_reg2)then
                k4=nz_end;k3=iz_sm-(nz_end-iz_sm);
            else
                k3=iz_sm-n_smz
                k4=iz_sm+n_smz
            endif

            n_total=(k4-k3+1)*(k2-k1+1)
            data_out(ir_sm,iz_sm)=sum(data_in(k1:k2,k3:k4))/n_total
        enddo
    enddo

    data_in=data_out;

end subroutine smooth_2d

subroutine escape_and_inject
    USE the_whole_varibles
    implicit none
    INTEGER*2::ip2,iloss
    real*8::zb0,rb1,zb1,zb2,rtp,ztp,v2_lim,vtp2!rb3,zb3,rb4,zb4,
    logical :: lost
    call find_func_cputime_1_of_2  !-------------------1/2

    v2_lim=(0.1*c)**2 !when ion velocity >0.1c (light velocity), make it escape.

    N_lost_particles=0
    ip=1
    do while(ip<=n_active)
        rtp=sqrt(x(ip,1)**2+x(ip,2)**2) !+1e-6;
        vtp2=v(ip,1)**2+v(ip,2)**2+v(ip,3)**2
        ztp=x(ip,3)
        lost=.false.
        
        if(mass_q_i_05*vtp2>life_and_ek(ip,2)) life_and_ek(ip,2)=mass_q_i_05*vtp2
        
        if     (  ztp<=zs )then
            call register_particle_loss(1)
            lost=.true.
        elseif (  ztp>=zd )then
            call register_particle_loss(2)
            lost=.true.
        elseif (  rtp>=rp .and. (ztp>zs .or. ztp<zd) )then
            call register_particle_loss(3)
            lost=.true.
        endif        
        if(.not. lost)then
            if(  isnan(x(ip,1)) .or. isnan(x(ip,2)) .or. isnan(x(ip,3)) .or.  vtp2>v2_lim .or. &
                & isnan(v(ip,1)) .or. isnan(v(ip,2)) .or. isnan(v(ip,3))  )then
                call register_particle_loss(4)
                lost=.true.
            endif
        endif

        if(lost)then
            call remove_particle(ip)
        else
            ip=ip+1
        endif
    enddo
    call ensure_particle_capacity(n_active+deltaN)
    do ip=n_active+1,n_active+deltaN
        call injection_inlet_sub(ip)
    enddo
    n_active=n_active+deltaN
    np_e=n_active
    num_inject=deltaN
    call find_func_cputime_2_of_2(func_time(2))  !-----2/2
end subroutine escape_and_inject


subroutine remove_particle(ip_remove)
    USE the_whole_varibles
    implicit none
    integer*4 :: ip_remove

    if(ip_remove<n_active)then
        x(ip_remove,1:3)=x(n_active,1:3)
        v(ip_remove,1:3)=v(n_active,1:3)
        v_e(ip_remove,1:3)=v_e(n_active,1:3)
        t_np(ip_remove)=t_np(n_active)
        life_and_ek(ip_remove,1:2)=life_and_ek(n_active,1:2)
        x_to_grid(ip_remove,1:2)=x_to_grid(n_active,1:2)
        ir1_iz1_grid(ip_remove,1:2)=ir1_iz1_grid(n_active,1:2)
    endif
    n_active=n_active-1
end subroutine remove_particle

subroutine register_particle_loss(iloss)
    USE the_whole_varibles
    implicit none
    INTEGER*2::iloss
    real*8::rtp,delat_t
    !when the ion escapes and then inject a new particle
    !replace the information of old particle (ip) with new born particle  
    !Particle energy loss at the boundary is counted
    
    if (rec_xv_loss == 1) then
        ip_loss_rec=ip_loss_rec+1
        if(ip_loss_rec>np_loss_rec)ip_loss_rec=1
        xv_loss(ip_loss_rec,1:3)=x(ip,1:3)
        xv_loss(ip_loss_rec,4:6)=v(ip,1:3)
        xv_loss(ip_loss_rec,7)=iloss
    end if
    
    ip_loss_for_tau_p=ip_loss_for_tau_p+1
    if(t-t1_taup>1e-6)then
        tau_p_ave=n_active*(t-t1_taup)/(ip_loss_for_tau_p+1e-3) 
        
        t1_taup=t
        ip_loss_for_tau_p=0
    endif
 
    Ek_loss_tol(iloss)=Ek_loss_tol(iloss)+mass_q_i_05*sum(v(ip,1:3)**2)
    num_loss(iloss)=num_loss(iloss)+1
    N_lost_particles=N_lost_particles+1

end subroutine register_particle_loss

subroutine ensure_particle_capacity(required_capacity)
    USE the_whole_varibles
    implicit none
    integer*4 :: required_capacity, new_capacity
    integer, allocatable :: np_tmp(:), ir_tmp(:,:)
    real*8, allocatable :: x_tmp(:,:), v_tmp(:,:), v_e_tmp(:,:), t_tmp(:), xg_tmp(:,:), life_tmp(:,:)

    if(required_capacity<=n_capacity)return

    new_capacity=max(required_capacity,n_capacity+max(deltaN,1)*100)

    allocate(np_tmp(1:new_capacity))
    allocate(ir_tmp(1:new_capacity,1:2))
    allocate(x_tmp(1:new_capacity,1:3))
    allocate(v_tmp(1:new_capacity,1:3))
    allocate(v_e_tmp(1:new_capacity,1:3))
    allocate(t_tmp(1:new_capacity))
    allocate(xg_tmp(1:new_capacity,1:2))
    allocate(life_tmp(1:new_capacity,1:2))

    np_tmp=0
    ir_tmp=0
    x_tmp=0.0d0
    v_tmp=0.0d0
    v_e_tmp=0.0d0
    t_tmp=0.0d0
    xg_tmp=0.0d0
    life_tmp=0.0d0

    np_tmp(1:n_active)=np(1:n_active)
    ir_tmp(1:n_active,1:2)=ir1_iz1_grid(1:n_active,1:2)
    x_tmp(1:n_active,1:3)=x(1:n_active,1:3)
    v_tmp(1:n_active,1:3)=v(1:n_active,1:3)
    v_e_tmp(1:n_active,1:3)=v_e(1:n_active,1:3)
    t_tmp(1:n_active)=t_np(1:n_active)
    xg_tmp(1:n_active,1:2)=x_to_grid(1:n_active,1:2)
    life_tmp(1:n_active,1:2)=life_and_ek(1:n_active,1:2)

    deallocate(np,ir1_iz1_grid,x,v,v_e,t_np,x_to_grid,life_and_ek)
    allocate(np(1:new_capacity))
    allocate(ir1_iz1_grid(1:new_capacity,1:2))
    allocate(x(1:new_capacity,1:3))
    allocate(v(1:new_capacity,1:3))
    allocate(v_e(1:new_capacity,1:3))
    allocate(t_np(1:new_capacity))
    allocate(x_to_grid(1:new_capacity,1:2))
    allocate(life_and_ek(1:new_capacity,1:2))

    np=np_tmp
    ir1_iz1_grid=ir_tmp
    x=x_tmp
    v=v_tmp
    v_e=v_e_tmp
    t_np=t_tmp
    x_to_grid=xg_tmp
    life_and_ek=life_tmp
    n_capacity=new_capacity
end subroutine ensure_particle_capacity

subroutine find_power
    use the_whole_varibles
    implicit none
    real*8 :: Ek_i_total, Ek_e_total, Ek_total_current_joules
    real*8 :: Ek_loss_interval_joules, energy_gain_interval

    if(mod(t,10.*trf)<dt)then
        Ek_loss_interval_joules = sum(Ek_loss_tol(1:4)) * n_macro * qe_abs

        power_loss_ave(1:4)=n_macro*qe_abs*Ek_loss_tol(1:4)/(10.*trf)
        Ek_loss_ave(1:4)=Ek_loss_tol(1:4)/num_loss(1:4)
        rec_t(30:33)=num_loss(1:4)
        Ek_loss_cumulative = Ek_loss_cumulative + sum(Ek_loss_tol(1:4))
        
        Ek_i_total = 0.0
        do ip=1,n_active
            Ek_i_total = Ek_i_total + mass_q_i_05*sum(v(ip,1:3)**2)
        enddo
        Ek_e_total = mass_q_e_05 * sum(v_e(1:n_active,1:3)**2)
        Ek_total_current_joules = (Ek_i_total + Ek_e_total) * n_macro * qe_abs

        energy_gain_interval = (Ek_total_current_joules - Ek_total_last_joules) &
                             + Ek_loss_interval_joules
        
        if (10.*trf > 1.0e-12) then
            absorbed_power = energy_gain_interval / (10.*trf)
        else
            absorbed_power = 0.0
        endif

        Ek_total_last_joules = Ek_total_current_joules

        num_loss=1e-2
        Ek_loss_tol=0.
        
        call find_irf_now
    endif
end subroutine find_power

!subroutine filter_1d(ndata,data_in)
    !implicit none
    !integer::ndata,idata
    !real*8::w,data_in(1:ndata),data_out(1:ndata)
    !w=0.5
    !
    !idata=1;      data_out(idata)=.5*(data_in(idata)+data_in(idata+1))
    !idata=ndata;  data_out(idata)=.5*(data_in(idata)+data_in(idata-1))
    !do idata=2,ndata-1
    !    data_out(idata)=.25*(data_in(idata-1)+2.*data_in(idata)+data_in(idata+1))
    !enddo
    !data_in=data_out
!end subroutine filter_1d


!subroutine filter_2d(nrdata,nzdata,data_in)
    !implicit none
    !integer::nrdata,nzdata,irdata,izdata
    !real*8::data_in(1:nrdata,1:nzdata),data_out(1:nrdata,1:nzdata)
    !
    !data_out=data_in
    !do irdata=2,nrdata-1
    !    do izdata=2,nzdata-1
    !        data_out(irdata,izdata)=1./6.*(2.*data_in(irdata,izdata)+data_in(irdata-1,izdata)+&
    !            & data_in(irdata+1,izdata)+data_in(irdata,izdata-1)+data_in(irdata,izdata+1))
    !    enddo
    !enddo
    !irdata=1 ;     data_out(irdata,2:nzdata-1)=0.25*data_in(irdata,1:nzdata-2)+&
    !    &        .5*data_in(irdata,2:nzdata-1)+0.25*data_in(irdata,3:nzdata)
    !irdata=nrdata; data_out(irdata,2:nzdata-1)=0.25*data_in(irdata,1:nzdata-2)+&
    !    &       0.5*data_in(irdata,2:nzdata-1)+0.25*data_in(irdata,3:nzdata)
    !izdata=1 ;     data_out(2:nrdata-1,izdata)=0.25*data_in(1:nrdata-2,izdata)+&
    !    &       0.5*data_in(2:nrdata-1,izdata)+0.25*data_in(3:nrdata,izdata)
    !izdata=nzdata; data_out(2:nrdata-1,izdata)=0.25*data_in(1:nrdata-2,izdata)+&
    !    &       0.5*data_in(2:nrdata-1,izdata)+0.25*data_in(3:nrdata,izdata)
    !
    !data_in(1:nrdata,1:nzdata)=data_out(1:nrdata,1:nzdata)
!end subroutine filter_2d

!subroutine interp_1d(n_bef,n_aft,x_bef,x_aft,data_bef,data_aft)
    !implicit none
    !integer::n_bef,n_aft
    !real*8::x_bef(1:n_bef),x_aft(1:n_aft),data_bef(1:n_bef),data_aft(1:n_aft)
    !integer::I_x,ix_aft,ix_bef1,ix_bef2
    !real*8::s1,s2,dx_bef,min_x,xtp,dx_tp1,dx_tp2
    !
    !dx_bef=x_bef(2)-x_bef(1)
    !do ix_aft=1,n_aft
    !    xtp=x_aft(ix_aft)
    !    I_x=minloc(abs(xtp-x_bef),1);
    !    min_x=xtp-x_bef(I_x);
    !    if (min_x<0)then
    !        ix_bef1=I_x-1;
    !        ix_bef2=I_x;
    !    else
    !        ix_bef1=I_x;
    !        ix_bef2=I_x+1;
    !    endif
    !    if(ix_bef1<1)ix_bef1=1
    !    if(ix_bef2>n_bef)ix_bef2=n_bef
    !    dx_tp1=xtp-x_bef(ix_bef1);
    !    dx_tp2=x_bef(ix_bef2)-xtp;
    !    s1=dx_tp1/dx_bef;
    !    s2=1.-s1;
    !    data_aft(ix_aft)=s2*data_bef(ix_bef1)+s1*data_bef(ix_bef2);
    !enddo
    !continue
!end subroutine interp_1d

!subroutine interp_0d(n_bef,x_bef,x_aft,data_bef,data_aft)
    !implicit none
    !integer::n_bef
    !real*8::x_bef(1:n_bef),x_aft,data_bef(1:n_bef),data_aft
    !integer::I_x,ix_bef1,ix_bef2
    !real*8::s1,s2,dx_bef,min_x,xtp,dx_tp1,dx_tp2
    !dx_bef=x_bef(2)-x_bef(1)
    !
    !xtp=x_aft
    !I_x=minloc(abs(xtp-x_bef),1);
    !min_x=xtp-x_bef(I_x);
    !if (min_x<0)then
    !    ix_bef1=I_x-1;
    !    ix_bef2=I_x;
    !else
    !    ix_bef1=I_x;
    !    ix_bef2=I_x+1;
    !endif
    !if(ix_bef1<1)ix_bef1=1
    !if(ix_bef2>n_bef)ix_bef2=n_bef
    !dx_tp1=xtp-x_bef(ix_bef1);
    !dx_tp2=x_bef(ix_bef2)-xtp;
    !s1=dx_tp1/dx_bef;
    !s2=1.-s1;
    !data_aft=s2*data_bef(ix_bef1)+s1*data_bef(ix_bef2);
    !continue
!end subroutine interp_0d
