!ifx hypic.f90 fun_b0.f90 fun_record_display.f90 fun_grid.f90 fun_ini.f90 fun_antenna_irf.f90 fun_particles.f90 fun_mcc.f90 fun_fdfd.f90 -qmkl -o zcal.x


module types
    Type sparse
        Integer :: nzmax
        Integer :: n
        Integer,allocatable :: row(:)
        Integer,allocatable :: col(:)
        real*8,allocatable :: val(:)
    end Type sparse
    Type sparse_complex
        Integer :: nzmax
        Integer :: n
        Integer,allocatable :: row(:)
        Integer,allocatable :: col(:)
        Complex*16,allocatable :: val(:)
    end Type sparse_complex
end module types

module the_whole_varibles
    use types
    INTEGER*4 ::np_max
    ! real*8,parameter :: frequency=600e3 !@namelist
    real*8 :: frequency
    ! real*8,parameter :: trf=1/frequency; !@namelist indirect
    real*8 :: trf
    real*8,parameter :: frequency2=20e6

    real*8,parameter:: qe0=1.6022e-19,mp=1.6726231D-27;
    ! real*8,parameter:: mass_number   =40.0; !@namelist
    real*8 :: mass_number
    real*8,parameter:: charger_number=1.0;
    ! real*8,parameter:: mi=mass_number*mp !@namelist indirect
    real*8 :: mi
    real*8,parameter:: qi=charger_number*qe0
    ! real*8,parameter:: mass_q_i_05=0.5*mi/qi !@namelist indirect
    real*8 :: mass_q_i_05
    ! real*8,parameter:: q_mass_i=qi/mi !@namelist indirect
    real*8 :: q_mass_i
    INTEGER,allocatable :: np(:)
    INTEGER :: i_s


    !---------------constants---------------------!
    real*8,parameter::c=3e8,me=9.1e-31, qe_abs=qe0,qe=-qe_abs,  pi = 3.1415926,kB=1.3806505e-23;
    real*8,parameter:: epson0=8.8542e-12 ,mu0=pi*(4e-7),cm3=1e6
    real*8,parameter:: q_mass_e=qe/me,mass_q_e_05=0.5*me/qe_abs;
    real*8,parameter::sig_en=15e-20,sig_in=1e-18;
    complex*16,parameter::i=(0.0D0,1.0D0)
    !ln_A0 Coulomb logarithm;
    real*8,parameter::ln_A0=15.;
    real*8,parameter::qe_4=qe_abs**4,ep_2=epson0**2
    real*8,parameter::coeff_Te_1=(qe_abs**4*ln_A0)/(8.*sqrt(2.)*pi*ep_2*sqrt(me)), coeff_Te_4=sqrt(8.*qe_abs/(pi*me))
    real*8 :: coeff_Te_2, coeff_Te_3

    !maxwell_FDFD
    integer*4::maxwellcount=0
    type(sparse_complex) :: A_helicon
    real*8  :: dt_run_fdfd,time_call_fdfd=0.
    integer :: i3,ith,m,m_start,m_end,m_delta,iswitch_antenna_type,switch_z_boundary
    Integer :: switch_source_type,switch_plasma_vacuum_boundary,switch_divE
    Integer :: switch_boundary_half_variables,switch_r_0_boundary
    Integer :: Region,jieshu,iswitch_output_Region,nzmax_pt=0
    Integer, allocatable :: FindRegion(:,:)
    Integer(kind=8) :: count1,count2,count0,count,count_rate,count_max
    real*8 :: om,omc2,ds,time_Maxwell,ptotal,ptotal_inside,ptotal_boundary,Res_p,m_real
    real*8, external :: DBesselJZero
    real*8, allocatable :: ni(:,:),b0_DC(:,:,:),te_in_FDFD(:,:)
    real*8, allocatable ::fvm(:,:),ptotm(:,:,:),power_depo_rthz(:,:,:),th(:),ptotal_m(:),max_eth(:)
    real*8, allocatable :: power_depo_species_m(:,:,:,:)
    complex*16:: im,imp,iom,m_max
    complex*16, allocatable ::eq_b(:),eq_xlast(:,:),ja_m_nor(:,:,:), ja_record(:,:,:,:),ja_AC(:,:,:,:)
    complex*16, allocatable :: e_m(:,:,:),e_record(:,:,:,:),e_AC(:,:,:,:),e_int(:,:,:)
    complex*16, allocatable :: jp_m(:,:,:), jp_record(:,:,:,:),jp_AC(:,:,:,:)
    complex*16, allocatable :: b_m(:,:,:), b_record(:,:,:,:),b_AC(:,:,:,:)
    complex*16, allocatable ::ep(:,:,:),si(:,:,:,:),power_depo_2D(:,:), Xkz_e(:,:,:)
    complex*16, allocatable :: si_species(:,:,:,:,:)
    Complex*16, allocatable :: e_output(:,:,:,:),Erf6(:,:,:)

    !initial
    integer:: iswitch_display,i_switch_B0,iswitch_dielectric,i_switch_power_mode,i_switch_MHD_Te,use_file_defined_density
    integer:: export_real_B0
    integer :: iswitch_v_closure_type
    integer :: k_closure_type = 1
    real*8 :: bz_max,ti_ini,te_ini,pn,B0_correction_factor,B_resonance,t_power_on,nn,max_density_set,ni_ave_ratio
    real*8 :: power,i_now,power_Joule,resistance,irf_set,irf2_set,t_power_end,Te_set,k_constant=20.0d0
    real*8 :: width_ne_r,width_ne_z,z_inject_cener ,r_inject_cener,width_ne_source_r,width_ne_source_z

    !!antenna and grid
    integer,parameter :: ns=3,nth=200,n3=3
    integer iswitch_RF2
    integer::nr,nz,n_vac,nra,nza,i_vac,nz6
    integer*4 :: nzx2,nrd,nzs,nzd,nrp,itp,ir,iz,ierror
    integer:: Nrz_shield_region1(1:4),Nrz_shield_region2(1:4)
    integer, allocatable :: nr_met(:),nz_met(:),nr_vac(:),nz_vac(:)
    integer, allocatable::i_plasma_region(:,:)
    real*8 :: irf_mean,n_mean=1
    real*8 :: rp,zp,d2r,d2z
    real*8 :: dr,dz,rs,rd,rl,zs,zd,zl,ra,za,r1,drz,dr22,dz22,dr_2,dz_2
    real*8 :: tp1,tp2,tp3
    real*8, allocatable :: r(:),z(:),r2(:)
    real*8, allocatable :: r_vac(:),r_met(:) ,z_vac(:)

    !display and record
    character*10:: time_char(3,2),timefunc_char(3,2)
    integer :: date_time(8,2), timefunc_int(8,2),func_time_count(8,2),i_count=0;
    real*8 :: timefunc_int_real(8,2),time_run
    real*8 :: t=0,dt,td,run_time_0=0,run_time_1
    Real*8 :: func_time(1:10)=0, time_diff(1:2)=0,fun_te_time(1:10)=0,time_te_diff(1:2)=0
    Real*8 :: tau_E_ave=1e-6,tau_p_ave=1e-6
    integer*4:: n_pic,ip_loss_rec=0 ,ip_loss_for_tau_p=0
    integer:: np_loss_rec
    !real*8:: xv_loss(1:np_loss_rec,1:7)=0.
    real*8,allocatable :: xv_loss(:,:)
    real*8:: t1_taup=0.
    real*8, allocatable :: life_and_ek(:,:)
    real*8 :: Ek_total_last_joules=0.0d0, Ek_total_current_joules=0.0d0
    real*8 :: thrust_source=0.0d0, thrust_total=0.0d0
    real*8 :: ISP_global=0.0d0, ISP_plume=0.0d0
    real*8 :: perf_exit_momentum_z=0.0d0, perf_inject_momentum_z=0.0d0
    real*8 :: perf_exit_count=0.0d0, perf_inject_count=0.0d0, perf_t0=0.0d0
    real*8 :: loading_stability_data(1:3,1:20)=0.0d0
    real*8 :: loading_cv_percent=0.0d0, loading_s_percent=0.0d0
    integer*4 :: loading_stability_count=0

    !particles
    integer*4,parameter:: iseed=1234567,nrec_t=42,nrec_p=20,te_update_interval=1000
    integer,parameter::nz_length=120
    integer,allocatable :: ir1_iz1_grid(:,:)
    integer*4::ifig=0
    integer*4::num_inject,it=0,ip,np_e,N_lost_particles
    real*8:: vtp(1:3),dseed=123457.D0, s_drz
    real*8:: coeff1,t_cyclotron,vi_ex,ve_ex,dt_trf,n_macro

    !v: ion velocity; ve: electron velocity; 1:3 -- x y z direction
    !real*8:: x(1:np_max,1:3)
    real*8, allocatable :: x(:,:)
    !real*8:: v(1:np_max,1:3)
    real*8, allocatable :: v(:,:)
    !real*8:: v_e(1:np_max,1:3)
    real*8, allocatable :: v_e(:,:)
    !real*8:: t_np(1:np_max)
    real*8, allocatable :: t_np(:)
    !real*8:: x_to_grid(1:np_max,1:2)
    real*8, allocatable :: x_to_grid(:,:)

    real*8:: rec_t(nrec_t)=0.,rec_p(nrec_p)=0.,vz0,state_power_on_off=0.
    real*8:: Ek_loss_ave(1:4)=0,power_loss_ave(1:4)=0,Ek_loss_tol(1:4)=0.,num_loss(1:4)=1e-2
    real*8 :: Ek_inject_tol=0.0d0
    real*8, allocatable :: pdf_ne_r(:),pdf_ne_source_r(:),pdf_ne_z(:),pdf_ne_source_z(:),r_particle_inj(:)
    Real*8:: t_rec_evolution,t_rec_trajectory,rz_ant_region(1:4),t_rec_profiles
    real*8, allocatable :: density_2D(:,:),Es_2D(:,:,:),Ek_ion_2D(:,:),Ek_ion_2D_r(:,:)
    real*8, allocatable :: te_2D(:,:),te_mhd(:,:),u_mhd(:,:,:),Q_ie(:,:),u_pic(:,:,:)
    complex*16::expt,expt2

    !fun_Te
    integer nrun
    real*8 ::Da_perp,Da,mue_i_brz,coeff_1_br2,coeff_1_bz2,coeff_1_mue2_br2
    real*8 ::coeff_1_mue2_bz2,coeff_1_mue2_b2,de

    !mcc
    INTEGER*2::iswitch_mcc
    Real*8:: s_tab(1:22),a_tab(1:22)
    Real*8:: ln_A=10.9,dt_mcc(1:4),time_to_run_mcc(1:4)=0.,dt_over_dt_mcc(1:4)=1e-2,ne_mcc,mu_ab,sab_part
    Real*8:: density_ave,te_ave,ti_ave,Ek_i_ave,Ek_e_ave,tau_ii,tau_ie,tau_ee,tau_ei
    Real*8:: Ek_i_ave_ini,Ek_e_ave_ini,absorbed_power=0.0d0,Ek_total_initial_joules=0.0d0,Ek_loss_cumulative=0.0d0

    ! change
    integer*4 :: iswitch_analytic_profile
    integer*4 :: iswitch_T_closure_type
    real*8, allocatable :: ni_midVal(:,:)
    real*8, allocatable :: ni_read(:,:)
    real*8, allocatable :: kap_ll_read_raw(:,:)
    complex*16, allocatable :: kap_ll_read(:,:)
    real*8 :: till_uni
    real*8 :: tipp_uni
    real*8 :: tell_uni
    real*8 :: t2pp_uni
    real*8 :: TminLim
    ! change

    character(len=256) :: nmlfile
    character(len=256) :: outputDir
    character(len=256) :: b0Dir
    character(len=256) :: kapDir
    character(len=256) :: plasmaLoadDir

    integer*4 :: mTruncation

    integer :: iswitch_specific_mode
    integer :: m_specific

    real*8 :: vl_in_FDFD, alpha_v

    integer*4 :: iswitch_log
    integer*4 :: iswitch_play_trajectory

    integer*4 :: dev_Only4CalMag

    integer*4 :: runHowManyPeriods
    integer*4 :: runRFonWhichPeriod

    real*8 :: Mig
    real*8 :: deltaN_real
    integer*4 :: deltaNtype
    integer*4 :: deltaN
    integer*4 :: CrCzType
    real*8 :: Cr_guess,Cz_guess,Cr,Cz
    real*8 :: n0_set
    real*8 :: ntar_set

    integer*4 :: n_active
    integer*4 :: n_capacity

    integer*4 :: rec_xv_loss
    integer*4 :: rec_plasma
    integer*4 :: rec_rf_field
    integer*4 :: rec_Ek_ave_z
    integer*4 :: rec_ions_rz
    integer*4 :: rec_density_Es_2D
    integer*4 :: rec_Erf_all_m
    integer*4 :: rec_dielectric_stat
    integer*4 :: rec_deposition


    !integer*4 :: recDensity
    !integer*4 :: recDeposition
    !integer*4 :: recTemperature
    !integer*4 :: recRFfield
    !integer*4 :: recEs
    !integer*4 :: recParticleEBF
    !integer*4 :: recEkDensity_PeriodAveraged
    !integer*4 :: recIons
    !integer*4 :: recDensityEs2D
    !integer*4 :: recErf

    integer*4 :: dev_Only4CalE

    integer*4 :: iswitch_vacum_dielectric_type
    
    real*8,allocatable :: nu_depo(:,:,:), u2_depo(:,:,:), t_depo(:,:,:), n_depo(:,:), vl_closure_prev(:,:)
    real*8,allocatable :: nu_depo_in_FDFD(:,:,:), t_in_dielectric(:,:,:)
    integer*4 :: statistic_trigger=0, statistic_acc_number=0, statistic_ready=0
end module the_whole_varibles

Program hypic
    use the_whole_varibles
    implicit none
    integer :: argc

    argc = command_argument_count()
    if (argc < 1) then
        print *, "Usage: app.exe input.nml"
        stop
    end if
    call get_command_argument(1, nmlfile)

    call start_ftime

    call initialize
    
    if(iswitch_display/=0) write(*,*)'Fun_ini finished. Fun_fdfd is running.'
    call maxwell_FDFD
    if(iswitch_RF2==1 )call maxwell_FDFD
    call particles_initialization
    call power_ini

    it=0;t=0;
    if(iswitch_display/=0) write(*,*)'Main loop is running.'
    do while(1)
        it=it+1; t=t+dt

        if(mod(t,dt_run_fdfd)<dt  .and. t>t_power_on .and. t<t_power_end )then
            call prepareStatisticForFDFD
            call set_ne_Te
            call maxwell_FDFD
            if(iswitch_RF2==1 )call maxwell_FDFD
        endif

        call escape_and_inject
        call accumulateStatisticBeforeFDFD
        call mover
        if(iswitch_mcc/=0)call mcc_main

        if(i_switch_MHD_Te==1 .and. mod(it, te_update_interval) == 0) then
        call find_te
        endif

        call find_power      
        call rec_time_ave     
        call rec_time_trajectory        
        call record_profiles
        call display_main
        if( t>td ) exit
    enddo

    ! call resetParameters
    ! call maxwell_FDFD

    close(42)
    close(43)
    call display_main
    if(iswitch_play_trajectory==1) call play_trajectory
    write(*,*) 'code running has finished !!! '
End Program hypic
