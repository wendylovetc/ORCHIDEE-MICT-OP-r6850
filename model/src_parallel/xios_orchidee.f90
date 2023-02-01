! ================================================================================================================================
!  MODULE       : xios_orchidee
!
!  CONTACT      : orchidee-help _at_ listes.ipsl.fr
!
!  LICENCE      : IPSL (2006)
!  This software is governed by the CeCILL licence see ORCHIDEE/ORCHIDEE_CeCILL.LIC
!
!>\BRIEF   This module contains the initialization and interface to the XIOS code.
!!
!!\n DESCRIPTION: This module contains the interface for the use of the XIOS code. All call to XIOS are done in this module.
!!
!!                Summury of subroutines
!!                      xios_orchidee_comm_init       : First call to XIOS to get the MPI communicator 
!!                      xios_orchidee_init            : Initialize variables needed for use of XIOS 
!!                                                      Deactivation of fields not calculated due specific run options
!!                      xios_orchidee_update_calendar : Update the calandar in XIOS
!!                      xios_orchidee_finalize        : Last call to XIOS for finalization
!!                      xios_orchidee_send_field      : Interface to send fields with 1, 2 or 3 dimensions to XIOS
!!                      xios_orchidee_send_field_r1d  : Internal subroutine for 1D(array) fields
!!                      xios_orchidee_send_field_r2d  : Internal subroutine for 2D fields
!!                      xios_orchidee_send_field_r3d  : Internal subroutine for 3D fields
!!
!!                It is only possible to use XIOS2. Note that compilation must be done with the preprocessing key XIOS 
!!                and CPP_PARA. Compiling without these keys makes it impossible to activate XIOS. 
!!                To activate running using XIOS, the flag XIOS_ORCHIDEE_OK=y must be set in run.def and the file iodef.xml must exist.  
!!
!! RECENT CHANGE(S): Created by Arnaud Caubel(LSCE), Josefine Ghattas (IPSL) 2013
!!                   Removed possibility to use XIOS1, 21/10/2016
!!
!! REFERENCE(S) : None
!!
!! SVN          :
!! $HeadURL$
!! $Date$
!! $Revision$
!! \n
!_ ================================================================================================================================

MODULE xios_orchidee

#ifdef XIOS
  USE xios
#endif
  USE defprec
  USE pft_parameters_var, ONLY : nvm, nvmap
  USE constantes_var
  USE time, ONLY : dt_sechiba
  USE constantes_soil_var, ONLY : nstm, check_waterbal, diaglev, check_cwrr2, ok_freeze_cwrr, ndeep
  USE vertical_soil_var, ONLY : ngrnd, nslm
  USE IOIPSL, ONLY : ioget_calendar, ju2ymds
  USE mod_orchidee_para_var
  USE mod_orchidee_transfert_para
  USE ioipsl_para
  USE grid_var, ONLY : GridType

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xios_orchidee_comm_init, xios_orchidee_init, xios_orchidee_change_context, &
            xios_orchidee_update_calendar, xios_orchidee_context_finalize, xios_orchidee_finalize, &
            xios_orchidee_send_field, xios_orchidee_recv_field

  !
  !! Declaration of public variables
  !
  LOGICAL, PUBLIC, SAVE           :: xios_orchidee_ok=.TRUE.     !! Use XIOS for diagnostic files
  !$OMP THREADPRIVATE(xios_orchidee_ok)
  REAL(r_std), PUBLIC, SAVE       :: xios_default_val=0          !! Default value (missing value) used in XIOS. The value 0 will be overwritten with the value taken from XIOS.
  !$OMP THREADPRIVATE(xios_default_val)

  !
  !! Declaration of internal variables
  !
#ifdef XIOS
  TYPE(xios_context)              :: ctx_hdl_orchidee      !! Handel for ORCHIDEE
  !$OMP THREADPRIVATE(ctx_hdl_orchidee)
#endif
  CHARACTER(len=*),PARAMETER      :: id="client"           !! Id for initialization of ORCHIDEE in XIOS



  !! ==============================================================================================================================
  !! INTERFACE   : xios_orchidee_send_field
  !!
  !>\BRIEF         Send a field to XIOS.
  !!
  !! DESCRIPTION  :\n Send a field to XIOS. The field can have 1, 2 or 3 dimensions.
  !!                  This interface should be called at each time-step for each output varaiables.
  !!
  !! \n
  !_ ================================================================================================================================
  INTERFACE xios_orchidee_send_field
     MODULE PROCEDURE xios_orchidee_send_field_r1d, xios_orchidee_send_field_r2d, xios_orchidee_send_field_r3d, &
                      xios_orchidee_send_field_r4d, xios_orchidee_send_field_r5d
  END INTERFACE

  INTERFACE xios_orchidee_recv_field
     MODULE PROCEDURE xios_orchidee_recv_field_r1d, xios_orchidee_recv_field_r2d, xios_orchidee_recv_field_r3d
  END INTERFACE


CONTAINS
  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_comm_init 
  !!
  !>\BRIEF         Get the MPI communicator.
  !!
  !! DESCRIPTION  :\n First call to XIOS to get the MPI communicator. 
  !!                  Note that it is XIOS that initialize the MPI communicator.
  !!                  This subroutine is only called in ORCHIDEE offline mode. When running in coupled mode, the 
  !!                  atmospheric model must initlialize XIOS at the same time as initializing MPI. 
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_comm_init(comm_local)
    !
    !! 0. Variable and parameter declaration
    !
    !!    Output variables
    INTEGER, INTENT(OUT) :: comm_local

    !_ ================================================================================================================================

    IF (is_omp_root) THEN
#ifdef XIOS
       CALL xios_initialize(id,return_comm=comm_local)
#else
       CALL ipslerr_p(3, 'xios_orchidee_comm_init', 'Preprocessing key XIOS is missing to run ORCHIDEE with XIOS', &
            'Recompile with preprocessing flag XIOS or set XIOS_ORCHIDEE_OK=n in run.def','')
#endif
    END IF
  END SUBROUTINE xios_orchidee_comm_init


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_init 
  !!
  !>\BRIEF         Initialize variables needed for use of XIOS.
  !!
  !! DESCRIPTION  :\n Initialization of specific varaiables needed to use XIOS such as model domain and time step. 
  !!
  !!                  In this subroutine also a section containg deactivation of some fields is found. The variables are 
  !!                  deactivated of not according to the corresponding control flag. For exemple the variables cacluated by the 
  !!                  routing scheme will be deactivated if the routing is deactivated. This is done to be able to keep the same 
  !!                  iodef.xml input file for several options without geting empty fields in the output file. Note that a field that
  !!                  is activated in the code can always be deactivated from the iodef.xml external file. 
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_init(MPI_COMM_ORCH,                   &
       date0,    year,      month,             day, julian_diff, &
       lon_mpi,  lat_mpi,   soilth_lev )

    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(in)                            :: MPI_COMM_ORCH    !! Orchidee MPI communicator (from module mod_orchidee_mpi_data)
    REAL(r_std), INTENT(in)                               :: date0            !! Julian day at first time step
    INTEGER(i_std), INTENT(in)                            :: year, month, day !! Current date information
    REAL(r_std), INTENT(in)                               :: julian_diff      !! Current day in the year [1,365(366)]
    REAL(r_std),DIMENSION (iim_g,jj_nb), INTENT(in)       :: lon_mpi, lat_mpi !! Longitudes and latitudes on MPI local domain 2D domain
    REAL(r_std),DIMENSION (ngrnd), INTENT(in)             :: soilth_lev       !! Vertical soil levels for thermal scheme (m)
    !
    !! 0.2 Local variables
    !
#ifdef XIOS

    TYPE(xios_duration)            :: dtime_xios
    TYPE(xios_date)                :: start_date
    TYPE(xios_date)                :: time_origin
    TYPE(xios_fieldgroup)          :: fieldgroup_handle
    TYPE(xios_field)               :: field_handle
    TYPE(xios_file)                :: file_handle
#endif
    INTEGER(i_std)                 :: i
    INTEGER(i_std)                 :: year0, month0, day0 !! Time origin date information
    REAL(r_std)                    :: sec0                !! Time origin date information
    CHARACTER(LEN=20)              :: calendar_str        !! Name of current calendar
    CHARACTER(LEN=30)              :: start_str           !! Current date as character string
    CHARACTER(LEN=30)              :: startorig_str       !! Time origin date as character string
    !_ ================================================================================================================================
    
    
    IF (printlev>=3) WRITE(numout,*) 'Entering xios_orchidee_init'

    !Config Key   = XIOS_ORCHIDEE_OK
    !Config Desc  = Use XIOS for writing diagnostics file
    !Config If    = 
    !Config Def   = y 
    !Config Help  = Compiling and linking with XIOS library is necessary. 
    !Config Units = [FLAG]
    CALL getin_p('XIOS_ORCHIDEE_OK',xios_orchidee_ok)
    IF (printlev>=1) WRITE(numout,*)'In xios_orchidee_init, xios_orchidee_ok=',xios_orchidee_ok

    ! Coherence test between flag and preprocessing key
#ifndef XIOS
    IF (xios_orchidee_ok) THEN
       CALL ipslerr_p(3,'xios_orchidee_init', 'Preprocessing key XIOS is missing to run ORCHIDEE with XIOS',&
            'Recompile with preprocessing flag XIOS or set XIOS_ORCHIDEE_OK=n in run.def', '')
    END IF
#endif


    !
    !! 1. Set date and calendar information on the format needed by XIOS
    !

    ! Get the calendar from IOIPSL and modify the string to correspond to what XIOS expects
    CALL ioget_calendar(calendar_str)

    IF (calendar_str == 'gregorian') THEN
       calendar_str='gregorian'
    ELSE IF (calendar_str == 'noleap') THEN
       calendar_str='noleap'
    ELSE IF (calendar_str == '360d') THEN
       calendar_str='d360'
    END IF

    ! Transform the time origin from julian days into year, month, day and seconds
    CALL ju2ymds(date0, year0, month0, day0, sec0)



    IF (xios_orchidee_ok .AND. is_omp_root) THEN
#ifdef XIOS
       !
       !! 2. Context initialization
       !
       CALL xios_context_initialize("orchidee",MPI_COMM_ORCH)
       CALL xios_get_handle("orchidee",ctx_hdl_orchidee)
       CALL xios_set_current_context(ctx_hdl_orchidee)

       !
       !! 2. Calendar, timstep and date definition
       !
       dtime_xios%second=dt_sechiba

       CALL xios_define_calendar(type=calendar_str, start_date=xios_date(year,month,day,0,0,0), &
            time_origin=xios_date(year0,month0,day0,0,0,0), timestep=dtime_xios)

       !
       !! 3. Domain definition
       !
       ! Global domain
       CALL xios_set_domain_attr("domain_landpoints", ni_glo=iim_g, nj_glo=jjm_g)

       ! Local MPI domain
       IF ( GridType == "RegLonLat" ) THEN
          CALL xios_set_domain_attr("domain_landpoints",type="rectilinear", ibegin=0, ni=iim_g, jbegin=jj_begin-1, nj=jj_nb)
       ELSE IF ( GridType == "RegXY" ) THEN
          CALL xios_set_domain_attr("domain_landpoints",type="curvilinear", ibegin=0, ni=iim_g, jbegin=jj_begin-1, nj=jj_nb)
       ELSE
          WRITE(numout,*) 'Following GridType is not supported: GridType =', GridType
          CALL ipslerr_p(3, 'xios_orchidee_init', 'GridType not yet supported.','Problem for defining local MPI domain','')
       ENDIF

       ! Define how data is stored on memory : 1D array for only continental points
       CALL xios_set_domain_attr("domain_landpoints",data_dim=1, data_ibegin=0, data_ni=nbp_mpi)
       CALL xios_set_domain_attr("domain_landpoints",data_ni=nbp_mpi, data_i_index=kindex_mpi-1)     

       ! Define longitudes and latitudes on local MPI domain depending on GridType
       IF ( GridType == "RegLonLat" ) THEN
          CALL xios_set_domain_attr("domain_landpoints",lonvalue_1d=lon_mpi(:,1),latvalue_1d=lat_mpi(1,:))
       ELSE IF ( GridType == "RegXY" ) THEN
          CALL xios_set_domain_attr("domain_landpoints",lonvalue_2d=lon_mpi,latvalue_2d=lat_mpi)
       ELSE
          WRITE(numout,*) 'Following GridType is not supported: GridType =', GridType
          CALL ipslerr_p(3, 'xios_orchidee_init', 'GridType not yet supported. ','Problem for defining longitudes and latitudes','')
       ENDIF

       !
       !! 4. Axis definition
       !
       CALL xios_set_axis_attr("nvm",n_glo=nvm ,VALUE=(/(REAL(i,r_std),i=1,nvm)/))
       CALL xios_set_axis_attr("nvmap",n_glo=nvmap ,VALUE=(/(REAL(i,r_std),i=1,nvmap)/))
       CALL xios_set_axis_attr("nlut",n_glo=nlut ,VALUE=(/(REAL(i,r_std),i=1,nlut)/))
       CALL xios_set_axis_attr("ncarb",n_glo=ncarb ,VALUE=(/(REAL(i,r_std),i=1,ncarb)/))
       CALL xios_set_axis_attr("nlaip1", n_glo=nlai+1,VALUE=(/(REAL(i,r_std),i=1,nlai+1)/))
       CALL xios_set_axis_attr("ngrnd",n_glo=ngrnd ,VALUE=soilth_lev(:))
       CALL xios_set_axis_attr("nstm", n_glo=nstm,VALUE=(/(REAL(i,r_std),i=1,nstm)/))
       CALL xios_set_axis_attr("nnobio", n_glo=nnobio,VALUE=(/(REAL(i,r_std),i=1,nnobio)/))
       CALL xios_set_axis_attr("albtyp", n_glo=2,VALUE=(/(REAL(i,r_std),i=1,2)/))
       CALL xios_set_axis_attr("nslm", n_glo=nslm,VALUE=(/(REAL(i,r_std),i=1,nslm)/))
       CALL xios_set_axis_attr("nwp", n_glo=nwp,VALUE=(/(REAL(i,r_std),i=1,nwp)/))
       CALL xios_set_axis_attr("ndeep", n_glo=ndeep,VALUE=(/(REAL(i,r_std),i=1,ndeep)/))
       CALL xios_set_axis_attr("P10", n_glo=10,VALUE=(/(REAL(i,r_std), i=1,10)/))
       CALL xios_set_axis_attr("P100", n_glo=100,VALUE=(/(REAL(i,r_std), i=1,100)/))
       CALL xios_set_axis_attr("P11", n_glo=11,VALUE=(/(REAL(i,r_std), i=1,11)/))
       CALL xios_set_axis_attr("P101", n_glo=101,VALUE=(/(REAL(i,r_std), i=1,101)/))
       IF (ok_explicitsnow) THEN
          CALL xios_set_axis_attr("nsnow", n_glo=nsnow,VALUE=(/(REAL(i,r_std),i=1,nsnow)/))
       ELSE
          CALL xios_set_axis_attr("nsnow", n_glo=1,VALUE=(/(REAL(i,r_std),i=1,1)/))
       END IF

       
       !
       !! 5. Get the default value (missing value) used by XIOS. This value is set in field_def_orchidee.xml
       !
       CALL xios_get_fieldgroup_attr("field_definition", default_value=xios_default_val)
       IF (printlev>=2) WRITE(numout,*) 'Default value read from XIOS, xios_default_val=',xios_default_val

       !
       !! 5. Deactivation of some fields if they are not calculated
       !
       IF ( OFF_LINE_MODE ) THEN
          CALL xios_set_field_attr("q2m",enabled=.FALSE.)
          CALL xios_set_field_attr("t2m",enabled=.FALSE.)
          CALL xios_set_field_attr("riverflow_cpl",enabled=.FALSE.)
          CALL xios_set_field_attr("coastalflow_cpl",enabled=.FALSE.)
       END IF

       IF ( .NOT. river_routing ) THEN
          CALL xios_set_field_attr("basinmap",enabled=.FALSE.)
          CALL xios_set_field_attr("nbrivers",enabled=.FALSE.)
          CALL xios_set_field_attr("riversret",enabled=.FALSE.)
          CALL xios_set_field_attr("hydrographs",enabled=.FALSE.)
          CALL xios_set_field_attr("fastr",enabled=.FALSE.)
          CALL xios_set_field_attr("slowr",enabled=.FALSE.)
          CALL xios_set_field_attr("streamr",enabled=.FALSE.)
          CALL xios_set_field_attr("laker",enabled=.FALSE.)
          CALL xios_set_field_attr("lake_overflow",enabled=.FALSE.)
          CALL xios_set_field_attr("mask_coast",enabled=.FALSE.)
          CALL xios_set_field_attr("pondr",enabled=.FALSE.)
          CALL xios_set_field_attr("floodr",enabled=.FALSE.)
          CALL xios_set_field_attr("slowflow",enabled=.FALSE.)
          CALL xios_set_field_attr("delfastr",enabled=.FALSE.)
          CALL xios_set_field_attr("delslowr",enabled=.FALSE.)
          CALL xios_set_field_attr("delstreamr",enabled=.FALSE.)
          CALL xios_set_field_attr("dellaker",enabled=.FALSE.)
          CALL xios_set_field_attr("delpondr",enabled=.FALSE.)
          CALL xios_set_field_attr("delfloodr",enabled=.FALSE.)
          CALL xios_set_field_attr("irrigmap",enabled=.FALSE.)
          CALL xios_set_field_attr("swampmap",enabled=.FALSE.)
          CALL xios_set_field_attr("wbr_stream",enabled=.FALSE.)
          CALL xios_set_field_attr("wbr_fast",enabled=.FALSE.)
          CALL xios_set_field_attr("wbr_slow",enabled=.FALSE.)
          CALL xios_set_field_attr("wbr_lake",enabled=.FALSE.)
          CALL xios_set_field_attr("reinfiltration",enabled=.FALSE.)
          CALL xios_set_field_attr("irrigation",enabled=.FALSE.)
          CALL xios_set_field_attr("netirrig",enabled=.FALSE.)
          CALL xios_set_field_attr("SurfStor",enabled=.FALSE.)
       END IF


       IF (hydrol_cwrr ) THEN
          CALL xios_set_field_attr("dss",enabled=.FALSE.)
          CALL xios_set_field_attr("gqsb",enabled=.FALSE.)
          CALL xios_set_field_attr("bqsb",enabled=.FALSE.)
          CALL xios_set_field_attr("rsol",enabled=.FALSE.)
       ELSE
          CALL xios_set_field_attr("frac_bare",enabled=.FALSE.)
          CALL xios_set_field_attr("twbr",enabled=.FALSE.)
          CALL xios_set_field_attr("nroot",enabled=.FALSE.)
          CALL xios_set_field_attr("dlh",enabled=.FALSE.)
          CALL xios_set_field_attr("mcs",enabled=.FALSE.)
          CALL xios_set_field_attr("water2infilt",enabled=.FALSE.)
          CALL xios_set_field_attr("reinf_slope",enabled=.FALSE.)
          CALL xios_set_field_attr("evapnu_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("drainage_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("transpir_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("runoff_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("tmc",enabled=.FALSE.)
          CALL xios_set_field_attr("njsc",enabled=.FALSE.)
          CALL xios_set_field_attr("k_litt",enabled=.FALSE.)
          CALL xios_set_field_attr("soilmoist",enabled=.FALSE.)
          CALL xios_set_field_attr("mc",enabled=.FALSE.)
          CALL xios_set_field_attr("kfact_root",enabled=.FALSE.)
          CALL xios_set_field_attr("vegetmax_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("undermcr",enabled=.FALSE.)
          CALL xios_set_field_attr("wtd",enabled=.FALSE.)
          CALL xios_set_field_attr("ru_corr",enabled=.FALSE.)
          CALL xios_set_field_attr("ru_corr2",enabled=.FALSE.)
          CALL xios_set_field_attr("dr_corr",enabled=.FALSE.)
          CALL xios_set_field_attr("dr_force",enabled=.FALSE.)
          CALL xios_set_field_attr("qinfilt",enabled=.FALSE.)
          CALL xios_set_field_attr("ru_infilt",enabled=.FALSE.)
          ! tws is defined in field_def.xml as a sum of several variables calculated only for cwrr
          CALL xios_set_field_attr("tws",enabled=.FALSE.)
       END IF

       IF (.NOT. ok_freeze_cwrr) THEN
          CALL xios_set_field_attr("profil_froz_hydro",enabled=.FALSE.)
          CALL xios_set_field_attr("temp_hydro",enabled=.FALSE.)
       END IF

       
       IF (.NOT. check_cwrr2) THEN
          CALL xios_set_field_attr("check_infilt",enabled=.FALSE.)
          CALL xios_set_field_attr("check_tr",enabled=.FALSE.)
          CALL xios_set_field_attr("check_over",enabled=.FALSE.)
          CALL xios_set_field_attr("check_under",enabled=.FALSE.)
       END IF

       IF ( .NOT. do_floodplains ) THEN
          CALL xios_set_field_attr("floodmap",enabled=.FALSE.)
          CALL xios_set_field_attr("floodh",enabled=.FALSE.)       
          CALL xios_set_field_attr("floodout",enabled=.FALSE.)       
       END IF

       ! Deactivate some stomate fields. 
       ! These fields were traditionally added in sechiba_history.nc output file.
       IF ( .NOT. ok_stomate ) THEN
          CALL xios_set_field_attr("nee",enabled=.FALSE.)
          CALL xios_set_field_attr("maint_resp",enabled=.FALSE.)
          CALL xios_set_field_attr("hetero_resp",enabled=.FALSE.)
          CALL xios_set_field_attr("growth_resp",enabled=.FALSE.)
          CALL xios_set_field_attr("npp",enabled=.FALSE.)
       END IF

       IF ( .NOT. do_irrigation ) THEN
          CALL xios_set_field_attr("irrigation",enabled=.FALSE.)
          CALL xios_set_field_attr("netirrig",enabled=.FALSE.)
          CALL xios_set_field_attr("irrigmap",enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_co2)THEN
          CALL xios_set_field_attr("vbetaco2",enabled=.FALSE.)
!          CALL xios_set_field_attr("cimean",enabled=.FALSE.)
!          CALL xios_set_field_attr("cim",enabled=.FALSE.)
          CALL xios_set_field_attr("gpp",enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_bvoc)THEN
          CALL xios_set_field_attr("PAR",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_fertil_no",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_iso",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_mono",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_ORVOC",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_MBO",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_methanol",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_acetone",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_acetal",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_formal",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_acetic",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_formic",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_no_soil",enabled=.FALSE.)
          CALL xios_set_field_attr("flx_no",enabled=.FALSE.)
          CALL xios_set_field_attr('flx_apinen'   ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_bpinen'   ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_limonen'  ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_myrcen'   ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_sabinen'  ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_camphen'  ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_3caren'   ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_tbocimen' ,enabled=.FALSE.)
          CALL xios_set_field_attr('flx_othermono',enabled=.FALSE.)
          CALL xios_set_field_attr('flx_sesquiter',enabled=.FALSE.)
          CALL xios_set_field_attr("CRF",enabled=.FALSE.)
          CALL xios_set_field_attr("fco2",enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_bvoc .OR. .NOT. ok_radcanopy ) THEN
          CALL xios_set_field_attr("PARdf",enabled=.FALSE.)
          CALL xios_set_field_attr("PARdr",enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_bvoc .OR. .NOT. ok_radcanopy .OR. .NOT. ok_multilayer ) THEN
          CALL xios_set_field_attr( 'PARsuntab',enabled=.FALSE.)
          CALL xios_set_field_attr( 'PARshtab' ,enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_bvoc .OR. .NOT. ok_radcanopy .OR. ok_multilayer ) THEN
          CALL xios_set_field_attr("PARsun",enabled=.FALSE.)
          CALL xios_set_field_attr("PARsh",enabled=.FALSE.)
          CALL xios_set_field_attr("laisun",enabled=.FALSE.)
          CALL xios_set_field_attr("laish",enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_bvoc .OR. .NOT. ok_bbgfertil_Nox) THEN
          CALL xios_set_field_attr("flx_co2_bbg_year",enabled=.FALSE.)
       END IF

       IF ( .NOT. ok_bvoc .OR. .NOT. ok_cropsfertil_Nox) THEN
          CALL xios_set_field_attr("N_qt_WRICE_year",enabled=.FALSE.)
          CALL xios_set_field_attr("N_qt_OTHER_year",enabled=.FALSE.)
       END IF

       ! Set record_offset for enable start in the middle of the year.
       ! julian_diff is the day of the year where the current run start
       IF (printlev>=3) WRITE(numout,*) 'In xios_orchidee_init, julian_diff, INT(julian_diff) =', &
            julian_diff, INT(julian_diff)

       IF (ok_nudge_mc .AND. nudge_interpol_with_xios) THEN
          ! Activate the input file with id="nudge_moistc" specified in file_def_orchidee.xml. 
          ! The nudging file should be called nudge_moistc.nc (see name in the xml file) and is 
          ! supposed to contain daily values for the full year for the variable moistc.
          CALL xios_set_file_attr("nudge_moistc",enabled=.TRUE.)
          ! Set record_offset to start read at correct day in the nudging file. 
          CALL xios_set_file_attr("nudge_moistc",record_offset=INT(julian_diff))
       ELSE
          ! Deactivate input file for nudging of soil moisture
          CALL xios_set_file_attr("nudge_moistc",enabled=.FALSE.)
          ! Deactivate variables related to soil moisture nudgnig
          CALL xios_set_field_attr("mask_moistc_interp",enabled=.FALSE.)
          CALL xios_set_field_attr("moistc_interp",enabled=.FALSE.)

          ! Deactivate output variables related to soil moisture nudging
          CALL xios_set_field_attr("mc_read_current",enabled=.FALSE.)
          CALL xios_set_field_attr("mc_read_prev",enabled=.FALSE.)
          CALL xios_set_field_attr("mc_read_next",enabled=.FALSE.)
          CALL xios_set_field_attr("mask_mc_interp_out",enabled=.FALSE.)
          CALL xios_set_field_attr("nudgincsm",enabled=.FALSE.)
       END IF

       IF (ok_nudge_snow .AND. nudge_interpol_with_xios) THEN
          ! Activate the input file with id="nudge_snow" specified in file_def_orchidee.xml. 
          ! The nudging file should be called nudge_snow.nc (see name in the xml file) and is 
          ! supposed to contain daily values for the full year for the variables snowdz, snowtemp and snowrho.
          CALL xios_set_file_attr("nudge_snow",enabled=.TRUE.)
          ! Set record_offset to start read at correct day in the nudging file. 
          CALL xios_set_file_attr("nudge_snow",record_offset=INT(julian_diff))
       ELSE
          ! Deactivate input file for nudging of snow variables
          CALL xios_set_file_attr("nudge_snow",enabled=.FALSE.)

          ! Deactivate input variables related to snow nudging
          CALL xios_set_field_attr("mask_snow_interp",enabled=.FALSE.)
          CALL xios_set_field_attr("snowdz_interp",enabled=.FALSE.)
          CALL xios_set_field_attr("snowrho_interp",enabled=.FALSE.)
          CALL xios_set_field_attr("snowtemp_interp",enabled=.FALSE.)

          ! Deactivate output variables related to snow nudging
          CALL xios_set_field_attr("snowdz_read_current",enabled=.FALSE.)
          CALL xios_set_field_attr("snowdz_read_prev",enabled=.FALSE.)
          CALL xios_set_field_attr("snowdz_read_next",enabled=.FALSE.)
          CALL xios_set_field_attr("snowrho_read_current",enabled=.FALSE.)
          CALL xios_set_field_attr("snowrho_read_prev",enabled=.FALSE.)
          CALL xios_set_field_attr("snowrho_read_next",enabled=.FALSE.)
          CALL xios_set_field_attr("snowtemp_read_current",enabled=.FALSE.)
          CALL xios_set_field_attr("snowtemp_read_prev",enabled=.FALSE.)
          CALL xios_set_field_attr("snowtemp_read_next",enabled=.FALSE.)
          CALL xios_set_field_attr("mask_snow_interp_out",enabled=.FALSE.)
          CALL xios_set_field_attr("nudgincswe",enabled=.FALSE.)
       END IF

       IF (impaze) THEN
          CALL xios_set_field_attr("soilalb_vis",enabled=.FALSE.)
          CALL xios_set_field_attr("soilalb_nir",enabled=.FALSE.)
          CALL xios_set_field_attr("vegalb_vis",enabled=.FALSE.)
          CALL xios_set_field_attr("vegalb_nir",enabled=.FALSE.)
       END IF

       IF (ok_explicitsnow) THEN
          ! The variable fusion is not calculated for ok_explicitsnow
          CALL xios_set_field_attr("Qf",enabled=.FALSE.)
       ELSE
          CALL xios_set_field_attr("pkappa_snow",enabled=.FALSE.)
          CALL xios_set_field_attr("pcapa_snow",enabled=.FALSE.)
          CALL xios_set_field_attr("snowliq",enabled=.FALSE.)
          CALL xios_set_field_attr("snowrho",enabled=.FALSE.)
          CALL xios_set_field_attr("snowheat",enabled=.FALSE.)
          CALL xios_set_field_attr("snowgrain",enabled=.FALSE.)
          CALL xios_set_field_attr("snowtemp",enabled=.FALSE.)
          CALL xios_set_field_attr("snowtemp_weighted",enabled=.FALSE.)
       END IF

       !
       !! 6. Close context
       !
     !yidi  CALL xios_close_context_definition()      
       CALL xios_close_context_definition()      


       !
       !! 7. Activate almaoutput if needed 
       !! Some extra calculations have to be done for the variables  
       !! delsoilmoist, delintercept, delswe and soilwet.
       !! Set almaoutput=true if at least one of these variables are defined in an output file. 
       !! If not, keep the initial value of almaoutput. 
       IF ( xios_field_is_active("delsoilmoist") .OR. xios_field_is_active("delintercept") .OR. &
            xios_field_is_active("delswe")       .OR. xios_field_is_active("soilwet")      .OR. &
            xios_field_is_active("twbr")) THEN

          almaoutput=.TRUE.
          IF (printlev >=3) WRITE(numout,*) 'The flag almaoutput has been activated in xios_orchidee_init'
       END IF
#endif
    END IF

    IF (xios_orchidee_ok) THEN
       ! Send variables to all OMP thredds
       CALL bcast(xios_default_val)
       CALL bcast(almaoutput)
    END IF

    IF (printlev>=3) WRITE(numout,*) 'End xios_orchidee_init'
  END SUBROUTINE xios_orchidee_init


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_change_context
  !!
  !>\BRIEF         Use this subroutine to switch between different context.
  !!               This subroutine must be called when running in coupled mode at each time ORCHIDEE is called, in the
  !!               begining and end of intersurf_gathered. First call is done after xios_orchidee_init is done. 
  !!
  !! DESCRIPTION  :\n 
  !!                  
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_change_context(new_context)
    !
    !! 0. Variable and parameter declaration
    !
    !!    Input variable
    CHARACTER(LEN=*),INTENT(IN)              :: new_context

    !! Local variables
#ifdef XIOS
    TYPE(xios_context) :: ctx_hdl
#endif
    !_ ================================================================================================================================

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
#ifdef XIOS
       CALL xios_get_handle(new_context,ctx_hdl)
       CALL xios_set_current_context(ctx_hdl)
#endif
    END IF
    
  END SUBROUTINE xios_orchidee_change_context

  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_update_calendar
  !!
  !>\BRIEF          Update the calandar in XIOS.
  !!
  !! DESCRIPTION  :\n Update the calendar in XIOS : let XIOS know that ORCHIDEE avanced one time-step.
  !!                  This subroutine should be called in the beginning of each time-step. The first 
  !!                  time-step in a new execution should always start at 1. Therefore, first calculate
  !!                  an offset that is substracted to the current time step in sechiba. 
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_update_calendar(itau_sechiba)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    INTEGER(i_std), INTENT(IN) :: itau_sechiba    !! Current time step of the model
    !
    !! 0.2 Local variables
    !
    LOGICAL, SAVE         :: first=.TRUE.         !! Flag for first entering in subroutine
    INTEGER(i_std), SAVE  :: offset               !! Offset to substract from itau_sechiba
    INTEGER(i_std)        :: itau_xios            !! Current time step for XIOS

    !_ ================================================================================================================================

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
#ifdef XIOS
       ! Calculate the offset
       IF (first) THEN
          offset=itau_sechiba-1
          first=.FALSE.
       END IF

       ! Substract the offset to the current time step in sechiba
       itau_xios=itau_sechiba-offset

       ! Send the new time step to XIOS
       IF (printlev>=3) WRITE(numout,*) 'xios_orchidee_update_calendar: itau_sechiba, itau_xios=',itau_sechiba,itau_xios
       CALL xios_update_calendar(itau_xios)
#endif
    END IF
  END SUBROUTINE xios_orchidee_update_calendar
  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_context_finalize
  !!
  !>\BRIEF         Finalize orchidee context.
  !!
  !! DESCRIPTION  :\n This subroutine finalizes the orchidee context without finalizing XIOS. In coupled mode, the atmospheric
  !!                  modele must finalize XIOS. This subroutine is called in the end of the execution of ORCHIDEE only in 
  !!                  coupeld mode.
  !!                  
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_context_finalize

    !_ ================================================================================================================================

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
       IF (printlev>=3) WRITE(numout,*) 'Entering xios_orchidee_context_finalize'
#ifdef XIOS
       CALL xios_context_finalize()
#endif
    END IF
  END SUBROUTINE xios_orchidee_context_finalize


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_finalize
  !!
  !>\BRIEF         Last call to XIOS for finalization.
  !!
  !! DESCRIPTION  :\n Last call to XIOS for finalization of the orchidee context and XIOS.
  !!                  This subroutine is called only when ORCHIDEE is run in offline mode. In coupled mode it is the atmospheric
  !!                  model that finalizes XIOS. In that case, the context orchidee must be finalized using the 
  !!                  subroutine xios_orchidee_context_finalize
  !!                  
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_finalize

    !_ ================================================================================================================================

    IF (xios_orchidee_ok .AND. is_omp_root) THEN
       IF (printlev>=3) WRITE(numout,*) 'Entering xios_orchidee_finalize'
#ifdef XIOS
       CALL xios_context_finalize()
       CALL xios_finalize()
#endif
    END IF
  END SUBROUTINE xios_orchidee_finalize


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r1d
  !!
  !>\BRIEF          Subroutine for sending 1D (array) fields to XIOS.
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 1D fields (array).
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r1d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)          :: field_id
    REAL(r_std), DIMENSION(:), INTENT(IN) :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_send_field_r1d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r1d


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r2d
  !!
  !>\BRIEF          Subroutine for sending 2D fields to XIOS.
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 2D fields.
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r2d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)            :: field_id
    REAL(r_std), DIMENSION(:,:), INTENT(IN) :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_send_field_r2d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r2d


  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r3d
  !!
  !>\BRIEF          Subroutine for sending 3D fields to XIOS. 
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 3D fields.
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r3d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    REAL(r_std), DIMENSION(:,:,:), INTENT(IN) :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2),size(field,3)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_send_field_r3d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r3d

  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r4d
  !!
  !>\BRIEF          Subroutine for sending 4D fields to XIOS. 
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 4D fields.
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r4d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    REAL(r_std), DIMENSION(:,:,:,:), INTENT(IN) :: field

    !! 0.2 Local variables
    INTEGER :: jv
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2),size(field,3),size(field,4)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_send_field_r4d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r4d

  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_send_field_r5d
  !!
  !>\BRIEF          Subroutine for sending 5D fields to XIOS. 
  !!
  !! DESCRIPTION  :\n Send one field to XIOS. This is the interface for 5D fields.
  !!                  NB! This subroutine should not be called directly. Use interface xios_orchidee_send_field.
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_send_field_r5d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    REAL(r_std), DIMENSION(:,:,:,:,:), INTENT(IN) :: field

    !! 0.2 Local variables
    INTEGER :: jv
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2),size(field,3),size(field,4),size(field,5)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_send_field_r5d, field_id=',field_id

       ! Gather all omp domains on the mpi domains
       CALL gather_omp(field, field_mpi)

       ! All master threads send the field to XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_send_field(field_id,field_mpi)
#endif
       END IF
    END IF
  END SUBROUTINE xios_orchidee_send_field_r5d
 
  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_recv_field_r2d
  !!
  !>\BRIEF          Subroutine for receiving 1D (kjpindex) fields to XIOS. 
  !!
  !! DESCRIPTION  :\n 
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_recv_field_r1d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    
    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:), INTENT(OUT)    :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi)           :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_recv_field_r1d, field_id=',field_id

       ! All master threads receive the field from XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_recv_field(field_id,field_mpi)
          IF (printlev>=5) WRITE(numout,*) 'Recieve done with xios_orchidee_recv_field_r1d, field_id=',field_id
#endif
       END IF

       ! Scatter the mpi domains on local omp domains
       CALL scatter_omp(field_mpi, field)

    END IF
  END SUBROUTINE xios_orchidee_recv_field_r1d

  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_recv_field_r2d
  !!
  !>\BRIEF          Subroutine for receiving 2D(kjpindex and 1 vertical axe) fields to XIOS. 
  !!
  !! DESCRIPTION  :\n 
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_recv_field_r2d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    
    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:,:), INTENT(OUT)  :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_recv_field_r2d, field_id=',field_id

       ! All master threads recieve the field from XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_recv_field(field_id,field_mpi)
          IF (printlev>=5) WRITE(numout,*) 'Recieve done with xios_orchidee_recv_field_r2d, field_id=',field_id
#endif
       END IF

       ! Scatter the mpi domains on local omp domains
       CALL scatter_omp(field_mpi, field)

    END IF
  END SUBROUTINE xios_orchidee_recv_field_r2d

  !! ==============================================================================================================================
  !! SUBROUTINE   : xios_orchidee_recv_field_r3d
  !!
  !>\BRIEF          Subroutine for receiving 3D(kjpindex and 2 vertical axes) fields to XIOS. 
  !!
  !! DESCRIPTION  :\n 
  !!
  !! \n
  !_ ================================================================================================================================
  SUBROUTINE xios_orchidee_recv_field_r3d(field_id,field)
    !
    !! 0. Variable and parameter declaration
    !
    !! 0.1 Input variables
    !
    CHARACTER(len=*), INTENT(IN)              :: field_id
    
    !! 0.2 Output variables
    REAL(r_std), DIMENSION(:,:,:), INTENT(OUT) :: field

    !! 0.2 Local variables
    REAL(r_std), DIMENSION(nbp_mpi,size(field,2),size(field,3)) :: field_mpi

    !_ ================================================================================================================================
    IF (xios_orchidee_ok) THEN
       IF (printlev>=4) WRITE(numout,*) 'Entering xios_orchidee_recv_field_r3d, field_id=',field_id

       ! All master threads receive the field from XIOS
       IF (is_omp_root) THEN
#ifdef XIOS
          CALL xios_recv_field(field_id,field_mpi)
          IF (printlev>=5) WRITE(numout,*) 'Recieve done with xios_orchidee_recv_field_r3d, field_id=',field_id
#endif
       END IF

       ! Scatter the mpi domains on local omp domains
       CALL scatter_omp(field_mpi, field)

    END IF
  END SUBROUTINE xios_orchidee_recv_field_r3d

END MODULE xios_orchidee

