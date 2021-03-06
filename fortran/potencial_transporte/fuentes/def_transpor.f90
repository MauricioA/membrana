module def_transpor

         
   double precision :: Faraday= 96485.34  !C/mol
   double precision :: R_cte= 8.314      ! J/K/mol
   double precision :: T_cte= 350.0      ! K                 
   double precision :: D_h = 12500.0     ! micro**2/sec
   double precision :: D_oh= 7050.0      ! micro**2/sec
   double precision :: D_Na = 1780.0     ! micro**2/sec
   double precision :: D_Cl= 2720.0      ! micro**2/sec
   double precision :: D_eff= 0.5      ! micro**2/sec                 
    
   double precision :: Kbolz= 8.61733e-5 ! ev/K/mol

   double precision :: kwf= 249.16  ! 249.16  ! micro**3/at/sec      1.5e+11  ! 1/M/sec
   double precision :: kwb= 2.7E-5  ! 1/sec



   double precision :: E_eq=1.23    ! Volt
   double precision :: I_eq=1.0E-18   ! A/micro**2

   double precision :: E2_eq=1.36    ! Volt
   double precision :: I2_eq=1.0E-12   ! A/micro**2


   
   double precision ::  epsilon = 78.5         !  Constante dielectrica del agua 
   double precision ::  epsilon_0 =8.85E-12    !  Constante de permitividad C**2 / (N m**2) 
   double precision ::  zh = 1.0               !	          Carga de los protones (e) */
   double precision ::  zoh = -1.0              !	          Carga de los OH (-e) */
   double precision ::  zna = 1.0               !	          Carga de los NA (e) */
   double precision ::  zcl = -1.0              !	          Carga de los CL (-e) */
   double precision ::  tabs = 293.0           !           Temperatura (Kelvin) */
!   double precision ::  e = 4.8E-10            !	 Carga del proton (uc)*/
   double precision ::  Ch2o=  3.34E+10        !  at/microm**3   55.5 M
   double precision ::  ch_inicial=  60.20        !  at/microm**3   1e-7M
   double precision ::  Coh_inicial= 60.20        !  at/microm**3   1e-7M

   double precision ::  cna_inicial=  96320000.0        !  at/microm**3   0.16M
   double precision ::  Ccl_inicial= 96320000.0        !  at/microm**3   0.16M


   double precision ::  Ch_anodo= 60.20  !1.5E+7        !  at/microm**3   
   double precision ::  Coh_anodo=60.20         !  at/microm**3   
   double precision ::  Cna_anodo=96320000.0 !1.0E+12        !  at/microm**3   
   double precision ::  Ccl_anodo=96320000.0 ! 0.0        !  at/microm**3   

   double precision ::  Ch_catodo= 60.20    !  at/microm**3   
   double precision ::  Coh_catodo= 60.20  ! 1.806E+7  !  at/microm**3   
   double precision ::  Cna_catodo=96320000.0 ! 0.0        !  at/microm**3   
   double precision ::  Ccl_catodo=96320000.0 !  0.0  !  at/microm**3   


   double precision,allocatable :: ch(:),grad_ch_x(:),grad_ch_y(:)
   double precision,allocatable :: coh(:),grad_coh_x(:),grad_coh_y(:)
   double precision,allocatable :: cna(:),grad_cna_x(:),grad_cna_y(:)
   double precision,allocatable :: ccl(:),grad_ccl_x(:),grad_ccl_y(:)

   
   double precision,allocatable ::rhs_h(:),rhs_oh(:),rhs_na(:),rhs_cl(:)  ! ***cambio***  ! agrego estos vectoes

   double precision :: Number_clave    


   double precision :: ccurrent,rresist,currentCl=0,currentH=0
   
end module def_transpor