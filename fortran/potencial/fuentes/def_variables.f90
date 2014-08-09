module def_variables

character(10):: problema
character(20):: archi_malla='    '
character(20):: archi_sistema='    '

double precision,allocatable :: coor_x(:),coor_y(:),coor_z(:),grad_x(:),grad_y(:),grad_z(:),cer(:),gradxel_x(:),gradxel_y(:),gradxel_z(:), corr_x(:), corr_y(:), corrxel_x(:), corrxel_y(:)

double precision,allocatable :: tension(:,:),deforma(:,:)
integer :: ntension=4

integer, allocatable :: conect(:,:), material(:),nnodtierra(:),nnodpoten(:),matelem(:),vec_tierra(:),vec_poten(:),vec_quieto(:)
integer, allocatable :: nod_dirichV(:),nod_dirichT(:),nod_quieto(:)
integer :: ndirichV,ndirichT,nquieto

double precision,allocatable ::  Zgrid(:)
integer, allocatable :: grid(:,:),grid_z(:,:,:)
integer :: nzgrid,Ngrid_z,Ngrid_y,Ngrid_x

integer :: nelements,nnodes,nodpel,nod_tierra,nod_poten,mat_externo,mat_membrana,mat_interno
integer :: num_mat

integer :: nopcion

double precision :: sigmaint,sigmaext,sigmamem,permit
double precision :: potencial,frecuencia


integer :: unit_data=1,  &! data de entrada
   &       unit_malla=2,   & ! malla de salida
   &       unit_cc=3,   &  ! condiciones de contorno    
   &       unit_cont=10,   & ! para controlas
   &       unit_sal=11,   &  ! resultados
   &       unit_2d=12,   &  
   &       unit_gra=13,   & 
   &       unit_grid=14,   & 
   &       unit_camp=15,   &    ! camp elemental
   &       unit_sal2d=16,   &   ! desplazameinos
   &       unit_el2d=17,    &   ! tensiones elementsles
   &       unit_corr=18			! corrientes
   
   
   
character*(30) :: filedata='input.in',  & ! entrada
     &          filemalla='mallado.fem' , &  ! malla
     &          filecc = 'contorno.fem',  &    ! CC
     &          file_aux='control.dat',  & ! controla
     &          file_sal='results.3D',   &   ! salidas de resultados
     &          file_2D='saleplano.dat',  & ! alguna salida especial     
     &          file_gra='gradiente.dat', & ! gradientes     
     &          file_grid='grid2D.dat',  &   ! alguna salida especial     
     &          file_camp='campo.dat',    & 
     &          file_sal2d='desplaz.csv' ,   & 
     &          file_el2d='elemen.csv', &
     &          file_corr='corriente.csv'	! corriente
     
end module def_variables
