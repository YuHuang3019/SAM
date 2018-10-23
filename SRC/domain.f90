! Set the domain dimensionality, size and number of subdomains.

module domain

       integer, parameter :: YES3D = 1  ! Domain dimensionality: 1 - 3D, 0 - 2D
       integer, parameter :: nx_gl = 384 ! Number of grid points in X
       integer, parameter :: ny_gl = 512 ! Number of grid points in Y
       integer, parameter :: nz_gl = 64 ! Number of pressure (scalar) levels
       integer, parameter :: nsubdomains_x  = 4 ! No of subdomains in x
       integer, parameter :: nsubdomains_y  = 6 ! No of subdomains in y


       ! define # of points in x and y direction to average for 
       !   output relating to statistical moments.
       ! For example, navgmom_x = 8 means the output will be   
       !  8 times coarser grid than the original.
       ! If don't wanna such output, just set them to -1 in both directions. 
       ! See Changes_log/README.UUmods for more details.
       integer, parameter :: navgmom_x = 96
       integer, parameter :: navgmom_y = 96

       integer, parameter :: ntracers = 0 ! number of transported tracers (dotracers=.true.)
       
! Note:
!  * nx_gl and ny_gl should be a factor of 2,3, or 5 (see User's Guide)
!  * if 2D case, ny_gl = nsubdomains_y = 1 ;
!  * nsubdomains_x*nsubdomains_y = total number of processors
!  * if one processor is used, than  nsubdomains_x = nsubdomains_y = 1;
!  * if ntracers is > 0, don't forget to set dotracers to .true. in namelist 

end module domain
