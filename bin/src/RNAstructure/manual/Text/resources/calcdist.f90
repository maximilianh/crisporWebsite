! make shared library calcdist.so by typing:
! gfortran -O3 -shared -fPIC calcdist.f90 -o calcdist.so
subroutine calcdist(no_nucl,no_struct,pairing_data,distance)

implicit none

integer, intent(in) :: no_struct, no_nucl, pairing_data(no_struct,no_nucl,no_nucl)
integer, intent(out) :: distance(no_struct,no_struct)
integer :: i,j

distance = 0
do i = 1, no_struct
   do j = 1, no_struct
      distance(i,j) = sum(abs(pairing_data(i,:,:)-pairing_data(j,:,:)))
   end do
end do


end subroutine calcdist
