program translate
	implicit none
	character*269 :: rec, outdir, infile, outfile
	character*11 :: ds
	character*4 :: element,year
	character*5 :: v(31)
	character*2 :: month
	integer :: nmonth, i, j, narg

	narg = command_argument_count()
	if (narg.ne.2) then
	  write(*,*) 'Syntax:  "translate <ds> <outdir>" where <ds> indicates the GHCN-DAILY data file and <dir> is the output directory'
	  call exit(0)
	end if
	call get_command_argument(1,ds)
	call get_command_argument(2,outdir)

	infile = trim(ds)//'.dly'
	outfile = trim(outdir)//trim(ds)//'.gms'
	write(*,'(a)') 'Input file: '//trim(infile)
	write(*,'(a)') 'Output file: '//trim(outfile)

	open(10,file=trim(ds)//'.dly')
	open(11,file=trim(outfile))
	nmonth = 0
	do
	  read(10,'(a)',end=100) rec
	  element = rec(18:21)
	  if (element.eq.'TMAX') then
	    nmonth = nmonth + 1
	    year = rec(12:15)
	    month = rec(16:17)
	    do i=1,31
	      j = 1 + 8*(i-1) + 21
	      v(i) = rec(j:j+4)
	      if (v(i).eq.'-9999') then
	        write(11,'(a,i2.2)') '* No data for '//year//'.'//month//'.',i
	      else
	        write(11,'(a,i2.2,2x,a)') year//'.'//month//'.',i,v(i)
	      endif
	    end do
	  end if
	end do
100	close(10)
	write(*,'(a,i0)') 'nmonth = ',nmonth
end program translate
