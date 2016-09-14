
	character*15 fin,fout,dset*5
	integer nlon,nlat
	real,allocatable:: lon(:), lat(:),var(:,:)

	character*80 head
	
	dset="cfsv2"
	fin="zs."//dset(:len_trim(dset))//".cmp"
	fout="zs."//dset(:len_trim(dset))//".dat"

	open(10,file=fin, action="read"
     &,form="unformatted")


	read(10) nlat;allocate(lat(nlat))
	read(10) lat
	print*, lat
	
	read(10) nlon;allocate(lon(nlon),var(nlon,nlat))
	read(10) lon
	print*, lon
		
	read(10) head
	print*, head
	read(10) ((var(i,j),i=1,nlon),j=1,nlat)
	print*, var(1,:)

	open(20,file=fout,action='write')
	write(20,'(3f10.2)') ((lon(i),lat(j),var(i,j),i=1,nlon),j=1,nlat)
	close(20)
	
	read(10) nlon
	print*, nlon

	end
