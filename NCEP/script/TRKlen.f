       program reformtrk
	implicit none
!rewrites trk file
	character*80 fin,fout,optarg
	character ch5*5,fline*97,nline*116,tline(200)*116

	integer i,nt,nnt,nit,narg,MinN,y1,y2

       narg=iargc()
	i=1
	do while(i.le.narg)
	call getarg(i,optarg)
	if(optarg.eq.'-i')then
	 i=i+1
	 call getarg(i,fin)
	elseif(optarg=='-o')then
	 i=i+1
	 call getarg(i,fout)
	elseif(optarg=='-n')then
	 i=i+1
	 call getarg(i,optarg)
	 read(optarg,*) MinN
	endif
	i=i+1
	enddo


	open(10,file=fin,action='read')
	open(20,file=fout,action='write')

	do i=1,66
	read(10,*)
	enddo

	nnt=0
	do 
	read(10,'(x,a5)',end=999)ch5
	if(ch5=='Track')then
	 backspace(10)
	 read(10,'(a)') fline
	 read(fline(56:59),'(i4)')nit
	 if(nit>200)then
	  write(*,'("ERROR(reformtrk: Check track length:",a)')
     &          fline(1:11)
	  stop
	 endif
 	 read(10,*)
	 read(10,'(a)')nline
	 do i=1,nit
	  read(10,'(a)')tline(i)
	 enddo
	 read(10,*)

	 if(nit>=MinN)then
	  nnt=nnt+1
	  write(20,*)
	  read(fline(61:62),'(i2)') y1;call yr(y1);
	  read(fline(74:75),'(i2)') y2;call yr(y2)
!for 2-digit year
!	  write(20,'(a,i5,"(",a,")",a)')fline(1:6),nnt,fline(7:11)
!     &                                 ,fline(12:)
!for 4-digit year
	  write(20,'(a,i5,"(",a,")",a,i4,a,i4,a)')fline(1:6),nnt
     &          ,fline(7:11),fline(12:60),y1,fline(63:73),y2,fline(76:) 
	  write(20,*)
!	  write(20,'(a)')nline  !for 2-digit yr
	  write(20,'(a,"  ",a)')nline(1:10),nline(11:) !for 4-digit yr
	  do i=1,nit
	  read(tline(i)(13:14),'(i2)') y1;call yr(y1)
!	   write(20,'(a,i4,a)')tline(i)  !for 2-digit yr
	   write(20,'(a,i4,a)')tline(i)(1:12),y1,tline(i)(15:)  !for 4-digit yr
	  enddo
	 endif
	
	else
	 write(*,'("ERROR(reformtrk): Check input file!")')

	endif
	enddo

999	close(10)
	close(20)

	stop
        end


	subroutine yr(y)

	integer y,ny


	if (y>25)then
	 ny=y+1900
	elseif (y<18)then
	 ny=y+2000
	elseif(y==18)then
	 ny=1998
	elseif(y==19)then
	 ny=1999
	elseif(y==20)then
	 ny=2000
	elseif(y==21)then
	 ny=2001
	elseif(y==22)then
	 ny=2002
	endif

	y=ny
	return
        end
	
