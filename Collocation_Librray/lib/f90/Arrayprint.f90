Module Array_Print
   include 'defs.fi'
   ! ---------------------------------------------------
   ! I finally got tired of doing write statements, so did this module to
   ! handle most situations. This is more difficult than it ought to be
   ! because of stupid purist FORTRAN restrictions. No wonder everyone moved
   ! to other languages.

   ! A FORTRAN commandment:
   !     "THOU SHALT NOT CHANGE TKR (Type, Kind or Rank)"
   ! Why can't A(5) be the same as A(5,1) or A(10,10) be the same as A(100)?
   ! Nope not allowed in "Clean" code.
   ! It can't be done if Interfaces are used everywhere, i.e. with "clean" code
   ! Maybe one of the newer 201x versions allow it, if so, I'd like to hear of it.
   ! All because of this silly rule there is no way to change rank without lots
   ! of useless copying. Imagine arrays with megabytes of data.
   ! It would be so easy if casting were allowed (like C++)
   !        End of Flame!
   ! ---------------------------------------------------
   Private
   Public OpenFile
   Public VectorPrint
   Interface VectorPrint
      Module Procedure VectorPrintR     ! horizontal print reals
      Module Procedure VectorPrintI     ! horizontal print integers
      Module Procedure VectorPrintQ     ! horizontal print quad
   End Interface
   Public ArrayPrint
   Interface ArrayPrint    ! tabular print, overload for various ranks
      Module Procedure ArrayPrint_aaaa    ! (v - vector, a - array)
      Module Procedure ArrayPrint_vaaa
      Module Procedure ArrayPrint_avaa
      Module Procedure ArrayPrint_vvaa
      Module Procedure ArrayPrint_aava
      Module Procedure ArrayPrint_vava
      Module Procedure ArrayPrint_avva
      Module Procedure ArrayPrint_vvva
      Module Procedure ArrayPrint_vvvv
      Module Procedure ArrayPrint_avvv
      Module Procedure ArrayPrint_vavv
      Module Procedure ArrayPrint_avav
      Module Procedure ArrayPrint_vaav
      Module Procedure ArrayPrint_vvav
      Module Procedure ArrayPrint_aavv
      Module Procedure ArrayPrint_aaav
      Module Procedure ArrayPrintI         ! integer arrays (all)
      Module Procedure ArrayPrintQ_aaaa    ! (v - vector, a - array)
      Module Procedure ArrayPrintQ_vaaa
      Module Procedure ArrayPrintQ_vava
      Module Procedure ArrayPrintQ_vvaa
      Module Procedure ArrayPrintQ_vvva
      Module Procedure ArrayPrintQ_vvvv
      Module Procedure ArrayPrintQ_avaa
      Module Procedure ArrayPrintQ_avva
      Module Procedure ArrayPrintQ_aava
   End Interface
   Character, Parameter :: fmti*(*) = '(20(I5,a1))'
   Character, Parameter :: fmtd*(*) = '(1x,20(g14.7,a1))'
  Contains
   ! -------------------------------------------------------
   Subroutine OpenFile(fname,ext,append)
      Character :: fname*(*)
      Character, optional :: ext*(3)
      Character, optional :: append*(*)
      Character :: fullname*(128),etyp*(3)
      integer :: iounit,iok
      iounit = iout2;  etyp = 'dat'
      if(present(ext))etyp = ext
      if(etyp == 'dat')iounit = iout
      if(etyp == 'log')iounit = ilog
      if(present(append))then
         iounit = iout2
         fullname = trim(fname)//trim(append)//'.'//etyp
      else
         fullname = trim(fname)//'.'//etyp
      endif
!      write(*,'(3a,i2)')' Open: "',trim(fullname),'", no.',iounit
      Open(iounit,file=trim(fullname),iostat=iok)
   End Subroutine OpenFile
   ! -------------------------------------------------------
   Subroutine VectorPrintR(title,v1,v2,v3,v4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: v1(:)
      Real(float), optional :: v2(:),v3(:),v4(:)
      Logical, optional :: flog
      Character :: fmtt*(80)
      integer :: i,iunit
      ! ---------------------------------------------------
      ! Print up to 4 vectors (strung together)
      !  title is an optional title
      !  v1,v2,v3,v4 - vectors, only 1 required
      !  fmtx is optional format
      !  flog optional file designator as follows
      !     = .true. use ilog file
      !     = .false. use iout2 file
      !     = .not.present use iout file
      ! ---------------------------------------------------
      if(Present(fmtx))then
         fmtt = fmtx
      else
         fmtt = fmtd
      endif
      iunit = iout
      if(present(flog))then
         iunit = iout2
         if(flog)iunit = ilog
      endif
      if(present(title))then
         if(title .ne. ' ')write(iunit,'(a)')title
      endif
      if(Present(v4).and.present(v3).and.present(v2))then
         write(iunit,fmtt)(v1(i),tb,i=1,size(v1)),(v2(i),tb,i=1,size(v2)), &
            (v3(i),tb,i=1,size(v3)),(v4(i),tb,i=1,size(v4))
      elseif(present(v3).and.present(v2))then
         write(iunit,fmtt)(v1(i),tb,i=1,size(v1)),(v2(i),tb,i=1,size(v2)),(v3(i),tb,i=1,size(v3))
      elseif(present(v2))then
         write(iunit,fmtt)(v1(i),tb,i=1,size(v1)),(v2(i),tb,i=1,size(v2))
      else
         write(iunit,fmtt)(v1(i),tb,i=1,size(v1))
      endif
   End Subroutine VectorPrintR
   ! -------------------------------------------------------
   Subroutine VectorPrintQ(title,v1,v2,v3,v4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(floatq), intent(in) :: v1(:)
      Real(floatq), optional :: v2(:),v3(:),v4(:)
      Logical, optional :: flog
      Character :: fmtt*(80)
      integer :: i,iunit
      ! ---------------------------------------------------
      ! Print up to 4 vectors (strung together)
      !  title is an optional title
      !  v1,v2,v3,v4 - vectors, only 1 required
      !  fmtx is optional format
      !  flog optional file designator as follows
      !     = .true. use ilog file
      !     = .false. use iout2 file
      !     = .not.present use iout file
      ! ---------------------------------------------------
      if(Present(fmtx))then
         fmtt = fmtx
      else
         fmtt = fmtd
      endif
      iunit = iout
      if(present(flog))then
         iunit = iout2
         if(flog)iunit = ilog
      endif
      if(present(title))then
         if(title .ne. ' ')write(iunit,'(a)')title
      endif
      if(Present(v4).and.present(v3).and.present(v2))then
         write(iunit,fmtt)(v1(i),tb,i=1,size(v1)),(v2(i),tb,i=1,size(v2)), &
            (v3(i),tb,i=1,size(v3)),(v4(i),tb,i=1,size(v4))
      elseif(present(v3).and.present(v2))then
         write(iunit,fmtt)(v1(i),tb,i=1,size(v1)),(v2(i),tb,i=1,size(v2)),(v3(i),tb,i=1,size(v3))
      elseif(present(v2))then
         write(iunit,fmtt)(v1(i),tb,i=1,size(v1)),(v2(i),tb,i=1,size(v2))
      else
         write(iunit,fmtt)(v1(i),tb,i=1,size(v1))
      endif
   End Subroutine VectorPrintQ
   ! -------------------------------------------------------
   Subroutine VectorPrintI(title,v1,v2,v3,v4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      integer, intent(in) :: v1(:)
      integer, optional :: v2(:),v3(:),v4(:)
      Logical, optional :: flog
      Character :: fmtt*(80)
      integer :: i,iunit
      ! ---------------------------------------------------
      ! Print up to 4 vectors (strung together)
      !  title is an optional title
      !  v1,v2,v3,v4 - vectors, only 1 required
      !  fmtx is optional format
      !  flog optional file designator as follows
      !     = .true. use ilog file
      !     = .false. use iout2 file
      !     = .not.present use iout file
      ! ---------------------------------------------------
      if(Present(fmtx))then
         fmtt = fmtx
      else
         fmtt = fmti
      endif
      iunit = iout
      if(present(flog))then
         iunit = iout2
         if(flog)iunit = ilog
      endif
      if(present(title))then
         if(title .ne. ' ')write(iunit,'(a)')title
      endif
      if(Present(v4).and.present(v3).and.present(v2))then
         write(iunit,fmtt)(v1(i),tb,i=1,size(v1)),(v2(i),tb,i=1,size(v2)), &
            (v3(i),tb,i=1,size(v3)),(v4(i),tb,i=1,size(v4))
      elseif(present(v3).and.present(v2))then
         write(iunit,fmtt)(v1(i),tb,i=1,size(v1)),(v2(i),tb,i=1,size(v2)),(v3(i),tb,i=1,size(v3))
      elseif(present(v2))then
         write(iunit,fmtt)(v1(i),tb,i=1,size(v1)),(v2(i),tb,i=1,size(v2))
      else
         write(iunit,fmtt)(v1(i),tb,i=1,size(v1))
      endif
   End Subroutine VectorPrintI
   ! -------------------------------------------------------
   Subroutine ArrayPrint_aaaa(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:,:)
      Real(float), optional :: A2(:,:),A3(:,:),A4(:,:)
      Logical, optional :: flog
      integer :: iunit
      Character :: fmtt*(80)
      integer :: i,j,n,n1(2),n2(2),n3(2),n4(2)
      ! ---------------------------------------------------
      ! Print up to 4 arrays (strung together side by side)
      !  title - is an optional title
      !     A1,A2,A3,A4 - arrays, only 1 required
      !  fmtx is optional format
      !  flog optional file designator as follows
      !     = .true. use ilog file
      !     = .false. use iout2 file
      !     = .not.present use iout file
      ! ---------------------------------------------------
      if(Present(fmtx))then
         fmtt = fmtx
      else
         fmtt = fmtd
      endif
      iunit = iout
      if(present(flog))then
         iunit = iout2
         if(flog)iunit = ilog
      endif
      if(present(title))then
         if(title .ne. ' ')write(iunit,'(a)')title
      endif
      n1 = shape(A1)
      if(Present(A2).and.Present(A3).and.Present(A4))then
         n2 = shape(A2); n3 = shape(A3); n4 = shape(A4); n = min(n1(1),n2(1),n3(1),n4(1))
         do i=1,n
            write(iunit,fmtt)(A1(i,j),tb,j=1,n1(2)),(A2(i,j),tb,j=1,n2(2)), &
               (A3(i,j),tb,j=1,n3(2)),(A4(i,j),tb,j=1,n4(2))
         enddo
      elseif(Present(A2).and.Present(A3))then
         n2 = shape(A2);  n3 = shape(A3);  n = min(n1(1),n2(1),n3(1))
         do i=1,n
            write(iunit,fmtt)(A1(i,j),tb,j=1,n1(2)),(A2(i,j),tb,j=1,n2(2)), &
               (A3(i,j),tb,j=1,n3(2))
         enddo
      elseif(Present(A2))then
         n2 = shape(A2); n = min(n2(1),n1(1))
         do i=1,n
            write(iunit,fmtt)(A1(i,j),tb,j=1,n1(2)),(A2(i,j),tb,j=1,n2(2))
         enddo
      else
         do i=1,n1(1)
            write(iunit,fmtt)(A1(i,j),tb,j=1,n1(2))
         enddo
      endif
   End Subroutine ArrayPrint_aaaa
   ! -------------------------------------------------------
   Subroutine ArrayPrint_vaaa(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:)
      Real(float), optional :: A2(:,:),A3(:,:),A4(:,:)
      Logical, optional :: flog
      Real(float) :: dum(size(A1),1)
      dum(:,1) = A1
      call ArrayPrint_aaaa(title,dum,A2,A3,A4,fmtx,flog)
   End Subroutine ArrayPrint_vaaa
   ! -------------------------------------------------------
   Subroutine ArrayPrint_avaa(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:,:),A2(:)
      Real(float), optional :: A3(:,:),A4(:,:)
      Logical, optional :: flog
      Real(float) :: dum2(size(A2),1)
      dum2(:,1) = A2
      call ArrayPrint_aaaa(title,A1,dum2,A3,A4,fmtx,flog)
   End Subroutine ArrayPrint_avaa
   ! -------------------------------------------------------
   Subroutine ArrayPrint_vvaa(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:),A2(:)
      Real(float), optional :: A3(:,:),A4(:,:)
      Logical, optional :: flog
      Real(float) :: dum1(size(A1),1),dum2(size(A2),1)
      dum1(:,1) = A1;  dum2(:,1) = A2
      call ArrayPrint_aaaa(title,dum1,dum2,A3,A4,fmtx,flog)
   End Subroutine ArrayPrint_vvaa
   ! -------------------------------------------------------
   Subroutine ArrayPrint_vava(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:),A2(:,:),A3(:)
      Real(float), optional :: A4(:,:)
      Logical, optional :: flog
      Real(float) :: dum(size(A1),1),dum3(size(A3),1)
      dum(:,1) = A1; dum3(:,1) = A3
      call ArrayPrint_aaaa(title,dum,A2,dum3,A4,fmtx,flog)
   End Subroutine ArrayPrint_vava
   ! -------------------------------------------------------
   Subroutine ArrayPrint_vvva(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:),A2(:),A3(:)
      Real(float), optional :: A4(:,:)
      Logical, optional :: flog
      Real(float) :: dum1(size(A1),1),dum2(size(A2),1),dum3(size(A3),1)
      dum1(:,1) = A1;  dum2(:,1) = A2;  dum3(:,1) = A3
      call ArrayPrint_aaaa(title,dum1,dum2,dum3,A4,fmtx,flog)
   End Subroutine ArrayPrint_vvva
   ! -------------------------------------------------------
   Subroutine ArrayPrint_avva(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:,:),A2(:),A3(:)
      Real(float), optional :: A4(:,:)
      Logical, optional :: flog
      Real(float) :: dum2(size(A2),1),dum3(size(A3),1)
      dum2(:,1) = A2;  dum3(:,1) = A3
      call ArrayPrint_aaaa(title,A1,dum2,dum3,A4,fmtx,flog)
   End Subroutine ArrayPrint_avva
   ! -------------------------------------------------------
   Subroutine ArrayPrint_aava(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:,:),A2(:,:),A3(:)
      Real(float), optional :: A4(:,:)
      Logical, optional :: flog
      Real(float) :: dum3(size(A3),1)
      dum3(:,1) = A3
      call ArrayPrint_aaaa(title,A1,A2,dum3,A4,fmtx,flog)
   End Subroutine ArrayPrint_aava
   ! -------------------------------------------------------
   Subroutine ArrayPrint_vvvv(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:),A2(:),A3(:),A4(:)
      Logical, optional :: flog
      Real(float) :: dum1(size(A1),1),dum2(size(A2),1),dum3(size(A3),1),dum4(size(A4),1)
      dum1(:,1) = A1;  dum2(:,1) = A2;  dum3(:,1) = A3;  dum4(:,1) = A4
      call ArrayPrint_aaaa(title,dum1,dum2,dum3,dum4,fmtx,flog)
   End Subroutine ArrayPrint_vvvv
   ! -------------------------------------------------------
   Subroutine ArrayPrint_aaav(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:,:),A2(:,:),A3(:,:),A4(:)
      Logical, optional :: flog
      Real(float) :: dum4(size(A4),1)
      dum4(:,1) = A4
      call ArrayPrint_aaaa(title,A1,A2,A3,dum4,fmtx,flog)
   End Subroutine ArrayPrint_aaav
   ! -------------------------------------------------------
   Subroutine ArrayPrint_aavv(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:,:),A2(:,:),A3(:),A4(:)
      Logical, optional :: flog
      Real(float) :: dum3(size(A3),1),dum4(size(A4),1)
      dum3(:,1) = A3;  dum4(:,1) = A4
      call ArrayPrint_aaaa(title,A1,A2,dum3,dum4,fmtx,flog)
   End Subroutine ArrayPrint_aavv
   ! -------------------------------------------------------
   Subroutine ArrayPrint_avav(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:,:),A2(:),A3(:,:),A4(:)
      Logical, optional :: flog
      Real(float) :: dum2(size(A2),1),dum4(size(A4),1)
      dum2(:,1) = A2;  dum4(:,1) = A4
      call ArrayPrint_aaaa(title,A1,dum2,A3,dum4,fmtx,flog)
   End Subroutine ArrayPrint_avav
   ! -------------------------------------------------------
   Subroutine ArrayPrint_vvav(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:),A2(:),A3(:,:),A4(:)
      Logical, optional :: flog
      Real(float) :: dum1(size(A1),1),dum2(size(A2),1),dum4(size(A4),1)
      dum1(:,1) = A1;  dum2(:,1) = A2;  dum4(:,1) = A4
      call ArrayPrint_aaaa(title,dum1,dum2,A3,dum4,fmtx,flog)
   End Subroutine ArrayPrint_vvav
   ! -------------------------------------------------------
   Subroutine ArrayPrint_vaav(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:),A2(:,:),A3(:,:),A4(:)
      Logical, optional :: flog
      Real(float) :: dum1(size(A1),1),dum4(size(A4),1)
      dum1(:,1) = A1;  dum4(:,1) = A4
      call ArrayPrint_aaaa(title,dum1,A2,A3,dum4,fmtx,flog)
   End Subroutine ArrayPrint_vaav
   ! -------------------------------------------------------
   Subroutine ArrayPrint_vavv(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:),A2(:,:),A3(:),A4(:)
      Logical, optional :: flog
      Real(float) :: dum1(size(A1),1),dum3(size(A3),1),dum4(size(A4),1)
      dum1(:,1) = A1;  dum3(:,1) = A3;  dum4(:,1) = A4
      call ArrayPrint_aaaa(title,dum1,A2,dum3,dum4,fmtx,flog)
   End Subroutine ArrayPrint_vavv
   ! -------------------------------------------------------
   Subroutine ArrayPrint_avvv(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(float), intent(in) :: A1(:,:),A2(:),A3(:),A4(:)
      Logical, optional :: flog
      Real(float) :: dum2(size(A2),1),dum3(size(A3),1),dum4(size(A4),1)
      dum2(:,1) = A2;  dum3(:,1) = A3;  dum4(:,1) = A4
      call ArrayPrint_aaaa(title,A1,dum2,dum3,dum4,fmtx,flog)
   End Subroutine ArrayPrint_avvv
   ! -------------------------------------------------------
   Subroutine ArrayPrintI(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Integer, intent(in) :: A1(:,:)
      Integer, optional :: A2(:,:),A3(:,:),A4(:,:)
      Logical, optional :: flog
      integer :: iunit
      Character :: fmtt*(80)
      integer :: i,j,n,n1(2),n2(2),n3(2),n4(2)
      if(Present(fmtx))then
         fmtt = fmtx
      else
         fmtt = fmti
      endif
      iunit = iout
      if(present(flog))then
         iunit = iout2
         if(flog)iunit = ilog
      endif
      if(present(title))then
         if(title .ne. ' ')write(iunit,'(a)')title
      endif
      n1 = shape(A1)
      if(Present(A2).and.Present(A3).and.Present(A4))then
         n2 = shape(A2); n3 = shape(A3); n4 = shape(A4); n = min(n1(1),n2(1),n3(1),n4(1))
         do i=1,n
            write(iunit,fmtt)(A1(i,j),tb,j=1,n1(2)),(A2(i,j),tb,j=1,n2(2)), &
               (A3(i,j),tb,j=1,n3(2)),(A4(i,j),tb,j=1,n4(2))
         enddo
      elseif(Present(A2).and.Present(A3))then
         n2 = shape(A2);  n3 = shape(A3);  n = min(n1(1),n2(1),n3(1))
         do i=1,n
            write(iunit,fmtt)(A1(i,j),tb,j=1,n1(2)),(A2(i,j),tb,j=1,n2(2)), &
               (A3(i,j),tb,j=1,n3(2))
         enddo
      elseif(Present(A2))then
         n2 = shape(A2); n = min(n2(1),n1(1))
         do i=1,n
            write(iunit,fmtt)(A1(i,j),tb,j=1,n1(2)),(A2(i,j),tb,j=1,n2(2))
         enddo
      else
         do i=1,n1(1)
            write(iunit,fmtt)(A1(i,j),tb,j=1,n1(2))
         enddo
      endif
   End Subroutine ArrayPrintI
   ! -------------------------------------------------------
   Subroutine ArrayPrintQ_aaaa(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(floatq), intent(in) :: A1(:,:)
      Real(floatq), optional :: A2(:,:),A3(:,:),A4(:,:)
      Logical, optional :: flog
      integer :: iunit
      Character :: fmtt*(80)
      integer :: i,j,n,n1(2),n2(2),n3(2),n4(2)
      ! ---------------------------------------------------
      ! Print up to 4 arrays (strung together side by side)
      !  title - is an optional title
      !     A1,A2,A3,A4 - arrays, only 1 required
      !  fmtx is optional format
      !  flog optional file designator as follows
      !     = .true. use ilog file
      !     = .false. use iout2 file
      !     = .not.present use iout file
      ! ---------------------------------------------------
      if(Present(fmtx))then
         fmtt = fmtx
      else
         fmtt = fmtd
      endif
      iunit = iout
      if(present(flog))then
         iunit = iout2
         if(flog)iunit = ilog
      endif
      if(present(title))then
         if(title .ne. ' ')write(iunit,'(a)')title
      endif
      n1 = shape(A1)
      if(Present(A2).and.Present(A3).and.Present(A4))then
         n2 = shape(A2); n3 = shape(A3); n4 = shape(A4); n = min(n1(1),n2(1),n3(1),n4(1))
         do i=1,n
            write(iunit,fmtt)(A1(i,j),tb,j=1,n1(2)),(A2(i,j),tb,j=1,n2(2)), &
               (A3(i,j),tb,j=1,n3(2)),(A4(i,j),tb,j=1,n4(2))
         enddo
      elseif(Present(A2).and.Present(A3))then
         n2 = shape(A2);  n3 = shape(A3);  n = min(n1(1),n2(1),n3(1))
         do i=1,n
            write(iunit,fmtt)(A1(i,j),tb,j=1,n1(2)),(A2(i,j),tb,j=1,n2(2)), &
               (A3(i,j),tb,j=1,n3(2))
         enddo
      elseif(Present(A2))then
         n2 = shape(A2); n = min(n2(1),n1(1))
         do i=1,n
            write(iunit,fmtt)(A1(i,j),tb,j=1,n1(2)),(A2(i,j),tb,j=1,n2(2))
         enddo
      else
         do i=1,n1(1)
            write(iunit,fmtt)(A1(i,j),tb,j=1,n1(2))
         enddo
      endif
   End Subroutine ArrayPrintQ_aaaa
   ! -------------------------------------------------------
   Subroutine ArrayPrintQ_vaaa(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(floatq), intent(in) :: A1(:)
      Real(floatq), optional :: A2(:,:),A3(:,:),A4(:,:)
      Logical, optional :: flog
      Real(floatq) :: dum(size(A1),1)
      dum(:,1) = A1
      call ArrayPrintQ_aaaa(title,dum,A2,A3,A4,fmtx,flog)
   End Subroutine ArrayPrintQ_vaaa
   ! -------------------------------------------------------
   Subroutine ArrayPrintQ_vava(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(floatq), intent(in) :: A1(:),A2(:,:),A3(:)
      Real(floatq), optional :: A4(:,:)
      Logical, optional :: flog
      Real(floatq) :: dum(size(A1),1),dum3(size(A3),1)
      dum(:,1) = A1; dum3(:,1) = A3
      call ArrayPrintQ_aaaa(title,dum,A2,dum3,A4,fmtx,flog)
   End Subroutine ArrayPrintQ_vava
   ! -------------------------------------------------------
   Subroutine ArrayPrintQ_vvaa(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(floatq), intent(in) :: A1(:),A2(:)
      Real(floatq), optional :: A3(:,:),A4(:,:)
      Logical, optional :: flog
      Real(floatq) :: dum1(size(A1),1),dum2(size(A2),1)
      dum1(:,1) = A1;  dum2(:,1) = A2
      call ArrayPrintQ_aaaa(title,dum1,dum2,A3,A4,fmtx,flog)
   End Subroutine ArrayPrintQ_vvaa
   ! -------------------------------------------------------
   Subroutine ArrayPrintQ_vvva(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(floatq), intent(in) :: A1(:),A2(:),A3(:)
      Real(floatq), optional :: A4(:,:)
      Logical, optional :: flog
      Real(floatq) :: dum1(size(A1),1),dum2(size(A2),1),dum3(size(A3),1)
      dum1(:,1) = A1;  dum2(:,1) = A2;  dum3(:,1) = A3
      call ArrayPrintQ_aaaa(title,dum1,dum2,dum3,A4,fmtx,flog)
   End Subroutine ArrayPrintQ_vvva
   ! -------------------------------------------------------
   Subroutine ArrayPrintQ_vvvv(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(floatq), intent(in) :: A1(:),A2(:),A3(:),A4(:)
      Logical, optional :: flog
      Real(floatq) :: dum1(size(A1),1),dum2(size(A2),1),dum3(size(A3),1),dum4(size(A4),1)
      dum1(:,1) = A1;  dum2(:,1) = A2;  dum3(:,1) = A3;  dum4(:,1) = A4
      call ArrayPrintQ_aaaa(title,dum1,dum2,dum3,dum4,fmtx,flog)
   End Subroutine ArrayPrintQ_vvvv
   ! -------------------------------------------------------
   Subroutine ArrayPrintQ_avaa(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(floatq), intent(in) :: A1(:,:),A2(:)
      Real(floatq), optional :: A3(:,:),A4(:,:)
      Logical, optional :: flog
      Real(floatq) :: dum2(size(A2),1)
      dum2(:,1) = A2
      call ArrayPrintQ_aaaa(title,A1,dum2,A3,A4,fmtx,flog)
   End Subroutine ArrayPrintQ_avaa
   ! -------------------------------------------------------
   Subroutine ArrayPrintQ_avva(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(floatq), intent(in) :: A1(:,:),A2(:),A3(:)
      Real(floatq), optional :: A4(:,:)
      Logical, optional :: flog
      Real(floatq) :: dum2(size(A2),1),dum3(size(A3),1)
      dum2(:,1) = A2;  dum3(:,1) = A3
      call ArrayPrintQ_aaaa(title,A1,dum2,dum3,A4,fmtx,flog)
   End Subroutine ArrayPrintQ_avva
   ! -------------------------------------------------------
   Subroutine ArrayPrintQ_aava(title,A1,A2,A3,A4,fmtx,flog)
      Character, optional :: title*(*),fmtx*(*)
      Real(floatq), intent(in) :: A1(:,:),A2(:,:),A3(:)
      Real(floatq), optional :: A4(:,:)
      Logical, optional :: flog
      Real(floatq) :: dum3(size(A3),1)
      dum3(:,1) = A3
      call ArrayPrintQ_aaaa(title,A1,A2,dum3,A4,fmtx,flog)
   End Subroutine ArrayPrintQ_aava
   ! -------------------------------------------------------
End Module Array_Print

