      SUBROUTINE FGNMOD ( name, namlen, xdata, xin, xout, xvar,
     1                    iniflg, ierflg)
      IMPLICIT REAL*8 (A-H, O-Z),  INTEGER*4 (I-N)
      DIMENSION xdata(*), xin(*), xout(*), xvar(*)
      CHARACTER*1 name(*)
      PARAMETER ( namcnt = 8 )
      CHARACTER*32 refnam(namcnt)
      CONTINUE !  --------------------------------------------------
      CONTINUE !  You may increase namcnt above to allow more names:
      CONTINUE !  --------------------------------------------------
      CONTINUE !  In the following lines, register your foreign model
      CONTINUE !  names as they are declared in your models:
      CONTINUE !   - use only uppercase characters for the name here
      CONTINUE !   - you can use any case for the name in the models
      CONTINUE !   - make a copy of the modifications you make to this
      CONTINUE !     file so that you don't lose them when installing
      CONTINUE !     a newer version of ATP later
      DATA refnam(1) / 'SAMPLE_MODEL' /  ! Do not modify this line
      DATA refnam(2) / 'ULM_LINE' /
      DATA refnam(3) / ' ' /
      DATA refnam(4) / ' ' /
      DATA refnam(5) / ' ' /
      DATA refnam(6) / ' ' /
      DATA refnam(7) / ' ' /
      DATA refnam(8) / ' ' /
      CONTINUE !  --------------------------------------------------
      CONTINUE !  Name identification loop
      CONTINUE !  -- no need to change anything here
      iname = 1
      lpflg = 1
      DO WHILE (iname.LE.namcnt .AND. lpflg.GT.0)
       ICHAR = 1
       DO WHILE (ichar.LE.namlen
     1           .AND. name(ichar).EQ.refnam(iname)(ichar:ichar))
        ichar = ichar + 1
       ENDDO
       IF (ichar.GT.namlen) THEN
        lpflg = 0
       ELSE
        iname = iname + 1
       ENDIF
      ENDDO
      IF (iname.GT.namcnt) THEN
       ierflg = 1
       RETURN
      ENDIF
      CONTINUE !  --------------------------------------------------
      CONTINUE !  In the following lines, this is where you call the
      CONTINUE !  actual foreign subroutines/procedures:
      CONTINUE !   - actual names may be different from the foreign
      CONTINUE !     names used in the models
      CONTINUE !   - notice how each one uses both an
      CONTINUE !     initialization routine and an execution routine
      IF ( iname.EQ.1 ) THEN
       IF (iniflg.EQ.1) THEN
        CALL sampli(xdata, xin, xout, xvar)
       ELSE
        CALL samplm(xdata, xin, xout, xvar)
       ENDIF
      CONTINUE !      -------------------------------------------
      ELSE IF ( iname.EQ.2 ) THEN
       IF (iniflg.EQ.1) THEN
        CALL ulm_i(xdata, xin, xout, xvar)
       ELSE
        CALL ulm_m(xdata, xin, xout, xvar)
       ENDIF
      CONTINUE !      -------------------------------------------
      ELSE IF ( iname.EQ.3 ) THEN
       IF (iniflg.EQ.1) THEN
       ELSE
       ENDIF
      CONTINUE !      -------------------------------------------
      ELSE IF ( iname.EQ.4 ) THEN
       IF (iniflg.EQ.1) THEN
       ELSE
       ENDIF
      CONTINUE !      -------------------------------------------
      ELSE IF ( iname.EQ.5 ) THEN
       IF (iniflg.EQ.1) THEN
       ELSE
       ENDIF
      CONTINUE !      -------------------------------------------
      ELSE IF ( iname.EQ.6 ) THEN
       IF (iniflg.EQ.1) THEN
       ELSE
       ENDIF
      CONTINUE !      -------------------------------------------
      ELSE IF ( iname.EQ.7 ) THEN
       IF (iniflg.EQ.1) THEN
       ELSE
       ENDIF
      CONTINUE !      -------------------------------------------
      ELSE IF ( iname.EQ.8 ) THEN
       IF (iniflg.EQ.1) THEN
       ELSE
       ENDIF
      ENDIF
      CONTINUE !      -------------------------------------------
      RETURN
      END
      SUBROUTINE FGNFUN ( name, namlen, xarg, nval, ierflg)
      IMPLICIT REAL*8 (A-H, O-Z),  INTEGER*4 (I-N)
      DIMENSION xarg(*)
      CHARACTER*1 name(*)
      PARAMETER ( namcnt = 8 )
      CHARACTER*32 refnam(namcnt)
      CONTINUE !  --------------------------------------------------
      CONTINUE !  You may increase namcnt above to allow more names:
      CONTINUE !  --------------------------------------------------
      CONTINUE !  In the following lines, register your foreign function
      CONTINUE !  names as they are declared in your models:
      CONTINUE !   - use only uppercase characters for the name here
      CONTINUE !   - you can use any case for the name in the models
      CONTINUE !   - make a copy of the modifications you make to this
      CONTINUE !     file so that you don't lose them when installing
      CONTINUE !     a newer version of ATP later
      DATA refnam(1) / 'SAMPLE_FUNCTION' /  ! Do not modify this line
      DATA refnam(2) / ' ' /
      DATA refnam(3) / ' ' /
      DATA refnam(4) / ' ' /
      DATA refnam(5) / ' ' /
      DATA refnam(6) / ' ' /
      DATA refnam(7) / ' ' /
      DATA refnam(8) / ' ' /
      CONTINUE !  --------------------------------------------------
      CONTINUE !  Name identification loop
      CONTINUE !  -- no need to change anything here
      iname = 1
      lpflg = 1
      DO WHILE (iname.LE.namcnt .AND. lpflg.GT.0)
       ichar = 1
       DO WHILE (ichar.LE.namlen
     1           .AND. name(ichar).EQ.refnam(iname)(ichar:ichar))
        ichar = ichar + 1
       ENDDO
       IF (ichar.GT.namlen) THEN
        lpflg = 0
       ELSE
        iname = iname + 1
       ENDIF
      ENDDO
      IF (iname.GT.namcnt) THEN
       ierflg = 1
       RETURN
      ENDIF
      CONTINUE !  --------------------------------------------------
      CONTINUE !  In the following lines, this is where you call the
      CONTINUE !  actual foreign functions:
      CONTINUE !   - actual names may be different from the foreign
      CONTINUE !     names used in the models
      IF      ( iname.EQ.1 ) THEN
                                  nval = samplf(xarg)
      ELSE IF ( iname.EQ.2 ) THEN
                                  nval = 0
      ELSE IF ( iname.EQ.3 ) THEN
                                  nval = 0
      ELSE IF ( iname.EQ.4 ) THEN
                                  nval = 0
      ELSE IF ( iname.EQ.5 ) THEN
                                  nval = 0
      ELSE IF ( iname.EQ.6 ) THEN
                                  nval = 0
      ELSE IF ( iname.EQ.7 ) THEN
                                  nval = 0
      ELSE IF ( iname.EQ.8 ) THEN
                                  nval = 0
      ENDIF
      CONTINUE !      -------------------------------------------
      RETURN
      END
      SUBROUTINE SAMPLM (xdata, xin, xout, xvar)
      IMPLICIT REAL*8 (A-H, O-Z),  INTEGER*4 (I-N)
      DIMENSION xdata(*), xin(*), xout(*), xvar(*)
      CHARACTER*80 text80  ! Buffer used to assemble messages for output
      CONTINUE !-------------------------------------------------------!
      CONTINUE ! This is a sample user-defined foreign model           !
      CONTINUE !-------------------------------------------------------!
      CONTINUE ! For this particular model:                            !
      CONTINUE !  - the dimension of xin,xout,xvar is given in xdata(1)!
      CONTINUE !  - xvar is to be initialized with HISTDEF             !
      CONTINUE !    in the USE statement                               !
      CONTINUE !-------------------------------------------------------!
      CONTINUE ! For any foreign model:                                !
      CONTINUE !  - the user can assign history to xin and xvar        !
      CONTINUE !    using HISTDEF in the USE statement                 !
      CONTINUE !  - MODELS always saves the values that are in xvar    !
      CONTINUE !    between uses, so that the memory of multiple       !
      CONTINUE !    instances of using the code can be managed         !
      CONTINUE !    automatically by MODELS                            !
      CONTINUE !-------------------------------------------------------!
      CONTINUE ! To introduce a foreign model in a simulation:         !
      CONTINUE !  - create an initialization and an execution routine  !
      CONTINUE !    for that model                                     !
      CONTINUE !  - provide access to your routines by modifying       !
      CONTINUE !    the routine fgnmod at the beginning of this file   !
      CONTINUE !  - declare and use the foreign model in a model       !
      CONTINUE !    of the simulation                                  !
      CONTINUE !-------------------------------------------------------!
      CONTINUE ! This sample foreign model could be used in a model    !
      CONTINUE ! as follows:                                           !
      CONTINUE !                                                       !
      CONTINUE !   MODEL anyname                                       !
      CONTINUE !   INPUT ...                                           !
      CONTINUE !   VAR a[1..3], ...                                    !
      CONTINUE !   ...                                                 !
      CONTINUE !   MODEL mymodel FOREIGN sample_model                  !
      CONTINUE !    followed by values for ixdata, ixin, ixout,        !
      CONTINUE !    and ixvar a pair of inside brackets                !
      CONTINUE !   -- the preceding line declares foreign model "sub"  !
      CONTINUE !   -- "mymodel" is the local name of the foreign model !
      CONTINUE !   -- "sample_model" is the identifier used for that   !
      CONTINUE !   --  foreign model in the user-modifiable            !
      CONTINUE !   --  MODELS interface routine "fgnmod"               !
      CONTINUE !   ...                                                 !
      CONTINUE !   EXEC                                                !
      CONTINUE !    ...                                                !
      CONTINUE !    USE mymodel AS mymodel                             !
      CONTINUE !     DATA xdata[1]:=3                                  !
      CONTINUE !     INPUT xin[1..3]:=[ list of 3 expressions... ]     !
      CONTINUE !     OUTPUT a[1..3]:=xout[1..3]                        !
      CONTINUE !     HISTORY histdef(xvar[1..3]):=0                    !
      CONTINUE !      -- can be any expression                         !
      CONTINUE !      -- note the fixed names xdata, xin, xout, xvar   !
      CONTINUE !    ENDUSE                                             !
      CONTINUE !    ...                                                !
      CONTINUE !   ENDEXEC                                             !
      CONTINUE !   ENDMODEL                                            !
      CONTINUE !-------------------------------------------------------!
      text80 = 'Executing model "sample_model".'   ! Assemble output mes
      CALL OUTSIX( text80, 80 )     ! Send text80 (arg 1) to ATP listing
      N8 = xdata(1)
      DO 1000 i=1, N8
       xvar(i)=(2*xvar(i)+xin(i))/2
       xout(i)=xvar(i)+100.0
 1000 CONTINUE
      RETURN
      ENTRY sampli(xdata, xin, xout, xvar)
      text80 = 'Initializing model "sample_model".'
      CALL OUTSIX( text80, 80 )
      RETURN
      END
      FUNCTION SAMPLF (arg)
      IMPLICIT REAL*8 (A-H, O-Z),  INTEGER*4 (I-N)
      DIMENSION arg(2)
      CONTINUE !-------------------------------------------------------!
      CONTINUE ! This is a sample user-defined foreign function        !
      CONTINUE !-------------------------------------------------------!
      CONTINUE ! For any foreign function:                             !
      CONTINUE !  - the function can use as many arguments as it       !
      CONTINUE !    expects to receive                                 !
      CONTINUE !  - the function receives its argument in an array     !
      CONTINUE !  - the function loads its output values in            !
      CONTINUE !    the same array                                     !
      CONTINUE !  - the function returns the number of output values   !
      CONTINUE !    as its value                                       !
      CONTINUE !-------------------------------------------------------!
      CONTINUE ! To introduce a foreign function in a simulation:      !
      CONTINUE !   - create a function procedure                       !
      CONTINUE !   - provide access to your functions by modifying     !
      CONTINUE !      the routine fgnfun at the beginning of this file !
      CONTINUE !   - declare and use the foreign function in a model   !
      CONTINUE !      of the simulation                                !
      CONTINUE !-------------------------------------------------------!
      CONTINUE ! This function could be declared in a model as follows:!
      CONTINUE !    FUNCTION myfunc FOREIGN sample_function            !
      CONTINUE !     followed by ixarg:2 inside a pair of brackets     !
      CONTINUE ! which defines the function "myfunc" to be a foreign   !
      CONTINUE ! function where:                                       !
      CONTINUE !  - "myfunc" is the local name of the foreign function !
      CONTINUE !  - "sample_function" is the identifier used for that  !
      CONTINUE !     foreign funcion in the user-modifiable            !
      CONTINUE !     MODELS interface routine "fgnfun"                 !
      CONTINUE !  - "ixarg" is a fixed keyword indicating              !
      CONTINUE !     the argument count                                !
      CONTINUE !-------------------------------------------------------!
      CONTINUE ! This function could be used in a model as follows:    !
      CONTINUE !    a[1..2] := 3*myfunc( 36.3, sqrt(c) )               !
      CONTINUE !-------------------------------------------------------!
      arg(1)=arg(1) +arg(2)
      arg(2)=10*arg(1)
      samplf=2
      RETURN
      END
