!***********************************************************************
! This is a header file to define macros for static and dynamic memory *
! allocation.  Define STATIC_MEMORY_ in SIS_memory.h for static memory *
! allocation.  Otherwise dynamic memory allocation will be assumed.    *
!***********************************************************************

#ifdef STATIC_MEMORY_
#  define NCATMEM_      NCAT_ICE_
#  define NCATMEM0_     0:NCAT_ICE_
#  define NKICEMEM_     NK_ICE_
#  define NKSNOWMEM_    NK_SNOW_

#  define SZCAT_(G)     NCAT_ICE_
#  define SZCAT0_(G)    0:NCAT_ICE_
#  define SZK_ICE_(G)   NK_ICE_
#  define SZK_SNOW_(G)   NK_SNOW_

#else
! Dynamic memory allocation

#  define NCATMEM_    :
#  define NCATMEM0_  0:
#  define NKICEMEM_   :
#  define NKSNOWMEM_  :

#  define SZCAT_(G)     G%CatIce
#  define SZCAT0_(G)    0:G%CatIce
#  define SZK_ICE_(G)   G%NkIce
#  define SZK_SNOW_(G)  G%NkSnow

#endif
