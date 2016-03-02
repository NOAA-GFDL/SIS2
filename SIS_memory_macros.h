!***********************************************************************
! This is a header file to define macros for static and dynamic memory *
! allocation.  Define STATIC_MEMORY_ in SIS_memory.h for static memory *
! allocation.  Otherwise dynamic memory allocation will be assumed.    *
!***********************************************************************

#ifdef STATIC_MEMORY_

#  define SZCAT_(IG)     NCAT_ICE_
#  define SZCAT0_(IG)    0:NCAT_ICE_
#  define SZK_ICE_(IG)   NK_ICE_
#  define SZK_SNOW_(IG)  NK_SNOW_

#else
! Dynamic memory allocation

#  define SZCAT_(IG)     IG%CatIce
#  define SZCAT0_(IG)    0:IG%CatIce
#  define SZK_ICE_(IG)   IG%NkIce
#  define SZK_SNOW_(IG)  IG%NkSnow

#endif
