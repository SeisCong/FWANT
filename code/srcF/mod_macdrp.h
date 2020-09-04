/****************************************************************************
 * This file is used for defining fd operators.                             *
 *                                                                          *
 * Author: Wei ZHANG     Email: zhangwei.zw@gmail.com                       *
 * Copyright (C) 2006 Wei ZHANG                                             *
 *                                                                          *
 * $Date: 2009-01-14 16:30:51 -0500 (Wed, 14 Jan 2009) $                    *
 * $Revision: 510 $                                                         *
 * $LastChangedBy: zhangw $                                                 *
*****************************************************************************/

/* ! 8th order finite difference
 ! #define FDWET real,parameter :: C0=9.0/8.0, C1=1.0/24.0,C2=0.0,C3=0.0
 ! #define FDWET real,parameter :: C0=1225.0/1024.0,C1=245.0/3072.0,C2=49.0/5120.0,C3=5.0/7168.0 */
/* #define FDOPLEN integer,parameter :: LenFD = 3 */
/* #define DECLAREVmod real,dimension(2*LenFD) :: Vmod */

#define DEFCFDWET real(SP),parameter,private :: C0=0.0_SP,C1=0.7709_SP,C2=-0.1667_SP,C3=0.02084_SP
#define DEFFDWET real(SP),parameter,private :: MCA_1=-0.30874_SP,MCA0=-0.6326_SP,MCA1=1.233_SP,MCA2=-0.3334_SP,MCA3=0.04168_SP
/* #define DEFFDWET real(SP),parameter,private :: MCA_1=-9.0_SP/30.0_SP,MCA0=-19.0_SP/30.0_SP,MCA1=36.0_SP/30.0_SP,MCA2=-9.0/30.0_SP,MCA3=1.0_SP/30.0_SP */
#define LenFD 3
#define LenFDS 1
#define LenFDL 3

#ifdef DataTypeDouble
#define SEISNC_DATATYPE NF90_DOUBLE
#define SEISMPI_DATATYPE MPI_DOUBLE_PRECISION
#else
#define SEISNC_DATATYPE NF90_FLOAT
#define SEISMPI_DATATYPE MPI_REAL
#endif

/**************************
 * 4-6 LDDRK coeffecients *
 **************************/
#define DEFLDDRK2A real(SP),parameter :: RK2a2=1.0_SP
#define DEFLDDRK2B real(SP),parameter :: RK2b1=1.0_SP/2.0_SP,RK2b2=1.0_SP/2.0_SP

#define DEFLDDRK4A real(SP),parameter :: RK4a2=0.5_SP,RK4a3=0.5_SP,RK4a4=1.0_SP
#define DEFLDDRK4B real(SP),parameter :: RK4b1=1.0_SP/6.0_SP,RK4b2=1.0_SP/3.0_SP,RK4b3=1.0_SP/3.0_SP,RK4b4=1.0_SP/6.0_SP

#define DEFLDDRK6A1 real(SP),parameter :: RK6a2=0.353323_SP,RK6a3=0.999597_SP,RK6a4=0.152188_SP
#define DEFLDDRK6A2 real(SP),parameter :: RK6a5=0.534216_SP,RK6a6=0.603907_SP
#define DEFLDDRK6B1 real(SP),parameter :: RK6b1=0.0467621_SP,RK6b2=0.137286_SP,RK6b3=0.170975_SP
#define DEFLDDRK6B2 real(SP),parameter :: RK6b4=0.197572_SP,RK6b5=0.282263_SP,RK6b6=0.165142_SP

/****************************************************************************************
 * #define DEFLDDRK6A1 real,parameter :: RK6a2=0.353323,RK6a3=0.353323,RK6a4=0.240823   *
 * #define DEFLDDRK6A2 real,parameter :: RK6a5=0.240823,RK6a6=0.341148                  *
 * #define DEFLDDRK6B1 real,parameter :: RK6b1=-0.766927,RK6b2=-0.519328,RK6b3=0.147469 *
 * #define DEFLDDRK6B2 real,parameter :: RK6b4=-0.140084,RK6b5=1.11946,RK6b6=1.15941    *
 ****************************************************************************************/

/*********************************************************
 * ! DRP FDx and FDy for 3d variables     *
 *********************************************************/
#define m3d_FDxC1(var,i,j,k)  (-C3*var(i-3,j,k)-C2*var(i-2,j,k)-C1*var(i-1,j,k)
#define m3d_FDxC2(var,i,j,k)   +C1*var(i+1,j,k)+C2*var(i+2,j,k)+C3*var(i+3,j,k))

#define m3d_FDyC1(var,i,j,k)  (-C3*var(i,j-3,k)-C2*var(i,j-2,k)-C1*var(i,j-1,k)
#define m3d_FDyC2(var,i,j,k)   +C1*var(i,j+1,k)+C2*var(i,j+2,k)+C3*var(i,j+3,k))

#define m3d_FDzC1(var,i,j,k)  (-C3*var(i,j,k-3)-C2*var(i,j,k-2)-C1*var(i,j,k-1)
#define m3d_FDzC2(var,i,j,k)   +C1*var(i,j,k+1)+C2*var(i,j,k+2)+C3*var(i,j,k+3))

/*********************************************************
 * ! DRP/opt MacCormack FDx and FDy for 3d variables     *
 *********************************************************/
#define m3d_FDxF1(var,i,j,k)  (MCA_1*var(i-1,j,k)+MCA0*var(i,j,k)
#define m3d_FDxF2(var,i,j,k)   +MCA1*var(i+1,j,k)+MCA2*var(i+2,j,k)+MCA3*var(i+3,j,k))

#define m3d_FDxB1(var,i,j,k)  (-MCA3*var(i-3,j,k)-MCA2*var(i-2,j,k)
#define m3d_FDxB2(var,i,j,k)  -MCA1*var(i-1,j,k)-MCA0*var(i,j,k)-MCA_1*var(i+1,j,k))

#define m3d_FDyF1(var,i,j,k)  (MCA_1*var(i,j-1,k)+MCA0*var(i,j,k)
#define m3d_FDyF2(var,i,j,k)   +MCA1*var(i,j+1,k)+MCA2*var(i,j+2,k)+MCA3*var(i,j+3,k))

#define m3d_FDyB1(var,i,j,k)  (-MCA3*var(i,j-3,k)-MCA2*var(i,j-2,k)
#define m3d_FDyB2(var,i,j,k)   -MCA1*var(i,j-1,k)-MCA0*var(i,j,k)-MCA_1*var(i,j+1,k))

#define m3d_FDzF1(var,i,j,k)  (MCA_1*var(i,j,k-1)+MCA0*var(i,j,k)
#define m3d_FDzF2(var,i,j,k)   +MCA1*var(i,j,k+1)+MCA2*var(i,j,k+2)+MCA3*var(i,j,k+3))

#define m3d_FDzB1(var,i,j,k)  (-MCA3*var(i,j,k-3)-MCA2*var(i,j,k-2)
#define m3d_FDzB2(var,i,j,k)   -MCA1*var(i,j,k-1)-MCA0*var(i,j,k)-MCA_1*var(i,j,k+1))

#define m4d_FDxF1(var,i,j,k)  (MCA_1*var(m,i-1,j,k)+MCA0*var(m,i,j,k)
#define m4d_FDxF2(var,i,j,k)   +MCA1*var(m,i+1,j,k)+MCA2*var(m,i+2,j,k)+MCA3*var(m,i+3,j,k))

#define m4d_FDxB1(var,i,j,k)  (-MCA3*var(m,i-3,j,k)-MCA2*var(m,i-2,j,k)
#define m4d_FDxB2(var,i,j,k)  -MCA1*var(m,i-1,j,k)-MCA0*var(m,i,j,k)-MCA_1*var(m,i+1,j,k))

#define m4d_FDyF1(var,i,j,k)  (MCA_1*var(m,i,j-1,k)+MCA0*var(m,i,j,k)
#define m4d_FDyF2(var,i,j,k)   +MCA1*var(m,i,j+1,k)+MCA2*var(m,i,j+2,k)+MCA3*var(m,i,j+3,k))

#define m4d_FDyB1(var,i,j,k)  (-MCA3*var(m,i,j-3,k)-MCA2*var(m,i,j-2,k)
#define m4d_FDyB2(var,i,j,k)   -MCA1*var(m,i,j-1,k)-MCA0*var(m,i,j,k)-MCA_1*var(m,i,j+1,k))

#define m4d_FDzF1(var,i,j,k)  (MCA_1*var(m,i,j,k-1)+MCA0*var(m,i,j,k)
#define m4d_FDzF2(var,i,j,k)   +MCA1*var(m,i,j,k+1)+MCA2*var(m,i,j,k+2)+MCA3*var(m,i,j,k+3))

#define m4d_FDzB1(var,i,j,k)  (-MCA3*var(m,i,j,k-3)-MCA2*var(m,i,j,k-2)
#define m4d_FDzB2(var,i,j,k)   -MCA1*var(m,i,j,k-1)-MCA0*var(m,i,j,k)-MCA_1*var(m,i,j,k+1))

/*********************************************************
 * ! DRP/opt MacCormack FDx and FDy for matrix variables *
 *********************************************************/
#define mat_FDxF1(var,i,j)  (MCA_1*var(i-1,j)+MCA0*var(i,j)
#define mat_FDxF2(var,i,j)   +MCA1*var(i+1,j)+MCA2*var(i+2,j)+MCA3*var(i+3,j))

#define mat_FDxB1(var,i,j)  (-MCA_1*var(i+1,j)-MCA0*var(i,j)
#define mat_FDxB2(var,i,j)   -MCA1*var(i-1,j)-MCA2*var(i-2,j)-MCA3*var(i-3,j))

#define mat_FDyF1(var,i,j)  (MCA_1*var(i,j-1)+MCA0*var(i,j)
#define mat_FDyF2(var,i,j)   +MCA1*var(i,j+1)+MCA2*var(i,j+2)+MCA3*var(i,j+3))

#define mat_FDyB1(var,i,j)  (-MCA_1*var(i,j+1)-MCA0*var(i,j)
#define mat_FDyB2(var,i,j)   -MCA1*var(i,j-1)-MCA2*var(i,j-2)-MCA3*var(i,j-3))

/**************************************
 * !DRP/opt MacCormack  FD for vecotr *
 **************************************/

#define vec_FD_F1(var,i)  (MCA_1*var(i-1)+MCA0*var(i)
#define vec_FD_F2(var,i)   +MCA1*var(i+1)+MCA2*var(i+2)+MCA3*var(i+3))

#define vec_FD_B1(var,i)  (-MCA_1*var(i+1)-MCA0*var(i)
#define vec_FD_B2(var,i)   -MCA1*var(i-1)-MCA2*var(i-2)-MCA3*var(i-3))


/********************************
 * ! 4/4 Compact MacCormack FDx *
 * all coeffecients should be nomaled by LA0
 ********************************/
/* #define HOCWETL real,parameter :: LA0=2.0, LA1=1.0
   #define HOCWETR real,parameter :: RA_1=-0.5, RA0=2.0, RA1=2.5 */
#define HOCWETL real(SP),parameter :: LA0=2.0_SP/3.0_SP, LA1=1.0_SP/3.0_SP
#define HOCWETR real(SP),parameter :: RA_1=-1.0_SP/6.0_SP, RA0=-2.0_SP/3.0_SP, RA1=5.0_SP/6.0_SP

#define mat_HOCxF1_RHS(var,i,j)  (RA_1*var(i-1,j)+RA0*var(i,j)
#define mat_HOCxF2_RHS(var,i,j)   +RA1*var(i+1,j))/LA0

#define mat_HOCxB1_RHS(var,i,j)  (-RA1*var(i-1,j)-RA0*var(i,j)
#define mat_HOCxB2_RHS(var,i,j)   -RA_1*var(i+1,j))/LA0

#define mat_HOCyF1_RHS(var,i,j)  (RA_1*var(i,j-1)+RA0*var(i,j)
#define mat_HOCyF2_RHS(var,i,j)   +RA1*var(i,j+1))/LA0

#define mat_HOCyB1_RHS(var,i,j)  (-RA1*var(i,j-1)-RA0*var(i,j)
#define mat_HOCyB2_RHS(var,i,j)   -RA_1*var(i,j+1))/LA0

#define m3d_HOCzF1_RHS(var,i,j,k)  (RA_1*var(i,j,k-1)+RA0*var(i,j,k)
#define m3d_HOCzF2_RHS(var,i,j,k)   +RA1*var(i,j,k+1))/LA0
#define m3d_HOCzB1_RHS(var,i,j,k)  (-RA1*var(i,j,k-1)-RA0*var(i,j,k)
#define m3d_HOCzB2_RHS(var,i,j,k)   -RA_1*var(i,j,k+1))/LA0

#define vec_HOC_F_LHS(var,i)  (LA1*var(i+1))/LA0
#define vec_HOC_B_LHS(var,i)  (LA1*var(i-1))/LA0

/**************************************
 * filter *
 **************************************/
#define DefFilterWet real(SP),parameter :: MCF0=164.0_SP/256.0_SP,MCF1=63.0_SP/256.0_SP,MCF2=-18.0_SP/256.0_SP,MCF3=1.0_SP/256.0_SP
#define m3d_FILTx1(var,i,j,k)  (MCF3*var(i-3,j,k)+MCF2*var(i-2,j,k)+MCF1*var(i-1,j,k)
#define m3d_FILTx2(var,i,j,k)   +MCF0*var(i,j,k)+MCF1*var(i+1,j,k)+MCF2*var(i+2,j,k)+MCF3*var(i+3,j,k))

#define m3d_FILTy1(var,i,j,k)  (MCF3*var(i,j-3,k)+MCF2*var(i,j-2,k)+MCF1*var(i,j-1,k)
#define m3d_FILTy2(var,i,j,k)   +MCF0*var(i,j,k)+MCF1*var(i,j+1,k)+MCF2*var(i,j+2,k)+MCF3*var(i,j+3,k))

#define m3d_FILTz1(var,i,j,k)  (MCF3*var(i,j,k-3)+MCF2*var(i,j,k-2)+MCF1*var(i,j,k-1)
#define m3d_FILTz2(var,i,j,k)   +MCF0*var(i,j,k)+MCF1*var(i,j,k+1)+MCF2*var(i,j,k+2)+MCF3*var(i,j,k+3))

#define m4d_FILTx1(var,i,j,k)  (MCF3*var(m,i-3,j,k)+MCF2*var(m,i-2,j,k)+MCF1*var(m,i-1,j,k)
#define m4d_FILTx2(var,i,j,k)   +MCF0*var(m,i,j,k)+MCF1*var(m,i+1,j,k)+MCF2*var(m,i+2,j,k)+MCF3*var(m,i+3,j,k))

#define m4d_FILTy1(var,i,j,k)  (MCF3*var(m,i,j-3,k)+MCF2*var(m,i,j-2,k)+MCF1*var(m,i,j-1,k)
#define m4d_FILTy2(var,i,j,k)   +MCF0*var(m,i,j,k)+MCF1*var(m,i,j+1,k)+MCF2*var(m,i,j+2,k)+MCF3*var(m,i,j+3,k))

#define m4d_FILTz1(var,i,j,k)  (MCF3*var(m,i,j,k-3)+MCF2*var(m,i,j,k-2)+MCF1*var(m,i,j,k-1)
#define m4d_FILTz2(var,i,j,k)   +MCF0*var(m,i,j,k)+MCF1*var(m,i,j,k+1)+MCF2*var(m,i,j,k+2)+MCF3*var(m,i,j,k+3))
/**************************************
 * filter *
 **************************************/
/*
#define DefFilterWet real(SP),parameter :: MCF0=20.0_SP/64.0_SP,MCF1=-15.0_SP/64.0_SP,MCF2=6.0_SP/64.0_SP,MCF3=-1.0_SP/64.0_SP
#define m3d_FILTx1(var,i,j,k)  (var(i,j,k)+(MCF3*var(i-3,j,k)+MCF2*var(i-2,j,k)+MCF1*var(i-1,j,k)
#define m3d_FILTx2(var,i,j,k)   +MCF0*var(i,j,k)+MCF1*var(i+1,j,k)+MCF2*var(i+2,j,k)+MCF3*var(i+3,j,k))/steph)

#define m3d_FILTy1(var,i,j,k)  (var(i,j,k)+(MCF3*var(i,j-3,k)+MCF2*var(i,j-2,k)+MCF1*var(i,j-1,k)
#define m3d_FILTy2(var,i,j,k)   +MCF0*var(i,j,k)+MCF1*var(i,j+1,k)+MCF2*var(i,j+2,k)+MCF3*var(i,j+3,k))/steph)

#define m3d_FILTz1(var,i,j,k)  (var(i,j,k)+(MCF3*var(i,j,k-3)+MCF2*var(i,j,k-2)+MCF1*var(i,j,k-1)
#define m3d_FILTz2(var,i,j,k)   +MCF0*var(i,j,k)+MCF1*var(i,j,k+1)+MCF2*var(i,j,k+2)+MCF3*var(i,j,k+3))/steph)

#define m4d_FILTx1(var,i,j,k)  (var(m,i,j,k)+(MCF3*var(m,i-3,j,k)+MCF2*var(m,i-2,j,k)+MCF1*var(m,i-1,j,k)
#define m4d_FILTx2(var,i,j,k)   +MCF0*var(m,i,j,k)+MCF1*var(m,i+1,j,k)+MCF2*var(m,i+2,j,k)+MCF3*var(m,i+3,j,k))/steph)

#define m4d_FILTy1(var,i,j,k)  (var(m,i,j,k)+(MCF3*var(m,i,j-3,k)+MCF2*var(m,i,j-2,k)+MCF1*var(m,i,j-1,k)
#define m4d_FILTy2(var,i,j,k)   +MCF0*var(m,i,j,k)+MCF1*var(m,i,j+1,k)+MCF2*var(m,i,j+2,k)+MCF3*var(m,i,j+3,k))/steph)

#define m4d_FILTz1(var,i,j,k)  (var(m,i,j,k)+(MCF3*var(m,i,j,k-3)+MCF2*var(m,i,j,k-2)+MCF1*var(m,i,j,k-1)
#define m4d_FILTz2(var,i,j,k)   +MCF0*var(m,i,j,k)+MCF1*var(m,i,j,k+1)+MCF2*var(m,i,j,k+2)+MCF3*var(m,i,j,k+3))/steph)
*/
/**************************************
 * damp *
 **************************************/
#define DEFDAMPWET real(SP),parameter :: MCD0=0.370563_SP,MCD1=-0.2411788_SP,MCD2=0.0647185_SP,MCD3=-0.0088212_SP

#define m3d_SMxF1(var,i,j,k)  (MCD3*var(i-3,j,k)+MCD2*var(i-2,j,k)+MCD1*var(i-1,j,k)
#define m3d_SMxF2(var,i,j,k)   +MCD0*var(i,j,k)+MCD1*var(i+1,j,k)+MCD2*var(i+2,j,k)+MCD3*var(i+3,j,k))

#define m3d_SMyF1(var,i,j,k)  (MCD3*var(i,j-3,k)+MCD2*var(i,j-2,k)+MCD1*var(i,j-1,k)
#define m3d_SMyF2(var,i,j,k)   +MCD0*var(i,j,k)+MCD1*var(i,j+1,k)+MCD2*var(i,j+2,k)+MCD3*var(i,j+3,k))

#define m3d_SMzF1(var,i,j,k)  (MCD3*var(i,j,k-3)+MCD2*var(i,j,k-2)+MCD1*var(i,j,k-1)
#define m3d_SMzF2(var,i,j,k)   +MCD0*var(i,j,k)+MCD1*var(i,j,k+1)+MCD2*var(i,j,k+2)+MCD3*var(i,j,k+3))

#define m4d_SMxF1(var,i,j,k)  (MCD3*var(m,i-3,j,k)+MCD2*var(m,i-2,j,k)+MCD1*var(m,i-1,j,k)
#define m4d_SMxF2(var,i,j,k)   +MCD0*var(m,i,j,k)+MCD1*var(m,i+1,j,k)+MCD2*var(m,i+2,j,k)+MCD3*var(m,i+3,j,k))

#define m4d_SMyF1(var,i,j,k)  (MCD3*var(m,i,j-3,k)+MCD2*var(m,i,j-2,k)+MCD1*var(m,i,j-1,k)
#define m4d_SMyF2(var,i,j,k)   +MCD0*var(m,i,j,k)+MCD1*var(m,i,j+1,k)+MCD2*var(m,i,j+2,k)+MCD3*var(m,i,j+3,k))

#define m4d_SMzF1(var,i,j,k)  (MCD3*var(m,i,j,k-3)+MCD2*var(m,i,j,k-2)+MCD1*var(m,i,j,k-1)
#define m4d_SMzF2(var,i,j,k)   +MCD0*var(m,i,j,k)+MCD1*var(m,i,j,k+1)+MCD2*var(m,i,j,k+2)+MCD3*var(m,i,j,k+3))

/*********************************************************
 * ! 2-4 MacCormack FDx and FDy for 3d variables     *
 *********************************************************/
#define DEFFDWET24 real(SP),parameter,private :: MAC24A0=-7.0_SP/6.0_SP,MAC24A1=8.0_SP/6.0_SP,MAC24A2=-1.0_SP/6.0_SP

#define m24_FDxF1(var,i,j,k)  (MAC24A0*var(i,j,k)
#define m24_FDxF2(var,i,j,k)   +MAC24A1*var(i+1,j,k)+MAC24A2*var(i+2,j,k))

#define m24_FDxB1(var,i,j,k)  (-MAC24A2*var(i-2,j,k)
#define m24_FDxB2(var,i,j,k)  -MAC24A1*var(i-1,j,k)-MAC24A0*var(i,j,k))

#define m24_FDyF1(var,i,j,k)  (MAC24A0*var(i,j,k)
#define m24_FDyF2(var,i,j,k)   +MAC24A1*var(i,j+1,k)+MAC24A2*var(i,j+2,k))

#define m24_FDyB1(var,i,j,k)  (-MAC24A2*var(i,j-2,k)
#define m24_FDyB2(var,i,j,k)   -MAC24A1*var(i,j-1,k)-MAC24A0*var(i,j,k))

#define m24_FDzF1(var,i,j,k)  (MAC24A0*var(i,j,k)
#define m24_FDzF2(var,i,j,k)   +MAC24A1*var(i,j,k+1)+MAC24A2*var(i,j,k+2))

#define m24_FDzB1(var,i,j,k)  (-MAC24A2*var(i,j,k-2)
#define m24_FDzB2(var,i,j,k)   -MAC24A1*var(i,j,k-1)-MAC24A0*var(i,j,k))

#define DEFFDWET22 real(SP),parameter,private :: MAC22A0=-1.0_SP,MAC22A1=1.0_SP

#define m22_FDxF1(var,i,j,k)  (MAC22A0*var(i,j,k)
#define m22_FDxF2(var,i,j,k)   +MAC22A1*var(i+1,j,k))

#define m22_FDxB1(var,i,j,k)  (-MAC22A1*var(i-1,j,k)
#define m22_FDxB2(var,i,j,k)  -MAC22A0*var(i,j,k))

#define m22_FDyF1(var,i,j,k)  (MAC22A0*var(i,j,k)
#define m22_FDyF2(var,i,j,k)   +MAC22A1*var(i,j+1,k))

#define m22_FDyB1(var,i,j,k)  (-MAC22A1*var(i,j-1,k)
#define m22_FDyB2(var,i,j,k)   -MAC22A0*var(i,j,k))

#define m22_FDzF1(var,i,j,k)  (MAC22A0*var(i,j,k)
#define m22_FDzF2(var,i,j,k)   +MAC22A1*var(i,j,k+1))

#define m22_FDzB1(var,i,j,k)  (-MAC22A1*var(i,j,k-1)
#define m22_FDzB2(var,i,j,k)   -MAC22A0*var(i,j,k))

/**************************************
 * parabolic *
 **************************************/

#define m3d_PAxF1(var,i,j,k)  (var(i-1,j,k)
#define m3d_PAxF2(var,i,j,k)   -2.0*var(i,j,k)+var(i+1,j,k))

#define m3d_PAyF1(var,i,j,k)  (var(i,j-1,k)
#define m3d_PAyF2(var,i,j,k)   -2.0*var(i,j,k)+var(i,j+1,k))

#define m3d_PAzF1(var,i,j,k)  (var(i,j,k-1)
#define m3d_PAzF2(var,i,j,k)   -2.0*var(i,j,k)+var(i,j,k+1))

#define m4d_PAxF1(var,i,j,k)  (var(m,i-1,j,k)
#define m4d_PAxF2(var,i,j,k)   -2.0*var(m,i,j,k)+var(m,i+1,j,k))

#define m4d_PAyF1(var,i,j,k)  (var(m,i,j-1,k)
#define m4d_PAyF2(var,i,j,k)   -2.0*var(m,i,j,k)+var(m,i,j+1,k))

#define m4d_PAzF1(var,i,j,k)  (var(m,i,j,k-1)
#define m4d_PAzF2(var,i,j,k)   -2.0*var(m,i,j,k)+var(m,i,j,k+1))

