//note 1: the friction coefficients, stiffness and damping coefficients are user defined at this point, if they are related, calculate that before function call for testing purpose, what about the scenario where stiffness are different for two bodies?
//note 2: dimension of the arrays can be the total number of bodies, total number of contacts, multiples of total number of contacts. Total number of potential contacts will not be considered at this point, since history is assumed to be handled already.
//note 3: assume each thread is handling one contact
//note 4: curvature is assumed to be known already, this is for testing function purpose. For ellipsoids, a bit more work needs to go in there.
//note 5: note that if no spinning is needed, contact frame does not need to be tracked, only contact point
//note 6: double precision is needed for tracking contact point and contact frame, due to small changes
//note 7: how to handle normal force? Treat as an input value? If included in the function, this needs to be an output array

__device__ void calcFrictionForces(int*   nc,                      // number of current contact  (this to be handled by thread index)
		                          int*    bodyIndices,             // size: nc*2, [nc1_i, nc1_j, nc2_i, nc2_j, ...], read only
                                  double* pos,                     // size: nb*3
                                  double* quaternion,              // size: nb*4
                                  float   dt,                      // time step size
                                  float*  normalForceMag,          // size: nc,   [N1_mag, N2_mag, N3_mag, ... ],   read only
                                  float*  curvature,               // size: nc*2, [curv1_i, curv1_j, curv2_i, ...], read only
                                  float*  fricParameters,          // size: nc*8, [mu_s, mu_k, K_sl, K_r, K_sp, D_sl, D_r, D_sp, ...]
                                  double* contactPointPrevLocal,   // size: nc*6, [CP1_i_x, CP1_i_y, CP1_i_z, CP1_j_x, CP1_j_y, CP1_j_z, CP2_i_x, ...],       read and updated
								  double* contactFrameUPrevLocal,  // size: nc*6, [CF1_i_Ux, CF1_i_Uy, CF1_i_Uz, CF1_j_Ux, CF1_j_Uy, CF1_j_Uz, CF2_i_Ux,...], read and updated
								  double* contactFrameWPrevLocal,  // size: nc*6, [CF1_i_Wx, CF1_i_Wy, CF1_i_Wz, CF1_j_Wx, CF1_j_Wy, CF1_j_Wz, CF2_i_Wx,...], read and updated
								  double* contactFrameNPrevLocal,  // size: nc*6, [CF1_i_Nx, CF1_i_Ny, CF1_i_Nz, CF1_j_Nx, CF1_j_Ny, CF1_j_Nz, CF2_i_Nx,...], read and updated
                                  double* contactFrameNCurrGlobal, // size: nc*3, [CF1_Nx, CF1_Ny, CF1_Nz, CF2_Nx, CF2_Ny, CF2_Nz, CF3_Nx,...],               read only
                                  double* slidingHistory,          // size: nc*3, [S1_x, S1_y, S1_z, S2_x, ...],                                              read and updated
								  double* rollingHistory,          // size: nc*6, [Theta1_i_x, Theta1_i_y, Theta1_i_z, Theta1_j_x, Theta1_j_y, Theta1_j_z, Theta2_i_x, ...], read and updated
								  double* spinningHistory,         // size: nc*3, [Psi1_x, Psi1_y, Psi1_z, Psi2_x, ...],                                      read and updated
                                  char*   slidingMode,             // size: nc,   [sl_mode_1, sl_mode_2, ...],                                                read and updated
                                  char*   rollingMode,             // size: nc*2, [r_mode_1i, r_mode_1j, r_mode_2i, ...],                                     read and updated
                                  char*   spinningMode,            // size: nc,   [sp_mode_1, sp_mode_2, sp_mode_3, ...],                                     read and updated
								  float*  out_slidingFr,           // size: nc*3, [Fr1_x, Fr1_y, Fr1_z, Fr2_x, ...],                                          updated only
                                  float*  out_rollingTr,           // size: nc*6, [Tr1_i_x, Tr1_i_y, Tr1_i_z, Tr1_j_x, Tr1_j_y, Tr1_j_z, Tr2_i_x, ...],       updated only
                                  float*  out_spinningTr           // size: nc*3, [Ts1_x, Ts1_y, Ts1_z, Ts2_x, ...],                                          udpated only
)
