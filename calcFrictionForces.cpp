//note 1: the friction coefficients, stiffness and damping coefficients are user defined at this point, if they are related, calculate that before function call for testing purpose, what about the scenario where stiffness are different for two bodies?
//note 2: dimension of the arrays can be the total number of bodies, total number of contacts, multiples of total number of contacts. Total number of potential contacts will not be considered at this point, since history is assumed to be handled already.
//note 3: assume each thread is handling one contact? or each thread is for a particle? (see note 11)
//note 4: curvature is assumed to be known already, this is for testing function purpose. For ellipsoids, a bit more work needs to go in there.
//note 5: note that if no spinning is needed, contact frame does not need to be tracked, only contact point
//note 6: double precision is needed for tracking contact point and contact frame, due to small changes
//note 7: how to handle normal force? Treat as an input value? If included in the function, this needs to be an output array


//note 8: function call probably is needed, a lot of rotation matrix multiplied with vector, mathUtil only has cross and dot product of vectors, and rotation matrix for vector a to vector b
//note 9: quaternion or orientation for bodies?
//note 10: 3*1 vector, or three doubles, or depending on the problem?
//note 11: what do you want the output to look like? each sphere has a sum of friction force and torque, or each contact out put friction force and torque
#define TOL        0.00000000001
#define TIGHT_TOL  0.0000000000000001
#define PI         3.1415926535897932

void calcFrictionForces(int*   nc,                      // number of current contact  (this to be handled by thread index)
		                    int*    bodyIndices,             // size: nc*2, [nc1_i, nc1_j, nc2_i, nc2_j, ...], read only
                        double* pos,                     // size: nb*3
                        double* quaternion,              // size: nb*4
                        float   dt,                      // time step size
                        float*  normalForceMag,          // size: nc,   [N1_mag, N2_mag, N3_mag, ... ],   read only
                        float*  curvature,               // size: nc*2, [curv1_i, curv1_j, curv2_i, ...], read only
                        float*  fricParameters,          // size: nc*8, [mu_s, mu_k, K_sl, K_r, K_sp, D_sl, D_r, D_sp, ...]
                        double* contactPointPrevLocal,   // size: nc*6, [CP1_i_x, CP1_i_y, CP1_i_z, CP1_j_x, CP1_j_y, CP1_j_z, CP2_i_x, ...],       read and updated
												double* contactPointCurrGlobal,   // size: nc*3, [CP1_x, CP1_y, CP1_z, CP2_x, CP2_y, CP2_z, CP3_x, ...],       read and updated
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
){
  int body_i = bodyIndices[2*nc];
	int body_j = bodyIndices[2*nc+1];
	double3 pos_i = make_float3(pos[3*body_i], pos[3*body_i+1], pos[3*body_i+2]);
	double3 pos_j = make_float3(pos[3*body_j], pos[3*body_j+1], pos[3*body_j+2]);
	// how is quaternion handled in chrono granular?

	float magN = normalForceMag[nc];
	float curv_i = curvature[2*nc];
	float curv_j = curvature[2*nc+1]
	float mu_s = fricParameters[8*nc];
	float mu_k = fricParameters[8*nc+1];
	float K_sl = fricParameters[8*nc+2];
	float K_r  = fricParameters[8*nc+3];
	float K_sp = fricParameters[8*nc+4];
	float D_sl = fricParameters[8*nc+5];
	float D_r  = fricParameters[8*nc+6];
	float D_sp = fricParameters[8*nc+7];

  float3 CP_prev_local_i = make_float3(contactPointPrevLocal[6*nc  ], contactPointPrevLocal[6*nc+1], contactPointPrevLocal[6*nc+2]);
	float3 CP_prev_local_j = make_float3(contactPointPrevLocal[6*nc+3], contactPointPrevLocal[6*nc+4], contactPointPrevLocal[6*nc+5]);
  float3 CP_curr_global  = make_float3(contactPointCurrGlobal[3*nc], contactPointCurrGlobal[3*nc+1], contactPointCurrGlobal[3*nc+2]);

	float3 CF_u_prev_local_i = make_float3(contactFrameUPrevLocal[6*nc  ], contactFrameUPrevLocal[6*nc+1], contactFrameUPrevLocal[6*nc+2]);
	float3 CF_u_prev_local_j = make_float3(contactFrameUPrevLocal[6*nc+3], contactFrameUPrevLocal[6*nc+4], contactFrameUPrevLocal[6*nc+5]);

	float3 CF_w_prev_local_i = make_float3(contactFrameWPrevLocal[6*nc  ], contactFrameWPrevLocal[6*nc+1], contactFrameWPrevLocal[6*nc+2]);
	float3 CF_w_prev_local_j = make_float3(contactFrameWPrevLocal[6*nc+3], contactFrameWPrevLocal[6*nc+4], contactFrameWPrevLocal[6*nc+5]);

	float3 CF_n_prev_local_i = make_float3(contactFrameNPrevLocal[6*nc  ], contactFrameNPrevLocal[6*nc+1], contactFrameNPrevLocal[6*nc+2]);
	float3 CF_n_prev_local_j = make_float3(contactFrameNPrevLocal[6*nc+3], contactFrameNPrevLocal[6*nc+4], contactFrameNPrevLocal[6*nc+5]);

	float3 CF_n_curr_global  = make_float3(contactFrameNCurrGlobal[3*nc], contactFrameNCurrGlobal[3*nc+1], contactFrameNCurrGlobal[3*nc+2]);


  // express previous contact frame in current global coordinate (what to do with rotation matrix Ai and Aj?)
	float3 CF_u_prev_global_curr_i = Ai*CF_u_prev_local_i;
	float3 CF_u_prev_global_curr_j = Aj*CF_u_prev_local_j;
	float3 CF_w_prev_global_curr_i = Ai*CF_w_prev_local_i;
	float3 CF_w_prev_global_curr_j = Aj*CF_w_prev_local_j;
	float3 CF_n_prev_global_curr_i = Ai*CF_n_prev_local_i;
	float3 CF_n_prev_global_curr_j = Aj*CF_n_prev_local_j;

	// express previous contact point in current global coordinate
	float3 CP_prev_global_curr_i = pos_i + Ai * CP_prev_local_i;
	float3 CP_prev_global_curr_j = pos_j + Aj * CP_prev_local_j;

  // use optimization method to find current global contact frame


	contactFrame contactFrame::getContactFrameSmallestRotation(){
		double a1, a2, b1, b2, theta_star, costFunc;
		a1 = (this->u).x; a2 = (this->u).y;
		b1 = (this->w).x; b2 = (this->w).y;

		if (fabs(a1+b2) < TIGHT_TOL){
			theta_star = PI/2.0;
		}
		else{
			theta_star = atan((a2-b1)/(a1+b2));
		}

		//***************************************************//
		//figure out a way to optimize this, seems expenseive//
		//***************************************************//
		costFunc = (a1+b2)*cos(theta_star) + (a2-b1)*sin(theta_star);

		if (costFunc < 0)
			theta_star = theta_star + PI;

		contactFrame newFrame;
		newFrame.u = vec3( cos(theta_star), sin(theta_star), 0);
		newFrame.w = vec3(-sin(theta_star), cos(theta_star), 0);
		newFrame.n = vec3( 0, 0, 1);
		return newFrame;
	}

	contactFrame contactFrame::getContactFrameSmallestRotation(vec3 n_curr){
		contactFrame rotationA;
		vec3 global_normal = vec3(0,0,1);
		// rotation matrix A*n_curr = [0;0;1];
		rotationA.getRotationMatrixFromAtoB(n_curr, global_normal);
		contactFrame rotatedFrame_prev;
		rotatedFrame_prev.u = rotationA.multiplyVec(this->u);
		rotatedFrame_prev.w = rotationA.multiplyVec(this->w);

		contactFrame CF_curr, CF_rotated;
		CF_curr = rotatedFrame_prev.getContactFrameSmallestRotation();
		CF_rotated.u = rotationA.transposeMultiplyVec(CF_curr.u);
		CF_rotated.w = rotationA.transposeMultiplyVec(CF_curr.w);
		CF_rotated.n = rotationA.transposeMultiplyVec(CF_curr.n);

		return CF_rotated;
	}







	// find increments for slide/roll/spin




	// enter force-displacement relation


	void spherePlaneContactKinematics::updateContactKinematics(){
		// express previous contact frame globally
		CF_prev_global_curr = CF_prev_local.expressContactFrameInGlobalRF(sphere_ptr);
		// express previous contact point globally
		CP_prev_global_curr = sphere_ptr->expressLocalPtInGlobalRF(CP_prev_local);
		// get current contact point
		CP_curr_global = plane_ptr->getProjectedPoint(sphere_ptr->getPos());
		// use optimization method to find current contact frame with minimum rotation
		CF_curr_global = CF_prev_global_curr.getContactFrameSmallestRotation(plane_ptr->normal);
		// find current contact point in local reference void spherePlaneContactKinematics::frame
		CP_curr_local = sphere_ptr->expressGlobalPtInLocalRF(CP_curr_global);
		// find current contac
		CF_curr_local = CF_curr_global.expressContactFrameInLocalRF(sphere_ptr);
	}

	// replace Previous Contact Frame


		this->CF_prev_global = this->CF_curr_global;
		this->CF_prev_local  = this->CF_curr_local;

	void spherePlaneContactKinematics::replacePreviousContactPoint(){
		this->CP_prev_global = this->CP_curr_global;
		this->CP_prev_local  = this->CP_curr_local;
	}

	vec3 spherePlaneContactKinematics::sphereSurfaceGeodesicProjection(){
		vec3 s_proj, chord, Ra, Rb, proj_chord_n;
		double ratio, angle, arc_length;
		chord = CP_curr_global - CP_prev_global_curr;
		Ra = CP_prev_global_curr - sphere_ptr->getPos();
		Rb = CP_curr_global - sphere_ptr->getPos();
		// use |Ra|*|Rb| instead of R^2
		// because when deformation is allowed in normal
		// direction, |Ra|!=R
		ratio = Ra.innerProduct(Rb)/(Ra.getNorm() * Rb.getNorm());
		assert(ratio <= 1.0);
		angle = acos(ratio);
		arc_length = Ra.getNorm() * angle;
		assert(fabs(CF_curr_global.n.getNorm() - 1) < TOL);
		proj_chord_n = chord - chord.innerProduct(CF_curr_global.n)*CF_curr_global.n;

		if (fabs(proj_chord_n.getNorm()) < TIGHT_TOL){
			s_proj = vec3(0,0,0);
		}
		else{
			s_proj = proj_chord_n/proj_chord_n.getNorm() * arc_length;
		}
		return s_proj;
	}

	void spherePlaneContactKinematics::evaluateIncrementFromKinematics(){
		pi = this->sphereSurfaceGeodesicProjection();
		pj_bar = CP_curr_global - CP_prev_global;
		// find relative rotation of the tangetial plane (psi)
		// body i CF at current time (u1, w1, n1)
		// body j CF at current time (u1_bar, w1_bar, n1)
		// share same n1
		// when j is the ground [u1_bar, w1_bar] = body i global CF at previous time
		psi = CF_curr_global.rotationFrom(&CF_prev_global);
		pj = pj_bar.rotationAboutAxis(CF_curr_global.n, -psi);
		delta  = pj - pi;
		excursion = pi/sphere_ptr->getRadius();
	}














}
