GeodesicFBORenderer = makeClass({
	//1024x1024 is just too big for my current card to handle
	lightTexWidth : 512,
	lightTexHeight : 512,
	lightPosVelChannels : [
		{pos_or_vel:'pos'},
		{pos_or_vel:'vel'}
	],
	
	init : function(glutil) {
		this.glutil = glutil;
		if (!this.glutil.context.getExtension('OES_texture_float')) {
			throw 'This requires OES_texture_float';
		}
	},
	
	initScene : function(skyTex) {
		var thiz = this;
		$.each(thiz.lightPosVelChannels, function(_,channel) {
			//working on how to organize this
			channel.uniforms = [
				'blackHoleMass', 
				'blackHoleCharge', 
				'blackHoleAngularVelocity', 
				'warpBubbleThickness', 
				'warpBubbleRadius', 
				'warpBubbleVelocity', 
				'objectDist', 
				'objectAngle', 
				'deltaLambda'
			];
			
			channel.texs = [];
			channel.fbos = [];
			for (var history = 0; history < 2; ++history) {
				var texs = [];
				channel.texs[history] = texs;
				var fbos = [];
				channel.fbos[history] = fbos;
				for (var side = 0; side < 6; ++side) {
					var tex = new thiz.glutil.Texture2D({
						internalFormat : gl.RGBA,
						format : gl.RGBA,
						type : gl.FLOAT,
						width : thiz.lightTexWidth,
						height : thiz.lightTexHeight,
						magFilter : gl.NEAREST,
						minFilter : gl.NEAREST,
						wrap : {
							s : gl.CLAMP_TO_EDGE,
							t : gl.CLAMP_TO_EDGE
						}
					});
					texs[side] = tex;
					
					var fbo = new thiz.glutil.Framebuffer();
					gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.obj);
					gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex.obj, 0);
					gl.bindFramebuffer(gl.FRAMEBUFFER, null);
					fbos[side] = fbo;
				}
			}
			
			channel.shaders = {};

			var metricCodes = {
				['Schwarzschild Black Hole'] : mlstr(function(){/*
uniform float blackHoleMass;

// g_tt = -(1 - m/(2r))^2 / (1 + m/(2r))^2
float g_tt(vec3 pos) {
	float r = length(pos);
	float _m_2r = blackHoleMass / (2. * r);
	float _1_plus_m_2r = 1. + _m_2r;
	float _1_minus_m_2r = 1. - _m_2r;
	float ratio = _1_minus_m_2r / _1_plus_m_2r;
	return -ratio * ratio; 
}

// g_ti = 0
vec3 g_ti(vec3 pos) { 
	return vec3(0., 0., 0.); 
}

// g_ij = delta_ij (1 + m/(2r))^4 
mat3 g_ij(vec3 pos) {
	float r = length(pos);
	float _m_2r = blackHoleMass / (2. * r);
	float _1_plus_m_2r = 1. + _m_2r;
	float _1_plus_m_2r_sq = _1_plus_m_2r * _1_plus_m_2r;
	float _1_plus_m_2r_toTheFourth = _1_plus_m_2r_sq * _1_plus_m_2r_sq;
	return mat3(_1_plus_m_2r_toTheFourth);
}

vec4 accel(vec4 pos, vec4 vel) {
	float r = length(pos.xyz);
	float posDotVel = dot(pos.xyz, vel.xyz);
	float posDotVelSq = posDotVel * posDotVel;
	float velSq = dot(vel.xyz, vel.xyz);
	float R = 2. * blackHoleMass;
	float r2 = r * r;
	float r3 = r * r2;


	float _m_2r = blackHoleMass / (2. * r);
	float _1_plus_m_2r = 1. + _m_2r;
	float _1_minus_m_2r = 1. - _m_2r;
	
	float _1_plus_m_2r_sq = _1_plus_m_2r * _1_plus_m_2r;
	float _1_plus_m_2r_tothe3 = _1_plus_m_2r_sq * _1_plus_m_2r;
	float _1_plus_m_2r_tothe6 = _1_plus_m_2r_tothe3 * _1_plus_m_2r_tothe3;
	
	vec4 result;
	result.w = -blackHoleMass * 2. * vel.w * posDotVel / (_1_plus_m_2r * _1_minus_m_2r * r3);
	result.xyz = -blackHoleMass * (
		vel.w * vel.w * pos.xyz * _1_minus_m_2r / _1_plus_m_2r_tothe6 
		- 2. * vel.xyz * posDotVel
		+ velSq * pos.xyz 
	) / (_1_plus_m_2r * r3);
	return result;
}

*/}),
				
				// substitution of a=0, Q=0, simlpified
				// This somewhat match the Schwarzschild isotropic geodesic above
				// The full Kerr geodesic below isn't working.
				['Kerr Black Hole degeneracy'] : mlstr(function(){/*
uniform float blackHoleMass;
uniform float blackHoleCharge;
uniform float blackHoleAngularMomentum;

float g_tt(vec3 pos) {
	return -1. + 2. * blackHoleMass / length(pos.xyz);
}

vec3 g_ti(vec3 pos) { 
	float R = 2. * blackHoleMass;
	return R * pos.xyz / dot(pos.xyz, pos.xyz);
}

mat3 g_ij(vec3 pos) {
	float R = 2. * blackHoleMass;
	float inv_r = 1. / length(pos.xyz);
	float R_r = R * inv_r;
	vec3 npos = pos.xyz * inv_r;	
	return mat3(1.) + R_r * outerProduct(npos, npos);
}

vec4 accel(vec4 pos, vec4 vel) {
	float R = 2. * blackHoleMass;
	
	float inv_r = 1. / length(pos.xyz);
	vec4 npos = pos * inv_r;	
	vec4 nvel = vel * inv_r;

	float R_r = R * inv_r;
	float npos_nvel = dot(npos.xyz, nvel.xyz);
	float nvel_nvel = dot(nvel.xyz, nvel.xyz);
	
	vec4 neg_accel;
	neg_accel.w = R * (
		.5 * R_r * nvel.w * nvel.w
		+ npos_nvel * (1. + R_r) * nvel.w
		+ npos_nvel * npos_nvel * (2. + .5 * R_r) 
		- nvel_nvel
	);
	neg_accel.xyz = (R * (
		.5 * (1. - R_r) * nvel.w * nvel.w
		- R_r * npos_nvel * nvel.w
		+ nvel_nvel
		- .5 * inv_r * npos_nvel * npos_nvel * (3. - R_r)
	)) * npos.xyz;
	return -neg_accel;
}
*/}),


/*
1 = (x^2 + y^2) / (r^2 + a^2) + z^2 / r^2
r^4 + r^2 (a^2 - x^2 - y^2 - z^2) - a^2 z^2 = 0

b = x^2 + y^2 + z^2 - a^2
r^2 = 1/2 ( -b +- sqrt( b*b + 4 * a^2 * z^2))

4 r^3 dr + 2 r dr (a^2 - x^2 - y^2 - z^2) + r^2 (-2 x dx - 2 y dy - 2 z dz) - 2 a^2 z dz = 0
(4 r^2 + 2 (a^2 - x^2 - y^2 - z^2)) r dr - 2 r^2 x dx - 2 r^2 y dy - 2 (r^2 + a^2) z dz = 0

dr/dx = r x / (2 r^2 + a^2 - x^2 - y^2 - z^2)
dr/dy = r y / (2 r^2 + a^2 - x^2 - y^2 - z^2)
dr/dz = (r + a^2/r) z / (2 r^2 + a^2 - x^2 - y^2 - z^2)
r_,i = (r + delta_iz a^2 / r) / (2 r^2 + a^2 - x^2 - y^2 - z^2) x^i
*/


				['Kerr Black Hole'] : mlstr(function(){/*
uniform float blackHoleMass;
uniform float blackHoleCharge;
uniform float blackHoleAngularMomentum;

float calc_Kerr_rSq(vec3 pos) {
	float a = blackHoleAngularMomentum / blackHoleMass;
	float z = pos.z;
	float b = a*a - dot(pos,pos);
	float rSq = .5*(-b + sqrt(b*b + 4.*a*a*z*z));
	return rSq;
}

float calc_Kerr_H(vec3 pos, float r) {
	float rs = 2. * blackHoleMass;
	float Q = blackHoleCharge;
	float a = blackHoleAngularMomentum / blackHoleMass;
	float z = pos.z;
	float H = .5 * (r * rs - Q * Q) / (r*r + a*a*z*z/(r*r));
	return H;
}

vec3 calc_Kerr_l(vec3 pos, float r) {
	float a = blackHoleAngularMomentum / blackHoleMass;
	float x = pos.x;
	float y = pos.y;
	float z = pos.z;
	vec3 l = vec3(
		(r*x + a*y)/(r*r + a*a),
		(r*y - a*x)/(r*r + a*a),
		z/r);
	return l;
}

float g_tt(vec3 pos) {
	float r = sqrt(calc_Kerr_rSq(pos));
	float H = calc_Kerr_H(pos, r);
	return -1. + 2. * H;
}

vec3 g_ti(vec3 pos) { 
	float r = sqrt(calc_Kerr_rSq(pos));
	float H = calc_Kerr_H(pos, r);
	vec3 l = calc_Kerr_l(pos, r);
	return 2. * H * l;
}

mat3 g_ij(vec3 pos) {
	float r = sqrt(calc_Kerr_rSq(pos));
	float H = calc_Kerr_H(pos, r);
	vec3 l = calc_Kerr_l(pos, r);
	return mat3(1.) + 2. * H * outerProduct(l,l);
}

//
//l_t = 1
//l_x = (r x + a y) / (r^2 + a^2)
//l_y = (r y - a x) / (r^2 + a^2)
//l_z = z / r
//
//H = .5 (r rs - Q^2) / (r^2 + a^2 z^2/r^2)
//H,i = 
//	(
//		.5 r_,i rs (r^2 + a^2 z^2 / r^2) 
//		- (r rs - Q^2) (r r_,i + a^2 (delta_zi r - z r_,i) z / r^3)
//	) / (r^2 + a^2 z^2 / r^2)^2
//
//g_uv = eta_uv + 2 H l_u l_v
//g_tt = -1 + 2 H
//g_ti = 2 H l_i
//g_ij = delta_ij + 2 H l_i l_j
//
//g^uv = eta^uv - 2 H (eta^ua l_a) (eta^vb l_b)
//g^tt = -1 - 2 H
//g^ti = 2 H l_i
//g^ij = delta^ij - 2 H l_i l_j
//
//g_uv,w = 2 (H,w l_u l_v + H (l_u,w l_v + l_u l_v,w))

vec4 accel(vec4 pos, vec4 vel) {
#if 1
	float t = pos.w;
	float x = pos.x;
	float y = pos.y;
	float z = pos.z;

	float zSq = z * z;
	float pos_pos = dot(pos.xyz, pos.xyz);

	float rs = 2. * blackHoleMass;
	float Q = blackHoleCharge;
	float QSq = Q * Q;
	float a = blackHoleAngularMomentum / blackHoleMass;
	float aSq = a*a;
	float b = aSq - pos_pos;
	float sqrtdiscr = sqrt(b*b + 4.*aSq*zSq);
	float rSq = .5*(-b + sqrtdiscr);
	float r = sqrt(rSq);

	vec3 zHat = vec3(0., 0., 1.);
	vec3 dr_dx = (pos.xyz * r + (pos.xyz + zHat) * aSq / r) / (2. * rSq + aSq - pos_pos);
	
	float H_denom = rSq + aSq * zSq / rSq;
	float H = .5 * (r * rs - QSq) / H_denom;

	vec4 dH_dx;
	dH_dx.xyz = (
			.5 * dr_dx * rs * H_denom
			- (r * rs - QSq) * (
				r * dr_dx + aSq * z / (r * rSq) * (zHat * r - z * dr_dx)
			)
		) / (H_denom * H_denom);

	float aSq_plus_rSq = rSq + aSq;
	
	vec4 l = vec4(
		(r*x + a*y)/aSq_plus_rSq,
		(r*y - a*x)/aSq_plus_rSq,
		z/r,
		1.
	);

	vec4 lU = l;
	lU.w = -lU.w;

	mat4 gInv = mat4(
		vec4(1., 0., 0., 0.),
		vec4(0., 1., 0., 0.),
		vec4(0., 0., 1., 0.),
		vec4(0., 0., 0., -1.)
	) - 2. * H * outerProduct(lU, lU);


	// l_a,b == dl_dx[a][b]
	mat4 dl_dx = mat4(
		vec4(
			(dr_dx[0] * (x * (aSq - rSq) - 2. * a * y * r) / aSq_plus_rSq + r) / aSq_plus_rSq,
			(dr_dx[1]  * (x * (aSq - rSq) - 2. * a * y * r) / aSq_plus_rSq + a) / aSq_plus_rSq,
			-dr_dx[2] * (x * (aSq - rSq) - 2. * a * y * r) / (aSq_plus_rSq * aSq_plus_rSq),
			0.),
		vec4(
			(dr_dx[0] * (y * (aSq - rSq) + 2. * a * x * r) / aSq_plus_rSq - a) / aSq_plus_rSq,
			(dr_dx[1] * (y * (aSq - rSq) + 2. * a * x * r) / aSq_plus_rSq + r) / aSq_plus_rSq,
			dr_dx[2] * (y * (aSq - rSq) + 2. * a * x * r) / (aSq_plus_rSq * aSq_plus_rSq),
			0.),
		vec4(
			-z * dr_dx[0] / rSq,
			-z * dr_dx[1] / rSq,
			1./r - z * dr_dx[2] / rSq,
			0.),
		vec4(0., 0., 0., 0.));
	
	// conn_abc = H,a l_b l_c + H,b l_a l_c + H,c l_a l_b
	//		+ H (l_a (l_b,c + l_c,b) + l_b (l_a,c + l_c,a) + l_c (l_a,b + l_b,a))
	vec4 conn_vel_vel;
	for (int a = 0; a < 4; ++a) {
		float sum = 0.;
		for (int b = 0; b < 4; ++b) {
			for (int c = 0; c < 4; ++c) {
				float conn_lll = dH_dx[a] * l[b] * l[c]
					+ dH_dx[b] * l[a] * l[c]
					+ dH_dx[c] * l[a] * l[b]
					+ H * (
						l[a] * (dl_dx[b][c] + dl_dx[c][b])
						+ l[b] * (dl_dx[a][c] + dl_dx[c][a])
						+ l[c] * (dl_dx[a][b] + dl_dx[b][a])
					);
				sum += conn_lll * vel[b] * vel[c];
			}
		}
		conn_vel_vel[a] = sum;
	}
	
	// conn^a_bc = g^ad conn_dbc
	return -gInv * conn_vel_vel;
#endif
	
#if 0	//degenerate case:
	float R = 2. * blackHoleMass;
	
	float inv_r = 1. / length(pos.xyz);
	vec4 npos = pos * inv_r;	
	vec4 nvel = vel * inv_r;

	float R_r = R * inv_r;
	float npos_nvel = dot(npos.xyz, nvel.xyz);
	float nvel_nvel = dot(nvel.xyz, nvel.xyz);
	
	vec4 neg_accel;
	neg_accel.w = R * (
		.5 * R_r * nvel.w * nvel.w
		+ npos_nvel * (1. + R_r) * nvel.w
		+ npos_nvel * npos_nvel * (2. + .5 * R_r) 
		- nvel_nvel
	);
	neg_accel.xyz = (R * (
		.5 * (1. - R_r) * nvel.w * nvel.w
		- R_r * npos_nvel * nvel.w
		+ nvel_nvel
		- .5 * inv_r * npos_nvel * npos_nvel * (3. - R_r)
	)) * npos.xyz;
	return -neg_accel;
#endif
}
*/}),
		
				
				['Kerr Black Hole full'] : mlstr(function(){/*
uniform float blackHoleMass;
uniform float blackHoleCharge;
uniform float blackHoleAngularMomentum;

float g_tt(vec3 pos) {
	float rs = 2. * blackHoleMass;
	float Q = blackHoleCharge;
	float a = blackHoleAngularMomentum / blackHoleMass;
	float z = pos.z;
	float b = a*a - dot(pos,pos);
	float rSq = .5*(-b + sqrt(b*b + 4.*a*a*z*z));
	float r = sqrt(rSq);
	float H = .5 * (r * rs - Q * Q) / (rSq + a*a*z*z/rSq);
	return -1. + 2. * H;
}

vec3 g_ti(vec3 pos) { 
	float rs = 2. * blackHoleMass;
	float Q = blackHoleCharge;
	float a = blackHoleAngularMomentum / blackHoleMass;
	float x = pos.x;
	float y = pos.y;
	float z = pos.z;
	float b = a*a - dot(pos,pos);
	float rSq = .5*(-b + sqrt(b*b + 4.*a*a*z*z));
	float r = sqrt(rSq);
	float H = .5 * (r * rs - Q * Q) / (rSq + a*a*z*z/rSq);
	vec3 l = vec3(
		(r*x + a*y)/(rSq + a*a),
		(r*y - a*x)/(rSq + a*a),
		z/r);
	return 2. * H * l;
}

mat3 g_ij(vec3 pos) {
	float rs = 2. * blackHoleMass;
	float Q = blackHoleCharge;
	float a = blackHoleAngularMomentum / blackHoleMass;
	float x = pos.x;
	float y = pos.y;
	float z = pos.z;
	float b = a*a - dot(pos,pos);
	float rSq = .5*(-b + sqrt(b*b + 4.*a*a*z*z));
	float r = sqrt(rSq);
	float H = .5 * (r * rs - Q * Q) / (rSq + a*a*z*z/rSq);
	vec3 l = vec3(
		(r*x + a*y)/(rSq + a*a),
		(r*y - a*x)/(rSq + a*a),
		z/r);
	return mat3(1.) + 2. * H * outerProduct(l,l);
}

vec4 accel(vec4 pos, vec4 vel) {
	float t = pos.w;
	float x = pos.x;
	float y = pos.y;
	float z = pos.z;

	float rs = 2. * blackHoleMass;
	float Q = blackHoleCharge;
	float a = blackHoleAngularMomentum / blackHoleMass;
	float aSq = a*a;
	float b = aSq - dot(pos.xyz, pos.xyz);
	float sqrtdiscr = sqrt(b*b + 4.*aSq*z*z);
	float rSq = .5*(-b + sqrtdiscr);
	float r = sqrt(rSq);
	float H = .5 * (r * rs - Q * Q) / (rSq + aSq*z*z/rSq);

	vec3 db_dxis = -2. * pos.xyz;
	vec3 dr_dxis;
	for (int i = 0; i < 3; ++i) {
		dr_dxis[i] = (
			-.25 * db_dxis[i] 
			+ .5 * (
				b * db_dxis[i] 
				+ (i == 2 ? (4.*aSq*z) : 0.)
			) / sqrtdiscr
		) / r;
	}
	vec4 dH_dx;	//stored xyzt
	dH_dx.w = 0.;
	for (int i = 0; i < 3; ++i) {
		dH_dx[i] = .5 * (
			dr_dxis[i] * rs * (rSq + aSq*z*z / rSq)
			- (r * rs - Q*Q) * (
				2. * r * dr_dxis[i]
				+ aSq*(
					(i == 2 ? (2. * z / rSq) : 0.)
					- 2. * z*z*r*dr_dxis[i]
				)
			)
		);
	}

	float aSq_plus_rSq = rSq + aSq;
	
	vec4 l = vec4(
		(r*x + a*y)/aSq_plus_rSq,
		(r*y - a*x)/aSq_plus_rSq,
		z/r,
		1.
	);

	vec4 lU = l;
	lU.w = -lU.w;

	mat4 gInv = mat4(
		vec4(1., 0., 0., 0.),
		vec4(0., 1., 0., 0.),
		vec4(0., 0., 1., 0.),
		vec4(0., 0., 0., -1.)
	) - 2. * H * outerProduct(lU, lU);


	// l_a,b == dl_dx[a][b]
	mat4 dl_dx = mat4(
		vec4(
			(dr_dxis[0] * (x * (aSq - rSq) - 2. * a * y * r) / aSq_plus_rSq + r) / aSq_plus_rSq,
			(dr_dxis[1]  * (x * (aSq - rSq) - 2. * a * y * r) / aSq_plus_rSq + a) / aSq_plus_rSq,
			-dr_dxis[2] * (x * (aSq - rSq) - 2. * a * y * r) / (aSq_plus_rSq * aSq_plus_rSq),
			0.),
		vec4(
			(dr_dxis[0] * (y * (aSq - rSq) + 2. * a * x * r) / aSq_plus_rSq - a) / aSq_plus_rSq,
			(dr_dxis[1] * (y * (aSq - rSq) + 2. * a * x * r) / aSq_plus_rSq + r) / aSq_plus_rSq,
			dr_dxis[2] * (y * (aSq - rSq) + 2. * a * x * r) / (aSq_plus_rSq * aSq_plus_rSq),
			0.),
		vec4(
			-z * dr_dxis[0] / rSq,
			-z * dr_dxis[1] / rSq,
			1./r - z * dr_dxis[2] / rSq,
			0.),
		vec4(0., 0., 0., 0.));
	
	// conn_abc = H,a l_b l_c + H,b l_a l_c + H,c l_a l_b
	//		+ H (l_a (l_b,c + l_c,b) + l_b (l_a,c + l_c,a) + l_c (l_a,b + l_b,a))
	vec4 conn_vel_vel;
	for (int a = 0; a < 4; ++a) {
		float sum = 0.;
		for (int b = 0; b < 4; ++b) {
			for (int c = 0; c < 4; ++c) {
				float conn_lll = dH_dx[a] * l[b] * l[c]
					+ dH_dx[b] * l[a] * l[c]
					+ dH_dx[c] * l[a] * l[b]
					+ H * (
						l[a] * (dl_dx[b][c] + dl_dx[c][b])
						+ l[b] * (dl_dx[a][c] + dl_dx[c][a])
						+ l[c] * (dl_dx[a][b] + dl_dx[b][a])
					);
				sum += conn_lll * vel[b] * vel[c];
			}
		}
		conn_vel_vel[a] = sum;
	}
	
	// conn^a_bc = g^ad conn_dbc
	return -gInv * conn_vel_vel;
}
*/}),
				
				['Alcubierre Warp Drive Bubble'] : mlstr(function(){/*
uniform float warpBubbleThickness;
uniform float warpBubbleVelocity;
uniform float warpBubbleRadius;

float Alcubierre_f(vec3 pos) {
	float rs = length(pos);
	float sigmaFront = warpBubbleThickness * (rs + warpBubbleRadius);
	float sigmaCenter = warpBubbleThickness * warpBubbleRadius;
	float sigmaBack = warpBubbleThickness * (rs - warpBubbleRadius);
	float tanhSigmaCenter = tanh(sigmaCenter);
	float f = (tanh(sigmaFront) - tanh(sigmaBack)) / (2. * tanhSigmaCenter);
	return f;
}

float g_tt(vec3 pos) {
	float f = Alcubierre_f(pos);
	float vf = f * warpBubbleVelocity;
	float vf2 = vf * vf;
	return -(1. - vf2); 
}

vec3 g_ti(vec3 pos) {
	float f = Alcubierre_f(pos);
	float vf = f * warpBubbleVelocity;
	float g_tx = -vf;
	return vec3(g_tx, 0., 0.);
}

mat3 g_ij(vec3 pos) {
	return mat3(1.);
}

vec4 accel(vec4 pos, vec4 vel) {
	float r = length(pos.xyz);
	float sigmaFront = warpBubbleThickness * (r + warpBubbleRadius);
	float sigmaCenter = warpBubbleThickness * warpBubbleRadius;
	float sigmaBack = warpBubbleThickness * (r - warpBubbleRadius);
	float tanhSigmaCenter = tanh(sigmaCenter);
	float f = (tanh(sigmaFront) - tanh(sigmaBack)) / (2. * tanhSigmaCenter);
	float sechDiff = sechSq(sigmaFront) - sechSq(sigmaBack);
	float dfScalar = sechDiff / (2. * r * tanhSigmaCenter);
	vec4 df;
	df.xyz = warpBubbleThickness * pos.xyz * dfScalar;
	df.w = -warpBubbleVelocity * warpBubbleThickness * pos.x * dfScalar;

	float u = f * warpBubbleVelocity;
	vec4 du = df * warpBubbleVelocity;
	
	vec4 result;
	result.w = 
		- du.x * vel.x * vel.x
		- du.y * vel.x * vel.y
		- du.z * vel.x * vel.z
		
		+ 2. * du.x * vel.w * vel.x * u
		- du.x * vel.w * vel.w * u * u
		+ du.y * vel.w * vel.y * u
		+ du.z * vel.w * vel.z * u
	;
	result.x = 
		du.w * vel.w * vel.w
		+ du.y * vel.w * vel.y
		+ du.y * vel.w * vel.y * u * u
		- du.y * vel.x * vel.y * u
		+ du.z * vel.w * vel.z
		+ du.z * vel.w * vel.z * u * u
		- du.z * vel.x * vel.z * u
		+ du.x * vel.w * vel.w * u
		- du.x * vel.w * vel.w * u * u * u
		+ 2. * du.x * vel.w * vel.x * u * u
		- du.x * vel.x * vel.x * u
	;
	result.y = -du.y * vel.w * (vel.x - vel.w * u);
	result.z = -du.z * vel.w * (vel.x - vel.w * u);
	return result;
}

*/})		
			};
					
			$.each(objectTypes, function(_,objType) {
				var pos_or_vel = channel.pos_or_vel;
				
				channel.shaders[objType] = {};

				var shaderTypes = ['reset', 'iterate'];
				
				$.each(shaderTypes, function(_,shaderType) {
					
					var vsh = [];
					var fsh = [shaderCommonCode];

					fsh.push(mlstr(function(){/*
uniform float objectDist;
uniform vec4 objectAngle;
*/}));

					if (shaderType == 'reset') {
						
						//new
						vsh.push(mlstr(function(){/*
attribute vec2 vertex;
varying vec3 vtxv;
uniform mat3 rotation;
void main() {
	vtxv = rotation * vec3(vertex * 2. - 1., 1.);
	gl_Position = vec4(vertex * 2. - 1., 0., 1.);
}
*/}));
						if (pos_or_vel == 'pos') {
							fsh.push(mlstr(function(){/*
varying vec3 vtxv;
void main() {
	gl_FragColor = vec4(-objectDist, 0., 0., 0.);
	gl_FragColor.xyz = quatRotate(objectAngle, gl_FragColor.xyz);
}
*/}));
						} else if (pos_or_vel == 'vel') {
							fsh.push(metricCodes[objType]);
							fsh.push(mlstr(function(){/*

//when initializing our metric:
//g_ab v^a v^b = 0 for our metric g
//g_tt (v^t)^2 + 2 g_ti v^t v^i + g_ij v^i v^j = 0
//(v^t)^2 + (2 v^i g_ti / g_tt) v^t + g_ij v^i v^j / g_tt = 0
//v_t = v^i g_ti / (-g_tt) +- sqrt( (v^i g_ti / g_tt)^2 + g_ij v^i v^j / (-g_tt) )
// I'll go with the plus
float reset_w(vec3 pos, vec3 vel) {
	float _g_tt = g_tt(pos);
	float negInv_g_tt = -1. / _g_tt;
	vec3 _g_ti = g_ti(pos);
	float v_t = dot(vel, _g_ti);	//v_t = g_ti v^i
	float a = v_t * negInv_g_tt;
	float bsq = a * a + dot(vel, g_ij(pos) * vel) * negInv_g_tt;
	if (bsq < 0.) bsq = 0.;
	return a + sqrt(bsq);
}

varying vec3 vtxv;
void main() {
	vec3 pos = vec3(-objectDist, 0., 0.);
	pos = quatRotate(objectAngle, pos);
	
	vec3 vel = normalize(vtxv);
	vel = quatRotate(objectAngle, vel);
	
	gl_FragColor.xyz = vel;
	gl_FragColor.w = reset_w(pos, vel);
}
*/}));
						}
					} else if (shaderType == 'iterate') {
						vsh.push(mlstr(function(){/*
varying vec2 uv;
varying vec3 vtxv;
attribute vec2 vertex;
uniform mat3 rotation;
void main() {
	uv = vertex;
	vtxv = rotation * vec3(vertex * 2. - 1., 1.);
	gl_Position = vec4(vertex * 2. - 1., 0., 1.);
}
*/}));
						fsh.push(mlstr(function(){/*
varying vec2 uv;
varying vec3 vtxv;
uniform float deltaLambda;
uniform sampler2D lightPosTex;
uniform sampler2D lightVelTex;
*/}));
						if (pos_or_vel == 'vel') {
							fsh.push(metricCodes[objType]);
						}
					
						fsh.push(mlstr(function(){/*
#define maxRadius 1e+6

void main() {
	vec4 pos = texture2D(lightPosTex, uv);
	vec4 vel = texture2D(lightVelTex, uv);

*/}));
						if (pos_or_vel == 'pos') {
							fsh.push(mlstr(function(){/*
	float distSq = dot(pos.xyz, pos.xyz);
	gl_FragColor = pos;
	if (distSq < maxRadius) {
		gl_FragColor += deltaLambda * vel;
		if (gl_FragColor != gl_FragColor) {
			gl_FragColor = pos;
		}
	}
}
*/}));
						} else if (pos_or_vel == 'vel') {
							fsh.push(mlstr(function(){/*
	float distSq = dot(pos.xyz, pos.xyz);
	gl_FragColor = vel;
	if (distSq < maxRadius) {
		gl_FragColor += deltaLambda * accel(pos, vel);
	}
}
*/}));
						}
					}

					var args = {
						vertexCode : vsh.join('\n'),
						fragmentCode : fsh.join('\n'),
					};

					args.context = thiz.glutil.context;
					args.uniforms = {
						lightPosTex : 0,
						lightVelTex : 1
					};
					args.vertexPrecision = 'best';
					args.fragmentPrecision = 'best';
					channel.shaders[objType][shaderType] = new thiz.glutil.ShaderProgram(args);
				});
			});
		});

		var cubeShader = new this.glutil.ShaderProgram({
			vertexPrecision : 'best',
			vertexCode : shaderCommonCode + mlstr(function(){/*
attribute vec2 vertex;
varying vec2 uv;
uniform mat4 projMat;
uniform vec4 angle;

const mat3 viewMatrix = mat3(
	0., 0., -1.,
	1., 0., 0., 
	0., -1., 0.);

void main() {
	uv = vertex;
	vec3 vtx3 = vec3(vertex * 2. - 1., 1.);
	vtx3 = quatRotate(vec4(angle.xyz, -angle.w), vtx3);
	vtx3 = viewMatrix * vtx3;
	vec4 vtx4 = vec4(vtx3, 1.);
	gl_Position = projMat * vtx4;
}
*/}),
			fragmentPrecision : 'best',
			fragmentCode : shaderCommonCode + mlstr(function(){/*
varying vec2 uv;
uniform samplerCube skyTex;
uniform sampler2D lightVelTex;
uniform vec4 viewAngle;

const mat3 viewMatrixInv = mat3(
	0., 1., 0.,
	0., 0., -1.,
	-1., 0., 0.);

void main() {
	vec3 dir = texture2D(lightVelTex, uv).xyz;
	dir = viewMatrixInv * viewMatrixInv * dir;
	dir = quatRotate(viewAngle, dir);
	gl_FragColor.xyz = textureCube(skyTex, dir).xyz;
	gl_FragColor.w = 1.; 
}
*/}),
			uniforms : {
				skyTex : 0,
				lightVelTex : 1
			}
		});

		//scene graph, 6 quads oriented in a cube
		//I would use a cubemap but the Mali-400 doesn't seem to want to use them as 2D FBO targets ...
		this.cubeSides = [];
		for (var side = 0; side < 6; ++side) {
			this.cubeSides[side] = new this.glutil.SceneObject({
				//geometry : glutil.unitQuad.geometry,
				mode : gl.TRIANGLE_STRIP,
				attrs : {
					vertex : glutil.unitQuadVertexBuffer
				},
				uniforms : {
					viewAngle : this.glutil.view.angle,
					angle : angleForSide[side]
				},
				shader : cubeShader,
				texs : [skyTex, this.lightPosVelChannels[1].texs[0][side]],
			});
		}
	},

	resetField : function() {
		gl.viewport(0, 0, this.lightTexWidth, this.lightTexHeight);
		$.each(this.lightPosVelChannels, function(_,channel) {
			for (var side = 0; side < 6; ++side) {
				var shader = channel.shaders[objectType].reset;
				
				var uniforms = {};
				if (channel.uniforms !== undefined) {
					for (var i = 0; i < channel.uniforms.length; ++i) {
						var uniformName = channel.uniforms[i];
						uniforms[uniformName] = window[uniformName];
					}
				}
				var rotationMat3 = mat3.create();
				mat3.fromQuat(rotationMat3, angleForSide[side]);
				uniforms.rotation = rotationMat3;
				
				var fbo = channel.fbos[0][side];
				fbo.draw({
					callback:function(){
						glutil.unitQuad.draw({
							shader:shader,
							uniforms:uniforms
						});
					}
				});
			}
		});
		gl.viewport(0, 0, this.glutil.canvas.width, this.glutil.canvas.height);
	},
	
	updateLightPosTex : function() {	
		var thiz = this;
		gl.viewport(0, 0, this.lightTexWidth, this.lightTexHeight);
		$.each(this.lightPosVelChannels, function(_,channel) {
			for (var side = 0; side < 6; ++side) {
				var shader = channel.shaders[objectType].iterate;
				
				var uniforms = {};
				if (channel.uniforms !== undefined) {
					for (var i = 0; i < channel.uniforms.length; ++i) {
						var uniformName = channel.uniforms[i];
						uniforms[uniformName] = window[uniformName];
					}
				}
				var rotationMat3 = mat3.create();
				mat3.fromQuat(rotationMat3, angleForSide[side]);
				uniforms.rotation = rotationMat3;
				uniforms.lightPosTex = 0;
				uniforms.lightVelTex = 1;
				
				var fbo = channel.fbos[1][side];
				fbo.draw({
					callback:function(){
						glutil.unitQuad.draw({
							shader:shader,
							uniforms:uniforms,
							texs:[
								thiz.lightPosVelChannels[0].texs[0][side],
								thiz.lightPosVelChannels[1].texs[0][side]
							]
						});
					}
				});
			}
		});
		gl.viewport(0, 0, this.glutil.canvas.width, this.glutil.canvas.height);
		$.each(this.lightPosVelChannels, function(_,channel) {
			var tmp;
			tmp = channel.texs[0];
			channel.texs[0] = channel.texs[1];
			channel.texs[1] = tmp;
			tmp = channel.fbos[0];
			channel.fbos[0] = channel.fbos[1];
			channel.fbos[1] = tmp;
		});
	},

	update : function() {
		this.updateLightPosTex();
		
		//turn on magnification filter
		for (var side = 0; side < 6; ++side) {
			var tex = this.lightPosVelChannels[1].texs[0][side].obj;
			gl.bindTexture(gl.TEXTURE_2D, tex);
			gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
			//gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
			//gl.generateMipmap(gl.TEXTURE_2D);
			gl.bindTexture(gl.TEXTURE_2D, null);
		}
		
		for (var side = 0; side < 6; ++side) {
			//texs[0] is skyTex
			this.cubeSides[side].texs[1] = this.lightPosVelChannels[1].texs[0][side];
		}

		this.glutil.draw();	
		
		//turn off magnification filter
		for (var side = 0; side < 6; ++side) {
			var tex = this.lightPosVelChannels[1].texs[0][side].obj;
			gl.bindTexture(gl.TEXTURE_2D, tex);
			gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
			//gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
			gl.bindTexture(gl.TEXTURE_2D, null);
		}	
	}
});

