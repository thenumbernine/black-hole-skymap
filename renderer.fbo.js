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
				'deltaLambda',
				'simTime'
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
/*
eh wikipedia ...
https://en.wikipedia.org/wiki/Derivation_of_the_Schwarzschild_solution
ds^2 = - (1-m/(2rho))^2 / (1+m/(2rho))^2 dt^2 + (1 + m/(2 rho))^4 (dx^2 + dy^2 + dz^2)
x = rho sin(theta) cos(phi)
y = rho sin(theta) cos(phi)
z = rho cos(theta)
... but then rho serves as the same definition as the radial coordinate of a spherical coordinate system
... so rho = sqrt(x^2 + y^2 + z^2)
... so we don't even need to use the Schwarzschild 'r' coordinate?
sure enough:
https://physics.stackexchange.com/questions/319133/isotropic-schwarzschild-coordinates?rq=1
https://descanso.jpl.nasa.gov/monograph/series2/Descanso2_all.pdf
yes, the isotropic 'rho' (not equal to the Schwarzschild-spherical 'r')
is what NASA interprets as ordinary euclidian distance
*/
				['Schwarzschild Black Hole'] : mlstr(function(){/*
uniform float blackHoleMass;

//common variables used by 
//* initialization of g_ab (in the vertex shader reset)
//* update of x''^a = -Gamma^a_bc x'^b x'^c (in the fragment shader update)
struct metricInfo_t {
	float inv_rho;
	float _1_plus_m_2rho;
	float _1_minus_m_2rho;
	float _1_plus_m_2rho_sq;
	float _1_plus_m_2rho_toTheFourth;
};

float sqr(float x) { return x * x; }

metricInfo_t init_metricInfo(vec4 pos) {
	float rho = length(pos.xyz);
	
	metricInfo_t m;
	m.inv_rho = 1. / rho;
	float _m_2rho = .5 * blackHoleMass * m.inv_rho;
	m._1_minus_m_2rho = 1. - _m_2rho;
	m._1_plus_m_2rho = 1. + _m_2rho;
	m._1_plus_m_2rho_sq = sqr(m._1_plus_m_2rho);
	m._1_plus_m_2rho_toTheFourth = sqr(m._1_plus_m_2rho_sq);
	return m;
}

// g_tt = -(1 - m/(2r))^2 / (1 + m/(2r))^2
float g_tt(metricInfo_t m) {
	float ratio = m._1_minus_m_2rho / m._1_plus_m_2rho;
	return -ratio * ratio;
}

// g_ti = 0
vec3 g_ti(metricInfo_t m) { 
	return vec3(0., 0., 0.); 
}

// g_ij = delta_ij (1 + m/(2r))^4
mat3 g_ij(metricInfo_t m) {
	return mat3(m._1_plus_m_2rho_toTheFourth);
}

vec4 accel(metricInfo_t m, vec4 pos, vec4 vel) {
	float inv_rho = m.inv_rho;
	float x_dot_v = dot(pos.xyz, vel.xyz);
	float velSq = dot(vel.xyz, vel.xyz);
	float inv_rhoSq = sqr(inv_rho);
	float inv_rho3 = inv_rho * inv_rhoSq;
	
	float _1_plus_m_2rho_tothe6 = m._1_plus_m_2rho_sq * m._1_plus_m_2rho_toTheFourth;
	
	vec4 result;
	result.w = 0.;
	result.xyz = -blackHoleMass * (
		vel.w * vel.w * pos.xyz * m._1_minus_m_2rho / _1_plus_m_2rho_tothe6
		+ velSq * pos.xyz
	) * inv_rho3 / m._1_plus_m_2rho;
	return result;
}

*/}),
				
				// substitution of a=0, Q=0, simlpified
				// This somewhat match the Schwarzschild isotropic geodesic above
				// The full Kerr geodesic below isn't working.
				['Kerr Black Hole degeneracy'] : mlstr(function(){/*
uniform float blackHoleMass;
uniform float blackHoleCharge;
uniform float blackHoleAngularVelocity;

struct metricInfo_t {
	vec3 pos;
	float inv_r;
};

metricInfo_t init_metricInfo(vec4 pos) {
	metricInfo_t m;
	m.pos = pos.xyz;
	float r = length(pos.xyz);
	m.inv_r = 1. / r;
	return m;
}

float g_tt(metricInfo_t m) {
	return -1. + 2. * blackHoleMass * m.inv_r;
}

vec3 g_ti(metricInfo_t m) { 
	return 2. * blackHoleMass * m.inv_r * m.inv_r * m.pos;
}

mat3 g_ij(metricInfo_t m) {
	float R = 2. * blackHoleMass;
	float R_r = R * m.inv_r;
	vec3 npos = m.pos.xyz * m.inv_r;	
	return mat3(1.) + R_r * outerProduct(npos, npos);
}

vec4 accel(metricInfo_t m, vec4 pos, vec4 vel) {
	float R = 2. * blackHoleMass;
	
	float inv_r = m.inv_r;
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


				['Kerr Black Hole'] : mlstr(function(){/*
uniform float blackHoleMass;
uniform float blackHoleCharge;
uniform float blackHoleAngularVelocity;

struct metricInfo_t {
	vec4 pos;
	float rSq;
	float r;
	float H;
	vec3 l;
};

metricInfo_t init_metricInfo(vec4 pos) {
	metricInfo_t m;
	m.pos = pos;
	float R = 2. * blackHoleMass;
	float a = blackHoleAngularVelocity;
	float aSq = a * a;
	float Q = blackHoleCharge;
	float QSq = Q * Q;
	float zSq = pos.z * pos.z;
	float b = aSq - dot(pos.xyz, pos.xyz);
	m.rSq = .5 * (-b + sqrt(b * b + 4. * aSq * zSq));
	m.r = sqrt(m.rSq);
	m.H = .5 * (m.r * R - QSq) / (m.rSq + aSq * zSq / m.rSq);
	m.l = vec3(
		(m.r * pos.x + a * pos.y) / (m.rSq + aSq),
		(m.r * pos.y - a * pos.x) / (m.rSq + aSq),
		pos.z / m.r);
	return m;
}

float g_tt(metricInfo_t m) {
	return -1. + 2. * m.H;
}

vec3 g_ti(metricInfo_t m) { 
	return 2. * m.H * m.l;
}

mat3 g_ij(metricInfo_t m) {
	return mat3(1.) + 2. * m.H * outerProduct(m.l, m.l);
}

vec4 accel(metricInfo_t m, vec4 pos, vec4 vel) {
	float zSq = pos.z * pos.z;
	float pos_pos = dot(pos.xyz, pos.xyz);

	float R = 2. * blackHoleMass;
	float Q = blackHoleCharge;
	float QSq = Q * Q;
	float a = blackHoleAngularVelocity;
	float aSq = a * a;
	
	float rSq = m.rSq;
	float r = m.r;
	
	float aSq_plus_rSq = rSq + aSq;
	
	vec3 zHat = vec3(0., 0., 1.);
	vec3 dr_dx = (pos.xyz * aSq_plus_rSq + zHat * aSq) / (r * (2. * rSq + aSq - pos_pos));
	
	float H_denom = rSq + aSq * zSq / rSq;
	float H = m.H;

	vec4 dH_dx;
	dH_dx.xyz = (
			.5 * dr_dx * R * H_denom
			- (r * R - QSq) * (
				r * dr_dx + aSq * pos.z / (r * rSq) * (zHat * r - pos.z * dr_dx)
			)
		) / (H_denom * H_denom);

	float aSq_plus_rSq_sq = aSq_plus_rSq * aSq_plus_rSq;  
	

	// l_a,b == dl_dx[a][b]
	mat4 dl_dx;
	dl_dx[0][0] = ((dr_dx.x * pos.x + r) * aSq_plus_rSq - (r * pos.x + a * pos.y) * (2. * r * dr_dx.x)) / aSq_plus_rSq_sq;
	dl_dx[0][1] = ((dr_dx.y * pos.x + a) * aSq_plus_rSq - (r * pos.x + a * pos.y) * (2. * r * dr_dx.y)) / aSq_plus_rSq_sq;
	dl_dx[0][2] = ((dr_dx.z * pos.x    ) * aSq_plus_rSq - (r * pos.x + a * pos.y) * (2. * r * dr_dx.z)) / aSq_plus_rSq_sq;
	dl_dx[1][0] = ((dr_dx.x * pos.y - a) * aSq_plus_rSq - (r * pos.y - a * pos.x) * (2. * r * dr_dx.x)) / aSq_plus_rSq_sq;
	dl_dx[1][1] = ((dr_dx.y * pos.y + r) * aSq_plus_rSq - (r * pos.y - a * pos.x) * (2. * r * dr_dx.y)) / aSq_plus_rSq_sq;
	dl_dx[1][2] = ((dr_dx.z * pos.y    ) * aSq_plus_rSq - (r * pos.y - a * pos.x) * (2. * r * dr_dx.z)) / aSq_plus_rSq_sq;
	dl_dx[2][0] = -pos.z * dr_dx.x / rSq;
	dl_dx[2][1] = -pos.z * dr_dx.y / rSq;
	dl_dx[2][2] = (r - pos.z * dr_dx.z) / rSq;

	vec4 l = vec4(m.l, 1.);
	float l_dot_vel = dot(l, vel);
	float dH_dx_dot_vel = dot(dH_dx, vel);
	float vel_dl_dx_vel = dot(vel, dl_dx * vel);
	vec4 conn_vel_vel = 
		l * 2. * (
			l_dot_vel * dH_dx_dot_vel
			+ H * vel_dl_dx_vel
		)
		- dH_dx * l_dot_vel * l_dot_vel
		+ 2. * H * l_dot_vel * (dl_dx - transpose(dl_dx)) * vel;

	vec4 lU = vec4(m.l, -1.);
#if 0
	mat4 gInv = mat4(
		vec4(1., 0., 0., 0.),
		vec4(0., 1., 0., 0.),
		vec4(0., 0., 1., 0.),
		vec4(0., 0., 0., -1.)
	) - 2. * H * outerProduct(lU, lU);

	// conn^a_bc = g^ad conn_dbc
	return -gInv * conn_vel_vel;
#endif
#if 1
	//g^ab = eta^ab - 2 H l*^a l*^b
	//c^a = g^ab c_b = (eta^ab - 2 H l*^a l*^b) c_b 
	//	= eta^ab c_b - 2 H l*^a l*^b c_b 
	//
	//c^t = eta^tb c_b - 2 H l*^t l*^b c_b 
	//c^t = -c_t + 2 H l*^b c_b 
	//
	//c^i = eta^ib c_b - 2 H l*^i l*^b c_b 
	//c^i = c_i - 2 H l*^i l*^b c_b 
	float conn_vel_vel_lU = dot(conn_vel_vel, lU);
	return -vec4(
		conn_vel_vel.xyz - lU.xyz * 2. * H * conn_vel_vel_lU,
		-conn_vel_vel.w + 2. * H * conn_vel_vel_lU);
#endif
}
*/}),
		
				
				['Kerr Black Hole full'] : mlstr(function(){/*
uniform float blackHoleMass;
uniform float blackHoleCharge;
uniform float blackHoleAngularVelocity;

struct metricInfo_t {
	vec3 pos;
};

metricInfo_t init_metricInfo(vec4 pos) {
	metricInfo_t m;
	m.pos = pos.xyz;
	return m;
}

float g_tt(metricInfo_t m) {
	vec3 pos = m.pos;
	float R = 2. * blackHoleMass;
	float Q = blackHoleCharge;
	float a = blackHoleAngularVelocity;
	float z = pos.z;
	float b = a*a - dot(pos,pos);
	float rSq = .5*(-b + sqrt(b*b + 4.*a*a*z*z));
	float r = sqrt(rSq);
	float H = .5 * (r * R - Q * Q) / (rSq + a*a*z*z/rSq);
	return -1. + 2. * H;
}

vec3 g_ti(metricInfo_t m) { 
	vec3 pos = m.pos;
	float R = 2. * blackHoleMass;
	float Q = blackHoleCharge;
	float a = blackHoleAngularVelocity;
	float x = pos.x;
	float y = pos.y;
	float z = pos.z;
	float b = a*a - dot(pos,pos);
	float rSq = .5*(-b + sqrt(b*b + 4.*a*a*z*z));
	float r = sqrt(rSq);
	float H = .5 * (r * R - Q * Q) / (rSq + a*a*z*z/rSq);
	vec3 l = vec3(
		(r*x + a*y)/(rSq + a*a),
		(r*y - a*x)/(rSq + a*a),
		z/r);
	return 2. * H * l;
}

mat3 g_ij(metricInfo_t m) {
	vec3 pos = m.pos;
	float R = 2. * blackHoleMass;
	float Q = blackHoleCharge;
	float a = blackHoleAngularVelocity;
	float x = pos.x;
	float y = pos.y;
	float z = pos.z;
	float b = a*a - dot(pos,pos);
	float rSq = .5*(-b + sqrt(b*b + 4.*a*a*z*z));
	float r = sqrt(rSq);
	float H = .5 * (r * R - Q * Q) / (rSq + a*a*z*z/rSq);
	vec3 l = vec3(
		(r*x + a*y)/(rSq + a*a),
		(r*y - a*x)/(rSq + a*a),
		z/r);
	return mat3(1.) + 2. * H * outerProduct(l,l);
}

vec4 accel(metricInfo_t m, vec4 pos, vec4 vel) {
	float x = pos.x;
	float y = pos.y;
	float z = pos.z;

	float R = 2. * blackHoleMass;
	float Q = blackHoleCharge;
	float a = blackHoleAngularVelocity;
	float aSq = a*a;
	float b = aSq - dot(pos.xyz, pos.xyz);
	float sqrtdiscr = sqrt(b*b + 4.*aSq*z*z);
	float rSq = .5*(-b + sqrtdiscr);
	float r = sqrt(rSq);
	float H = .5 * (r * R - Q * Q) / (rSq + aSq*z*z/rSq);

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
			dr_dxis[i] * R * (rSq + aSq*z*z / rSq)
			- (r * R - Q*Q) * (
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

struct metricInfo_t {
	float r;
	float sigmaFront;
	float sigmaCenter;
	float sigmaBack;
	float tanhSigmaCenter;
	float f;
};

metricInfo_t init_metricInfo(vec4 pos) {
	metricInfo_t m;
	m.r = length(pos.xyz);
	m.sigmaFront = warpBubbleThickness * (m.r + warpBubbleRadius);
	m.sigmaCenter = warpBubbleThickness * warpBubbleRadius;
	m.sigmaBack = warpBubbleThickness * (m.r - warpBubbleRadius);
	m.tanhSigmaCenter = tanh(m.sigmaCenter);
	m.f = (tanh(m.sigmaFront) - tanh(m.sigmaBack)) / (2. * m.tanhSigmaCenter);
	return m;
}

float g_tt(metricInfo_t m) {
	float vf = m.f * warpBubbleVelocity;
	float vf2 = vf * vf;
	return -1. + vf2; 
}

vec3 g_ti(metricInfo_t m) {
	return vec3(-m.f * warpBubbleVelocity, 0., 0.);
}

mat3 g_ij(metricInfo_t m) {
	return mat3(1.);
}

vec4 accel(metricInfo_t m, vec4 pos, vec4 vel) {
	float sechDiff = sechSq(m.sigmaFront) - sechSq(m.sigmaBack);
	float dfScalar = sechDiff / (2. * m.r * m.tanhSigmaCenter);
	
	vec4 df;
	df.xyz = warpBubbleThickness * pos.xyz * dfScalar;
	df.w = -warpBubbleVelocity * warpBubbleThickness * pos.x * dfScalar;

	float u = m.f * warpBubbleVelocity;
	vec4 du = df * warpBubbleVelocity;

	float du_v_dot3 = dot(du.xyz, vel.xyz);
	float du_v_dot4 = dot(du, vel);

	vec4 result;
	result.w = 
		- vel.x * du_v_dot3
		+ u * vel.w * (du.x * vel.x + du_v_dot3)
		- u * u * vel.w * vel.w * du.x
	;
	result.xyz = -du.xyz * vel.w * (vel.x - vel.w * u);
	result.x += 
		- vel.x * u * du_v_dot3
		+ vel.w * (
			du_v_dot4
			+ u * u * (du.x * vel.x + du_v_dot3)
			- vel.w * du.x * u * u * u
		);
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

vec4 calc_pos0() {
	vec4 pos = vec4(-objectDist, 0., 0., 0.);
	pos.xyz = quatRotate(quatConj(objectAngle), pos.xyz);
	return pos;
}

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
	gl_FragColor = calc_pos0();
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
float reset_w(vec4 pos, vec4 vel) {
	metricInfo_t m = init_metricInfo(pos);
	float _g_tt = g_tt(m);
	float negInv_g_tt = -1. / _g_tt;
	vec3 _g_ti = g_ti(m);
	float v_t = dot(vel.xyz, _g_ti);	//v_t = g_ti v^i
	float a = v_t * negInv_g_tt;
	float bsq = a * a + dot(vel.xyz, g_ij(m) * vel.xyz) * negInv_g_tt;
	if (bsq < 0.) bsq = 0.;
	return a + sqrt(bsq);
}

//relies on objectAngle
vec4 calc_vel0(vec4 pos, vec3 vtxv) {
	vec4 vel = vec4(normalize(vtxv.xyz), 0.);
	vel.xyz = quatRotate(quatConj(objectAngle), vel.xyz);
	vel.w = reset_w(pos, vel);
	return vel;
}

varying vec3 vtxv;
void main() {
	vec4 pos = calc_pos0();
	gl_FragColor = calc_vel0(pos, vtxv.xyz);	
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
		gl_FragColor += deltaLambda * accel(init_metricInfo(pos), pos, vel);
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

		var cubeVertexCode = shaderCommonCode + mlstr(function(){/*
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
*/});

		var cubeShaderUniforms = {
			skyTex : 0,
			lightPosTex : 1,
			lightVelTex : 2,
			hsvTex : 3
		};

		var cubeFragmentHeader = shaderCommonCode + mlstr(function(){/*
varying vec2 uv;
uniform samplerCube skyTex;
uniform sampler2D lightPosTex;
uniform sampler2D lightVelTex;
uniform vec4 viewAngle;

const mat3 viewMatrixInv = mat3(
	0., 1., 0.,
	0., 0., -1.,
	-1., 0., 0.);
*/});

		cubeBackgroundShader = new this.glutil.ShaderProgram({
			vertexPrecision : 'best',
			fragmentPrecision : 'best',
			uniforms : cubeShaderUniforms,
			vertexCode : cubeVertexCode,
			fragmentCode : cubeFragmentHeader + mlstr(function(){/*
void main() {
	vec3 dir = texture2D(lightVelTex, uv).xyz;
	dir = viewMatrixInv * viewMatrixInv * dir;
	dir = quatRotate(viewAngle, dir);
	gl_FragColor = vec4(textureCube(skyTex, dir).xyz, 1.);
}
*/})
		});
		
		cubePosXYZShader = new this.glutil.ShaderProgram({
			vertexPrecision : 'best',
			fragmentPrecision : 'best',
			uniforms : cubeShaderUniforms,
			vertexCode : cubeVertexCode,
			fragmentCode : cubeFragmentHeader + mlstr(function(){/*
void main() {
	vec3 pos = texture2D(lightPosTex, uv).xyz;
	pos = normalize(pos);
	pos = viewMatrixInv * viewMatrixInv * pos;
	pos = quatRotate(viewAngle, pos);
	gl_FragColor = vec4(.5 + .5 * pos.xyz, 1.);
}
*/})
		});

		cubeVelXYZShader = new this.glutil.ShaderProgram({
			vertexPrecision : 'best',
			fragmentPrecision : 'best',
			uniforms : cubeShaderUniforms,
			vertexCode : cubeVertexCode,
			fragmentCode : cubeFragmentHeader + mlstr(function(){/*
void main() {
	vec3 dir = texture2D(lightVelTex, uv).xyz;
	dir = normalize(dir);
	dir = viewMatrixInv * viewMatrixInv * dir;
	dir = quatRotate(viewAngle, dir);
	gl_FragColor = vec4(.5 + .5 * dir, 1.);
}
*/})
		});


		cubePosTShader = new this.glutil.ShaderProgram({
			vertexPrecision : 'best',
			fragmentPrecision : 'best',
			uniforms : cubeShaderUniforms,
			vertexCode : cubeVertexCode,
			fragmentCode : cubeFragmentHeader + mlstr(function(){/*
uniform float simTime;	//lambda
uniform sampler2D hsvTex;
void main() {
	vec4 pos = texture2D(lightPosTex, uv);
	float t = pos.w;
	float ratio = t / simTime;	//global time / local time
	float r = .5 * log(ratio);
	gl_FragColor = texture2D(hsvTex, vec2(r, .5));
}
*/})
		});

		cubeVelTShader = new this.glutil.ShaderProgram({
			vertexPrecision : 'best',
			fragmentPrecision : 'best',
			uniforms : cubeShaderUniforms,
			vertexCode : cubeVertexCode,
			fragmentCode : cubeFragmentHeader + mlstr(function(){/*
uniform sampler2D hsvTex;
void main() {
	vec4 vel = texture2D(lightVelTex, uv);
	// TODO divide by the original vel.w, which is g_tt of the initial position ...
	float vt = vel.w;
	float r = .5 * log(vt);
	gl_FragColor = texture2D(hsvTex, vec2(r, .5));
}
*/})
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
				shader : cubeBackgroundShader,
				texs : [
					skyTex,
					this.lightPosVelChannels[0].texs[0][side],
					this.lightPosVelChannels[1].texs[0][side],
					hsvTex
				],
			});
		}
	},

	resetField : function() {
		simTime = 0;
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

	runSimulation : true,

	update : function() {
		if (skyboxRenderer.runSimulation) {
			this.updateLightPosTex();
			simTime += deltaLambda;
		}
		
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
			this.cubeSides[side].texs[1] = this.lightPosVelChannels[0].texs[0][side];
			this.cubeSides[side].texs[2] = this.lightPosVelChannels[1].texs[0][side];
			//texs[3] is hsvTex
		}

		//if we are drawing the background (skytex lookup) ...
		var shader;
		if (drawMethod == 'background') {
			shader = cubeBackgroundShader;
		} else if (drawMethod == 'pos_xyz') {
			shader = cubePosXYZShader;
		} else if (drawMethod == 'vel_xyz') {
			shader = cubeVelXYZShader;
		} else if (drawMethod == 'pos_t') {
			shader = cubePosTShader;
		} else if (drawMethod == 'vel_t') {
			shader = cubeVelTShader;
		}
		for (var side = 0; side < 6; ++side) {
			this.cubeSides[side].shader = shader;
			
			//the only dynamic variable
			this.cubeSides[side].uniforms.simTime = simTime;
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

