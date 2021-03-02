#include "simulation.hpp"
#define SOLUTION

using namespace vcl;



#ifdef SOLUTION
void collision_sphere_plane(vcl::vec3& p, vcl::vec3& v, float r, vcl::vec3 const& n, vcl::vec3 const& p0)
{
    float const epsilon = 1e-5f;
    float const alpha_n = 0.95f;  // attenuation normal
    float const alpha_t = 0.90f;  // attenuation tangential

    float const s = dot(p-p0,n) - r;
    if( s<-epsilon )
    {
        vec3 const vn = dot(v,n) * n;
        vec3 const vt = v - vn;

        p = p - (s * n);
        v = -alpha_n * vn + alpha_t * vt;
    }

}
void collision_sphere_sphere(vcl::vec3& p1, vcl::vec3& v1, float r1, vcl::vec3& p2, vcl::vec3& v2, float r2)
{
	float const epsilon = 1e-5f;
    float const alpha = 0.95f;

    vec3 const p12 = p1-p2;
    float const d12 = norm(p12);

    if(d12 < r1+r2)
    {
        vec3 const u12 = p12/d12;
        float const collision_depth = r1+r2-d12;

        p1 += (collision_depth/2.0f+epsilon) * u12 ;
        p2 -= (collision_depth/2.0f+epsilon) * u12 ;

        if(norm(v1-v2)>0.2f){
            float const j = dot(v1-v2,u12);
            v1 = v1 - alpha*j*u12;
            v2 = v2 + alpha*j*u12;
        }
        else // Contact
        {
            v1 = v1/1.2f;
            v2 = v2/1.2f;
        }

    }
}
#endif

#ifdef SOLUTION
void simulate(std::vector<particle_structure>& particles, float dt_true)
{
	
	vec3 const g = {0,0,-9.81f};
	size_t const N_substep = 10;
	float const dt = dt_true/N_substep;
	for(size_t k_substep=0; k_substep<N_substep; ++k_substep)
	{
		size_t const N = particles.size();
		for (size_t k = 0; k < N; ++k)
		{
			particle_structure& particle = particles[k];

			vec3 const f = particle.m * g;

			particle.v = (1-0.9f*dt)*particle.v + dt*f;
			particle.p = particle.p + dt*particle.v;
		}

		// Collisions between spheres
		for(size_t k1=0; k1<N; ++k1)
		{
			for(size_t k2=k1+1; k2<N; ++k2)
			{
				particle_structure& p1 = particles[k1];
				particle_structure& p2 = particles[k2];

				collision_sphere_sphere(p1.p,p1.v,p1.r, p2.p,p2.v,p2.r);
			}
		}

		// Collisions with cube
		const std::vector<vec3> face_normal  = {{0, 1,0}, { 1,0,0}, {0,0, 1}, {0,-1,0}, {-1,0,0}, {0,0,-1}};
		const std::vector<vec3> face_position = {{0,-1,0}, {-1,0,0}, {0,0,-1}, {0, 1,0}, { 1,0,0}, {0,0, 1}};
		const size_t N_face = face_normal.size();
		for(size_t k=0; k<N; ++k){
			particle_structure& part = particles[k];
			for(size_t k_face=0; k_face<N_face; ++k_face)
				collision_sphere_plane(part.p, part.v, part.r, face_normal[k_face], face_position[k_face]);
		}
	}

}
#else
void simulate(std::vector<particle_structure>& particles, float dt)
{
	vec3 const g = {0,0,-9.81f};
	size_t const N = particles.size();
	for (size_t k = 0; k < N; ++k)
	{
		particle_structure& particle = particles[k];

		vec3 const f = particle.m * g;

		particle.v = (1-0.9f*dt)*particle.v + dt*f;
		particle.p = particle.p + dt*particle.v;
	}

	// To do :
	//  Handle collision ...
}
#endif