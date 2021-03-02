#include "simulation.hpp"
#define SOLUTION
using namespace vcl;

#ifdef SOLUTION
static vec3 spring_force(const vec3& pi, const vec3& pj, float L0, float K)
{
    vec3 const p = pi-pj;
    float const L = norm(p);
    vec3 const u = p/L;

    vec3 const F = -K * (L-L0) * u;
    return F;
}
#endif

// Fill value of force applied on each particle
// - Gravity
// - Drag
// - Spring force
// - Wind force
void compute_forces(grid_2D<vec3>& force, grid_2D<vec3> const& position, grid_2D<vec3> const& velocity, grid_2D<vec3>& normals, simulation_parameters const& parameters, float wind_magnitude)
{
    size_t const N = force.size();        // Total number of particles of the cloth Nu x Nv
    size_t const N_dim = force.dimension[0]; // Number of particles along one dimension (square dimension)

    float const K  = parameters.K;
    float const m  = parameters.mass_total/N;
    float const mu = parameters.mu;
    float const	L0 = 1.0f/(N_dim-1.0f);

    // Gravity
    const vec3 g = {0,0,-9.81f};
    for(size_t k=0; k<N; ++k)
        force[k] = m*g;

    // Drag
    for(size_t k=0; k<N; ++k)
        force[k] += -mu*m*velocity[k];


#ifdef SOLUTION
    // Springs
    const int N_neighbors = 2;

    for(int ku=0; ku<N_dim; ++ku) {
        for(int kv=0; kv<N_dim; ++kv) {
            vec3& f = force(ku,kv);

            // Loops over neighbors
            for(int du=-N_neighbors; du<=N_neighbors; ++du) {
                for(int dv=-N_neighbors; dv<=N_neighbors; ++dv) {
                    if(du!=0 || dv!=0)
                    {
                        // Neighbors indices
                        int const ku_n = ku+du;
                        int const kv_n = kv+dv;
                        float const alpha = float(std::sqrt(du*du+dv*dv)); // rest-length

                        if( ku_n>=0 && ku_n<N_dim && kv_n>=0 && kv_n<N_dim)
                            f += spring_force(position(ku,kv), position(ku_n,kv_n), alpha*L0, K/alpha) ;
                    }
                }
            }
        }
    }

    // Wind
    if(std::abs(wind_magnitude)>1e-5f){
        const vec3 wind = wind_magnitude * vec3(-1.0f,0,0);
        const vec3 wind_u = normalize(wind);
        for(size_t k=0; k<N; ++k)
        {
            const vec3& n = normals[k];
            const float w = std::abs(dot(wind,n));
            const vec3 f = w * wind_u * L0*L0;

            force[k] += 0*f;
        }
    }
#else
    // TO DO: Add spring forces ...
#endif


}

void numerical_integration(grid_2D<vec3>& position, grid_2D<vec3>& velocity, grid_2D<vec3> const& force, float mass, float dt)
{
    size_t const N = position.size();

    for(size_t k=0; k<N; ++k)
    {
        velocity[k] = velocity[k] + dt*force[k]/mass;
        position[k] = position[k] + dt*velocity[k];
    }

}

void apply_constraints(grid_2D<vec3>& position, grid_2D<vec3>& velocity, std::map<size_t, vec3> const& positional_constraints, obstacles_parameters const& obstacles)
{
    // Fixed positions of the cloth
    for(const auto& constraints : positional_constraints)
        position[constraints.first] = constraints.second;

#ifdef SOLUTION
    const size_t N = position.size();
    const float epsilon = 5e-3f;
    for(size_t k=0; k<N; ++k)
    {
        vec3& p = position[k];
        vec3& v = velocity[k];

        // Ground constraint
        {
            if( p.z<=obstacles.z_ground+epsilon)
            {
                p.z = obstacles.z_ground+epsilon;
                v.z = 0.0f;
            }
        }

        // Sphere constraint
        {
            vec3 const& p0 = obstacles.sphere_center;
            float const r = obstacles.sphere_radius;
            if( norm(p-p0)<(r+2*epsilon) )
            {
                const vec3 u = normalize(p-p0);
                p = (r+2*epsilon)*u+p0;
                v = v-dot(v,u)*u;
            }
        }
    }
#else
    // To do: apply external constraints
#endif
}


void initialize_simulation_parameters(simulation_parameters& parameters, float L_cloth, size_t N_cloth)
{
	parameters.mass_total = 0.8f;
	parameters.K  = 5.0f;
	parameters.mu = 10.0f;
}

bool detect_simulation_divergence(grid_2D<vec3> const& force, grid_2D<vec3> const& position)
{
    bool simulation_diverged = false;
    const size_t N = position.size();
    for(size_t k=0; simulation_diverged==false && k<N; ++k)
    {
        const float f = norm(force[k]);
        const vec3& p = position[k];

        if( std::isnan(f) ) // detect NaN in force
        {
            std::cout<<"NaN detected in forces"<<std::endl;
            simulation_diverged = true;
        }

        if( f>600.0f ) // detect strong force magnitude
        {
            std::cout<<" **** Warning : Strong force magnitude detected "<<f<<" at vertex "<<k<<" ****"<<std::endl;
            simulation_diverged = true;
        }

        if( std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.z) ) // detect NaN in position
        {
            std::cout<<"NaN detected in positions"<<std::endl;
            simulation_diverged = true;
        }
    }

    return simulation_diverged;

}