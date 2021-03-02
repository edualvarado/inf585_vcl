#include "simulation.hpp"

using namespace vcl;

#define SOLUTION

void divergence_free(grid_2D<vec2>& new_velocity, grid_2D<vec2> const& velocity, grid_2D<float>& divergence, grid_2D<float>& gradient_field)
{
    // v = projection of v0 on divergence free vector field
    //
    // v : Final vector field to be filled
    // v0: Initial vector field (non divergence free)
    // divergence: temporary buffer used to compute the divergence of v0
    // gradient_field: temporary buffer used to compute v = v0 - nabla(gradient_field)


#ifdef SOLUTION
    size_t const Nx = velocity.dimension[0];
    size_t const Ny = velocity.dimension[1];

    for(size_t x=1; x<Nx-1; ++x){
        for(size_t y=1; y<Ny-1; ++y){
            divergence(x,y) = 0.5f*(velocity(x+1,y).x-velocity(x-1,y).x + velocity(x,y+1).y-velocity(x,y-1).y);
        }
    }

    gradient_field.fill(0);
    set_boundary(divergence);
    set_boundary(gradient_field);

    // Gauss Seidel
    grid_2D<float>& q = gradient_field;
    for(size_t k_iter=0; k_iter<20; ++k_iter)
    {
        for(size_t x=1; x<Nx-1; ++x){
            for(size_t y=1; y<Ny-1; ++y){
                q(x,y) = (q(x+1,y)+q(x-1,y)+q(x,y+1)+q(x,y-1)-divergence(x,y))/4.0f;
            }
        }
        set_boundary(q);
    }

    // Velocity
    for(size_t x=1; x<Nx-1; ++x){
        for(size_t y=1; y<Ny-1; ++y){
            new_velocity(x,y) = velocity(x,y) - vec2(q(x+1,y)-q(x-1,y), q(x,y+1)-q(x,y-1))/2.0f;
        }
    }
    set_boundary_reflective(new_velocity);
#else
    // TO do:
    // 1. Compute divergence of v0
    // 2. Compute gradient_field such that nabla(gradient_field)^2 = div(v0)
    // 3. Compute v = v0 - nabla(gradient_field)
#endif

}