#include "ffd.hpp"
#define SOLUTION


using namespace vcl;

// Helper to compute permutation(n, k) - k among n
//   avoids the recursive formula (too costly)
int binomial_coeff(int n, int k)
{
    int res = 1;
    if(k>n-k)
        k = n - k;
    for(int i=0; i<k; ++i) {
        res *= (n-i);
        res /= (i+1);
    }
    return res;
}


#ifdef SOLUTION
float bezier_coeff(const float t, int k, int N)
{
    float const Nk = float(binomial_coeff(N,k));
    return Nk * float(std::pow(t,k)) * float(std::pow(1-t, N-k));
}


void ffd_deform(buffer<vec3>& position, grid_3D<vec3> const& grid, buffer<grid_3D<float> > const& weights)
{
	size_t const Nx = grid.dimension.x;
	size_t const Ny = grid.dimension.y;
	size_t const Nz = grid.dimension.z;

	size_t const N_vertex = position.size();
	for (size_t k = 0; k < N_vertex; ++k)
	{
		vec3 p;
		for (size_t kx = 0; kx < Nx; ++kx) {
			for (size_t ky = 0; ky < Ny; ++ky) {
				for (size_t kz = 0; kz < Nz; ++kz) {
					p += weights[k](kx,ky,kz) * grid(kx,ky,kz);
				}
			}
		}
		
		position[k] = p;
	}

}

buffer<grid_3D<float> > precompute_weights(buffer<vec3>& position, int Nx, int Ny, int Nz)
{
	buffer<grid_3D<float>> weights;
	weights.resize(position.size());

	for(size_t k=0; k<weights.size(); ++k){
		vec3 const& p = position[k];
		weights[k].resize(Nx,Ny,Nz);
		for (int kx = 0; kx < Nx; ++kx) {
			for (int ky = 0; ky < Ny; ++ky) {
				for (int kz = 0; kz < Nz; ++kz) {
					
					float const alpha_x = bezier_coeff(p.x, kx, Nx-1 );
					float const alpha_y = bezier_coeff(p.y, ky, Ny-1 );
					float const alpha_z = bezier_coeff(p.z, kz, Nz-1 );

					weights[k](kx,ky,kz) = alpha_x * alpha_y * alpha_z;
				}
			}
		}
	}

	return weights;
}
#else
// Computation of the FFD deformation on the position with respect to the grid
void ffd_deform(buffer<vec3>& position, grid_3D<vec3> const& grid /** You may want to add extra parameters*/)
{
	// Get dimension of the grid
	size_t const Nx = grid.dimension.x;
	size_t const Ny = grid.dimension.y;
	size_t const Nz = grid.dimension.z;

	// Number of position to deform
	size_t const N_vertex = position.size();
	for (size_t k = 0; k < N_vertex; ++k)
	{
		// Loop over all grid points
		for (size_t kx = 0; kx < Nx; ++kx) {
			for (size_t ky = 0; ky < Ny; ++ky) {
				for (size_t kz = 0; kz < Nz; ++kz) {

					// TO DO: Should do some computations depending on the relative coordinates of the vertex k (you may need extra parameters to handle this), and the grid position.
					// Note: A grid position can be accessed as grid(kx,ky,kz)
					// ...
				}
			}
		}
		
	}
}
#endif