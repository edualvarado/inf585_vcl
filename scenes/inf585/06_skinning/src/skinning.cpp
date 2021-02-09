#include "skinning.hpp"
#define SOLUTION

namespace vcl
{
	void normalize_weights(buffer<buffer<float>>& weights)
	{
		size_t const N = weights.size();
		for (size_t k = 0; k < N; ++k) {
			float s = 0.0f;
			for(float w : weights[k]) s += w;
			assert_vcl_no_msg(s>1e-5f);
			for(float& w : weights[k]) w /= s;
		}
	}


	// Linear Blend Skinning
	void skinning_LBS_compute(
		buffer<vec3>& position_skinned,  // position to deform
		buffer<vec3>& normal_skinned,    // normal to deform
		buffer<affine_rt> const& skeleton_current,    // rigid transforms for skeleton joints in current pose
		buffer<affine_rt> const& skeleton_rest_pose,  // rigid transforms of skeleton joints in rest pose
		buffer<vec3> const& position_rest_pose,       // vertex positions of the mesh in rest pose
		buffer<vec3> const& normal_rest_pose,         // normal coordinates of the mesh in rest pose
		rig_structure const& rig)                     // information of the skinning weights (joints and weights associated to a vertex)
	{
		size_t const N_vertex = position_rest_pose.size();
		size_t const N_joint = skeleton_current.size();

		// Sanity check on sizes of buffers
		assert_vcl_no_msg(position_skinned.size()==N_vertex);
		assert_vcl_no_msg(normal_skinned.size()==N_vertex);
		assert_vcl_no_msg(normal_rest_pose.size()==N_vertex);
		assert_vcl_no_msg(skeleton_rest_pose.size()==N_joint);
		assert_vcl_no_msg(rig.joint.size()==N_vertex);
		assert_vcl_no_msg(rig.weight.size()==N_vertex);

		// To do
		//   Compute the Linear Blend Skinning ...

#ifdef SOLUTION
		buffer<mat4> M;
		M.resize(N_joint);
		for(size_t kj=0; kj<N_joint; ++kj)
			M[kj] = skeleton_current[kj].matrix() * inverse(skeleton_rest_pose[kj]).matrix();

		for (size_t k_vertex = 0; k_vertex < N_vertex; ++k_vertex)
		{
			vec3& p = position_skinned[k_vertex];
			vec3& n = normal_skinned[k_vertex];
			vec3 const& p0 = position_rest_pose[k_vertex];
			vec3 const& n0 = normal_rest_pose[k_vertex];

			size_t const N_dep = rig.joint[k_vertex].size();
			mat4 M_cumulative;
			for (size_t kj = 0; kj < N_dep; ++kj)
			{
				size_t const joint = rig.joint[k_vertex][kj];
				float const w = rig.weight[k_vertex][kj];

				M_cumulative += w*M[joint];
			}
			p = (M_cumulative*vec4(p0,1.0f)).xyz();
			n = (M_cumulative*vec4(n0,0.0f)).xyz();
		}
#endif

	}

#ifdef SOLUTION
quaternion_dual::quaternion_dual()
    :q(), qe()
{}
quaternion_dual::quaternion_dual(quaternion const& q_arg, quaternion const qe_arg)
    :q(q_arg),qe(qe_arg)
{}

quaternion_dual::quaternion_dual(quaternion const& q_arg, vec3 const& t)
    :q(q_arg), qe()
{
    qe = 0.5f * quaternion(t.x, t.y, t.z, 0.0f) * q_arg;
//    qe.x = 0.5f * ( t.x*q.w + t.y*q.z - t.z*q.y);
//    qe.y = 0.5f * (-t.x*q.z + t.y*q.w + t.z*q.x);
//    qe.z = 0.5f * ( t.x*q.y - t.y*q.x + t.z*q.w);
//    qe.w = 0.5f * (-t.x*q.x - t.y*q.y - t.z*q.z);
}

vcl::vec3 quaternion_dual::translation() const
{
    quaternion q_t = 2 * qe * conjugate(q);
    return {q_t.x, q_t.y, q_t.z};

//    float const tx = 2.0f * ( -qe.w*q.x + qe.x*q.w - qe.y*q.z + qe.z*q.y );
//    float const ty = 2.0f * ( -qe.w*q.y + qe.x*q.z + qe.y*q.w - qe.z*q.x );
//    float const tz = 2.0f * ( -qe.w*q.z - qe.x*q.y + qe.y*q.x + qe.z*q.w );
//    return {tx,ty,tz};
}

quaternion_dual& operator+=(quaternion_dual& a, quaternion_dual const& b)
{
    a.q += b.q;
    a.qe += b.qe;
    return a;
}
quaternion_dual operator+(quaternion_dual const& a, quaternion_dual const& b)
{
    return {a.q+b.q, a.qe+b.qe};
}
quaternion_dual operator*(float s, quaternion_dual const& d)
{
    return {s*d.q, s*d.qe};
}
quaternion_dual operator/(quaternion_dual const& d, float s)
{
    return {d.q/s, d.qe/s};
}

void skinning_DQS_compute(
		buffer<vec3>& position_skinned,
		buffer<vec3>& normal_skinned,
		buffer<affine_rt> const& skeleton_current,
		buffer<affine_rt> const& skeleton_rest_pose,
		buffer<vec3> const& position_rest_pose,
		buffer<vec3> const& normal_rest_pose,
		rig_structure const& rig)
	{
		size_t const N_vertex = position_rest_pose.size();
		size_t const N_joint = skeleton_current.size();

		assert_vcl_no_msg(position_skinned.size()==N_vertex);
		assert_vcl_no_msg(normal_skinned.size()==N_vertex);
		assert_vcl_no_msg(normal_rest_pose.size()==N_vertex);
		assert_vcl_no_msg(skeleton_rest_pose.size()==N_joint);
		assert_vcl_no_msg(rig.joint.size()==N_vertex);
		assert_vcl_no_msg(rig.weight.size()==N_vertex);

		buffer<quaternion_dual> dq;
		dq.resize(N_joint);
		for (size_t kj = 0; kj < N_joint; ++kj) {
			rotation const& r  = skeleton_current[kj].rotate;
			rotation const& r0 = skeleton_rest_pose[kj].rotate;
			vec3 const& t = skeleton_current[kj].translate;
			vec3 const& t0 = skeleton_rest_pose[kj].translate;
			dq[kj] = quaternion_dual( (r*inverse(r0)).data, r*inverse(r0)*(-t0)+t );
		}

		
		for (size_t k_vertex = 0; k_vertex < N_vertex; ++k_vertex)
		{
			vec3& p = position_skinned[k_vertex];
			vec3& n = normal_skinned[k_vertex];
			vec3 const& p0 = position_rest_pose[k_vertex];
			vec3 const& n0 = normal_rest_pose[k_vertex];

			size_t const N_dep = rig.joint[k_vertex].size();
			quaternion_dual d = { {0,0,0,0}, {0,0,0,0} };
			for (size_t kj = 0; kj < N_dep; ++kj)
			{
				size_t const joint = rig.joint[k_vertex][kj];
				float const w = rig.weight[k_vertex][kj];
				d += w*dq[joint];
			}

			float const q_n = norm(d.q);
			if( std::abs(q_n)>1e-5f)
				d = d/q_n;

			p = rotation(d.q)*p0 + d.translation();
			n = rotation(d.q)*n0;
		}

	}
#endif
}