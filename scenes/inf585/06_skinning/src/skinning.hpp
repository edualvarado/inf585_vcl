#pragma once

#include "vcl/vcl.hpp"
#define SOLUTION


namespace vcl
{
	struct rig_structure
	{
		buffer<buffer<int>> joint;
		buffer<buffer<float>> weight;
	};

	void normalize_weights(buffer<buffer<float>>& weights);

	void skinning_LBS_compute(
		buffer<vec3>& position_skinned,
		buffer<vec3>& normal_skinned,
		buffer<affine_rt> const& skeleton_current,
		buffer<affine_rt> const& skeleton_rest_pose,
		buffer<vec3> const& position_rest_pose, 
		buffer<vec3> const& normal_rest_pose, 
		rig_structure const& rig);

#ifdef SOLUTION
struct quaternion_dual
{
    quaternion q;
    quaternion qe;

    quaternion_dual();
    quaternion_dual(quaternion const& q, quaternion const qe);
    quaternion_dual(quaternion const& q, vcl::vec3 const& tr);

    vcl::vec3 translation() const;
};

quaternion_dual& operator+=(quaternion_dual& a, quaternion_dual const& b);
quaternion_dual operator+(quaternion_dual const& a, quaternion_dual const& b);
quaternion_dual operator*(float s, quaternion_dual const& d);
quaternion_dual operator/(quaternion_dual const& d, float s);

void skinning_DQS_compute(
		buffer<vec3>& position_skinned,
		buffer<vec3>& normal_skinned,
		buffer<affine_rt> const& skeleton_current,
		buffer<affine_rt> const& skeleton_rest_pose,
		buffer<vec3> const& position_rest_pose, 
		buffer<vec3> const& normal_rest_pose, 
		rig_structure const& rig);
#endif

}