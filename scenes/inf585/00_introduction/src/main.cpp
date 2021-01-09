/**
	Introduction to use the VCL library
*/


#include "vcl/vcl.hpp"
#include <iostream>


// Add vcl namespace within the current one - Allows to use function from vcl library without explicitely preceeding their name with vcl::
using namespace vcl;


// ****************************************** //
// Structures associated to the current scene
//   In general these structures will be pre-coded 
//   for you in the exercises
// ****************************************** //

// Structure storing the variables used in the GUI interface
struct gui_parameters {
	bool display_frame = true; // Display a frame representing the coordinate system
	bool is_wireframe = false;
};


// Structure storing user-related interaction data and GUI parameter
struct user_interaction_parameters {
	vec2 mouse_prev;     // Current position of the mouse
	bool cursor_on_gui;  // Indicate if the cursor is on the GUI widget
	gui_parameters gui;  // The gui structure
};
user_interaction_parameters user; // (declaration of user as a global variable)


// Structure storing the global variable of the 3D scene
//   can be use to send uniform parameter when displaying a shape
struct scene_environment
{
	camera_around_center camera; // A camera looking at, and rotating around, a specific center position
	mat4 projection;             // The perspective projection matrix
	vec3 light;                  // Position of the light in the scene
};
scene_environment scene;  // (declaration of scene as a global variable)


// ****************************************** //
// Functions signatures
// ****************************************** //

// Callback functions
//   Functions called when a corresponding event is received by GLFW (mouse move, keyboard, etc).
void mouse_move_callback(GLFWwindow* window, double xpos, double ypos);
void window_size_callback(GLFWwindow* window, int width, int height);

// Specific functions:
//    Initialize the data of this scene - executed once before the start of the animation loop.
void initialize_data();                      
//    Display calls - executed a each frame
void display_scene(float current_time);
//    Display the GUI widgets
void display_interface();
// Deformation
void evolve_shape(float time);



// ****************************************** //
// Declaration of Global variables
// ****************************************** //

mesh_drawable global_frame;
mesh_drawable cube;
mesh_drawable ground;
mesh_drawable cylinder;
mesh_drawable sphere;
curve_drawable curve;

timer_basic timer;

// Deformation
mesh shape;
buffer <vec3> initial_position;
mesh_drawable shape_visual;


// ****************************************** //
// Functions definitions
// ****************************************** //


// Main function with creation of the scene and animation loop
int main(int, char* argv[])
{
	std::cout << "Run " << argv[0] << std::endl;

	// create GLFW window and initialize OpenGL
	GLFWwindow* window = create_window(1280,1024); 
	window_size_callback(window, 1280, 1024);
	std::cout << opengl_info_display() << std::endl;;

	imgui_init(window); // Initialize GUI library

	// Set GLFW callback functions
	glfwSetCursorPosCallback(window, mouse_move_callback); 
	glfwSetWindowSizeCallback(window, window_size_callback);
	
	std::cout<<"Initialize data ..."<<std::endl;
	initialize_data();

	std::cout<<"Start animation loop ..."<<std::endl;
	timer.start();
	glEnable(GL_DEPTH_TEST);
	while (!glfwWindowShouldClose(window))
	{
		scene.light = scene.camera.position();
		timer.update(); // update the time at this current frame

		// Clear screen
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_DEPTH_BUFFER_BIT);

		// Create GUI interface for the current frame
		imgui_create_frame();
		ImGui::Begin("GUI",NULL,ImGuiWindowFlags_AlwaysAutoResize);

		// Set the GUI interface (widgets: buttons, checkbox, sliders, etc)
		display_interface();

		// Display the objects of the scene
		display_scene(timer.t);


		// Display GUI
		ImGui::End();
		imgui_render_frame(window);

		// Swap buffer and handle GLFW events
		glfwSwapBuffers(window);
		glfwPollEvents();
	}


	imgui_cleanup();
	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}

void initialize_data()
{
	// Load and set the common shaders
	// *************************************** //

	//   - Shader used to display meshes
	GLuint const shader_mesh = opengl_create_shader_program(opengl_shader_preset("mesh_vertex"), opengl_shader_preset("mesh_fragment"));
	//   - Shader used to display constant color (ex. for curves)
	GLuint const shader_single_color = opengl_create_shader_program(opengl_shader_preset("single_color_vertex"), opengl_shader_preset("single_color_fragment"));
	//   - Default white texture
	GLuint const texture_white = opengl_texture_to_gpu(image_raw{1,1,image_color_type::rgba,{255,255,255,255}});

	// Set default shader and texture to drawable mesh
	mesh_drawable::default_shader = shader_mesh;
	mesh_drawable::default_texture = texture_white;
	curve_drawable::default_shader = shader_single_color;

	// Set the initial position of the camera
	// *************************************** //

	vec3 const camera_position = {2.0f, -3.5f, 2.0f};        // position of the camera in space
	vec3 const camera_target_position = {0,0,0}; // position the camera is looking at / point around which the camera rotates
	vec3 const up = {0,0,1};                     // approximated "up" vector of the camera
	scene.camera.look_at(camera_position, camera_target_position, up); 

	// Prepare the objects visible in the scene
	// *************************************** //

	// Create a visual frame representing the coordinate system
	global_frame = mesh_drawable(mesh_primitive_frame());

	// Create a cube as a mesh
	mesh const cube_mesh = mesh_primitive_cube(/*center*/{0,0,0}, /*edge length*/ 1.0f);
	// Create a mesh drawable from a mesh structure
	//   - mesh : store buffer of data (vertices, indices, etc) on the CPU. The mesh structure is convenient to manipulate in the C++ code but cannot be displayed (data is not on GPU).
	//   - mesh_drawable : store VBO associated to elements on the GPU + associated uniform parameters. A mesh_drawable can be displayed using the function draw(mesh_drawable, scene). It only stores the indices of the buffers on the GPU - the buffer of data cannot be directly accessed in the C++ code through a mesh_drawable.
	//   Note: a mesh_drawable can be created from a mesh structure in calling explicitely the constructor mesh_drawable(mesh)
	cube = mesh_drawable(cube_mesh);  // note: cube is a mesh_drawable declared as a global variable
	cube.shading.color = {1,1,0};     // set the color of the cube (R,G,B) - sent as uniform parameter to the shader when display is called

	// Create the ground plane
	ground  = mesh_drawable(mesh_primitive_quadrangle({-2,-2,-1},{ 2,-2,-1},{ 2, 2,-1},{-2, 2,-1}));

	// Create the cylinder
	cylinder = mesh_drawable(mesh_primitive_cylinder(/*radius*/ 0.2f, /*first extremity*/ {0,-1,0}, /*second extremity*/{0,1,0}));
	cylinder.shading.color = {0.8f,0.8f,1};

	// Create the sphere
	mesh const sphere_mesh = mesh_primitive_sphere(0.3f, {0,-1.2,0});
	sphere = mesh_drawable(sphere_mesh);
	sphere.shading.color = { 0,1,1 };

	// Create a parametric curve
    // **************************************** //
	buffer<vec3> curve_positions;   // the basic structure of a curve is a vector of vec3
	size_t const N_curve = 150;     // Number of samples of the curve
    for(size_t k=0; k<N_curve; ++k)
    {
        const float u = k/(N_curve-1.0f); // u \in [0,1]

        // curve oscillating as a cosine
        const float x = 0;
        const float y = 4.0f * (u - 0.5f);
        const float z = 0.1f * std::cos(u*16*3.14f);

        curve_positions.push_back({x,y,z});
    }
    // send data to GPU and store it into a curve_drawable structure
    curve = curve_drawable(curve_positions);
	curve.color = {0,1,0};

	// Deformable mesh
	size_t const N = 100;
	shape = mesh_primitive_grid({ 0,0,0 }, { 1,0,0 }, { 1,1,0 }, { 0,1,0 }, N, N);
	initial_position = shape.position;
	shape_visual = mesh_drawable(shape);
	shape_visual.shading.color = { 0.6f, 0.6f, 0.9f };
}



void display_scene(float time)
{
	// the general syntax to display a mesh is:
    //   draw(objectDrawableName, scene);
	//     Note: scene is used to set the uniform parameters associated to the camera, light, etc. to the shader
	draw(ground, scene);


	if(user.gui.display_frame) // conditional display of the global frame (set via the GUI)
		draw(global_frame, scene);

	if (user.gui.is_wireframe) // conditional display of the global frame (set via the GUI)
		draw_wireframe(sphere, scene, { 0,0,1 });

	// Display cylinder
    // ********************************************* //


	// Cylinder rotated around the axis (1,0,0), by an angle = time/2
	vec3 const axis_of_rotation = {1,0,0};
	float const angle_of_rotation = time/2.0f;
	rotation const rotation_cylinder = rotation(axis_of_rotation, angle_of_rotation);

	// Set translation and rotation parameters (send and used in shaders using uniform variables)
	cylinder.transform.rotate = rotation_cylinder;
	cylinder.transform.translate = {1.5f,0,0};
	draw(cylinder, scene); // Display of the cylinder

	// Meshes can also be displayed as wireframe using the specific draw_wireframe call
	draw_wireframe(cylinder, scene, /*color of the wireframe*/ {1,0.3f,0.3f} );
	
	// Display cube
    // ********************************************* //

	cube.transform.rotate = rotation({0,0,1},std::sin(3*time));
    cube.transform.translate = {-1,0,0};
    draw(cube, scene);

    curve.transform.translate = {1.9f,0,0};
    curve.transform.rotate = rotation({0,1,0},time);
	draw(curve, scene);
	
	// Display sphere
	vec3 const axis_of_sphere_rotation = { 0,0,1 };
	float const angle_of_sphere_rotation = time / 2.0f;
	rotation const rotation_sphere = rotation(axis_of_sphere_rotation, angle_of_sphere_rotation);
	sphere.transform.rotate = rotation_sphere;
	sphere.transform.translate.z = std::sin(time * 0.5);
	sphere.shading.color = vec3(1 + std::cos(time), 1 + std::sin(time), 2.0) / 2.0f;
	draw(sphere, scene);

	// Display plane to deform
	draw(shape_visual, scene);
	if (user.gui.is_wireframe)
		draw_wireframe(shape_visual, scene, { 0,0,0 });
	
	evolve_shape(time);
	shape_visual.update_position(shape.position);

	// Recompute normals on the CPU (given the position and the connectivity currently in the mesh structure)
	shape.compute_normal();
	// Send updated normals on the GPU
	shape_visual.update_normal(shape.normal);

	// Texture
	// Reset the color of the shape to white (only the texture image will be seen)
	shape_visual.shading.color = { 1,1,1 };

	// Load the image and associate the texture id to the structure
	shape_visual.texture = opengl_texture_to_gpu(image_load_png("assets/squirrel.png"));
}

void evolve_shape(float time)
{
	size_t const N = initial_position.size();
	for (size_t k = 0; k < N; ++k)
	{
		// Wave
		vec3 const& p0 = initial_position[k];
		vec3& p = shape.position[k];
		//p.z = p0.z + 0.1f * std::cos(10 * p.x + 4 * time);

		//Perlin
		float const dz = 0.3f * noise_perlin({ p0.x + 0.2f * time, p0.y, 0 }, 2) + 0.015f * noise_perlin({ 4 * p0.x, 4 * p0.y, time }, 2);
		p = p0 + vec3(0, 0, dz);
	}
}

// Display the GUI
void display_interface()
{
	ImGui::Checkbox("Display frame", &user.gui.display_frame);
	ImGui::SliderFloat("Time Scale", &timer.scale, 0.0f, 2.0f, "%.1f");
	ImGui::Checkbox("Wireframe", &user.gui.is_wireframe);
}

// Function called every time the screen is resized
void window_size_callback(GLFWwindow* , int width, int height)
{	
	glViewport(0, 0, width, height); // The image is displayed on the entire window
	float const aspect = width / static_cast<float>(height); // Aspect ratio of the window

	// Generate the perspective matrix for this aspect ratio
	float const field_of_view = 50.0f*pi/180.0f; // the angle of the field of view
	float const z_near = 0.1f;  // closest visible point
	float const z_far = 100.0f; // furthest visible point
	scene.projection = projection_perspective(field_of_view, aspect, z_near, z_far); 
}

// Function called every time the mouse is moved
void mouse_move_callback(GLFWwindow* window, double xpos, double ypos)
{
	vec2 const  p1 = glfw_get_mouse_cursor(window, xpos, ypos);
	vec2 const& p0 = user.mouse_prev;

	glfw_state state = glfw_current_state(window);
	user.cursor_on_gui = ImGui::IsAnyWindowFocused();

	// Handle camera rotation
	auto& camera = scene.camera;
	if(!user.cursor_on_gui){
		if(state.mouse_click_left && !state.key_ctrl)
			scene.camera.manipulator_rotate_trackball(p0, p1);
		if(state.mouse_click_left && state.key_ctrl)
			camera.manipulator_translate_in_plane(p1-p0);
		if(state.mouse_click_right)
			camera.manipulator_scale_distance_to_center( (p1-p0).y );
	}

	user.mouse_prev = p1;
}

// Uniform data used when displaying an object in this scene
void opengl_uniform(GLuint shader, scene_environment const& current_scene)
{
	opengl_uniform(shader, "projection", current_scene.projection);
	opengl_uniform(shader, "view", scene.camera.matrix_view());
	opengl_uniform(shader, "light", scene.light, false);
}



