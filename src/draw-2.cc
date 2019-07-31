#ifdef USE_PYTHON
#include <Python.h>
#endif // USE_PYTHON

#define GLM_ENABLE_EXPERIMENTAL // # for norm things
#include <glm/ext.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "globjects.h"
#include "trackball.h"
#include "graphics-info.h"

#include "draw.hh"
#include "draw-2.hh"

gint idle_contour_function(gpointer data);

void init_central_cube_shaders() {

   GLenum err = glGetError();
   std::cout << "start init_central_cube_shaders() err is " << err << std::endl;

   std::string  shader_file_name = "central-cube.shader";
   std::cout << "----------- parse and create shader " << shader_file_name << std::endl;
   shader_program_source sps = parse_shader(shader_file_name);
   err = glGetError();
   std::cout << "init_central_cube_shaders() parse_sjaderhader() err is " << err << std::endl;
   unsigned int programID = CreateShader(sps.VertexSource, sps.FragmentSource);
   err = glGetError();
   std::cout << "init_central_cube_shaders() CreateShader() returned programID " << programID
       << " for file " << shader_file_name << " with err " << err << std::endl;
   graphics_info_t::programID_for_central_cube = programID;
   std::cout << "----------- created shader program " << programID << std::endl;

   glBindAttribLocation(programID, 0, "position");
   err = glGetError();
   std::cout << "init_central_cube_shaders() glBindAttribLocation() err is " << err << std::endl;
   if (err == GL_INVALID_VALUE)
      std::cout << "ERROR:: programID " << programID
                << " is not a value generated by OpenGL" << std::endl;

   /* get the location of the "mvp" uniform */ // No! Not here - not in this program.
   // graphics_info_t::central_cube_mvp_location = glGetUniformLocation(programID, "mvp");

   err = glGetError();
   std::cout << "finish init_central_cube_shaders() err is " << err << std::endl;

}

void init_shaders() {
   graphics_info_t::shader_for_maps.init("map.shader", Shader::Entity_t::MAP);
}

void init_central_cube();

void init_buffers() {
   init_central_cube();
}

void init_central_cube() {

   {
      float positions[24] = {
         -0.5,  -0.5, -0.5,
         -0.5,  -0.5,  0.5,
         -0.5,   0.5, -0.5,
         -0.5,   0.5,  0.5,
         0.5,  -0.5, -0.5,
         0.5,  -0.5,  0.5,
         0.5,   0.5, -0.5,
         0.5,   0.5,  0.5
      };

      // number of lines * 2:
      unsigned int indices[24] { 0,1, 1,5, 5,4, 4,0, 2,3, 3,7, 7,6, 6,2, 0,2, 1,3, 5,7, 4,6 };

      GLenum err = glGetError();

      // GLuint VertexArrayID;
      glGenVertexArrays(1, &graphics_info_t::central_cube_vertexarray_id);
      glBindVertexArray(graphics_info_t::central_cube_vertexarray_id);

      // GLuint vertexbuffer;
      glGenBuffers(1, &graphics_info_t::central_cube_array_buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, graphics_info_t::central_cube_array_buffer_id);
      glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 24, &positions[0], GL_STATIC_DRAW);
      glEnableVertexAttribArray(0);
      glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

      // unsigned int ibo;
      glGenBuffers(1, &graphics_info_t::central_cube_index_buffer_id);
      err = glGetError();
      if (err) std::cout << "init_central_cube() index glGenBuffers() err is " << err << std::endl;
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, graphics_info_t::central_cube_index_buffer_id);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * 24, &indices[0], GL_STATIC_DRAW);
      err = glGetError();
      if (err) std::cout << "init_central_cube() glBufferData() err is " << err << std::endl;

   }
}

glm::mat4 get_molecule_mvp() {

   // presumes that we are in the correct programID

   float w = static_cast<float>(graphics_info_t::graphics_x_size);
   float h = static_cast<float>(graphics_info_t::graphics_y_size);
   float screen_ratio = w/h;

   // I don't think that the quaternion belongs to the model matrix, it should be
   // part of the view matrix I think.
   // Yes. That's right.
   glm::mat4 model_matrix = glm::mat4(1.0);

   float z = graphics_info_t::zoom * 0.04;
   glm::vec3 sc(z,z,z);
   float ortho_size = 90.0;

   // start: clipping front and back are 0
   // with large depth of field: clipping_front is -10, clipping_back is -9
   // with narrow depth of field: 5 and 6.

   GLfloat near_scale = 0.3;
   GLfloat far  =      -near_scale*graphics_info_t::zoom * (graphics_info_t::clipping_front*-0.3 + 1.0);
   GLfloat near =  0.30*near_scale*graphics_info_t::zoom * (graphics_info_t::clipping_back* -0.3 + 1.0);

   if (false)
      std::cout << "near: " << near << " far " << far
                << " clipping_front " << graphics_info_t::clipping_front
                << " clipping_back "  << graphics_info_t::clipping_back
                << std::endl;

   glm::mat4 projection_matrix = glm::ortho(-ortho_size * screen_ratio, ortho_size * screen_ratio,
                                            -ortho_size, ortho_size,
                                            near, far); // wrong way round?

   glm::vec3 rc = graphics_info_t::get_rotation_centre();
   // std::cout << "rotation centre " << glm::to_string(rc) << std::endl;
   glm::mat4 view_matrix = glm::toMat4(graphics_info_t::glm_quat);
   view_matrix = glm::scale(view_matrix, sc);
   view_matrix = glm::translate(view_matrix, -rc);

#if 0
   // for fun/testing
   // turn off view scaling when tinkering with this?
   // there should not be a concept of "zoom" with perspective view, just translation
   // along screen-Z.
   float fov = 60.0;
   std::cout << "fov " << fov << std::endl;
   glm::mat4 projection_matrix_persp = glm::perspective(glm::radians(fov), screen_ratio, 2.1f, 1000.0f);
#endif

   glm::mat4 mvp = projection_matrix * view_matrix * model_matrix;

   return mvp;
}

glm::mat4 get_view_rotation() {

   // need to be in the correct program

   glm::mat4 view_matrix = glm::toMat4(graphics_info_t::glm_quat);
   return view_matrix;
}

void draw_map_molecules() {
   glLineWidth(1.0f);
   GLenum err = glGetError();
   if (err) std::cout << "gtk3_draw_molecules() glLineWidth " << err << std::endl;

   GLuint pid = graphics_info_t::shader_for_maps.get_program_id();
   glUseProgram(pid);
   err = glGetError();
   if (err) std::cout << "   gtk3_draw_molecules() glUseProgram with GL err "
                      << err << std::endl;


   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 view_rotation = get_view_rotation(); // hhmm... naming

   glEnable(GL_DEPTH_TEST); // this needs to be in the draw loop!?
   glDepthFunc(GL_LESS);

   for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {
      if (! graphics_info_t::is_valid_map_molecule(ii)) continue;
      if (! graphics_info_t::molecules[ii].draw_it_for_map) continue;
      if (graphics_info_t::molecules[ii].n_vertices_for_VertexArray > 0) {

         bool draw_with_lines = true;
         if (draw_with_lines) {
            glBindVertexArray(graphics_info_t::molecules[ii].m_VertexArrayID_for_map);
            err = glGetError();
            if (err) std::cout << "   draw_map_molecules() glBindVertexArray() "
                               << graphics_info_t::molecules[ii].m_VertexArrayID_for_map
                               << " with GL err " << err << std::endl;

            glBindBuffer(GL_ARRAY_BUFFER,         graphics_info_t::molecules[ii].m_VertexBufferID);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, graphics_info_t::molecules[ii].m_IndexBufferID);

            glUniformMatrix4fv(graphics_info_t::shader_for_maps.mvp_uniform_location,           1, GL_FALSE, &mvp[0][0]);
            err = glGetError();
            if (err) std::cout << "   draw_map_molecules() glUniformMatrix4fv() mvp " << err << std::endl;
            glUniformMatrix4fv(graphics_info_t::shader_for_maps.view_rotation_uniform_location, 1, GL_FALSE, &view_rotation[0][0]);
            err = glGetError();
            if (err) std::cout << "   draw_map_molecules() glUniformMatrix4fv() vr  " << err << std::endl;

            GLuint background_colour_uniform_location = graphics_info_t::shader_for_maps.background_colour_uniform_location;
            glm::vec4 bgc(graphics_info_t::background_colour, 1.0);
            glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
            err = glGetError();
            if (err) std::cout << "   draw_map_molecules() glUniform4fv() for bg  " << err << std::endl;

            glDrawElements(GL_LINES, graphics_info_t::molecules[ii].n_vertices_for_VertexArray,
                           GL_UNSIGNED_INT, nullptr);
            err = glGetError();
            if (err) std::cout << "   draw_map_molecules() glDrawElements() n_vertices: "
                               << graphics_info_t::molecules[ii].n_vertices_for_VertexArray
                               << " with GL err " << err << std::endl;
         }

         if (!draw_with_lines) { // draw as a solid object
            if (true)
               std::cout << "   draw_map_molecules(): imol " << ii
                         << " array_id and n_vertices_for_VertexArray: "
                         << graphics_info_t::molecules[ii].m_VertexArrayID_for_map << " "
                         << graphics_info_t::molecules[ii].n_indices_for_triangles
                         << std::endl;

            glBindVertexArray(graphics_info_t::molecules[ii].m_VertexArrayID_for_map);
            err = glGetError();
            if (err) std::cout << "   draw_map_molecules() glBindVertexArray() "
                               << graphics_info_t::molecules[ii].m_VertexArrayID_for_map
                               << " with GL err " << err << std::endl;

            glBindBuffer(GL_ARRAY_BUFFER,         graphics_info_t::molecules[ii].m_VertexBufferID);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, graphics_info_t::molecules[ii].m_IndexBuffer_for_triangles_ID);

            glUniformMatrix4fv(graphics_info_t::mvp_location, 1, GL_FALSE, &mvp[0][0]);
            err = glGetError();
            if (err) std::cout << "   draw_map_molecules() glUniformMatrix4fv() " << err << std::endl;
            glUniformMatrix4fv(graphics_info_t::view_rotation_location, 1, GL_FALSE, &view_rotation[0][0]);
            err = glGetError();
            if (err) std::cout << "   draw_map_molecules() glUniformMatrix4fv() " << err << std::endl;

            // glDrawElements() uses a vertex count, not n indices, needs checking
            glDrawElements(GL_TRIANGLES, graphics_info_t::molecules[ii].n_indices_for_triangles,
                           GL_UNSIGNED_INT, nullptr);

            err = glGetError();
            if (err) std::cout << "   draw_map_molecules() glDrawElements() n_indices_for_triangles "
                               << graphics_info_t::molecules[ii].n_indices_for_triangles
                               << " with GL err " << err << std::endl;
         }
      }
   }
}

void
draw_model_molecules() {

   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 view_rotation = get_view_rotation(); // hhmm... naming

   for (int ii=graphics_info_t::n_molecules()-1; ii>=0; ii--) {
      if (! graphics_info_t::is_valid_model_molecule(ii)) continue;
      if (! graphics_info_t::molecules[ii].draw_it) continue;

      if (false)
         std::cout << "imol " << ii << " n_vertices_for_model_VertexArray "
                   << graphics_info_t::molecules[ii].n_vertices_for_model_VertexArray << std::endl;
      if (graphics_info_t::molecules[ii].n_vertices_for_model_VertexArray > 0) {
         // OOps - every model has its own shader - that's a mistake
         GLuint pid = graphics_info_t::molecules[ii].shader.get_program_id();
         glUseProgram(pid);
         GLuint err = glGetError(); if (err) std::cout << "   error draw_model_molecules() glUseProgram() "
                                                       << err << std::endl;

         glBindVertexArray(graphics_info_t::molecules[ii].m_VertexArray_for_model_ID);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glBindVertexArray() "
                               << graphics_info_t::molecules[ii].m_VertexArray_for_model_ID
                               << " with GL err " << err << std::endl;

         glBindBuffer(GL_ARRAY_BUFFER,         graphics_info_t::molecules[ii].m_VertexBuffer_for_model_ID);
         err = glGetError(); if (err) std::cout << "   error draw_model_molecules() glBindBuffer() v " << err << std::endl;
         glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, graphics_info_t::molecules[ii].m_IndexBuffer_for_model_ID);
         err = glGetError(); if (err) std::cout << "   error draw_model_molecules() glBindBuffer() i " << err << std::endl;

         GLuint mvp_location = graphics_info_t::molecules[ii].shader.mvp_uniform_location;
         GLuint view_rotation_location = graphics_info_t::molecules[ii].shader.view_rotation_uniform_location;
         // std::cout << "  draw_model_molecules() locations " << mvp_location << " " << view_rotation_location << std::endl;

         glUniformMatrix4fv(mvp_location, 1, GL_FALSE, &mvp[0][0]);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniformMatrix4fv() " << err << std::endl;
         glUniformMatrix4fv(view_rotation_location, 1, GL_FALSE, &view_rotation[0][0]);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniformMatrix4fv() for mvp " << err << std::endl;

         GLuint background_colour_uniform_location = graphics_info_t::molecules[ii].shader.background_colour_uniform_location;
         glm::vec4 bgc(graphics_info_t::background_colour, 1.0);
         glUniform4fv(background_colour_uniform_location, 1, glm::value_ptr(bgc));
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glUniform4fv() for background " << err << std::endl;

         // draw with the vertex count, not the index count.
         GLuint n_verts = graphics_info_t::molecules[ii].n_indices_for_model_triangles;
         //std::cout << "   Drawing " << n_verts << " model vertices" << std::endl;
         glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
         err = glGetError();
         if (err) std::cout << "   error draw_model_molecules() glDrawElements() "
                            << n_verts << " with GL err " << err << std::endl;

      }
   }
}

void
draw_molecules() {

   draw_map_molecules();
   draw_model_molecules();

}

void
draw_central_cube(GtkGLArea *glarea) {

   gtk_gl_area_make_current(glarea);
   glLineWidth(2.0);  // GLv4 antialiasing - OpenGL implementations are not required to support this
   GLenum err = glGetError();
   if (err) std::cout << "draw_triangle() A err " << err << std::endl;

   // To see the possible values of the line width in aliased mode:
   // GLfloat line_width_max_min[2] = {0.0f, 0.0f};
   // glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, lineWidthRange);
   // This may not be possible in GL_LINE_SMOOTH mode.

   if (true) {


      // Note to self: use the same shader as for models. You will need
      // to do move the cube to the rotattion centre.

      GtkAllocation allocation;
      gtk_widget_get_allocation(GTK_WIDGET(glarea), &allocation);
      int w = allocation.width;
      int h = allocation.height;
      float aspect_ratio = static_cast<float>(h)/static_cast<float>(w);

      glBindVertexArray(graphics_info_t::central_cube_vertexarray_id);
      err = glGetError();
      if (err) std::cout << "draw_central_cube() B err " << err << std::endl;
      glUseProgram(graphics_info_t::programID_for_central_cube);
      err = glGetError();
      if (err) std::cout << "draw_central_cube() C err " << err << std::endl;

      err = glGetError();
      if (err) std::cout << "draw_central_cube() D err " << err << std::endl;

      glm::mat4 view_orientation = glm::toMat4(graphics_info_t::glm_quat);
      float z = graphics_info_t::zoom * 0.0002;
      glm::vec3 sc(z,z,z);
      // std::cout << "z " << z << std::endl;
      // glm::vec3 sc(0.2f, 0.2f, 0.2f);
      glm::mat4 mvp  = glm::scale(view_orientation, sc);
      // glm::mat4 mvp = glm::scale(m1, glm::vec3(4.0*aspect_ratio, 1.0, 1.0));

      if (false) {
         std::cout << "debug:: draw_central_cube()       local mvp: ";
         for (unsigned int i=0; i<4; i++)
            for (unsigned int j=0; j<4; j++)
               std::cout << std::setw(8) << mvp[i][j] << " ";
         std::cout << std::endl;
      }
      glUniformMatrix4fv(graphics_info_t::mvp_location, 1, GL_FALSE, &mvp[0][0]);
      err = glGetError();
      if (err) std::cout << "draw_central_cube() glUniformMatrix4fv() " << err << std::endl;

      // glBindBuffer(GL_ARRAY_BUFFER, graphics_info_t::central_cube_array_buffer_id);
      // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, graphics_info_t::central_cube_index_buffer_id);
      glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, nullptr);
      err = glGetError();
      if (err) std::cout << "draw_central_cube() F glDrawElements() " << err << std::endl;
   }

   glBindVertexArray(0); // unbind
   glUseProgram(0);

}

GtkWidget *my_gtkglarea(GtkWidget *vbox) {

   GtkWidget *w = gtk_gl_area_new();
   gtk_widget_set_size_request(w, 900, 900);
   gtk_box_pack_start(GTK_BOX(vbox), w, TRUE, TRUE, 2);
   return w;
}

void
on_glarea_realize(GtkGLArea *glarea) {

   std::cout << "realize!" << std::endl;

   gtk_gl_area_make_current(glarea);
   gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(glarea), TRUE);

   GLenum err = glGetError();
   std::cout << "start on_glarea_realize() err is " << err << std::endl;

   init_shaders();
   err = glGetError();
   std::cout << "on_glarea_realize() post init_shaders() err is " << err << std::endl;

   init_central_cube_shaders();
   err = glGetError();
   std::cout << "on_glarea_realize() post init_central_cube_shaders() err is " << err << std::endl;

   init_buffers();
   err = glGetError();
   std::cout << "on_glarea_realize() post init_buffer(): err is " << err << std::endl;

   gtk_gl_area_set_has_depth_buffer(GTK_GL_AREA(glarea), TRUE);

   glEnable(GL_DEPTH_TEST);
   glDepthFunc(GL_GREATER); // what does this do?
   // glDisable(GL_DEPTH_TEST);
   // glDepthFunc(GL_ALWAYS);

   // Make antialised lines
   if (false) {
      glEnable (GL_BLEND);
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_LINE_SMOOTH);
   }

#if 0
   transform = Transform(glm::vec3(0.0, 0.0, 0.0),   // position
       glm::vec3(0.0, 0.0, 0.0),   // rotation
       glm::vec3(1.0, 1.0, 1.0));  // scales
#endif

}


gboolean
on_glarea_render(GtkGLArea *glarea) {

   auto tp_0 = std::chrono::high_resolution_clock::now();

   // is this needed?
   gtk_gl_area_make_current(glarea);
   GLenum err = glGetError();
   if (err) std::cout << "on_glarea_render() start " << err << std::endl;

   glClearColor (0.24, 0.24, 0.24, 1.0);
   const glm::vec3 &bg = graphics_info_t::background_colour;
   glClearColor (bg[0], bg[1], bg[2], 1.0);
   err = glGetError();
   if (err) std::cout << "on_glarea_render B err " << err << std::endl;
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   err = glGetError();
   if (err) std::cout << "on_glarea_render C err " << err << std::endl;

   //

   // get rid of this?
   // if (graphics_info_t::draw_the_other_things)
   // draw_other_triangle(glarea);


   // Has the start triangle been correctly init? It errors on draw now.
   // draw_triangle(glarea);

   draw_central_cube(glarea);
   draw_molecules();

   err = glGetError();
   if (err) std::cout << "on_glarea_render gtk3_draw_molecules() " << err << std::endl;

   glFlush ();
   err = glGetError();
   if (err) std::cout << "on_glarea_render E err " << err << std::endl;

   graphics_info_t::frame_counter++;
   // auto tp_1 = std::chrono::high_resolution_clock::now();
   // auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
   // std::cout << "INFO:: Timing for frame render " << d10 << " microseconds" << std::endl;

  return FALSE;
}


void
on_glarea_resize(GtkGLArea *glarea, gint width, gint height) {

   std::cout << "resize!" << std::endl;
   graphics_info_t g;
   g.graphics_x_size = width;
   g.graphics_y_size = height;

}

gboolean
on_glarea_scroll(GtkWidget *widget, GdkEventScroll *event) {

   int direction = 1;
   if (event->direction == GDK_SCROLL_UP)
      direction = -1;

   std::cout << "scroll " << direction << std::endl;

   graphics_info_t g;
   int imol_scroll = graphics_info_t::scroll_wheel_map;

   if (g.is_valid_map_molecule(imol_scroll)) {
      // use direction
      if (direction == 1)
         graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count--;
      if (direction == -1)
         graphics_info_t::molecules[imol_scroll].pending_contour_level_change_count++;
      int contour_idle_token = g_idle_add(idle_contour_function, g.glarea);
      std::cout << "####### Now contour level for map " << imol_scroll << "is "
                << g.molecules[imol_scroll].contour_level << std::endl;
      g.set_density_level_string(imol_scroll, g.molecules[imol_scroll].contour_level);
      g.display_density_level_this_image = 1;
      g.update_maps();
      gtk_widget_queue_draw(widget);
   } else {
      std::cout << "No map" << std::endl;
   }
   return TRUE;
}

gboolean
on_glarea_button_press(GtkWidget *widget, GdkEventButton *event) {

   std::cout << "button press!" << std::endl;
   graphics_info_t g;
   g.SetMouseBegin(event->x,event->y);
   return TRUE;
}

gboolean
on_glarea_button_release(GtkWidget *widget, GdkEventButton *event) {

   std::cout << "button release!" << std::endl;
   return TRUE;
}

gboolean
on_glarea_motion_notify(GtkWidget *widget, GdkEventMotion *event) {

   int r = 0;
   graphics_info_t g;

   // split this function up before it gets too big.

   g.mouse_current_x = event->x;
   g.mouse_current_y = event->y;

   if (event->state & GDK_BUTTON1_MASK) {

      GtkAllocation allocation;
      gtk_widget_get_allocation(widget, &allocation);
      int w = allocation.width;
      int h = allocation.height;

      glm::quat tb_quat =
         g.trackball_to_quaternion((2.0*g.GetMouseBeginX() - w)/w, (h - 2.0*g.GetMouseBeginY())/h,
           (2.0*g.mouse_current_x - w)/w,  (h - 2.0*g.mouse_current_y)/h,
           g.get_trackball_size());

      glm::mat4 mat_from_quat = glm::toMat4(tb_quat);

      glm::quat product = tb_quat * graphics_info_t::glm_quat;
      graphics_info_t::glm_quat = glm::normalize(product);
   }


   if (event->state & GDK_BUTTON2_MASK) {

      // View Panning

      GtkAllocation allocation;
      gtk_widget_get_allocation(widget, &allocation);
      int w = allocation.width;
      int h = allocation.height;

      glm::mat4 mvp = get_molecule_mvp(); // modeglml matrix includes orientation with the quaternion


      float mouseX_1 = g.GetMouseBeginX() / (w * 0.5f) - 1.0f;
      float mouseY_1 = g.GetMouseBeginY() / (h * 0.5f) - 1.0f;
      float mouseX_2 = g.mouse_current_x  / (w * 0.5f) - 1.0f;
      float mouseY_2 = g.mouse_current_y  / (h * 0.5f) - 1.0f;

      glm::mat4 vp_inv = glm::inverse(mvp);

      glm::vec4 screenPos_1 = glm::vec4(mouseX_1, -mouseY_1, 1.0f, 1.0f);
      glm::vec4 screenPos_2 = glm::vec4(mouseX_2, -mouseY_2, 1.0f, 1.0f);
      glm::vec4 worldPos_1 = vp_inv * screenPos_1;
      glm::vec4 worldPos_2 = vp_inv * screenPos_2;

      glm::vec4 delta(worldPos_1 - worldPos_2);

      g.add_to_rotation_centre(delta);
      g.update_maps();
      int contour_idle_token = g_idle_add(idle_contour_function, g.glarea);

   }

   if (event->state & GDK_BUTTON3_MASK) {

      // Zooming

      double delta_x = event->x - g.GetMouseBeginX();
      double delta_y = event->y - g.GetMouseBeginY();
      double fx = 1.0f +  delta_x/300.0;
      double fy = 1.0f +  delta_y/300.0;
      if (fx > 0.0) g.zoom *= fx;
      if (fy > 0.0) g.zoom *= fy;
      // std::cout << "now zoom: " << g.zoom << std::endl;
   }

   // for next motion
   g.SetMouseBegin(event->x,event->y);
   gtk_widget_queue_draw(widget);
   return TRUE;
}

gint
spin_func(gpointer data) {

   float delta = 0.02;
   glm::vec3 EulerAngles(0, delta, 0);
   glm::quat quat_delta(EulerAngles);
   glm::quat normalized_quat_delta(glm::normalize(quat_delta));
   glm::quat product = normalized_quat_delta * graphics_info_t::glm_quat;
   graphics_info_t::glm_quat = glm::normalize(product);
   gtk_widget_queue_draw(graphics_info_t::glarea);

   std::chrono::time_point<std::chrono::system_clock> tp_now = std::chrono::high_resolution_clock::now();
   std::chrono::time_point<std::chrono::system_clock> tp_nowa = std::chrono::high_resolution_clock::now();
   std::chrono::time_point<std::chrono::system_clock> tp_nowb = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed_seconds = tp_now - graphics_info_t::previous_frame_time;
   std::chrono::duration<double> elapsed_secondsA = tp_now - graphics_info_t::previous_frame_time;
   std::chrono::duration<double> elapsed_secondsB = tp_now - graphics_info_t::previous_frame_time;
   if (elapsed_seconds.count() > 1.0) {
      float nf = graphics_info_t::frame_counter - graphics_info_t::frame_counter_at_last_display;
      std::cout << "Frame/second: " << 1000 * elapsed_seconds.count()/nf << " milliseconds\n";
      graphics_info_t::previous_frame_time = tp_now;
      graphics_info_t::frame_counter_at_last_display = graphics_info_t::frame_counter;
   }

   // kludge/race condition?
   if (graphics_info_t::idle_function_spin_rock_token == -1)
      return FALSE;
   else
      return TRUE;
}

gboolean
on_glarea_key_press_notify(GtkWidget *widget, GdkEventKey *event) {

   std::cout << "on_glarea_key_press_notify() " << std::endl;
   graphics_info_t g;
   gboolean handled = false;

   if (event->keyval == GDK_KEY_n) {
      std::cout << "Zoom in " << std::endl;
      graphics_info_t::zoom *= 0.9;
   }
   if (event->keyval == GDK_KEY_m) {
      std::cout << "Zoom out " << std::endl;
      graphics_info_t::zoom *= 1.1;
   }

   // think about the more generic adjust_clipping()
   if (event->keyval == GDK_KEY_d) {
      adjust_clipping(1.0);
   }

   if (event->keyval == GDK_KEY_f) {
      adjust_clipping(-1.0);
   }

   if (event->keyval == GDK_KEY_i) {
      std::cout << "Debug idle_function_spin_rock_token " << graphics_info_t::idle_function_spin_rock_token
                << std::endl;
      if (graphics_info_t::idle_function_spin_rock_token != -1) {
         std::cout << "Removing the idle function\n";
         g_idle_remove_by_data(GINT_TO_POINTER(66)); // just a kludge for the moment
         graphics_info_t::idle_function_spin_rock_token = -1;
      } else {
         int toi = g_timeout_add(5, spin_func, GINT_TO_POINTER(66));
         graphics_info_t::idle_function_spin_rock_token = toi;
      }
   }

   if (event->keyval == GDK_KEY_minus || event->keyval == GDK_KEY_plus) {
      int s = graphics_info_t::scroll_wheel_map;
      if (graphics_info_t::is_valid_map_molecule(s)) {
         if (event->keyval == GDK_KEY_minus)
            graphics_info_t::molecules[s].pending_contour_level_change_count--;
         if (event->keyval == GDK_KEY_plus)
            graphics_info_t::molecules[s].pending_contour_level_change_count++;
         int contour_idle_token = g_idle_add(idle_contour_function, g.glarea);
         g.set_density_level_string(s, g.molecules[s].contour_level);
         g.display_density_level_this_image = 1;
      }
   }

   gtk_widget_queue_draw(widget);

   return handled;

}

gboolean
on_glarea_key_release_notify(GtkWidget *widget, GdkEventKey *event) {

   return TRUE;
}

void my_glarea_add_signals_and_events(GtkWidget *glarea) {

   gtk_widget_add_events(glarea, GDK_SCROLL_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON_PRESS_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON_RELEASE_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON1_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON2_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON3_MOTION_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON1_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON2_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON3_MASK);
   gtk_widget_add_events(glarea, GDK_KEY_PRESS_MASK);

   // key presses for the glarea:
   gtk_widget_set_can_focus(glarea, TRUE);
   gtk_widget_grab_focus(glarea);

   g_signal_connect(glarea, "realize", G_CALLBACK(on_glarea_realize), NULL);
   g_signal_connect(glarea, "render",  G_CALLBACK(on_glarea_render),  NULL);
   g_signal_connect(glarea, "resize",  G_CALLBACK(on_glarea_resize),  NULL);
   g_signal_connect(glarea, "scroll-event",          G_CALLBACK(on_glarea_scroll),             NULL);
   g_signal_connect(glarea, "button-press-event",    G_CALLBACK(on_glarea_button_press),       NULL);
   g_signal_connect(glarea, "button-release-event",  G_CALLBACK(on_glarea_button_release),     NULL);
   g_signal_connect(glarea, "motion-notify-event",   G_CALLBACK(on_glarea_motion_notify),      NULL);
   g_signal_connect(glarea, "key-press-event",       G_CALLBACK(on_glarea_key_press_notify),   NULL);
   g_signal_connect(glarea, "key-release-event",     G_CALLBACK(on_glarea_key_release_notify), NULL);

}
