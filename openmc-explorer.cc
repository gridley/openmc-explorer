#include <SDL2/SDL.h>

#include <openmc/cell.h>
#include <openmc/constants.h>
#include <openmc/error.h>
#include <openmc/geometry.h>
#include <openmc/material.h>
#include <openmc/math_functions.h>
#include <openmc/particle.h>
#include <openmc/plot.h>
#include <openmc/position.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

/* TODO
  - Use timers to steady the scrolling rate

  - Define an abstract base class that defines how MC code
    can interact with this program.

  - Make the distance_to_boundary_from_void work faster. Just need
    to build the list of surfaces that bound the root universe, rather
    than looping over cells in the root universe.

  - Make camera move in alignment with the current direction, not
    just with xyz when pressing wasd.

  - Add capability to press F5 to view current position, orientation, etc.

  - Add atomic operations where needed to remove weird fuzzy effects when
    using OpenMP

  - Currently causes a segfault if openmc loses a particle. This happens
    when particles perfectly intersect the corner of a Cartesian lattice
    at 45 deg to the corner. If the camera is slightly adjusted, the error
    does not happen.
*/

std::array<unsigned char, 3> hsv2rgb(float h, float s, float l) {
  std::array<unsigned char,3> result;
  std::array<float,3> resultf;
  float c = (1. - abs(2.*l - 1.)) * s;
  int hprime = h / 60.;
  float x = c * (1.-abs(fmod(hprime,2)-1));
  switch(hprime) {
    case 0:
      resultf[0]=c;
      resultf[1]=x;
      resultf[2]=0;
      break;
    case 1:
      resultf[0]=x;
      resultf[1]=c;
      resultf[2]=0;
      break;
    case 2:
      resultf[0]=0;
      resultf[1]=c;
      resultf[2]=x;
      break;
    case 3:
      resultf[0]=0;
      resultf[1]=x;
      resultf[2]=c;
      break;
    case 4:
      resultf[0]=x;
      resultf[1]=0;
      resultf[2]=c;
      break;
    case 5:
      resultf[0]=c;
      resultf[1]=0;
      resultf[2]=x;
      break;
  }
  float m = l - c/2;
  result[0] = static_cast<unsigned char>(255*(resultf[0]+m));
  result[1] = static_cast<unsigned char>(255*(resultf[1]+m));
  result[2] = static_cast<unsigned char>(255*(resultf[2]+m));
  return result;
}

// Simple struct that holds settings for the program, and has default values
struct ExplorerSettings
{
  static const int default_win_size_x = 800;
  static const int default_win_size_y = 600;

  bool running;
  bool use_fullscreen;
  double delta_angle; // degrees to rotate camera on keypress
  double delta_scoot; // distance to scoot cammera (cm)

  // How many lines to draw between frames. This dynamically adjusts to keep
  // a desirable framerate
  unsigned n_lines_to_draw;
  
  int window_size_x;
  int window_size_y;
  unsigned n_threads;

  ExplorerSettings();
};
ExplorerSettings::ExplorerSettings() : running(true), delta_angle(5.0), delta_scoot(1.0),
  n_lines_to_draw(5), window_size_x(ExplorerSettings::default_win_size_x),
  window_size_y(ExplorerSettings::default_win_size_y), n_threads(1),
  use_fullscreen(false) 
{
#ifdef _OPENMP
  n_threads = omp_get_max_threads();
  omp_set_num_threads(n_threads);
  n_lines_to_draw = n_threads;
  std::cout << "Ray tracing with " << n_threads << " threads." << std::endl;
#endif
}

class WindowManager
{
  int window_width, window_height;
  ExplorerSettings& settings;
  public:
    SDL_Window* window;

    WindowManager(ExplorerSettings& sett);
    ~WindowManager();
    void toggleFullscreen();
};

void WindowManager::toggleFullscreen() {
  int flag;
  if (settings.use_fullscreen) {
    settings.use_fullscreen = false;
    settings.window_size_x = ExplorerSettings::default_win_size_x;
    settings.window_size_y = ExplorerSettings::default_win_size_y;
    flag = SDL_SetWindowFullscreen(window, 0);
    if (flag) openmc::fatal_error("Unable to remove fullscreen.\n");
  }
  else {
    settings.use_fullscreen = true;
    flag = SDL_SetWindowFullscreen(window, SDL_WINDOW_FULLSCREEN);
    if (flag) openmc::fatal_error("Unable to go full-screen.\n");
    SDL_GetWindowSize(window, &settings.window_size_x, &settings.window_size_y);
  }
}
WindowManager::WindowManager(ExplorerSettings& sett) : window_width(sett.window_size_x),
  window_height(sett.window_size_y),
  window(nullptr),
  settings(sett)
{
  SDL_Init(SDL_INIT_VIDEO);
  window = SDL_CreateWindow("OpenMC Explorer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
      window_width, window_height, SDL_WINDOW_SHOWN);

  if (window == nullptr)
    openmc::fatal_error("Failed to create window instance.\n");

  if (sett.use_fullscreen) toggleFullscreen();
}
WindowManager::~WindowManager()
{
  SDL_DestroyWindow(window);
  SDL_Quit();
}

// Simple class to do RAII for SDL renderers
class Renderer
{
  SDL_Renderer* renderer_;
  public:
    SDL_Renderer* rend() { return renderer_; }
    Renderer(WindowManager& win);
    ~Renderer();
    void present();
};
Renderer::Renderer(WindowManager& win) : renderer_(nullptr) {
  renderer_ = SDL_CreateRenderer(win.window, -1, SDL_RENDERER_ACCELERATED || SDL_HINT_RENDER_VSYNC);
  if (renderer_ == nullptr)
  {
    openmc::fatal_error("unable to make renderer instance: \n"
        + std::string(SDL_GetError()));
  }
  // Initialize black screen
  SDL_RenderClear(renderer_);
  SDL_SetRenderDrawColor(renderer_, 0, 0, 0, 0);
  SDL_RenderClear(renderer_);
}
void Renderer::present() { SDL_RenderPresent(renderer_); }
Renderer::~Renderer() { SDL_DestroyRenderer(renderer_); }

// Track a particle from outside the OpenMC geometry up to the
// boundary. Returns true if an OpenMC boundary is intersected.
bool distance_to_boundary_from_void(openmc::Particle& p) {
  // Check intersections with surfaces all surfaces used in universe 0
  constexpr double scoot = 1e-5; // fairly big scoot, because whatever.
  double min_dist {INFINITY};

  // The root level coordinates are always being used in this case
  auto coord = p.coord_[0];
  openmc::Universe* uni = openmc::model::universes[openmc::model::root_universe].get();
  for (auto c_i : uni->cells_) {
    auto dist = openmc::model::cells.at(c_i)->distance(coord.r, coord.u, 0, &p);
    if (dist.first < min_dist) min_dist = dist.first;
  }

  if (min_dist > 1e300) return false;
  else {
    // advance
    // std::cout << min_dist << std::endl;
    for (int j = 0; j < p.n_coord_; ++j) {
      p.coord_[j].r += (min_dist+scoot) * p.coord_[j].u;
    }
    return true;
  }
}

class Camera
{
  openmc::Position pos;
  openmc::Direction dir;

  // whether to update the rendering
  bool viewChanged_;

  // Amount of roll to the camera (radians)
  double roll;

  // Field of view. Initially set to values close to the human eye (radians)
  double fov_x;
  double fov_y;

  // current view angle
  double phi;
  double mu;

  // How transpartent the geometry is
  double opacity;

  // Compute direction cosines
  void computeDir();

  // Colors of the materials. Randomized...
  std::vector<openmc::RGBColor> colors;

  bool isRendering_;
  
  unsigned line_index;

  // Needs a window manager object to know how many rays to create
  WindowManager& window;

  public:

    // Initialize at (1,1,1) pointing at origin
    Camera(WindowManager& win);

    // Phi and mu are given in degrees here:
    Camera(double x, double y, double z, double phi, double mu, WindowManager& win);

    // Start the rendering process
    void beginRendering();
    bool isRendering();
    void renderlines(Renderer& rend, ExplorerSettings& sett);

    void toggleFullscreen() {
      window.toggleFullscreen();
      viewChanged_ = true;
    }

    void rotatePhi(double angle_d);
    void rotateMu(double angle_d);
    void increment_roll(double angle_d);
    void increment_fov_horiz(double angle_d);
    void increment_fov_vert(double angle_d);
    void increment_opacity(double d_opac);

    void randomizeColors();

    // Scoot in a Cartesian direction
    template<unsigned Dir>
    void scoot(double distance) {
      pos[Dir] += distance;
      viewChanged_ = true;
    }

    // Check whether view has been changed (resets internal state after checking)
    bool viewChanged();
};
bool Camera::isRendering() { return isRendering_; }

void Camera::increment_opacity(double d_opac) {
  opacity += d_opac;
  if (opacity < 0) opacity = .00001;
  viewChanged_ = true;
}
void Camera::computeDir() {
  openmc::Direction result;
  result[0] = cos(phi) * sin(mu);
  result[1] = sin(phi) * sin(mu);
  result[2] = cos(mu);
  dir = result;
}

void Camera::randomizeColors() {
  for (unsigned n=0; n<openmc::model::materials.size(); ++n) {
    colors[n] = openmc::random_color();
  }
  viewChanged_ = true;
}

Camera::Camera(double x, double y, double z, double phid, double mud,
    WindowManager& win) : pos(x, y, z),
  fov_x(135.0*M_PI/180.),
  fov_y(75.0*M_PI/180.),
  phi(phid * M_PI/180.0),
  mu(mud * M_PI/180.0),
  roll(0.0),
  window(win),
  viewChanged_(true),
  opacity(.001)
{
  computeDir();

  unsigned n_materials = openmc::model::materials.size();
  if (n_materials == 0) openmc::fatal_error("No materials found in problem!");

  // Create random colors for materials (need to make sure these have lightness=0.5!
  for (unsigned m=0; m<n_materials; ++m) {
    colors.push_back(openmc::random_color());
  }
}
Camera::Camera(WindowManager& win) : Camera(-500, -501, 0.0, 45.0, 90.0, win) {}

void Camera::rotatePhi(double angle_d) {
  double angle_r = angle_d * M_PI/180.0;
  phi += angle_r;

  if (phi < 0.0) phi += 2.*M_PI;
  else if (phi > 2.*M_PI) phi -= 2.*M_PI;

  computeDir();
  viewChanged_ = true;
}

void Camera::rotateMu(double angle_d) {
  double phi_zero = 0.0;
  double angle_r = angle_d * M_PI/180.;
  mu += angle_r;

  if (mu < 0.0) mu = 0.;
  else if (mu > M_PI) mu = M_PI;

  computeDir();
  viewChanged_ = true;
}
bool Camera::viewChanged() {
  bool result = viewChanged_;
  return result;
}
void Camera::increment_roll(double angle_d) {
  roll += angle_d * M_PI/180.;
  viewChanged_ = true;
}
void Camera::increment_fov_horiz(double angle_d) {
  fov_x += angle_d * M_PI/180.;
  viewChanged_ = true;
}
void Camera::increment_fov_vert(double angle_d) {
  fov_y += angle_d * M_PI/180.;
  viewChanged_ = true;
}
void Camera::beginRendering() {
  isRendering_ = true;
  line_index = 0;
}
void Camera::renderlines(Renderer& rend, ExplorerSettings& sett) {
  constexpr unsigned max_intersections = 1000000;
  const double half_mu = fov_y/2.;
  const double half_phi = fov_x/2.;
  const double d_mu = fov_y / sett.window_size_y;
  const double d_phi = fov_x / sett.window_size_x;

  #pragma omp parallel for
  for (unsigned n=0; n<sett.n_threads; ++n) {

    openmc::Direction thisdir; // thread local direction
    openmc::Particle part;
    openmc::Particle::Bank thisbank;
    thisbank.E = 1;
    thisbank.wgt = 1;
    thisbank.delayed_group = 0;
    thisbank.particle = openmc::Particle::Type::photon; //why not?
    thisbank.parent_id = 1;
    thisbank.progeny_id = 2;
    thisbank.r = pos;

    // Loop over lines per thread to draw
    for (unsigned l=0; l<sett.n_lines_to_draw; ++l) {

      int vert_indx = n * sett.n_lines_to_draw + l + line_index;
      if (vert_indx >= sett.window_size_y) {
        isRendering_ = false;
        viewChanged_ = false;
        break;
      }

      // Loop over horizontal pixels
      // Need to track distance through each material
      std::vector<double> materials_tracks(colors.size(), 0.0);
      for (unsigned horiz_indx=0; horiz_indx<sett.window_size_x; horiz_indx++) {

        // Reset material tracks to zero
        for (int mi=0; mi<colors.size(); ++mi)
          materials_tracks[mi] = 0.0;

        // Roll is not currently handled... shouldn't be too hard to add?
        double phi_t = d_phi*horiz_indx - half_phi + phi;
        double mu_t = half_mu - d_mu*vert_indx + mu;
        thisdir[0] = cos(phi_t) * sin(mu_t);
        thisdir[1] = sin(phi_t) * sin(mu_t);
        thisdir[2] = cos(mu_t);

        // Trace ray through geometry
        thisbank.u = thisdir;
        part.from_source(&thisbank);

        #ifdef DEBUG
        std::cout << "New particle" << std::endl;
        #endif

        bool hitsomething = false;
        bool intersection_found = true;
        unsigned n_loops = 0;
        while (intersection_found) {
          bool inside_cell = openmc::find_cell(&part, false);
          if (inside_cell) {
            // Inside geometry, so normal tracing routines may be used
            hitsomething = true;
            intersection_found = true;
            auto dist = openmc::distance_to_boundary(&part);
            materials_tracks[part.material_] += dist.distance;

            // Advance particle
            for (int lev=0; lev<part.n_coord_; ++lev) {
              part.coord_[lev].r += dist.distance * part.coord_[lev].u;
            }

            #ifdef DEBUG
            std::cout << part.coord_[0].r << std::endl;
            #endif

            part.surface_ = dist.surface_index;
            part.n_coord_ = dist.coord_level;
            part.n_coord_last_ = part.n_coord_;
            if (dist.lattice_translation[0] != 0 ||
                dist.lattice_translation[1] != 0 ||
                dist.lattice_translation[2] != 0) {
              // Particle crosses lattice boundary
              cross_lattice(&part, dist);
            }

          } else {
            intersection_found = distance_to_boundary_from_void(part);
          }
          n_loops++;

          if (n_loops > max_intersections) openmc::fatal_error("Infinite loop in tracking!\n");
        }

        // TODO the color calculation here is a bit ad hoc. It just
        // linearly mixes colors based on pathlength through material,
        // which is clearly bad. More consideration should be given to
        // which material had been passed through first. Unable to come
        // up with a good way to handle this at the moment.

        // Convert track lengths to attenuation values
        double overall_product = 1.0;
        for (int mi=0; mi<colors.size(); ++mi) {
          materials_tracks[mi] = exp(-opacity*materials_tracks[mi]);
          overall_product *= materials_tracks[mi];
        }

        std::array<double, 3> final_color;
        for (int ci=0; ci<3; ++ci) final_color[ci] = 0.0; // zero out stuff

        for (int mi=0; mi<colors.size(); ++mi) { // materials
          // This operation just feels right. It maintains
          // the expected attenuation of color in limiting
          // cases, and maintains boundedness always. No indexing
          // into openmc::RGBColor right now so this is manually unrolled.
          final_color[0] += (1.0-materials_tracks[mi]) * colors[mi].red *
            overall_product / materials_tracks[mi];
          final_color[1] += (1.0-materials_tracks[mi]) * colors[mi].green *
            overall_product / materials_tracks[mi];
          final_color[2] += (1.0-materials_tracks[mi]) * colors[mi].blue *
            overall_product / materials_tracks[mi];
        }

        #pragma omp critical
        {
          SDL_SetRenderDrawColor(rend.rend(),
              static_cast<unsigned char>(final_color[0]),
              static_cast<unsigned char>(final_color[1]),
              static_cast<unsigned char>(final_color[2]),
              255);
          SDL_RenderDrawPoint(rend.rend(), horiz_indx, vert_indx);
        }
      }
    }
  }
  line_index += sett.n_threads * sett.n_lines_to_draw;
}

void handle_events(Camera& cam, ExplorerSettings& sett) {
  SDL_Event event;
  // Process keyboard input
  while (SDL_PollEvent(&event)) {
    switch (event.type) {
      case SDL_QUIT:
        sett.running=false;
        break;
      case SDL_KEYDOWN:
        switch (event.key.keysym.sym) {

          // Camera rotation
          case SDLK_UP:
            cam.rotateMu(sett.delta_angle);
            break;
          case SDLK_DOWN:
            cam.rotateMu(-sett.delta_angle);
            break;
          case SDLK_LEFT:
            cam.rotatePhi(-sett.delta_angle);
            break;
          case SDLK_RIGHT:
            cam.rotatePhi(sett.delta_angle);
            break;

          // Camera movement
          case SDLK_w:
            cam.scoot<0>(sett.delta_scoot);
            break;
          case SDLK_s:
            cam.scoot<0>(-sett.delta_scoot);
            break;
          case SDLK_a:
            cam.scoot<1>(-sett.delta_scoot);
            break;
          case SDLK_d:
            cam.scoot<1>(sett.delta_scoot);
            break;
          case SDLK_LSHIFT:
            cam.scoot<2>(-sett.delta_scoot);
            break;
          case SDLK_SPACE:
            cam.scoot<2>(sett.delta_scoot);
            break;

          // Camera roll
          case SDLK_e:
            cam.increment_roll(sett.delta_angle);
            break;
          case SDLK_q:
            cam.increment_roll(-sett.delta_angle);
            break;

          // Field of view adjustment
          case SDLK_r:
            cam.increment_fov_horiz(sett.delta_angle);
            break;
          case SDLK_f:
            cam.increment_fov_horiz(-sett.delta_angle);
            break;
          case SDLK_t:
            cam.increment_fov_vert(sett.delta_angle);
            break;
          case SDLK_g:
            cam.increment_fov_vert(-sett.delta_angle);
            break;

          // Material opacity adjustment
          case SDLK_j:
            cam.increment_opacity(-.0001);
            break;
          case SDLK_k:
            cam.increment_opacity(.0001);
            break;

          case SDLK_F11:
            cam.toggleFullscreen();
            break;

          case SDLK_c:
            cam.randomizeColors();
            break;
        }
        break;
    }
  }
}

int main(int argc, char** argv)
{
  ExplorerSettings settings;
  WindowManager window(settings);
  Renderer renderer(window);

  // Set up openmc
  int err;
  err = openmc_init(argc, argv, nullptr);
  if (err == -1) { return 0; }
  else if (err) { openmc::fatal_error(openmc_err_msg); }

  // Camera must be constructed after openmc_init
  Camera cam(window);
  cam.beginRendering();

  // Get geometry bounding box, and move camera so all is initially visible
  auto bBox = openmc::model::universes[openmc::model::root_universe]->bounding_box();
  std::cout << "Model bounding box is:" << std::endl;
  std::cout << "x: " << bBox.xmin << " " << bBox.xmax << std::endl;
  std::cout << "x: " << bBox.ymin << " " << bBox.ymax << std::endl;
  std::cout << "x: " << bBox.zmin << " " << bBox.zmax << std::endl;

  std::array<unsigned char, 3> loading_color;
  while (settings.running) {
    handle_events(cam, settings);

    if (cam.viewChanged() || cam.isRendering()) {
      if (not cam.isRendering()) cam.beginRendering();
      cam.renderlines(renderer, settings);

      loading_color[0] = 240;
      loading_color[1] = 0;
      loading_color[2] = 0;
    } else {
      loading_color[0] = 0;
      loading_color[1] = 240;
      loading_color[2] = 0;
    }

    // Draw red indicator rectangle for if it's rendering, otherwise green
    SDL_Rect rect;
    rect.x = 10;
    rect.y = 10;
    rect.w = 10;
    rect.h = 10;
    SDL_SetRenderDrawColor(renderer.rend(), loading_color[0], loading_color[1], loading_color[2], 240);
    SDL_RenderDrawRect(renderer.rend(), &rect);

    renderer.present();
  }

  err = openmc_finalize();
  if (err) openmc::fatal_error(openmc_err_msg);
}
