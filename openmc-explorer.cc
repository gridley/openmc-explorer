#include <SDL2/SDL.h>

#include <openmc/error.h>
#include <openmc/geometry.h>
#include <openmc/math_functions.h>
#include <openmc/particle.h>
#include <openmc/position.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

/* TODO
  - Use timers to steady the scrolling rate
  - Make TraceableGeometry into an abstract base class
    that defines how Monte Carlo codes should integrate
    with this program
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
  unsigned n_lines_to_draw = 1;
  
  int window_size_x;
  int window_size_y;
  unsigned n_threads;

  ExplorerSettings();
};
ExplorerSettings::ExplorerSettings() : running(true), delta_angle(5.0), delta_scoot(5.),
  n_lines_to_draw(5), window_size_x(ExplorerSettings::default_win_size_x),
  window_size_y(ExplorerSettings::default_win_size_y), n_threads(1),
  use_fullscreen(false) {}

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

class Camera
{
  openmc::Position pos;
  openmc::Direction dir;

  // whether to update the rendering
  bool viewChanged_;

  // Amount of roll to the camera (radians)
  double roll;

  // Field of view. Initially set to values close to the human eye
  double fov_x;
  double fov_y;

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

    void toggleFullscreen() { window.toggleFullscreen(); }

    void rotatePhi(double angle_d);
    void rotateMu(double angle_d);
    void increment_roll(double angle_d);
    void increment_fov_horiz(double angle_d);
    void increment_fov_vert(double angle_d);

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

Camera::Camera(double x, double y, double z, double phi, double mu,
    WindowManager& win) : dir(x, y, z),
  fov_x(135.0*M_PI/180.),
  fov_y(75.0*M_PI/180.),
  roll(0.0),
  window(win)
{
  double phi_rad = phi * M_PI/180.0;
  double mu_rad = mu * M_PI/180.0;
  dir.x = cos(phi_rad) * sin(mu_rad);
  dir.y = sin(phi_rad) * sin(mu_rad);
  dir.z = cos(mu_rad);
}
Camera::Camera(WindowManager& win) : Camera(100.0, 100.0, 100.0, 225.0, 135.0, win) {}

void Camera::rotatePhi(double angle_d) {
  double angle_r = angle_d * M_PI/180.0;
  dir = openmc::rotate_angle(dir, 0.0, &angle_r, nullptr);
  viewChanged_ = true;
}

void Camera::rotateMu(double angle_d) {
  double phi_zero = 0.0;
  double angle_r = angle_d * M_PI/180.;
  dir = openmc::rotate_angle(dir, angle_r, &phi_zero, nullptr);
  viewChanged_ = true;
}
bool Camera::viewChanged() {
  bool result = viewChanged_;
  viewChanged_ = false;
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
  double half_mu = fov_y/2.;
  double half_phi = fov_x/2.;
  double d_mu = fov_y / sett.window_size_y;
  double d_phi = fov_x / sett.window_size_x;

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

    // Loop over lines per thread to draw
    for (unsigned l=0; l<sett.n_lines_to_draw; ++l) {

      int vert_indx = n * sett.n_lines_to_draw + l + line_index;
      if (vert_indx >= sett.window_size_y) {
        isRendering_ = false;
        break;
      }

      // Loop over horizontal pixels
      for (unsigned horiz_indx=0; horiz_indx<sett.window_size_x; horiz_indx++) {

        // Roll is not currently handled... shouldn't be too hard to add?
        double phi = d_phi*horiz_indx - half_phi;
        double mu = half_mu - d_mu*vert_indx;
        thisdir = openmc::rotate_angle(dir, mu, &phi, nullptr);

        // Trace ray through geometry
        thisbank.r = pos;
        thisbank.u = thisdir;
        part.from_source(&thisbank);

        int cellid = openmc::find_cell(&part, false);
        auto dist = openmc::distance_to_boundary(&part);
        bool hitsomething;
        if (dist.distance == INFINITY) hitsomething = false;
        else hitsomething = true;

        double value = 240.*(1-exp(dist.distance));
        auto c = static_cast<unsigned char>(value);
        if (hitsomething) {
          SDL_SetRenderDrawColor(rend.rend(), c, c, c, 255);
        }
        else
          SDL_SetRenderDrawColor(rend.rend(), 0, 0, 0, 0);

        SDL_RenderDrawPoint(rend.rend(), horiz_indx, vert_indx);
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

          case SDLK_F11:
            cam.toggleFullscreen();
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
  Camera cam(window);

  // Set up openmc
  int err;
  err = openmc_init(argc, argv, nullptr);
  if (err == -1) { return 0; }
  else if (err) { openmc::fatal_error(openmc_err_msg); }

  while (settings.running) {
    handle_events(cam, settings);

    if (cam.viewChanged() || cam.isRendering())
    {
      if (not cam.isRendering()) cam.beginRendering();
      cam.renderlines(renderer, settings);
    }

    renderer.present();
  }

  err = openmc_finalize();
  if (err) openmc::fatal_error(openmc_err_msg);
}
